# =============================================================================
#  R/03_change_analysis.R — MODULE 3: LAND COVER CHANGE ANALYSIS  v2.3
#
#  FIXES IN v2.3 (on top of v2.2):
#  1. pix_km2 unit bug: GEE exports EPSG:4326, terra::res() returns degrees.
#     prod(res)/1e6 gave ~7e-14 instead of ~0.0009 km². Fixed by detecting
#     CRS type and converting degrees to metres via latitude correction.
#  2. Alluvial plot: terra::values() returns a 1-column matrix named after the
#     layer (e.g. "sm_2008"). as.vector() wrapper ensures columns are named
#     "first"/"mid"/"last" in the data.frame, preventing "object not found".
#
#  FIXES IN v2.2:
#  - STREAM B: class=0 (GEE NoData mask) excluded from change comparison.
#    Up to 86% false positives in historical year pairs before this fix.
#    New output columns: n_nodata_excluded, pct_valid_overlap, etc.
#  - CONF_THRESHOLD lowered to 50 in config.R (only meaningful with calibrated).
#
#  DESCRIPTION
#  Two change detection streams:
#
#  STREAM A (primary product — smoothed maps, 100% coverage):
#    Confidence-weighted temporal smoothing (window ±h years) produces
#    complete, NoData-free maps for every year. Area statistics, OLS trends,
#    transition matrix, and 7 standard plots are derived from these maps.
#
#  STREAM B (uncertainty diagnostic — raw maps, valid-overlap pixels only):
#    Pixel-level change detection with confidence gating. Reports what fraction
#    of apparent changes exceeds the confidence threshold in BOTH years.
#    Meaningful only with calibrated confidence from Module 1.
#
#  INPUT  (config: GEE_DIR, OUT_CALIBRATION, OUT_QUALITY, etc.)
#    land_cover_YEAR.tif         all processed years
#    confidence_YEAR.tif         all processed years
#    calibrated_confidence_YEAR.tif  (optional, from Module 1)
#    class_lookup_table.csv
#    quality_index.csv           (optional, from Module 2)
#
#  OUTPUT (config: OUT_CHANGE)
#    smoothed_YEAR.tif (per year)   confidence-weighted smoothed map (100% coverage)
#    area_by_class_year.csv         area (km²) per class per year
#    trend_by_class.csv             OLS slope, CI, R², p-value per class
#    transition_{Y1}_{Y2}.csv       transition matrix (km²) first-last
#    change_first_last.tif          encoded transition raster
#    stable_mask.tif                1 = stable across full series
#    uncertainty_reliability.csv    Stream B: overlap, fraction_certain, etc.
#    area_timeseries.png, stacked_area.png, trends.png, transition_heatmap.png
#    uncertainty_reliability.png, alluvial_transitions.png, pct_change_by_class.png
# =============================================================================
if (!exists("OUTPUT_ROOT")) source("config.R")
source("R/00_utils.R")
.check_packages(c("terra","dplyr","tidyr","ggplot2","scales","stringr","RColorBrewer"))

section("Module 3 — Land Cover Change Analysis (v2.3)")

lc_f  <- find_files(GEE_DIR, "^land_cover_\\d{4}\\.tif$")
cf_f  <- find_files(GEE_DIR, "^confidence_\\d{4}\\.tif$")
years <- sort(intersect(lc_f$years, cf_f$years))
if (length(years) == 0) stop("No years with both land_cover and confidence rasters")
n_yr <- length(years)
lc_f$files <- lc_f$files[lc_f$years %in% years]
cf_f$files <- cf_f$files[cf_f$years %in% years]
lookup <- read.csv(LOOKUP_CSV, stringsAsFactors = FALSE) |>
  dplyr::rename(class_code = final_code, class_name = class_name)
message(.ts(), sprintf("  %d years: %d\u2013%d", n_yr, min(years), max(years)))

cal_files <- if (USE_CALIBRATED_CONF) {
  sapply(years, function(yr) {
    f <- file.path(OUT_CALIBRATION, sprintf("calibrated_confidence_%d.tif", yr))
    if (file.exists(f)) f else NA_character_
  })
} else rep(NA_character_, n_yr)
message(.ts(), sprintf("  Calibrated confidence: %d/%d years", sum(!is.na(cal_files)), n_yr))

qi_path <- file.path(OUT_QUALITY, "quality_index.csv")
qi_df <- if (file.exists(qi_path)) {
  message(.ts(), "  Quality index loaded"); read.csv(qi_path, stringsAsFactors = FALSE)
} else {
  message(.ts(), "  quality_index.csv not found — run Module 2 first")
  data.frame(year = years, quality_index = NA_real_, low_quality = FALSE)
}

ref <- terra::rast(lc_f$files[1])

# FIX v2.3: pix_km2 — geographic (degrees) vs projected (metres) CRS detection.
# GEE exports EPSG:4326. terra::res() returns degrees, NOT metres.
# prod(res)/1e6 gives ~7e-14 instead of ~0.0009 km².
pix_km2 <- if (terra::is.lonlat(ref)) {
  lat_mid      <- (terra::ymax(ref) + terra::ymin(ref)) / 2
  m_per_deg_x  <- 111320 * cos(lat_mid * pi / 180)
  m_per_deg_y  <- 110540
  (terra::res(ref)[1] * m_per_deg_x) * (terra::res(ref)[2] * m_per_deg_y) / 1e6
} else {
  prod(terra::res(ref)) / 1e6
}
message(.ts(), sprintf("  CRS: %s | Pixel area: %.6f km²",
                        ifelse(terra::is.lonlat(ref), "geographic (degrees)", "projected (metres)"),
                        pix_km2))

lc_stk  <- load_stack(lc_f$files, ref, names_prefix = "lc",   years = years)
cf_stk  <- load_stack(cf_f$files, ref, names_prefix = "conf", years = years)
best_cf <- terra::rast(lapply(seq_along(years), function(i) {
  if (!is.na(cal_files[i])) terra::rast(cal_files[i]) * 100 else cf_stk[[i]]
}))
names(best_cf) <- paste0("bcf_", years)
all_cls <- sort(unique(lookup$class_code))

# =============================================================================
# STREAM A: Confidence-weighted temporal smoothing
# =============================================================================
section("Stream A — Confidence-weighted smoothing")
sm_list <- vector("list", n_yr)
for (i in seq_along(years)) {
  yr  <- years[i]
  win <- seq(max(1, i - SMOOTH_HALF_WINDOW), min(n_yr, i + SMOOTH_HALF_WINDOW))
  wt_stk <- terra::rast(lapply(all_cls, function(cls)
    Reduce("+", lapply(win, function(s)
      terra::ifel(lc_stk[[s]] == cls, best_cf[[s]] / 100, 0)))))
  sm_cls <- terra::classify(terra::which.max(wt_stk), cbind(seq_along(all_cls), all_cls))
  names(sm_cls) <- paste0("sm_", yr)
  sm_list[[i]]  <- sm_cls
  terra::writeRaster(sm_cls, file.path(OUT_CHANGE, sprintf("smoothed_%d.tif", yr)),
                     datatype = "INT2S", overwrite = TRUE)
  message(.ts(), sprintf("  %d smoothed [window %d\u2013%d | %s]",
    yr, years[min(win)], years[max(win)],
    ifelse(!is.na(cal_files[i]), "calibrated", "raw conf")))
}
sm_stk <- terra::rast(sm_list); names(sm_stk) <- paste0("sm_", years)

# Area statistics (terra::freq excludes NAs by default)
area_df <- dplyr::bind_rows(lapply(seq_along(years), function(i) {
  terra::freq(sm_stk[[i]]) |> as.data.frame() |>
    dplyr::rename(class_code = value, n_pixels = count) |>
    dplyr::mutate(year = years[i], area_km2 = n_pixels * pix_km2) |>
    dplyr::left_join(lookup, by = "class_code") |>
    dplyr::select(year, class_code, class_name, n_pixels, area_km2)
})) |>
  dplyr::group_by(year) |>
  dplyr::mutate(area_pct = area_km2 / sum(area_km2) * 100) |>
  dplyr::ungroup() |>
  dplyr::left_join(qi_df |> dplyr::select(year, quality_index, low_quality), by = "year") |>
  dplyr::arrange(year, class_code)
write.csv(area_df, file.path(OUT_CHANGE, "area_by_class_year.csv"), row.names = FALSE)
message(.ts(), sprintf("  Area table: %d rows | total AOI: %.1f km²",
                        nrow(area_df),
                        sum(area_df$area_km2[area_df$year == years[1]], na.rm = TRUE)))

# QI-weighted OLS trend
trend_df <- area_df |>
  dplyr::group_by(class_code, class_name) |>
  dplyr::filter(dplyr::n() >= 3) |>
  dplyr::do({
    d <- .; w <- ifelse(is.na(d$quality_index), 1, d$quality_index)
    m <- lm(area_km2 ~ year, data = d, weights = w); s <- summary(m)
    ci <- tryCatch(confint(m, "year"), error = function(e) matrix(c(NA, NA), 1))
    data.frame(slope_km2_yr = coef(m)["year"], slope_ci_lo = ci[1], slope_ci_hi = ci[2],
               r_squared = s$r.squared, p_value = s$coefficients["year", "Pr(>|t|)"],
               n_years = nrow(d), n_highquality = sum(!d$low_quality, na.rm = TRUE),
               area_first = d$area_km2[which.min(d$year)],
               area_last  = d$area_km2[which.max(d$year)])
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(significant = p_value < 0.05,
                pct_change_total = (area_last - area_first) / area_first * 100) |>
  dplyr::arrange(p_value)
write.csv(trend_df, file.path(OUT_CHANGE, "trend_by_class.csv"), row.names = FALSE)

cls_ord <- area_df |> dplyr::group_by(class_name) |>
  dplyr::summarise(m = mean(area_km2), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(m)) |> dplyr::pull(class_name)
col_map <- get_class_colours(cls_ord, CLASS_COLOURS)
area_df$class_name <- factor(area_df$class_name, levels = cls_ord)
lq_band <- lq_band_layer(qi_df)

# Plot 1: Area time series
p1 <- ggplot2::ggplot(area_df,
    ggplot2::aes(x = year, y = area_km2, colour = class_name, group = class_name)) +
  lq_band + ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_colour_manual(values = col_map, name = "Class") +
  ggplot2::scale_x_continuous(breaks = seq(min(years), max(years), 5)) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(title = "Land cover area time series",
    subtitle = sprintf("Smoothed | Grey = low-quality years | %d\u2013%d", min(years), max(years)),
    x = "Year", y = "Area (km\u00b2)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right", panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p1, file.path(OUT_CHANGE, "area_timeseries.png"), 14, 7)

# Plot 2: Stacked composition
p2 <- ggplot2::ggplot(area_df, ggplot2::aes(x = year, y = area_pct, fill = class_name)) +
  lq_band + ggplot2::geom_area(colour = "white", linewidth = 0.2) +
  ggplot2::scale_fill_manual(values = col_map, name = "Class") +
  ggplot2::scale_x_continuous(breaks = seq(min(years), max(years), 5)) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggplot2::labs(title = "Land cover composition", x = "Year", y = "% of AOI") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right", panel.grid = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p2, file.path(OUT_CHANGE, "stacked_area.png"), 14, 7)

# Plot 3: Trend slopes
p3 <- ggplot2::ggplot(
    trend_df |> dplyr::mutate(class_name = factor(class_name,
                                levels = class_name[order(slope_km2_yr)])),
    ggplot2::aes(x = slope_km2_yr, y = class_name, colour = significant,
                 xmin = slope_ci_lo, xmax = slope_ci_hi)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_errorbarh(height = 0.35, linewidth = 0.8) + ggplot2::geom_point(size = 4) +
  ggplot2::scale_colour_manual(values = c("TRUE" = "#C00000", "FALSE" = "#595959"),
    labels = c("TRUE" = "p<0.05 (*)", "FALSE" = "n.s."), name = "") +
  ggplot2::labs(title = "Linear trend \u2014 QI-weighted OLS",
    subtitle = "Bars = 95% CI | Weights = annual Quality Index",
    x = "Slope (km\u00b2/yr)", y = NULL) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p3, file.path(OUT_CHANGE, "trends.png"), 12, 6)

# Transition matrix + stable mask
sm1 <- sm_stk[[1]]; smN <- sm_stk[[n_yr]]
terra::writeRaster(sm1 * 10000L + smN, file.path(OUT_CHANGE, "change_first_last.tif"),
                   datatype = "INT4S", overwrite = TRUE)
terra::writeRaster(terra::ifel(sm1 == smN, 1L, 0L), file.path(OUT_CHANGE, "stable_mask.tif"),
                   datatype = "INT1U", overwrite = TRUE)
pct_stable <- mean(terra::values(terra::ifel(sm1 == smN, 1L, 0L)), na.rm = TRUE) * 100
message(.ts(), sprintf("  Stable pixels (smoothed, %d\u2013%d): %.1f%%", years[1], years[n_yr], pct_stable))

trans_df <- as.data.frame(
    table(from = terra::values(sm1, dataframe = FALSE),
          to   = terra::values(smN, dataframe = FALSE)) * pix_km2) |>
  dplyr::rename(from_code = from, to_code = to, area_km2 = Freq) |>
  dplyr::mutate(dplyr::across(c(from_code, to_code), ~as.integer(as.character(.)))) |>
  dplyr::left_join(lookup |> dplyr::rename(from_code = class_code, from_name = class_name),
                   by = "from_code") |>
  dplyr::left_join(lookup |> dplyr::rename(to_code = class_code, to_name = class_name),
                   by = "to_code") |>
  dplyr::filter(!is.na(from_name), !is.na(to_name))
write.csv(trans_df,
  file.path(OUT_CHANGE, sprintf("transition_%d_%d.csv", years[1], years[n_yr])),
  row.names = FALSE)

co <- lookup |>
  dplyr::filter(class_code %in% unique(c(trans_df$from_code, trans_df$to_code))) |>
  dplyr::arrange(class_code) |> dplyr::pull(class_name)
trans_df$from_name <- factor(trans_df$from_name, levels = rev(co))
trans_df$to_name   <- factor(trans_df$to_name,   levels = co)

p4 <- ggplot2::ggplot(trans_df, ggplot2::aes(x = to_name, y = from_name, fill = area_km2)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(area_km2 > 0.01, sprintf("%.1f", area_km2), "")),
    size = 3, colour = "white", fontface = "bold") +
  ggplot2::scale_fill_gradient(low = "#FFFDE7", high = "#1F4E79", name = "Area (km\u00b2)",
    trans = "sqrt", labels = scales::comma) +
  ggplot2::labs(title = sprintf("Transition matrix: %d \u2192 %d (smoothed)", years[1], years[n_yr]),
    subtitle = "Values in km\u00b2 | Diagonal = stable pixels",
    x = sprintf("Destination (%d)", years[n_yr]), y = sprintf("Origin (%d)", years[1])) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                 panel.grid  = ggplot2::element_blank(),
                 plot.title  = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p4, file.path(OUT_CHANGE, "transition_heatmap.png"), 10, 8)

# =============================================================================
# STREAM B: Uncertainty-gated change detection (v2.2 NoData fix)
# =============================================================================
# FIX v2.2: class=0 in raw land_cover = GEE NoData mask (pixels with no valid
# Landsat observations). Before this fix, 0->class transitions were counted as
# "change", causing up to 86% false positives in early year pairs.
# Fix: compare only pixels with class>0 in BOTH years (valid_both mask).
section("Stream B — Uncertainty-gated change (v2.2 NoData-masked)")
rel_df <- dplyr::bind_rows(lapply(seq_len(n_yr - 1), function(i) {
  lc_t  <- lc_stk[[i]];   lc_t1 <- lc_stk[[i + 1]]
  bc_t  <- best_cf[[i]];  bc_t1 <- best_cf[[i + 1]]
  valid_both <- (lc_t > 0) & (lc_t1 > 0)
  ct <- terra::ifel(
    !valid_both, NA_integer_,
    terra::ifel(lc_t == lc_t1, 0L,
      terra::ifel((bc_t >= CONF_THRESHOLD) & (bc_t1 >= CONF_THRESHOLD), 1L, 2L)))
  v   <- terra::values(ct, dataframe = FALSE)
  ns  <- sum(v == 0L, na.rm = TRUE); nc <- sum(v == 1L, na.rm = TRUE)
  nu  <- sum(v == 2L, na.rm = TRUE); nnd <- sum(is.na(v))
  n_valid <- ns + nc + nu
  data.frame(
    year_from = years[i], year_to = years[i + 1],
    conf_t  = ifelse(!is.na(cal_files[i]),     "calibrated", "raw"),
    conf_t1 = ifelse(!is.na(cal_files[i + 1]), "calibrated", "raw"),
    n_nodata_excluded = nnd, n_valid_overlap = n_valid,
    pct_valid_overlap = n_valid / (n_valid + nnd) * 100,
    n_stable = ns, n_certain = nc, n_uncertain = nu,
    area_nodata_km2    = nnd     * pix_km2,
    area_valid_km2     = n_valid * pix_km2,
    area_stable_km2    = ns      * pix_km2,
    area_certain_km2   = nc      * pix_km2,
    area_uncertain_km2 = nu      * pix_km2,
    fraction_certain   = ifelse(nc + nu > 0, nc / (nc + nu), NA_real_),
    pct_valid_changed  = ifelse(n_valid > 0, (nc + nu) / n_valid * 100, NA_real_))
}))
write.csv(rel_df, file.path(OUT_CHANGE, "uncertainty_reliability.csv"), row.names = FALSE)
message(.ts(), sprintf("  Mean fraction certain (valid overlap): %.1f%%",
                        mean(rel_df$fraction_certain, na.rm = TRUE) * 100))
message(.ts(), sprintf("  Mean valid overlap: %.1f%% of AOI",
                        mean(rel_df$pct_valid_overlap, na.rm = TRUE)))
message(.ts(), sprintf("  Year pairs with <50%% valid overlap: %s",
  paste(sprintf("%d\u2192%d", rel_df$year_from, rel_df$year_to)[rel_df$pct_valid_overlap < 50],
        collapse = ", ")))

rel_long <- rel_df |>
  dplyr::mutate(pair = paste0(year_from, "\u2192", year_to),
                both_cal = (conf_t == "calibrated") & (conf_t1 == "calibrated")) |>
  dplyr::select(pair, both_cal, area_stable_km2, area_certain_km2,
                area_uncertain_km2, area_nodata_km2) |>
  tidyr::pivot_longer(c(area_stable_km2, area_certain_km2, area_uncertain_km2, area_nodata_km2),
                      names_to = "cat", values_to = "area") |>
  dplyr::mutate(cat = factor(cat,
    levels = c("area_stable_km2","area_certain_km2","area_uncertain_km2","area_nodata_km2"),
    labels = c("Stable",
               sprintf("Certain change (\u2265%d%%)", CONF_THRESHOLD),
               sprintf("Uncertain change (<%d%%)",   CONF_THRESHOLD),
               "NoData in \u22651 year (excluded)")))

fill_vals_b <- setNames(c("#2c7bb6","#d7191c","#fdae61","#CCCCCC"),
  c("Stable",
    sprintf("Certain change (\u2265%d%%)", CONF_THRESHOLD),
    sprintf("Uncertain change (<%d%%)",   CONF_THRESHOLD),
    "NoData in \u22651 year (excluded)"))

p5 <- ggplot2::ggplot(rel_long, ggplot2::aes(x = pair, y = area, fill = cat)) +
  ggplot2::geom_bar(stat = "identity", width = 0.85) +
  ggplot2::scale_fill_manual(values = fill_vals_b, name = NULL) +
  ggplot2::geom_point(
    data = rel_df |> dplyr::mutate(pair = paste0(year_from, "\u2192", year_to),
      both_cal = (conf_t == "calibrated") & (conf_t1 == "calibrated")),
    ggplot2::aes(x = pair, y = -max(rel_long$area, na.rm = TRUE) * 0.05, shape = both_cal),
    colour = "#1F4E79", size = 2.5, inherit.aes = FALSE) +
  ggplot2::scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
    labels = c("TRUE" = "Both: calibrated", "FALSE" = "At least one: raw"),
    name = "Confidence source") +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(
    title   = sprintf("Land cover change reliability (threshold = %d%%)", CONF_THRESHOLD),
    subtitle = paste0("Grey = NoData in \u22651 year (excluded). \u25b2 = calibrated confidence. ",
                      "Percentages over valid overlap only."),
    x = "Year pair", y = "Area (km\u00b2)",
    caption = "v2.2 fix: class=0 excluded. v2.3 fix: area units corrected.") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom", legend.box = "vertical",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p5, file.path(OUT_CHANGE, "uncertainty_reliability.png"), 16, 8)

# =============================================================================
# PLOT 6: Alluvial (Sankey) — three time snapshots
# =============================================================================
yr_first <- years[1]; yr_mid <- years[which.min(abs(years - median(years)))]; yr_last <- years[n_yr]
message(.ts(), sprintf("  Snapshot years: %d / %d / %d", yr_first, yr_mid, yr_last))

if (requireNamespace("ggalluvial", quietly = TRUE)) {
  library(ggalluvial)
  # FIX v2.3: as.vector() required — terra::values() returns a 1-column matrix
  # named after the layer (e.g. "sm_2008"), causing "object 'mid' not found".
  vals_first <- as.vector(terra::values(sm_stk[[which(years == yr_first)]], dataframe = FALSE))
  vals_mid   <- as.vector(terra::values(sm_stk[[which(years == yr_mid)]],   dataframe = FALSE))
  vals_last  <- as.vector(terra::values(sm_stk[[which(years == yr_last)]],  dataframe = FALSE))
  alluvial_raw <- data.frame(first = vals_first, mid = vals_mid, last = vals_last) |>
    dplyr::filter(!is.na(first), !is.na(mid), !is.na(last)) |>
    dplyr::group_by(first, mid, last) |>
    dplyr::summarise(n_pixels = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(area_km2 = n_pixels * pix_km2) |>
    dplyr::left_join(lookup |> dplyr::rename(first = class_code, first_name = class_name), by = "first") |>
    dplyr::left_join(lookup |> dplyr::rename(mid   = class_code, mid_name   = class_name), by = "mid")   |>
    dplyr::left_join(lookup |> dplyr::rename(last  = class_code, last_name  = class_name), by = "last")  |>
    dplyr::filter(!is.na(first_name), !is.na(mid_name), !is.na(last_name))
  alluvial_long <- alluvial_raw |>
    dplyr::mutate(alluvium = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = c(first_name, mid_name, last_name),
                        names_to = "time_point", values_to = "class_name") |>
    dplyr::mutate(time_point = factor(time_point,
        levels = c("first_name","mid_name","last_name"),
        labels = c(as.character(yr_first), as.character(yr_mid), as.character(yr_last))),
      class_name = factor(class_name, levels = cls_ord))
  p6 <- ggplot2::ggplot(alluvial_long,
      ggplot2::aes(x = time_point, y = area_km2, alluvium = alluvium,
                   stratum = class_name, fill = class_name, label = class_name)) +
    ggalluvial::geom_flow(stat = "alluvium", aes.flow = "forward",
                          alpha = 0.55, colour = "white", linewidth = 0.2) +
    ggalluvial::geom_stratum(alpha = 0.85, colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(stat = "stratum", size = 3, fontface = "bold",
                       colour = "white", min.y = 5) +
    ggplot2::scale_fill_manual(values = col_map, name = "Land cover class") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title    = sprintf("Land cover transitions: %d \u2192 %d \u2192 %d", yr_first, yr_mid, yr_last),
      subtitle = "Band width proportional to area (km\u00b2) | Smoothed maps",
      x = "Year", y = "Area (km\u00b2)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "right", panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"),
                   axis.text.x = ggplot2::element_text(size = 12, face = "bold"))
  save_plot(p6, file.path(OUT_CHANGE, "alluvial_transitions.png"), 14, 9)
} else {
  message(.ts(), "  Alluvial skipped — install ggalluvial: install.packages('ggalluvial')")
}

# =============================================================================
# PLOT 7: Percentage change relative to first year
# =============================================================================
snapshot_years <- unique(c(yr_first, yr_mid, yr_last))
area_snap <- area_df |>
  dplyr::filter(year %in% snapshot_years) |>
  dplyr::select(year, class_name, area_km2) |>
  dplyr::mutate(snapshot = dplyr::case_when(
    year == yr_first ~ sprintf("Initial (%d)", yr_first),
    year == yr_mid   ~ sprintf("Mid (%d)",     yr_mid),
    year == yr_last  ~ sprintf("Final (%d)",   yr_last)))
baseline <- area_snap |>
  dplyr::filter(year == yr_first) |>
  dplyr::select(class_name, baseline_km2 = area_km2)
pct_change_df <- area_snap |>
  dplyr::left_join(baseline, by = "class_name") |>
  dplyr::filter(!is.na(baseline_km2), baseline_km2 > 0) |>
  dplyr::mutate(pct_change = (area_km2 - baseline_km2) / baseline_km2 * 100,
    snapshot = factor(snapshot, levels = c(sprintf("Initial (%d)", yr_first),
      sprintf("Mid (%d)", yr_mid), sprintf("Final (%d)", yr_last))),
    class_name = factor(class_name, levels = cls_ord))
write.csv(pct_change_df, file.path(OUT_CHANGE, "pct_change_by_class.csv"), row.names = FALSE)

snap_fill <- setNames(c("#74add1","#1F4E79"),
  c(sprintf("Mid (%d)", yr_mid), sprintf("Final (%d)", yr_last)))
p7 <- ggplot2::ggplot(pct_change_df |> dplyr::filter(year != yr_first),
    ggplot2::aes(x = class_name, y = pct_change, fill = snapshot)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.8) +
  ggplot2::geom_bar(stat = "identity",
                    position = ggplot2::position_dodge(width = 0.75), width = 0.65, alpha = 0.9) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f%%", pct_change),
                 vjust = ifelse(pct_change >= 0, -0.4, 1.2)),
    position = ggplot2::position_dodge(width = 0.75),
    size = 3, colour = "grey20", fontface = "bold") +
  ggplot2::scale_fill_manual(values = snap_fill, name = "Snapshot") +
  ggplot2::scale_y_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), x, "%"),
    expand = ggplot2::expansion(mult = c(0.15, 0.20))) +
  ggplot2::labs(
    title   = sprintf("Land cover change relative to initial area (%d)", yr_first),
    subtitle = sprintf("(area at T \u2212 area at %d) / area at %d \u00d7 100 | Smoothed maps",
                        yr_first, yr_first),
    x = NULL, y = sprintf("Change relative to %d (%%)", yr_first)) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    legend.position = "top", panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p7, file.path(OUT_CHANGE, "pct_change_by_class.png"), 13, 7)

message(.ts(), sprintf("  Module 3 complete — outputs in: %s", OUT_CHANGE))
