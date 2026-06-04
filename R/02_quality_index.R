# =============================================================================
#  R/02_quality_index.R — MODULE 2: ANNUAL QUALITY INDEX
#
#  Computes a per-year composite Quality Index (QI) from four components:
#   coverage   = fraction of AOI with valid Landsat observations
#   confidence = spatial mean of max RF probability
#   bal_acc    = Balanced Accuracy from accuracy CSV
#   s1_bonus   = binary: Sentinel-1 available (post-S1_START_YEAR)
#
#  QI in [0,1]. Years with QI < QI_LOW_THRESHOLD flagged as low-quality
#  and shown as grey bands in Module 3 plots and down-weighted in OLS.
#
#  KNOWN LIMITATION:
#  coverage is always 1.0 because GEE encodes NoData as class=0, not NA.
#  The QI is therefore driven mainly by BA and S1_bonus. Future fix:
#  use n_scenes_per_season from GEE metadata instead of fraction non-NA.
#
#  INPUT (config: GEE_DIR, QI_WEIGHTS, QI_LOW_THRESHOLD, S1_START_YEAR)
#    land_cover_YEAR.tif    all available years
#    confidence_YEAR.tif    all available years
#    accuracy_YEAR.csv      all available years
#
#  OUTPUT (config: OUT_QUALITY)
#    quality_index.csv      per-year QI + component values
#    quality_index.png      time series with stacked component contributions
# =============================================================================
if (!exists("OUTPUT_ROOT")) source("config.R")
source("R/00_utils.R")
.check_packages(c("terra","dplyr","tidyr","ggplot2","stringr","scales"))

section("Module 2 — Annual Quality Index")

lc  <- find_files(GEE_DIR, "^land_cover_\\d{4}\\.tif$")
cf  <- find_files(GEE_DIR, "^confidence_\\d{4}\\.tif$")
acc <- find_files(GEE_DIR, "^accuracy_\\d{4}\\.csv$")
years <- sort(Reduce(intersect, list(lc$years, cf$years, acc$years)))
if (length(years) == 0)
  stop("No years have all three required files (land_cover, confidence, accuracy)")
message(.ts(), sprintf("  %d years: %d\u2013%d", length(years), min(years), max(years)))

lc$files  <- lc$files[lc$years   %in% years]
cf$files  <- cf$files[cf$years   %in% years]
acc$files <- acc$files[acc$years %in% years]

qi_rows <- lapply(seq_along(years), function(i) {
  yr       <- years[i]
  coverage <- mean(!is.na(terra::values(terra::rast(lc$files[i]), dataframe = FALSE)))
  mean_conf <- mean(terra::values(terra::rast(cf$files[i]) / 100,
                                  dataframe = FALSE), na.rm = TRUE)
  acc_df   <- read.csv(acc$files[i], stringsAsFactors = FALSE)
  ba_row   <- acc_df[acc_df$final_code == -1, ]
  ba <- if (nrow(ba_row) > 0 && "balanced_accuracy" %in% names(ba_row))
          as.numeric(ba_row$balanced_accuracy[1]) else NA_real_
  data.frame(year = yr, coverage = coverage, mean_confidence = mean_conf,
             balanced_accuracy = ba, s1_available = as.integer(yr >= S1_START_YEAR))
})
qi_df <- dplyr::bind_rows(qi_rows)

norm01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1]) / diff(r)
}
qi_df <- qi_df |>
  dplyr::mutate(
    coverage_norm   = norm01(coverage),
    confidence_norm = norm01(mean_confidence),
    ba_norm         = norm01(balanced_accuracy),
    quality_index   = QI_WEIGHTS$coverage   * coverage_norm +
                      QI_WEIGHTS$confidence * confidence_norm +
                      QI_WEIGHTS$bal_acc    * ba_norm +
                      QI_WEIGHTS$s1_bonus   * s1_available,
    low_quality = quality_index < QI_LOW_THRESHOLD
  )

out_csv <- file.path(OUT_QUALITY, "quality_index.csv")
write.csv(qi_df, out_csv, row.names = FALSE)
n_lq <- sum(qi_df$low_quality, na.rm = TRUE)
message(.ts(), sprintf("  Mean QI = %.2f | Low-quality years: %d/%d",
                        mean(qi_df$quality_index, na.rm = TRUE), n_lq, nrow(qi_df)))
if (n_lq > 0)
  message(.ts(), "  Low-quality years: ",
          paste(qi_df$year[qi_df$low_quality], collapse = ", "))

comp_labs <- c(coverage_norm   = sprintf("Coverage (w=%.0f%%)", QI_WEIGHTS$coverage*100),
               confidence_norm = sprintf("Mean confidence (w=%.0f%%)", QI_WEIGHTS$confidence*100),
               ba_norm         = sprintf("Balanced Accuracy (w=%.0f%%)", QI_WEIGHTS$bal_acc*100),
               s1_available    = sprintf("S1 bonus (w=%.0f%%)", QI_WEIGHTS$s1_bonus*100))
comp_w <- setNames(c(QI_WEIGHTS$coverage, QI_WEIGHTS$confidence,
                     QI_WEIGHTS$bal_acc,  QI_WEIGHTS$s1_bonus),
                   names(comp_labs))

qi_long <- qi_df |>
  dplyr::select(year, coverage_norm, confidence_norm, ba_norm, s1_available) |>
  tidyr::pivot_longer(-year, names_to = "component", values_to = "value") |>
  dplyr::mutate(component    = comp_labs[component],
                contribution = value * comp_w[match(component, comp_labs)])

p <- ggplot2::ggplot() +
  ggplot2::geom_rect(data = qi_df[qi_df$low_quality, ],
    ggplot2::aes(xmin = year - 0.5, xmax = year + 0.5, ymin = 0, ymax = 1),
    fill = "#FFE0E0", alpha = 0.5) +
  ggplot2::geom_bar(data = qi_long,
    ggplot2::aes(x = year, y = contribution, fill = component),
    stat = "identity") +
  ggplot2::geom_line(data = qi_df,
    ggplot2::aes(x = year, y = quality_index), linewidth = 1.2, colour = "#1F4E79") +
  ggplot2::geom_point(data = qi_df,
    ggplot2::aes(x = year, y = quality_index), size = 2, colour = "#1F4E79") +
  ggplot2::geom_hline(yintercept = QI_LOW_THRESHOLD,
    linetype = "dashed", colour = "#C00000", linewidth = 0.7) +
  ggplot2::scale_fill_brewer(palette = "Set2", name = "Component") +
  ggplot2::scale_x_continuous(breaks = seq(min(years), max(years), 5)) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(title = "Annual Classification Quality Index",
    subtitle = sprintf("Red band = low quality (QI < %.1f) | Navy line = QI composite",
                        QI_LOW_THRESHOLD),
    x = "Year", y = "Quality Index [0, 1]") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right",
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p, file.path(OUT_QUALITY, "quality_index.png"), 14, 6)
message(.ts(), sprintf("  Module 2 complete — outputs in: %s", OUT_QUALITY))
