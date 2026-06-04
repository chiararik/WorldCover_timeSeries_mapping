# =============================================================================
#  R/04_vip_analysis.R — MODULE 4: VARIABLE IMPORTANCE OVER TIME
#
#  Aggregates Random Forest variable importance across all processed years
#  and produces comparative visualisations showing how the feature ranking
#  evolves before and after Sentinel-1 integration (2015).
#
#  SCIENTIFIC NOTE
#  Each year's RF is trained independently on that year's feature stack.
#  Variable importance reflects the discriminative value of each feature
#  for that year's spectral and climatic conditions.
#  Key finding: srtm_twi (Topographic Wetness Index) ranks in the top 3
#  for ALL 32 years, including pre-SAR years. S1_VV becomes top-ranked
#  from 2015 onward. srtm_slope is consistently last (flat floodplain).
#
#  INPUT (config: GEE_DIR, TOP_N_FEATURES, VIP_MIN_YEARS_PRESENT,
#                 FEATURE_GROUPS, GROUP_COLOURS, S1_START_YEAR)
#    variable_importance_YEAR.csv    one per processed year
#
#  OUTPUT (config: OUT_VIP)
#    vip_table.csv              long-format: year × feature × importance × rank
#    vip_top_features.csv       features ranked by mean importance across years
#    vip_heatmap.png            feature × year heatmap (magma palette)
#    vip_top10_timeseries.png   top-N importance time series per feature
#    vip_rank_stability.png     Spearman rank correlation between years
#    vip_pre_post_s1.png        delta importance before vs after S1
# =============================================================================
if (!exists("OUTPUT_ROOT")) source("config.R")
source("R/00_utils.R")
.check_packages(c("dplyr","tidyr","ggplot2","scales","stringr",
                  "RColorBrewer","viridis","purrr","forcats"))

section("Module 4 — Variable Importance Over Time")

vip_f <- find_files(GEE_DIR, "^variable_importance_\\d{4}\\.csv$")
if (length(vip_f$files) == 0)
  stop("No variable_importance_YEAR.csv files found in GEE_DIR")
message(.ts(), sprintf("  %d VIP files: %d\u2013%d",
                        length(vip_f$years), min(vip_f$years), max(vip_f$years)))

grp_map <- stack(FEATURE_GROUPS) |>
  as.data.frame() |>
  dplyr::rename(feature = values, group = ind) |>
  dplyr::mutate(feature = as.character(feature), group = as.character(group))

vip_all <- dplyr::bind_rows(lapply(seq_along(vip_f$files), function(i) {
  read.csv(vip_f$files[i], stringsAsFactors = FALSE) |>
    dplyr::select(rank, feature, importance, importance_norm) |>
    dplyr::mutate(year = vip_f$years[i], has_s1 = year >= S1_START_YEAR)
})) |>
  dplyr::left_join(grp_map, by = "feature") |>
  dplyr::mutate(group = ifelse(is.na(group), "Other", group))

write.csv(vip_all, file.path(OUT_VIP, "vip_table.csv"), row.names = FALSE)
message(.ts(), sprintf("  %d rows | %d unique features | %d years",
                        nrow(vip_all), dplyr::n_distinct(vip_all$feature),
                        dplyr::n_distinct(vip_all$year)))

vip_sum <- vip_all |>
  dplyr::group_by(feature, group) |>
  dplyr::summarise(n_years       = dplyr::n(),
                   mean_imp      = mean(importance,      na.rm = TRUE),
                   mean_imp_norm = mean(importance_norm, na.rm = TRUE),
                   mean_rank     = mean(rank,            na.rm = TRUE),
                   median_rank   = median(rank,          na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(mean_rank)
write.csv(vip_sum, file.path(OUT_VIP, "vip_top_features.csv"), row.names = FALSE)

message(.ts(), "  Top 10 features by mean rank:")
for (i in seq_len(min(10, nrow(vip_sum)))) {
  r <- vip_sum[i, ]
  message(sprintf("    %2d. %-25s [%-10s] mean_rank=%.1f present=%d yrs",
                  i, r$feature, r$group, r$mean_rank, r$n_years))
}

yrs         <- sort(unique(vip_all$year))
feats_keep  <- vip_sum |> dplyr::filter(n_years >= VIP_MIN_YEARS_PRESENT) |> dplyr::pull(feature)
vip_heat    <- vip_all |>
  dplyr::filter(feature %in% feats_keep) |>
  dplyr::mutate(feature = factor(feature, levels = rev(feats_keep)))

# Plot 1: Heatmap feature × year
p1 <- ggplot2::ggplot(vip_heat, ggplot2::aes(x = year, y = feature, fill = importance_norm)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = S1_START_YEAR - 0.5, linetype = "dashed",
                      colour = "white", linewidth = 1) +
  ggplot2::annotate("text", x = S1_START_YEAR, y = length(feats_keep) + 0.7,
    label = sprintf("S1 (%d)", S1_START_YEAR), colour = "white", size = 3, angle = 90, hjust = 0) +
  viridis::scale_fill_viridis(option = "magma", name = "Norm. importance", na.value = "grey90") +
  ggplot2::scale_x_continuous(breaks = seq(min(yrs), max(yrs), 5)) +
  ggplot2::labs(title = "RF variable importance over time",
    subtitle = sprintf("Sorted by mean rank | Dashed = S1 introduced (%d)", S1_START_YEAR),
    x = "Year", y = NULL) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8),
                 plot.title  = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p1, file.path(OUT_VIP, "vip_heatmap.png"), 14, 8)

# Plot 2: Top-N feature time series
top_feats <- head(feats_keep, TOP_N_FEATURES)
grp_col_map <- setNames(unlist(GROUP_COLOURS), names(GROUP_COLOURS))
vip_ts <- vip_all |>
  dplyr::filter(feature %in% top_feats) |>
  dplyr::left_join(dplyr::select(vip_sum, feature, group), by = "feature") |>
  dplyr::mutate(feature = factor(feature, levels = top_feats))

p2 <- ggplot2::ggplot(vip_ts, ggplot2::aes(x = year, y = importance_norm,
                                            colour = group, group = feature)) +
  ggplot2::geom_vline(xintercept = S1_START_YEAR, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 1.5) +
  ggplot2::facet_wrap(~feature, ncol = 5) +
  ggplot2::scale_colour_manual(values = grp_col_map, name = "Feature group") +
  ggplot2::scale_x_continuous(breaks = seq(min(yrs), max(yrs), 10)) +
  ggplot2::labs(title = sprintf("Top %d features: importance time series", TOP_N_FEATURES),
    subtitle = "Dashed = Sentinel-1 introduction (2015)", x = "Year", y = "Norm. importance") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position  = "bottom",
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p2, file.path(OUT_VIP, "vip_top10_timeseries.png"), 16, 9)

# Plot 3: Rank stability (Spearman rho between consecutive years)
stab_df <- dplyr::bind_rows(lapply(seq_len(length(yrs) - 1), function(i) {
  yr1_df <- vip_all |> dplyr::filter(year == yrs[i])
  yr2_df <- vip_all |> dplyr::filter(year == yrs[i+1])
  both   <- dplyr::inner_join(yr1_df, yr2_df, by = "feature", suffix = c("_1","_2"))
  rho    <- if (nrow(both) >= 3)
              cor(both$rank_1, both$rank_2, method = "spearman") else NA_real_
  data.frame(year_pair = sprintf("%d\u2192%d", yrs[i], yrs[i+1]),
             year_mid  = (yrs[i] + yrs[i+1]) / 2, rho = rho)
}))

p3 <- ggplot2::ggplot(stab_df, ggplot2::aes(x = year_mid, y = rho)) +
  ggplot2::geom_vline(xintercept = S1_START_YEAR, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_hline(yintercept = 0.7, linetype = "dotted", colour = "grey50") +
  ggplot2::geom_line(colour = "#1F4E79", linewidth = 0.9) +
  ggplot2::geom_point(colour = "#1F4E79", size = 2.5) +
  ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.01)) +
  ggplot2::labs(title = "Feature rank stability between consecutive years",
    subtitle = "Spearman \u03c1 of VIP rank vectors | Dashed = S1 (2015) | Dotted = \u03c1 = 0.7",
    x = "Year", y = "Spearman \u03c1") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p3, file.path(OUT_VIP, "vip_rank_stability.png"), 12, 5)

# Plot 4: Pre vs post S1 delta importance
pre_mean  <- vip_all |> dplyr::filter(!has_s1) |> dplyr::group_by(feature) |>
  dplyr::summarise(mean_pre = mean(importance_norm, na.rm = TRUE), .groups = "drop")
post_mean <- vip_all |> dplyr::filter(has_s1) |> dplyr::group_by(feature) |>
  dplyr::summarise(mean_post = mean(importance_norm, na.rm = TRUE), .groups = "drop")
delta_df <- dplyr::inner_join(pre_mean, post_mean, by = "feature") |>
  dplyr::left_join(grp_map, by = "feature") |>
  dplyr::mutate(group = ifelse(is.na(group), "Other", group),
                delta = mean_post - mean_pre) |>
  dplyr::arrange(delta) |>
  dplyr::mutate(feature = factor(feature, levels = feature))

p4 <- ggplot2::ggplot(delta_df, ggplot2::aes(x = delta, y = feature, fill = group)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.8) +
  ggplot2::geom_bar(stat = "identity", alpha = 0.85) +
  ggplot2::scale_fill_manual(values = grp_col_map, name = "Feature group") +
  ggplot2::labs(title = "Change in normalised importance: post-S1 vs pre-S1",
    subtitle = sprintf("Positive = feature gained importance after %d | Negative = displaced by S1",
                        S1_START_YEAR),
    x = "\u0394 (post-S1 mean \u2212 pre-S1 mean)", y = NULL) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1F4E79"))
save_plot(p4, file.path(OUT_VIP, "vip_pre_post_s1.png"), 12, 10)
message(.ts(), sprintf("  Module 4 complete — outputs in: %s", OUT_VIP))
