# =============================================================================
#  R/01_rf_calibration.R — MODULE 1: RF PROBABILITY CALIBRATION v2.1.2
#
#  Corrects systematic under-dispersion of Random Forest class probabilities
#  using isotonic regression (Pool-Adjacent-Violators algorithm).
#
#  BUG FIXES:
#  v2.1.1: Added norm_name() for space/underscore mismatch between prob_* band
#          names and class_name_actual column.
#  v2.1.2: norm_name() replaces ALL non-alphanumeric sequences with underscore
#          (handles /, -, . and space simultaneously). Fixes n_pos=0 for
#          "Grassland / Cropland" and "Non-vegetated / Bare".
#
#  INPUT  (config: CALIB_YEAR, GEE_DIR)
#    validation_points_YEAR.csv  per-point RF probabilities + correctness flag
#    multiprob_YEAR.tif          per-class probabilities raster (Byte 0-100)
#
#  OUTPUT (config: OUT_CALIBRATION)
#    calibrated_confidence_YEAR.tif    max calibrated probability (Float32 [0,1])
#    calibrated_uncertainty_YEAR.tif   normalised Shannon entropy (Float32 [0,1])
#    calibration_functions_YEAR.csv    isotonic lookup table per class
#    calibration_plot_YEAR.png         reliability diagrams + entropy violin
# =============================================================================
if (!exists("OUTPUT_ROOT")) source("config.R")
source("R/00_utils.R")
.check_packages(c("terra","dplyr","tidyr","ggplot2","jsonlite","scales","gridExtra"))

yr <- CALIB_YEAR
section(sprintf("Module 1 — RF Calibration | Year %d", yr))

# ── Load validation points ────────────────────────────────────────────────────
val_csv <- file.path(GEE_DIR, sprintf("validation_points_%d.csv", yr))
if (!file.exists(val_csv)) stop("validation_points not found: ", val_csv)

pts_raw   <- read.csv(val_csv, stringsAsFactors = FALSE)
prob_cols <- grep("^prob_", names(pts_raw), value = TRUE)
class_names <- sub("^prob_", "", prob_cols)
n_classes   <- length(prob_cols)
if (n_classes == 0)
  stop("No prob_* columns found — re-export validation_points with per-class probabilities from GEE")

coords <- do.call(rbind, lapply(pts_raw[[".geo"]], function(g) {
  j <- tryCatch(jsonlite::fromJSON(g), error = function(e) NULL)
  if (is.null(j)) return(c(NA_real_, NA_real_))
  as.numeric(j$coordinates[1:2])
}))
pts <- pts_raw |>
  dplyr::mutate(longitude = coords[,1], latitude = coords[,2],
                correct   = as.integer(correct),
                LC        = as.integer(LC),
                predicted = as.integer(predicted),
                dplyr::across(dplyr::all_of(prob_cols), as.numeric)) |>
  dplyr::filter(!is.na(correct), !is.na(longitude))

message(.ts(), sprintf("  %d validation points | %d classes | OA = %.1f%%",
                        nrow(pts), n_classes, mean(pts$correct) * 100))

# ── Class name normalisation (v2.1.2) ─────────────────────────────────────────
# GEE prob_* band names use underscores; class_name_actual uses spaces/slashes.
# gsub replaces ALL non-alphanumeric sequences with "_" for consistent matching.
norm_name <- function(x) gsub("[^A-Za-z0-9]+", "_", trimws(x))

# ── Isotonic regression per class ─────────────────────────────────────────────
calib_fns   <- list()
calib_table <- list()

for (cn in class_names) {
  col   <- paste0("prob_", cn)
  p_raw <- pts[[col]]
  y_c   <- as.integer(norm_name(pts$class_name_actual) == norm_name(cn))
  ord   <- order(p_raw)
  fit   <- isoreg(p_raw[ord], y_c[ord])
  p_cal <- approx(p_raw[ord], fit$yf, xout = p_raw, rule = 2, ties = "ordered")$y
  calib_fns[[cn]]   <- list(p_raw = p_raw[ord], p_cal = fit$yf)
  calib_table[[cn]] <- data.frame(class = cn, p_raw = p_raw, p_cal = p_cal)

  status <- if (sum(y_c) == 0) "(class absent in validation set — skipped)"
            else sprintf("mean_raw=%.3f  mean_cal=%.3f  n_pos=%d",
                         mean(p_raw), mean(p_cal, na.rm = TRUE), sum(y_c))
  message(.ts(), sprintf("  %-25s  %s", cn, status))
}

out_csv <- file.path(OUT_CALIBRATION, sprintf("calibration_functions_%d.csv", yr))
write.csv(do.call(rbind, calib_table), out_csv, row.names = FALSE)
message(.ts(), "  Saved: ", basename(out_csv))

# ── Apply calibration to validation points (entropy plot) ────────────────────
for (cn in class_names) {
  fn <- calib_fns[[cn]]
  pts[[paste0("pcal_", cn)]] <-
    approx(fn$p_raw, fn$p_cal, xout = pts[[paste0("prob_", cn)]],
           rule = 2, ties = "ordered")$y
}
pcal_cols <- paste0("pcal_", class_names)
rsums <- rowSums(pts[, pcal_cols], na.rm = TRUE)
rsums[rsums == 0] <- 1
for (cn in class_names)
  pts[[paste0("pcal_", cn)]] <- pts[[paste0("pcal_", cn)]] / rsums

eps   <- 1e-10
pts$H_raw <- apply(pts[, prob_cols], 1,
  function(p) { p <- pmax(p/100, eps); -sum(p * log(p)) / log(n_classes) })
pts$H_cal <- apply(pts[, pcal_cols], 1,
  function(p) { p <- pmax(p, eps);    -sum(p * log(p)) / log(n_classes) })
message(.ts(), sprintf("  Entropy: raw [%.3f, %.3f]  cal [%.3f, %.3f]",
                        min(pts$H_raw), max(pts$H_raw),
                        min(pts$H_cal), max(pts$H_cal)))

# ── Apply calibration to multiprob raster ─────────────────────────────────────
mp_path <- file.path(GEE_DIR, sprintf("multiprob_%d.tif", yr))
if (!file.exists(mp_path)) {
  message(.ts(), "  multiprob not found — skipping raster calibration")
} else {
  mp_stk <- terra::rast(mp_path) / 100  # Byte 0-100 → [0,1]
  cal_stk <- terra::rast(lapply(seq_along(class_names), function(j) {
    cn <- class_names[j]
    fn <- calib_fns[[cn]]
    p  <- mp_stk[[j]]
    terra::app(p, function(x) approx(fn$p_raw, fn$p_cal, xout = x, rule = 2, ties = "ordered")$y)
  }))
  # Renormalise
  row_sums <- sum(cal_stk, na.rm = TRUE)
  row_sums[row_sums == 0] <- 1
  cal_stk <- cal_stk / row_sums

  # Calibrated confidence = max calibrated class probability
  conf_cal <- max(cal_stk)
  terra::writeRaster(conf_cal,
    file.path(OUT_CALIBRATION, sprintf("calibrated_confidence_%d.tif", yr)),
    datatype = "FLT4S", overwrite = TRUE)

  # Calibrated uncertainty = normalised Shannon entropy
  unc_cal <- -sum(ifel(cal_stk > eps, cal_stk * log(cal_stk), 0)) / log(n_classes)
  terra::writeRaster(unc_cal,
    file.path(OUT_CALIBRATION, sprintf("calibrated_uncertainty_%d.tif", yr)),
    datatype = "FLT4S", overwrite = TRUE)

  message(.ts(), sprintf("  Calibrated rasters written for year %d", yr))
}

# ── Calibration plot ──────────────────────────────────────────────────────────
n_bins <- 10
rl_list <- lapply(class_names, function(cn) {
  p   <- pts[[paste0("prob_", cn)]] / 100
  y   <- as.integer(norm_name(pts$class_name_actual) == norm_name(cn))
  pc  <- pts[[paste0("pcal_", cn)]]
  cuts <- cut(p, breaks = seq(0, 1, 1/n_bins), include.lowest = TRUE)
  df <- tapply(y, cuts, mean, na.rm = TRUE)
  data.frame(class = cn, bin_mid = seq(0.05, 0.95, 0.1),
             obs_freq_raw = as.numeric(df),
             obs_freq_cal = as.numeric(tapply(y, cut(pc, breaks = seq(0, 1, 1/n_bins),
                                                     include.lowest = TRUE), mean, na.rm = TRUE)))
})
rl_df <- dplyr::bind_rows(rl_list)

p_rel <- ggplot2::ggplot(rl_df, ggplot2::aes(x = bin_mid)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_line(ggplot2::aes(y = obs_freq_raw), colour = "#d73027", linewidth = 0.8) +
  ggplot2::geom_point(ggplot2::aes(y = obs_freq_raw), colour = "#d73027", size = 2) +
  ggplot2::geom_line(ggplot2::aes(y = obs_freq_cal), colour = "#1a9641", linewidth = 0.8) +
  ggplot2::geom_point(ggplot2::aes(y = obs_freq_cal), colour = "#1a9641", size = 2) +
  ggplot2::facet_wrap(~class, nrow = 2) +
  ggplot2::labs(title = sprintf("Reliability diagrams — Year %d", yr),
    subtitle = "Red = raw RF | Green = calibrated (isotonic regression)",
    x = "Mean predicted probability", y = "Fraction correct") +
  ggplot2::coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  ggplot2::theme_bw(base_size = 10)

p_ent <- ggplot2::ggplot(
    tidyr::pivot_longer(pts[, c("H_raw","H_cal")],
                        everything(), names_to = "type", values_to = "H"),
    ggplot2::aes(x = type, y = H, fill = type)) +
  ggplot2::geom_violin(alpha = 0.7) +
  ggplot2::scale_fill_manual(values = c(H_raw = "#d73027", H_cal = "#1a9641"), guide = "none") +
  ggplot2::scale_x_discrete(labels = c(H_raw = "Raw", H_cal = "Calibrated")) +
  ggplot2::labs(title = "Entropy distribution", x = NULL,
                y = "Normalised Shannon entropy H") +
  ggplot2::theme_bw(base_size = 10)

g <- gridExtra::arrangeGrob(p_rel, p_ent, ncol = 2, widths = c(3, 1))
out_png <- file.path(OUT_CALIBRATION, sprintf("calibration_plot_%d.png", yr))
ggplot2::ggsave(out_png, g, width = 14, height = 7, dpi = PLOT_DPI)
message(.ts(), "  Saved: ", basename(out_png))
message(.ts(), sprintf("  Module 1 complete — outputs in: %s", OUT_CALIBRATION))
