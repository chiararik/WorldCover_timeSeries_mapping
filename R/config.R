# =============================================================================
#  config.R — CENTRALIZED USER CONFIGURATION  v2.2
#  Wetland Land Cover Classification — Ibaji Floodplain, Nigeria
#  ENEA — Division of Models and Technologies for Disaster Risk Reduction
# =============================================================================
#
#  INSTRUCTIONS
#  This is the ONLY file you need to edit before running any module.
#  All paths, thresholds and parameters are defined here.
#  Run any module independently: source("R/0X_module.R")
#  Run the full pipeline:       source("run_all.R")
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. INPUT PATHS
# -----------------------------------------------------------------------------
GEE_DIR    <- "path/to/your/gee_outputs/"
LOOKUP_CSV <- file.path(GEE_DIR, "class_lookup_table.csv")

# -----------------------------------------------------------------------------
# 2. OUTPUT PATHS
# -----------------------------------------------------------------------------
OUTPUT_ROOT     <- "path/to/your/outputs/"
OUT_CALIBRATION <- file.path(OUTPUT_ROOT, "01_calibration")
OUT_QUALITY     <- file.path(OUTPUT_ROOT, "02_quality")
OUT_CHANGE      <- file.path(OUTPUT_ROOT, "03_change")
OUT_VIP         <- file.path(OUTPUT_ROOT, "04_vip")

# -----------------------------------------------------------------------------
# 3. MODULE 1 — RF PROBABILITY CALIBRATION
# -----------------------------------------------------------------------------
# Year to calibrate. Run Module 1 once per year that has multiprob + validation.
# Change CALIB_YEAR and re-run source("R/01_rf_calibration.R") for each year.
CALIB_YEAR           <- 2020
CALIB_MIN_RELIABLE_N <- 10

# -----------------------------------------------------------------------------
# 4. MODULE 2 — ANNUAL QUALITY INDEX
# -----------------------------------------------------------------------------
QI_WEIGHTS <- list(
  coverage   = 0.30,
  confidence = 0.30,
  bal_acc    = 0.30,
  s1_bonus   = 0.10
)
QI_LOW_THRESHOLD <- 0.5
S1_START_YEAR    <- 2015

# -----------------------------------------------------------------------------
# 5. MODULE 3 — LAND COVER CHANGE ANALYSIS
# -----------------------------------------------------------------------------
SMOOTH_HALF_WINDOW <- 1   # temporal smoothing window: ±h years

# Confidence threshold for uncertainty-gated change detection (0-100 scale).
# A changed pixel is "certain" only if BOTH years have confidence >= this value.
#
# IMPORTANT — two very different regimes:
#
# With RAW confidence (confidence_YEAR.tif, Byte 0-100):
#   RF probabilities are under-dispersed. On pixels that actually changed class,
#   only 0.1-8% exceed 60 in both years -> Stream B uninformative at any threshold.
#
# With CALIBRATED confidence (calibrated_confidence_YEAR.tif, from Module 1):
#   Calibrated probabilities are well-dispersed. On changed pixels ~11-25% exceed
#   60% and ~27-56% exceed 50% in both years -> meaningful signal available.
#   Recommended: CONF_THRESHOLD = 50 when USE_CALIBRATED_CONF = TRUE.
#   Stream B only works well when Module 1 has been run for both years of each pair.
CONF_THRESHOLD <- 50

# If TRUE and calibrated_confidence_YEAR.tif exists in OUT_CALIBRATION,
# it is used instead of raw confidence in smoothing and uncertainty gating.
USE_CALIBRATED_CONF <- TRUE

# -----------------------------------------------------------------------------
# 6. MODULE 4 — VARIABLE IMPORTANCE OVER TIME
# -----------------------------------------------------------------------------
TOP_N_FEATURES        <- 10
VIP_MIN_YEARS_PRESENT <- 3

FEATURE_GROUPS <- list(
  Water      = c("MNDWI","AWEI","NDWI","LSWI","NDPI","WRI","S1_VV","S1_VH","S1_VVVH_ratio"),
  Vegetation = c("NDVI","EVI2","SAVI","GCVI","MSAVI","NDMI","NDAVI","NDRE"),
  Soil_Bare  = c("NDBI","BSI","NDTI","SBI","EBBI"),
  Climate    = c("tc_tmmn_mean","tc_tmmx_mean","tc_def_mean","tc_vap_mean","tc_soil_mean"),
  Terrain    = c("srtm_elevation","srtm_slope","srtm_roughness","srtm_twi")
)

# -----------------------------------------------------------------------------
# 7. VISUAL STYLE
# -----------------------------------------------------------------------------
CLASS_COLOURS <- c(
  "Woody vegetation"       = "#1a9641",
  "Tree cover"             = "#1a9641",
  "Shrubland"              = "#a6d96a",
  "Grassland / Cropland"   = "#ffffbf",
  "Grassland"              = "#d9ef8b",
  "Cropland"               = "#fdae61",
  "Non-vegetated / Bare"   = "#d7191c",
  "Permanent water bodies" = "#2c7bb6",
  "Herbaceous wetland"     = "#74add1",
  "Mangroves"              = "#006837"
)
GROUP_COLOURS <- c(Water      = "#2c7bb6",
                   Vegetation = "#1a9641",
                   Soil_Bare  = "#d7191c",
                   Climate    = "#fdae61",
                   Terrain    = "#9970ab",
                   Other      = "#cccccc")
PLOT_DPI <- 180

# =============================================================================
# INTERNAL — do not edit below this line
# =============================================================================
for (d in c(OUTPUT_ROOT, OUT_CALIBRATION, OUT_QUALITY, OUT_CHANGE, OUT_VIP))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

.check_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0)
    stop(sprintf("Install missing packages:\n  install.packages(c(%s))",
                 paste0('"', missing, '"', collapse = ", ")))
  invisible(lapply(pkgs, library, character.only = TRUE))
}
.ts <- function() format(Sys.time(), "[%H:%M:%S]")
message(.ts(), " config.R v2.2 loaded — output root: ", OUTPUT_ROOT)
