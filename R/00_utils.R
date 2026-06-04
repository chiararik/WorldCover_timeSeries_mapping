# =============================================================================
#  R/00_utils.R — SHARED UTILITY FUNCTIONS
#  Sourced automatically by each module — do not run directly.
# =============================================================================

# ── File discovery ─────────────────────────────────────────────────────────────

#' Find GEE output files matching a pattern and extract years
#' @param dir     Input directory
#' @param pattern Regex pattern with \d{4} for year
#' @return Named list: files = sorted paths, years = integer years
find_files <- function(dir, pattern) {
  files <- sort(list.files(dir, pattern = pattern, full.names = TRUE))
  years <- as.integer(stringr::str_extract(basename(files), "\\d{4}"))
  list(files = files, years = years)
}

#' Select years present in all provided year vectors (intersection)
intersect_years <- function(...) Reduce(intersect, list(...))

# ── Raster helpers ─────────────────────────────────────────────────────────────

#' Load a raster and reproject to a reference grid if needed
load_aligned <- function(path, ref, method = "near") {
  r <- terra::rast(path)
  if (!terra::compareGeom(r, ref, stopOnError = FALSE))
    r <- terra::project(r, ref, method = method)
  r
}

#' Build a multi-layer SpatRaster aligned to a reference grid
load_stack <- function(files, ref, method = "near",
                        names_prefix = "y", years = NULL) {
  stk <- terra::rast(lapply(files, load_aligned, ref = ref, method = method))
  if (!is.null(years)) names(stk) <- paste0(names_prefix, "_", years)
  stk
}

# ── Colour helpers ─────────────────────────────────────────────────────────────

#' Map class names to colours; fall back to RColorBrewer for unknowns
get_class_colours <- function(class_names, colour_map) {
  cols    <- colour_map[class_names]
  missing <- is.na(cols)
  if (any(missing)) {
    fb <- RColorBrewer::brewer.pal(max(3, sum(missing)), "Set2")[seq_len(sum(missing))]
    cols[missing] <- fb
  }
  setNames(cols, class_names)
}

# ── Low-quality year shading ──────────────────────────────────────────────────

#' ggplot layer: grey bands for low-quality years (requires quality_index.csv)
lq_band_layer <- function(qi_df, ymin = -Inf, ymax = Inf) {
  lq <- qi_df[!is.na(qi_df$low_quality) & qi_df$low_quality, "year", drop = TRUE]
  if (length(lq) == 0) return(ggplot2::geom_blank())
  ggplot2::geom_rect(
    data = data.frame(year = lq),
    ggplot2::aes(xmin = year - 0.5, xmax = year + 0.5, ymin = ymin, ymax = ymax),
    fill = "grey80", alpha = 0.4, inherit.aes = FALSE
  )
}

# ── Console helpers ───────────────────────────────────────────────────────────

#' Print a section separator
section <- function(msg, width = 70) {
  cat("\n", paste(rep("-", width), collapse = ""), "\n  ", msg,
      "\n", paste(rep("-", width), collapse = ""), "\n", sep = "")
}

#' Save ggplot with consistent settings and log the path
save_plot <- function(plot, path, width, height, dpi = PLOT_DPI) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  message(.ts(), "  Saved: ", basename(path))
}
