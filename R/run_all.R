# =============================================================================
#  run_all.R — MASTER PIPELINE SCRIPT
#  Wetland Land Cover Classification — Ibaji Floodplain, Nigeria
# =============================================================================
#
#  USAGE
#  -----
#  Open WetlandNigeria.Rproj in RStudio, edit config.R, then:
#    source("run_all.R")            # run all modules
#    source("run_all.R"); run(1)    # run only Module 1
#    source("run_all.R"); run(c(2,3))  # run Modules 2 and 3
#
#  PIPELINE ORDER
#  --------------
#  Module 1 — RF probability calibration (per year, run once per year)
#  Module 2 — Annual quality index (all years)
#  Module 3 — Land cover change analysis (all years)
#  Module 4 — Variable importance over time (all years)
#
#  MODULE DEPENDENCIES
#  -------------------
#  Module 1 → Module 3 (calibrated confidence used if available)
#  Module 2 → Module 3 (QI weights OLS trend; flags low-quality years)
#  Modules 1, 2, 4 are independent of each other.
#
# =============================================================================

source("config.R")

run <- function(modules = 1:4) {
  t0_all <- proc.time()
  results <- list()

  for (m in modules) {
    script <- switch(as.character(m),
      "1" = "R/01_rf_calibration.R",
      "2" = "R/02_quality_index.R",
      "3" = "R/03_change_analysis.R",
      "4" = "R/04_vip_analysis.R",
      stop(sprintf("Unknown module: %d (valid: 1-4)", m))
    )

    cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
    cat(.ts(), sprintf("  MODULE %d — %s\n", m, basename(script)))
    cat(paste(rep("=", 70), collapse = ""), "\n\n")

    t0 <- proc.time()
    result <- tryCatch({
      source(script, local = FALSE)
      list(status = "OK", elapsed = round((proc.time() - t0)["elapsed"]))
    }, error = function(e) {
      cat("\n  ERROR:", conditionMessage(e), "\n")
      list(status = "ERROR",
           message = conditionMessage(e),
           elapsed = round((proc.time() - t0)["elapsed"]))
    })
    results[[as.character(m)]] <- result
    cat(sprintf("\n  Module %d: %s (%.0fs)\n", m, result$status, result$elapsed))
  }

  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat(.ts(), " PIPELINE COMPLETE\n")
  cat(sprintf("  Total elapsed: %.0fs\n", (proc.time() - t0_all)["elapsed"]))
  for (m in names(results)) {
    r   <- results[[m]]
    sym <- if (r$status == "OK") "\u2713" else "\u2717"
    cat(sprintf("  %s Module %s: %s (%.0fs)\n", sym, m, r$status, r$elapsed))
  }
  cat(paste(rep("=", 70), collapse = ""), "\n")
  invisible(results)
}

if (!interactive()) run(1:4)
