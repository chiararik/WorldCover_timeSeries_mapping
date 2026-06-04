# Annual Wetland Land Cover Classification from Landsat Time Series

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![GEE](https://img.shields.io/badge/Google%20Earth%20Engine-JS%20API-4CAF50)](gee/wetland_classification.js)
[![R](https://img.shields.io/badge/R-%E2%89%A54.1-276DC3)](R/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXXX-orange)](https://doi.org/10.5281/zenodo.XXXXXXX)

**Annual Wetland Land Cover Classification from Landsat Time Series  
Ibaji Floodplain, Kogi State, Nigeria — 1985–2024**

Part of the UNOOSA [Space4Water](https://www.space4water.org) programme — Challenge 64:  
*Disappearing Wetlands in parts of North Central Nigeria*

---

## Overview

This repository provides a fully automated, cloud-based pipeline for annual wetland land cover mapping using the Landsat archive (1985–2024) in Google Earth Engine, with a four-module R post-processing pipeline for probability calibration, quality assessment, change analysis, and variable importance analysis.

**Five land cover classes** are mapped annually:
| Class | WorldCover source | Ecological role |
|---|---|---|
| Permanent water bodies | WC80 | Open water; fisheries habitat |
| Herbaceous wetland | WC90 | Seasonally flooded vegetation |
| Woody vegetation | WC10 + WC20 | Tree cover and shrubland |
| Grassland / Cropland | WC30 + WC40 | Smallholder agriculture and pasture |
| Non-vegetated / Bare | WC50 + WC60 | Bare ground and built-up |

### Key results (case study Ibaji LGA, 2,281.8 km², 1985–2024)
- **Permanent water bodies: −20.1%** (R² = 0.85, p < 10⁻¹³) — the strongest and most robust signal
- **Grassland/Cropland: +8.5%** (p = 0.016) — significant agricultural expansion
- **Herbaceous wetland: +60.4%** (p = 0.21, n.s.) — non-linear increase concentrated 1985–1999
- **Water → Bare ground transition: 20.6 km²** — ecologically critical signal of water body retreat

---

## Repository structure

```
wetland-ibaji-nigeria/
├── README.md               ← this file
├── LICENSE                 ← MIT licence
├── CITATION.cff            ← machine-readable citation (APA/BibTeX via GitHub)
├── .gitignore
│
├── gee/
│   ├── README.md           ← GEE-specific instructions
│   └── wetland_classification.js   ← main GEE script (run per year)
│
├── R/
│   ├── WetlandNigeria.Rproj
│   ├── config.R            ← EDIT THIS — all user settings
│   ├── run_all.R           ← master pipeline script
│   ├── 00_utils.R          ← shared utility functions
│   ├── 01_rf_calibration.R ← Module 1: isotonic probability calibration
│   ├── 02_quality_index.R  ← Module 2: annual quality index
│   ├── 03_change_analysis.R← Module 3: land cover change (v2.3)
│   └── 04_vip_analysis.R   ← Module 4: variable importance over time
│
├── docs/
│   ├── ATBD_v3.0.md        ← Algorithm Theoretical Basis Document (summary)
│   └── formulas/           ← publication-quality formula PNG images (300 DPI)
│       ├── F01_composite.png
│       ├── F09_shannon_entropy.png
│       └── ... (16 total)
│
└── data/
    └── README.md           ← data description and download instructions
```

---

## Quick start

### Step 1 — GEE classification (one year at a time)

1. Open the [GEE Code Editor](https://code.earthengine.google.com/)
2. Copy `gee/wetland_classification.js` into the editor
3. Set `YEAR`, `AOI_ASSET`, and `OUTPUT_FOLDER` in the **USER CONFIGURATION** block
4. Run → outputs exported to Google Drive
5. Repeat for each year (1985–2024)

See `gee/README.md` for the full list of outputs and configuration options.

### Step 2 — R post-processing

```r
# 1. Clone the repository and open the R project
# 2. Edit R/config.R — set GEE_DIR and OUTPUT_ROOT at minimum
# 3. Run the full pipeline:
source("run_all.R")

# Or run individual modules:
source("run_all.R")
run(2)        # Module 2 only
run(c(1, 3))  # Modules 1 and 3
```

---

## Module overview

| Module | Input (from GEE) | Key outputs |
|---|---|---|
| **1 — RF Calibration** | `validation_points_YEAR.csv`, `multiprob_YEAR.tif` | `calibrated_confidence_YEAR.tif`, `calibrated_uncertainty_YEAR.tif` |
| **2 — Quality Index** | `land_cover_*.tif`, `confidence_*.tif`, `accuracy_*.csv` | `quality_index.csv`, `quality_index.png` |
| **3 — Change Analysis** | `land_cover_*.tif`, `confidence_*.tif`, `class_lookup_table.csv` | `smoothed_YEAR.tif`, `area_by_class_year.csv`, `trend_by_class.csv`, 7 plots |
| **4 — VIP Analysis** | `variable_importance_YEAR.csv` | `vip_heatmap.png`, `vip_top10_timeseries.png`, `vip_pre_post_s1.png` |

Module dependencies:
```
Module 1 ──► Module 3 (calibrated_confidence → uncertainty gating + smoothing)
Module 2 ──► Module 3 (quality_index.csv → QI-weighted OLS, LQ year flags)
Module 4 is independent
```

Modules 2 and 3 can run without Module 1 (fallback to raw confidence).  
Module 3 can run without Module 2 (fallback to uniform OLS weights).

---

## Required R packages

```r
install.packages(c(
  "terra", "dplyr", "tidyr", "ggplot2", "scales",
  "stringr", "purrr", "jsonlite", "RColorBrewer",
  "viridis", "gridExtra", "forcats"
))
# Optional (for alluvial plot in Module 3):
install.packages("ggalluvial")
```

---

## Known limitations

| Limitation | Effect |
|---|---|
| **Circular validation** | Labels (WorldCover 2020) used for both training and validation. Metrics measure RF reproducibility of WorldCover, not true accuracy. Independent field validation required for publication. |
| **Sparse Landsat archive 1985–2003** | Coverage as low as 18.8% (1985). Temporal smoothing fills gaps but introduces information from adjacent years. |
| **RF probability under-dispersion** | Raw confidence maps are uninformative. Corrected by Module 1 (isotonic calibration). |
| **Sentinel-1 absent before 2014** | SAR discrimination of water/flooded vegetation unavailable for the first 30 years. |
| **WorldCover label temporal distance** | 2020 labels applied to 1985–2024. Areas with genuine change between 1985 and 2020 may be systematically misclassified in historical years. |

---

## AI assistance declaration

The GEE JavaScript pipeline and R post-processing modules in this repository were co-developed with AI assistance using [Claude](https://claude.ai) (Anthropic, models claude-sonnet-4-5 and claude-sonnet-4-6, accessed 2025). AI assistance included code generation, debugging, and identification of methodological limitations. All outputs were reviewed and validated by the authors, who take full responsibility for the content.

See the ATBD (Algorithm Theoretical Basis Document) for full documentation.

---

## Citation

If you use this code or methodology, please cite:

```bibtex
@software{richiardi_wetland_2025,
  author    = {Richiardi, Chiara and Awe-Peter, Helen and Steinbach, Stefanie},
  title     = {Annual Wetland Land Cover Classification from Landsat Time Series:
               Ibaji Floodplain, Nigeria (1985–2024)},
  year      = {2025},
  publisher = {GitHub},
  url       = {https://github.com/chiararichiardi/wetland-ibaji-nigeria},
  note      = {UNOOSA Space4Water Programme — Challenge 64}
}
```

---

## Licence

MIT — see [LICENSE](LICENSE).  
The code is freely reusable. Please cite this repository if you adapt the methodology.

---

## Contacts

**Chiara Richiardi** — ENEA, Laboratory Biodiversity and Ecosystems  
Division of Anthropic and Climate Change Impacts  
[chiara.richiardi@enea.it](mailto:chiara.richiardi@enea.it)

**Helen Awe-Peter** — NASRDA, Strategic Space Applications Department, Abuja, Nigeria
