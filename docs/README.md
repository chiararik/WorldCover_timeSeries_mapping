# Documentation

## ATBD — Algorithm Theoretical Basis Document

The full ATBD (v3.0) is available as a Word document in the project Google Drive.
A summary of the key algorithmic components is provided below.

### Algorithm overview

```
INPUT: Landsat 4/5/7/8/9 + Sentinel-1 (post-2014) + WorldCover 2020 labels
                    │
                    ▼
        ┌───────────────────────────┐
        │  GEE PIPELINE (per year) │
        │  1. Cloud masking         │
        │  2. Per-scene indices     │
        │  3. Dry-season median     │
        │  4. JM class schema       │
        │  5. Z-score filtering     │
        │  6. Balanced RF training  │
        │  7. MULTIPROBABILITY      │
        └───────────┬───────────────┘
                    │ Google Drive exports
                    ▼
        ┌───────────────────────────┐
        │  R PIPELINE (all years)  │
        │  M1: Isotonic calibration │
        │  M2: Quality Index        │
        │  M3: Change analysis      │
        │      Stream A (smoothed)  │
        │      Stream B (pixel LCC) │
        │  M4: VIP analysis         │
        └───────────────────────────┘
```

### Key formula reference

All 16 key formulas are available as publication-ready PNG images (300 DPI)
in the `formulas/` subdirectory:

| File | Formula | Section |
|---|---|---|
| F01_composite.png | Dry-season median composite | 4.3 |
| F02_zscorefilter.png | Z-score spectral filter | 7.2 |
| F03_bhattacharyya.png | Bhattacharyya divergence | 6.3 |
| F04_jm_distance.png | Jeffries-Matusita distance | 6.3 |
| F05_raw_confidence.png | Raw RF confidence | 9.2 |
| F06_isotonic.png | Isotonic (PAV) calibration | 13.1 |
| F07_renorm.png | Probability renormalisation | 13.1 |
| F08_cal_confidence.png | Calibrated confidence | 9.3.1 |
| F09_shannon_entropy.png | Shannon entropy (uncertainty) | 9.3.2 |
| F10_qi.png | Quality Index | 14.1 |
| F11_smoothing.png | Confidence-weighted smoothing | 15.2 |
| F12_pix_km2.png | Pixel area (geographic CRS fix) | 15.1 |
| F13_stable.png | Stable pixel definition | 15.3 |
| F14_certain_change.png | Certain change definition | 15.3 |
| F15_fraction_certain.png | Fraction certain | 15.3 |
| F16_ols_trend.png | QI-weighted OLS trend | Results |

### Version history

| Version | Date | Changes |
|---|---|---|
| 1.0 | 2024 | Initial release. BAP compositing. Multi-year loop. |
| 2.0 | 2025 | Median-of-indices compositing. Single-year mode. JM_YEAR fixed. |
| 3.0 | 2025 | Full R pipeline docs. Calibrated uncertainty (Shannon entropy). Stream B NoData fix (v2.2). pix_km2 CRS fix (v2.3). AI assistance documented. |

### Scientific references

- Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32.
- ESA (2021). WorldCover 10 m 2020 v100. https://doi.org/10.5281/zenodo.5571936
- Jeffries, H. & Matusita, K. (1951). Annals of the Institute of Statistical Mathematics, 10, 211–221.
- Niculescu-Mizil, A. & Caruana, R. (2005). Predicting good probabilities with supervised learning. *ICML*.
- Shannon, C.E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27(3), 379–423.
- Zadrozny, B. & Elkan, C. (2002). Transforming classifier scores into accurate multiclass probability estimates. *KDD*.
