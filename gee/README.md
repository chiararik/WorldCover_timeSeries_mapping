# GEE Script — wetland_classification.js

## Prerequisites

- A Google Earth Engine account ([earthengine.google.com](https://earthengine.google.com/))
- A GEE asset containing the AOI polygon (see below)
- A Google Drive folder for outputs

## Setup

1. Open the [GEE Code Editor](https://code.earthengine.google.com/)
2. Create a new script and paste the contents of `wetland_classification.js`
3. Edit the **USER CONFIGURATION** block at the top of the file

## Configuration parameters

| Parameter | Default | Description |
|---|---|---|
| `YEAR` | `2020` | Year to process (1985–2024) |
| `JM_YEAR` | `2020` | **Fixed** reference year for class schema — do NOT change between runs |
| `AOI_ASSET` | `(set me)` | GEE asset path to your study area polygon |
| `OUTPUT_FOLDER` | `(set me)` | Google Drive output folder name |
| `DRY_SEASON` | `{11,1,5,31}` | Nov 1 – May 31 |
| `JM_THRESHOLD` | `1.8` | JM separability merge threshold (0–2 scale) |
| `TARGET_SAMPLES_PER_CLASS` | `150` | Max balanced training samples per class |
| `S1_MIN_YEAR` | `2015` | First year with Sentinel-1 data |

## AOI asset

Load your study area polygon from a GEE FeatureCollection asset:

```javascript
var AOI_ASSET = 'projects/your-project/assets/your_polygon';
```

For Ibaji LGA, use the FAO GAUL Level-2 dataset:
```javascript
var aoi = ee.FeatureCollection("FAO/GAUL/2015/level2")
  .filter(ee.Filter.eq('ADM2_NAME', 'Ibaji')).geometry();
```

## Running for multiple years

Run the script once per year, changing only `YEAR`. JM_YEAR must remain fixed
at 2020 for all years to ensure a consistent class schema across the time series.

A GEE task is submitted automatically for each output file. Monitor progress in
the **Tasks** tab. Outputs land in the specified Google Drive folder.

## Outputs per year

| File | Format | Content |
|---|---|---|
| `land_cover_YEAR.tif` | Int16 GeoTIFF | Classified map; class=0 = NoData |
| `confidence_YEAR.tif` | Byte 0–100 | Raw RF max class probability × 100 |
| `multiprob_YEAR.tif` | Byte 0–100 | N-band per-class probabilities × 100 |
| `accuracy_YEAR.csv` | CSV | PA, UA, F1 per class; OA, Kappa, BA overall |
| `validation_points_YEAR.csv` | CSV | Per-point probabilities (for Module 1) |
| `variable_importance_YEAR.csv` | CSV | Feature importance (Gini, normalised, rank) |
| `class_lookup_table.csv` | CSV | Class code ↔ class name mapping |
| `jm_distances.csv` | CSV | Pairwise JM separability matrix |

## Troubleshooting

**"Cannot export array bands"**: The `arrayFlatten()` call is required for
MULTIPROBABILITY mode. Ensure the GEE API version is ≥ April 2024.

**"User memory limit exceeded"**: Reduce `TARGET_SAMPLES_PER_CLASS` or
increase `Z_THRESHOLD` to 4 (fewer samples retained per class).

**Sentinel-1 errors before 2015**: S1 availability is checked client-side
(`year >= S1_MIN_YEAR`). If you set `YEAR < 2015`, S1 features are
automatically excluded from the feature stack.
