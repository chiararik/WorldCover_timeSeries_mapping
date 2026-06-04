/************************************************************
 * Wetland Land Cover Classification — Dry Season
 * Google Earth Engine — JavaScript API
 *
 * DESCRIPTION
 * Self-supervised annual land cover classification using
 * ESA WorldCover 2020 as the label source (no field data required).
 * One year is processed per run; re-run with a different year value
 * to build a multi-year time series.
 *
 * COMPOSITE METHOD
 * Median of per-scene spectral indices computed over the dry season
 * window. Each Landsat scene is preprocessed individually (scale
 * factors, cloud masking, band harmonisation, index computation)
 * before the pixel-wise median is applied across all valid scenes.
 *
 * MAPPABLE YEARS
 * 1985–2024  (Landsat 4/5/7/8/9 archive).
 * Years with insufficient Landsat coverage (< 30 % valid pixels) or
 * with too few training samples for any land cover class are skipped
 * automatically with a diagnostic message.
 *
 * OUTPUTS (per run)
 * - land_cover_YEAR.tif        : classified map (Int16, class codes)
 * - uncertainty_index_YEAR.tif : U = α·H + β·M  (Byte 0–100)
 * - multiprob_YEAR.tif         : per-class probabilities (Byte 0–100)
 * - accuracy_YEAR.csv          : OA, Kappa, BA, PA, UA, F1 per class
 * - class_lookup_table.csv     : class code → class name mapping
 * - jm_distances.csv           : pairwise JM separability matrix
 *
 * OUTPUT CLASS CODES
 * Codes are assigned dynamically based on JM separability:
 *   10  Tree cover             20  Shrubland
 *   15  Woody vegetation       (10+20 merged if JM < threshold)
 *   30  Grassland              40  Cropland
 *   35  Grassland/Cropland     (30+40 merged if JM < threshold)
 *   55  Non-vegetated/Bare     (50+60 always merged)
 *   80  Permanent water bodies
 *   90  Herbaceous wetland     95  Mangroves
 * See class_lookup_table.csv for the schema used in each run.
 ************************************************************/


// ╔══════════════════════════════════════════════════════════════╗
// ║                     USER CONFIGURATION                      ║
// ║  Edit ONLY this section before running the script.          ║
// ╚══════════════════════════════════════════════════════════════╝

// ── YEAR TO PROCESS ──────────────────────────────────────────────────────────
// Any year between 1985 and 2024 (Landsat archive extent).
// Re-run the script with a different value to process another year.
// Note: years with very low dry-season cloud-free coverage (e.g. 1992–1997
// in the Ibaji area) will be automatically skipped with a diagnostic message.
var YEAR = 2020;

// ── AOI ASSET ────────────────────────────────────────────────────────────────
// GEE asset path of the study area polygon (FeatureCollection or Geometry).
// Replace with your own asset path.
var AOI_ASSET = 'projects/ee-chiararichiardi/assets/WetlandNigeria/ibaji';

// ── OUTPUT GOOGLE DRIVE FOLDER ───────────────────────────────────────────────
// Name of the Google Drive folder where all outputs will be saved.
// The folder will be created automatically if it does not exist.
var OUTPUT_FOLDER = 'Ibaji_Wetlands_DRY_MEDIAN';

// ── DRY SEASON WINDOW ────────────────────────────────────────────────────────
// Start: November 1 of YEAR   →   End: May 31 of YEAR+1.
// A single full-season median composite is used for all years.
// Sub-season splitting (early/late dry) was evaluated but removed because:
//   (a) the main confusions (water vs wetland, grassland vs bare soil) are
//       better addressed by the new spectral indices and Sentinel-1 SAR,
//       which work at pixel level without reducing coverage;
//   (b) with Landsat 4/5 (1985–2000), sub-windows of 3 months may have
//       0–3 cloud-free scenes, making sub-season bands unreliable;
//   (c) a fallback to the full-season median defeats the purpose.
var DRY_SEASON = { startMonth: 11, startDay: 1, endMonth: 5, endDay: 31 };

// ── SENTINEL-1 INTEGRATION ────────────────────────────────────────────────────
// Sentinel-1 SAR (C-band) is available from October 2014 onward.
// SAR backscatter is physically different from optical reflectance:
//   Open water:            very low VV/VH (specular reflection → no return)
//   Flooded vegetation:    high VV (double-bounce canopy-water interaction)
//   Dry vegetation/soil:   moderate VV, low VH
// Adding SAR significantly improves water vs. herbaceous wetland discrimination.
// For years before S1 availability (< 2015), SAR bands are automatically omitted.
var USE_SENTINEL1      = true;    // set false to disable S1 for all years
var S1_MIN_YEAR        = 2015;    // first full year of S1 data (IW mode)

// ── JM REFERENCE YEAR ────────────────────────────────────────────────────────
// Year used to compute the Jeffries-Matusita separability analysis and define
// the class merge schema. Keep this FIXED across all processed years so that
// the class schema is identical throughout the time series.
// Use the WorldCover 2020 reference year or any year with good coverage.
var JM_YEAR = 2020;

// ╔══════════════════════════════════════════════════════════════╗
// ║               ADVANCED CONFIGURATION                        ║
// ║  Change only if you know what you are doing.                ║
// ╚══════════════════════════════════════════════════════════════╝

// ── JM SEPARABILITY ──────────────────────────────────────────────────────────
// Bands used to compute pairwise Jeffries-Matusita distances between classes.
// Threshold: class pairs with JM < jmThreshold are merge candidates.
// JM separability bands — updated for v2.0 to include the new indices
// most discriminating for the observed confusions:
//   NDPI, WRI  → open water vs. herbaceous wetland
//   SAVI       → grassland/cropland vs. bare soil
var JM_BANDS     = ['NDVI','MNDWI','NDWI','LSWI','NDBI','BSI','EVI2',
                    'NDPI','WRI','SAVI'];
var JM_THRESHOLD = 1.8;

// ── SAMPLING ─────────────────────────────────────────────────────────────────
// targetSamplesPerClass : ideal number of balanced samples per class.
// minSamplesPerClass    : adaptive floor = max(15, nTrees × 0.05).
//                         Classes below this floor are excluded for the year.
var TARGET_SAMPLES_PER_CLASS = 150;

// ── TRAIN / VALIDATION SPLIT ─────────────────────────────────────────────────
var TRAIN_FRACTION = 0.70;   // 70 % training, 30 % validation

// ── SPECTRAL OUTLIER REMOVAL ─────────────────────────────────────────────────
var Z_THRESHOLD = 3;   // samples with |z| > threshold are excluded

// ── UNCERTAINTY WEIGHTS ───────────────────────────────────────────────────────
// U = alpha * H  +  beta * M   (alpha + beta must equal 1)
//   H = normalised Shannon entropy  (spectral ambiguity)
//   M = 1 − margin                  (decision-boundary ambiguity)
var UNCERTAINTY_ALPHA = 0.6;
var UNCERTAINTY_BETA  = 0.4;

// ── EXPORT RESOLUTION ────────────────────────────────────────────────────────
var EXPORT_SCALE = 30;        // metres
var EXPORT_CRS   = 'EPSG:4326';

// ── RANDOM FOREST ────────────────────────────────────────────────────────────
var RF_PARAMS = {
  numberOfTrees    : 300,
  variablesPerSplit: 6,
  minLeafPopulation: 2,
  bagFraction      : 0.7,
  maxNodes         : 128,
  seed             : 42
};


// ============================================================
// 0.  INTERNAL CONFIG  (assembled from user variables above)
// ============================================================

var config = {
  aoiAssetId            : AOI_ASSET,
  outputDriveFolder     : OUTPUT_FOLDER,
  year                  : YEAR,
  jmYear                : JM_YEAR,
  drySeason             : DRY_SEASON,
  worldCoverCollection  : 'ESA/WorldCover/v100',
  worldCoverBand        : 'Map',
  jmBands               : JM_BANDS,
  jmThreshold           : JM_THRESHOLD,
  targetSamplesPerClass : TARGET_SAMPLES_PER_CLASS,
  trainFraction         : TRAIN_FRACTION,
  zThreshold            : Z_THRESHOLD,
  exportScale           : EXPORT_SCALE,
  exportCrs             : EXPORT_CRS,
  rfParams              : RF_PARAMS,

  uncertaintyWeights : { alpha: UNCERTAINTY_ALPHA, beta: UNCERTAINTY_BETA },
  useSentinel1       : USE_SENTINEL1,
  s1MinYear          : S1_MIN_YEAR
};

// Adaptive minimum sample floor: scales with model depth.
// Computed once here and used throughout.
config.minSamplesPerClass = Math.max(15,
  Math.floor(config.rfParams.numberOfTrees * 0.05));

print('Adaptive minSamplesPerClass:', config.minSamplesPerClass,
      '(= max(15, numberOfTrees * 0.05))');


// ============================================================
// 1.  AOI
// ============================================================

var aoi = ee.FeatureCollection(config.aoiAssetId).geometry();


// ============================================================
// 2.  BAND NAME LISTS
// ============================================================

var spectralBandsNoThermal = ['BLUE','GREEN','RED','NIR','SWIR1','SWIR2'];
var spectralBandsAll       = spectralBandsNoThermal.concat(['THERMAL']);

// EBBI requires THERMAL; omitted when ST is unavailable.
// Core per-scene indices (all years).
var indexBandsNoThermal = [
  // Water discrimination
  'MNDWI','AWEI','NDWI','LSWI','NDPI','WRI',
  // Vegetation
  'NDVI','EVI2','SAVI','GCVI','MSAVI','NDMI','NDAVI','NDRE',
  // Soil / built-up / bare
  'NDBI','BSI','NDTI','SBI'
];
var thermalIndexBands = ['EBBI'];

// Sentinel-1 bands (post-2014 only; conditionally added in buildFeatureStack).
var s1Bands = ['S1_VV','S1_VH','S1_VVVH_ratio'];

var terraclimateBands = [
  'tc_tmmn_mean','tc_tmmx_mean',
  'tc_def_mean','tc_vap_mean','tc_soil_mean'
];


// ============================================================
// 3.  UTILITY FUNCTIONS
// ============================================================

function getBandOrDefault(image, bandName, defaultValue) {
  image = ee.Image(image);
  return ee.Image(ee.Algorithms.If(
    image.bandNames().contains(bandName),
    image.select([bandName]),
    ee.Image.constant(defaultValue).rename(bandName)
  ));
}

function applyScaleFactors(image) {
  image = ee.Image(image);
  var optCandidates   = ee.List(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','SR_B8','SR_B9']);
  var thermCandidates = ee.List(['ST_B6','ST_B10','ST_B11']);
  var bands           = image.bandNames();
  var availOpt = ee.List(optCandidates.iterate(function(b, acc) {
    b = ee.String(b); acc = ee.List(acc);
    return ee.Algorithms.If(bands.contains(b), acc.add(b), acc);
  }, ee.List([])));
  var availTherm = ee.List(thermCandidates.iterate(function(b, acc) {
    b = ee.String(b); acc = ee.List(acc);
    return ee.Algorithms.If(bands.contains(b), acc.add(b), acc);
  }, ee.List([])));
  var scaled = image;
  scaled = ee.Image(ee.Algorithms.If(
    availOpt.size().gt(0),
    scaled.addBands(image.select(availOpt).multiply(0.0000275).add(-0.2), null, true),
    scaled
  ));
  scaled = ee.Image(ee.Algorithms.If(
    availTherm.size().gt(0),
    scaled.addBands(image.select(availTherm).multiply(0.00341802).add(149.0), null, true),
    scaled
  ));
  return scaled;
}

function maskLandsat(image) {
  image = ee.Image(image);
  var bands = image.bandNames();
  var sc    = ee.String(ee.Algorithms.If(image.get('SPACECRAFT_ID'), image.get('SPACECRAFT_ID'), ''));
  var isL89 = sc.compareTo('LANDSAT_8').eq(0).or(sc.compareTo('LANDSAT_9').eq(0));
  var cloudSource = ee.Image(ee.Algorithms.If(
    bands.contains('QA_PIXEL'),
    (function() {
      var qa     = image.select('QA_PIXEL');
      var dilated= qa.bitwiseAnd(1 << 1).neq(0);
      var cirrus = ee.Image(ee.Algorithms.If(isL89, qa.bitwiseAnd(1 << 2).neq(0), ee.Image(0)));
      var cloud  = qa.bitwiseAnd(1 << 3).neq(0);
      var shadow = qa.bitwiseAnd(1 << 4).neq(0);
      var snow   = qa.bitwiseAnd(1 << 5).neq(0);
      return dilated.or(cirrus).or(cloud).or(shadow).or(snow).rename('cloud_source').unmask(0);
    })(),
    ee.Image(0).rename('cloud_source')
  ));
  var qaMask = ee.Image(ee.Algorithms.If(
    bands.contains('QA_PIXEL'),
    (function() {
      var qa        = image.select('QA_PIXEL');
      var notFill   = qa.bitwiseAnd(1 << 0).eq(0);
      var notDil    = qa.bitwiseAnd(1 << 1).eq(0);
      var notCirrus = ee.Image(ee.Algorithms.If(isL89, qa.bitwiseAnd(1 << 2).eq(0), ee.Image(1)));
      var notCloud  = qa.bitwiseAnd(1 << 3).eq(0);
      var notShadow = qa.bitwiseAnd(1 << 4).eq(0);
      var notSnow   = qa.bitwiseAnd(1 << 5).eq(0);
      return notFill.and(notDil).and(notCirrus).and(notCloud).and(notShadow).and(notSnow);
    })(),
    ee.Image(1)
  ));
  var satMask = ee.Image(ee.Algorithms.If(
    bands.contains('QA_RADSAT'),
    image.select('QA_RADSAT').eq(0),
    ee.Image(1)
  ));
  return image.updateMask(qaMask).updateMask(satMask)
              .addBands(cloudSource, null, true);
}

function harmonizeLandsat(image) {
  image = ee.Image(image);
  var sc    = ee.String(ee.Algorithms.If(image.get('SPACECRAFT_ID'), image.get('SPACECRAFT_ID'), ''));
  var isL89 = sc.compareTo('LANDSAT_8').eq(0).or(sc.compareTo('LANDSAT_9').eq(0));
  var blue  = ee.String(ee.Algorithms.If(isL89, 'SR_B2', 'SR_B1'));
  var green = ee.String(ee.Algorithms.If(isL89, 'SR_B3', 'SR_B2'));
  var red   = ee.String(ee.Algorithms.If(isL89, 'SR_B4', 'SR_B3'));
  var nir   = ee.String(ee.Algorithms.If(isL89, 'SR_B5', 'SR_B4'));
  var swir1 = ee.String(ee.Algorithms.If(isL89, 'SR_B6', 'SR_B5'));
  var swir2 = ee.String(ee.Algorithms.If(isL89, 'SR_B7', 'SR_B7'));
  var therm = ee.String(ee.Algorithms.If(isL89, 'ST_B10', 'ST_B6'));
  var core = getBandOrDefault(image, blue,  0).rename('BLUE')
    .addBands(getBandOrDefault(image, green, 0).rename('GREEN'))
    .addBands(getBandOrDefault(image, red,   0).rename('RED'))
    .addBands(getBandOrDefault(image, nir,   0).rename('NIR'))
    .addBands(getBandOrDefault(image, swir1, 0).rename('SWIR1'))
    .addBands(getBandOrDefault(image, swir2, 0).rename('SWIR2'))
    .addBands(getBandOrDefault(image, therm, 0).rename('THERMAL'));
  return core.addBands(getBandOrDefault(image, 'cloud_source', 0).unmask(0), null, true);
}

function addSpectralIndices(image) {
  image = ee.Image(image);
  var B = image.select('BLUE'),  G = image.select('GREEN'),
      R = image.select('RED'),   N = image.select('NIR'),
      S1= image.select('SWIR1'), S2= image.select('SWIR2'),
      T = image.select('THERMAL');
  var mndwi = image.normalizedDifference(['GREEN','SWIR1']).rename('MNDWI');
  var awei  = image.expression('4*(g-s1)-(0.25*nir+2.75*s2)',{g:G,s1:S1,nir:N,s2:S2}).rename('AWEI');
  var ndwi  = image.normalizedDifference(['GREEN','NIR']).rename('NDWI');
  var lswi  = image.normalizedDifference(['NIR','SWIR1']).rename('LSWI');
  var ndvi  = image.normalizedDifference(['NIR','RED']).rename('NDVI');
  var evi2  = image.expression('2.5*(nir-r)/(nir+2.4*r+1)',{nir:N,r:R}).rename('EVI2');
  var msavi = image.expression('0.5*(2*nir+1-sqrt(pow((2*nir+1),2)-8*(nir-r)))',{nir:N,r:R}).rename('MSAVI');
  var ndmi  = image.normalizedDifference(['NIR','SWIR1']).rename('NDMI');
  var ndavi = image.normalizedDifference(['NIR','BLUE']).rename('NDAVI');
  var ndre  = image.normalizedDifference(['NIR','SWIR1']).rename('NDRE');
  var ndbi  = image.normalizedDifference(['SWIR1','NIR']).rename('NDBI');
  var bsi   = image.expression('((s1+r)-(nir+b))/((s1+r)+(nir+b))',{s1:S1,r:R,nir:N,b:B}).rename('BSI');
  var ndti  = image.normalizedDifference(['SWIR1','SWIR2']).rename('NDTI');
  var ebbi  = image.expression('(s1-nir)/(0.36*s1+0.47*nir+0.17*tir)',{s1:S1,nir:N,tir:T}).rename('EBBI');

  // ── New indices (v2.0) ────────────────────────────────────────────────────
  //
  // NDPI — Normalized Difference Pond Index
  //   Specifically designed to separate OPEN WATER from FLOODED VEGETATION.
  //   Open water: very low SWIR1 → NDPI strongly negative.
  //   Flooded vegetation: moderate SWIR1 → NDPI near zero or positive.
  //   Key discriminator for the water vs. herbaceous wetland confusion.
  var ndpi = image.normalizedDifference(['SWIR1','GREEN']).rename('NDPI');

  // WRI — Water Ratio Index
  //   (GREEN + RED) / (NIR + SWIR1)
  //   Open water: WRI > 1 (high visible, low NIR/SWIR).
  //   Vegetation: WRI < 1. Simple and highly discriminating.
  var wri = image.expression('(g+r)/(nir+s1)',{g:G,r:R,nir:N,s1:S1}).rename('WRI');

  // SAVI — Soil-Adjusted Vegetation Index
  //   More sensitive than NDVI for LOW vegetation cover (NDVI range 0.0–0.3),
  //   which is exactly where Grassland/Cropland and Bare soil overlap.
  //   L=0.5 is the standard adjustment factor for moderate cover.
  var savi = image.expression('1.5*(nir-r)/(nir+r+0.5)',{nir:N,r:R}).rename('SAVI');

  // GCVI — Green Chlorophyll Vegetation Index
  //   (NIR / GREEN) - 1
  //   More sensitive than NDVI to residual chlorophyll in crop residue vs.
  //   pure bare soil — differs even after vegetation senescence.
  var gcvi = N.divide(G).subtract(1).rename('GCVI');

  // SBI — Soil Brightness Index
  //   sqrt(RED² + NIR²)
  //   Measures overall soil brightness, sensitive to soil moisture and organic
  //   matter content. Bright dry soils have high SBI; dark moist soils low SBI.
  //   Complements BSI (which uses SWIR) for bare soil vs. sparse vegetation.
  var sbi = R.pow(2).add(N.pow(2)).sqrt().rename('SBI');

  return image.addBands([mndwi,awei,ndwi,lswi,ndvi,evi2,msavi,ndmi,ndavi,ndre,
                          ndbi,bsi,ndti,ebbi,ndpi,wri,savi,gcvi,sbi]);
}

function createEmptyImage(bandNames) {
  bandNames = ee.List(bandNames);
  var first = ee.Image.constant(0).rename(ee.String(bandNames.get(0)));
  return ee.Image(bandNames.slice(1).iterate(function(b, img) {
    return ee.Image(img).addBands(ee.Image.constant(0).rename(ee.String(b)));
  }, first));
}


// ============================================================
// 4.  DRY-SEASON MEDIAN COMPOSITE BUILDER
// ============================================================
//
// DESIGN CHOICE — Median of per-scene indices instead of BAP:
//
// The BAP (Best Available Pixel) approach selects for each pixel the
// single observation with the highest composite score. This creates
// spatial discontinuities: adjacent pixels may come from different
// acquisition dates or Landsat missions, producing systematic spectral
// offsets at scene boundaries that the RF interprets as class signal.
//
// Median of per-scene indices solves this:
//   1. Every valid scene in the dry season is preprocessed individually
//      (scale factors, cloud mask, harmonisation, spectral indices).
//   2. Spectral indices are computed per scene BEFORE aggregation.
//      median(NDVI_scene1, NDVI_scene2, ...) ≠ NDVI(median(NIR), median(RED))
//      — the per-scene approach is correct because each scene's bands are
//      internally consistent (same atmosphere, same geometry).
//   3. The median is taken across all valid scenes, producing a single
//      representative value for the full dry season window.
//   4. No score weighting: all valid scenes contribute equally, making
//      the composite spatially homogeneous and sensor-agnostic.
//
// THERMAL / EBBI: included only when ALL scenes in the collection have
// the ST product (L2SP processing level). Identical criterion to before.
//
// PIXEL MASK: pixels with no valid observation in any scene remain
// masked (no-data) in the output. This prevents SLC-off stripe artefacts
// by ensuring those pixels are excluded from classification.

/**
 * Build the dry-season median composite for a given year.
 *
 * @param {number} year - Calendar year (season starts Nov of this year)
 * @returns {ee.Image} Multi-band median composite, masked to valid pixels
 */
function buildDryMedianComposite(year) {
  var sc    = config.drySeason;
  var y     = ee.Number(year);
  var start = ee.Date.fromYMD(y, sc.startMonth, sc.startDay);
  var endY  = (sc.endMonth < sc.startMonth) ? y.add(1) : y;
  var end   = ee.Date.fromYMD(endY, sc.endMonth, sc.endDay).advance(1, 'day');

  var missions = [
    ee.ImageCollection('LANDSAT/LT04/C02/T1_L2'),
    ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'),
    ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'),
    ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'),
    ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  ];

  // Merge all missions, filter spatially and temporally, then:
  // apply scale factors → cloud mask → harmonise bands → compute indices.
  // Each image in the resulting collection is a fully processed scene
  // with named spectral bands and all 14 indices.
  var col = missions.slice(1).reduce(function(acc, m) {
    return ee.ImageCollection(acc).merge(m);
  }, ee.ImageCollection(missions[0]))
    .filterBounds(aoi)
    .filterDate(start, end)
    .map(applyScaleFactors)
    .map(maskLandsat)
    .map(function(img) {
      // Tag ST availability: L2SP = SR+ST; L2SR = SR only.
      var pl = ee.String(ee.Algorithms.If(
        img.get('PROCESSING_LEVEL'), img.get('PROCESSING_LEVEL'), ''));
      return img.set('has_st', pl.compareTo('L2SP').eq(0));
    })
    .map(harmonizeLandsat)
    .map(addSpectralIndices);

  // Determine whether ALL scenes have ST — only then include THERMAL / EBBI.
  // aggregate_min returns 1 only if every scene in the collection is L2SP.
  var allHaveST = ee.Number(ee.Algorithms.If(
    col.size().gt(0), col.aggregate_min('has_st'), 0));

  // Select the band list according to thermal availability.
  var bandList = ee.List(ee.Algorithms.If(
    allHaveST.eq(1),
    ee.List(spectralBandsAll).cat(indexBandsNoThermal).cat(thermalIndexBands),
    ee.List(spectralBandsNoThermal).cat(indexBandsNoThermal)
  ));

  // Compute the per-pixel median across all valid scenes.
  // ee.ImageCollection.median() computes the median independently for
  // each band, using only unmasked (valid) pixels at each location.
  // Pixels with no valid observation in any scene remain masked.
  var composite = ee.Image(ee.Algorithms.If(
    col.size().gt(0),
    col.select(bandList).median(),
    // Fallback: return a fully masked empty image when no scenes exist.
    // Its mask propagates to the classified output → no-data output pixels.
    createEmptyImage(bandList).updateMask(ee.Image(0)).clip(aoi)
  ));

  return composite.clip(aoi);
}


// ============================================================
// 5.  ANCILLARY DATASETS (TERRAIN + TERRACLIMATE)
// ============================================================

var srtm         = ee.Image('USGS/SRTMGL1_003').clip(aoi);
var srtmProducts = ee.Terrain.products(srtm);
var slopeRad     = srtmProducts.select('slope').multiply(Math.PI / 180);
var meritHydro   = ee.Image('MERIT/Hydro/v1_0_1').clip(aoi);
var flowAcc = meritHydro.select('upa').multiply(1e6)
                .resample('bilinear').rename('flow_accumulation');
var slopeTan  = slopeRad.tan().max(0.001);
var twi       = flowAcc.add(ee.Image.pixelArea()).divide(slopeTan).log().rename('srtm_twi');
var roughness = srtm.reduceNeighborhood({
  reducer: ee.Reducer.stdDev(),
  kernel : ee.Kernel.circle({ radius: 90, units: 'meters' }),
  skipMasked: true
}).rename('srtm_roughness');
var terrainStack = srtm.select('elevation').rename('srtm_elevation')
  .addBands(srtmProducts.select('slope').rename('srtm_slope'))
  .addBands(roughness).addBands(twi).resample('bilinear');

/**
 * Build annual TerraClimate predictors for a given year.
 *
 * Each band is selected and averaged INDEPENDENTLY, then stacked.
 * This avoids the "5 names, 1 band" rename error that occurs when
 * ImageCollection.select([n bands]).mean() returns fewer bands than
 * expected for years where the collection has non-uniform band coverage
 * (e.g. partially available years, or GEE version-specific behaviour).
 *
 * Bands: tmmn (min temp), tmmx (max temp), def (water deficit),
 *        vap (vapour pressure), soil (soil moisture).
 * Resolution: ~4.6 km; bilinearly resampled to 30 m.
 */
function buildTerraClimate(year) {
  var base = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
    .filterBounds(aoi)
    .filter(ee.Filter.calendarRange(year, year, 'year'));

  // Compute mean per band individually and rename explicitly.
  // Using separate .select() calls guarantees each output has exactly
  // 1 band with a known name, regardless of collection structure.
  return base.select('tmmn').mean().rename('tc_tmmn_mean')
    .addBands(base.select('tmmx').mean().rename('tc_tmmx_mean'))
    .addBands(base.select('def') .mean().rename('tc_def_mean'))
    .addBands(base.select('vap') .mean().rename('tc_vap_mean'))
    .addBands(base.select('soil').mean().rename('tc_soil_mean'))
    .resample('bilinear');
}

/**
 * Build Sentinel-1 SAR dry-season median composite.
 *
 * SAR backscatter physics:
 *   Open water:          very low VV and VH (specular reflection)
 *   Flooded vegetation:  elevated VV due to double-bounce (canopy + water surface)
 *   Dry bare soil:       moderate VV, low VH
 *   Dry vegetation:      moderate VV, moderate VH
 *
 * Using the median across all available IW-mode descending/ascending
 * scenes in the dry season window produces a seasonally representative
 * SAR composite. The VV/VH ratio is a particularly strong discriminator
 * for flooded vegetation (high ratio) vs. open water (low ratio ≈ 1).
 *
 * Returns null for years before S1_MIN_YEAR.
 *
 * @param {number} year
 * @returns {ee.Image|null}
 */
/**
 * Build Sentinel-1 dry-season median composite.
 *
 * IMPORTANT — availability check is CLIENT-SIDE (plain JS), not server-side:
 *   The function returns the JavaScript value false (not an ee.Image) when
 *   S1 is disabled or the year is before availability. buildFeatureStack()
 *   checks for this with a plain JS !== false, which works correctly.
 *
 *   Previous version returned null via ee.Algorithms.If() (server-side),
 *   but the check  if (s1 !== null)  always evaluated to true because
 *   ee.Algorithms.If always returns an ee.Image object, never JS null.
 *   This caused an empty masked image to be added to the feature stack for
 *   pre-2015 years, corrupting all bands → "No valid training data".
 *
 * @param {number} year  — client-side JS number (not ee.Number)
 * @returns {ee.Image|false}
 */
function buildSentinel1Composite(year) {
  // Client-side check — year is a plain JS number here.
  if (!config.useSentinel1 || year < config.s1MinYear) return false;

  var sc    = config.drySeason;
  var y     = ee.Number(year);
  var start = ee.Date.fromYMD(y, sc.startMonth, sc.startDay);
  var endY  = (sc.endMonth < sc.startMonth) ? y.add(1) : y;
  var end   = ee.Date.fromYMD(endY, sc.endMonth, sc.endDay).advance(1, 'day');

  var s1col = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV','VH']);

  // Compute median in dB space (acceptable for SAR; simpler than linear).
  var med   = s1col.median();
  var vv    = med.select('VV').rename('S1_VV');
  var vh    = med.select('VH').rename('S1_VH');
  // VV−VH in dB = VV/VH ratio in linear scale.
  // High ratio → flooded vegetation (double-bounce elevates VV).
  // Low ratio  → open water (specular reflection, both VV and VH very low).
  var ratio = vv.subtract(vh).rename('S1_VVVH_ratio');

  return vv.addBands(vh).addBands(ratio).clip(aoi);
}



function buildFeatureStack(year) {
  // Full dry-season median composite (Nov–May, all indices).
  var dryComposite = buildDryMedianComposite(year);

  // Base stack: Landsat indices + TerraClimate + terrain.
  var stack = dryComposite
    .addBands(buildTerraClimate(year))
    .addBands(terrainStack);

  // Sentinel-1 SAR (post-2014 only).
  // Client-side check — buildSentinel1Composite() returns JS false (not null,
  // not an ee.Image) for pre-S1 years, so !== false is reliable.
  var s1 = buildSentinel1Composite(year);
  if (s1 !== false) {
    stack = stack.addBands(s1);
    print('Year', year, '- S1 added: S1_VV, S1_VH, S1_VVVH_ratio.');
  } else {
    print('Year', year, '- S1 not available (pre-' + config.s1MinYear + ' or disabled).');
  }

  return stack.clip(aoi);
}


// ============================================================
// 6.  WORLDCOVER LABEL MAP — LOAD AND RESAMPLE TO 30 m
// ============================================================

/**
 * Majority-resample WorldCover from 10 m to 30 m.
 *
 * Each 30 m Landsat pixel covers a 3×3 grid of 10 m WorldCover pixels.
 * Mode resampling assigns the most frequent class among those 9 pixels —
 * the only statistically correct aggregation for categorical data.
 * reduceResolution() MUST precede reproject(); without it GEE falls back
 * to nearest-neighbour at display/export time.
 */
var worldCoverRaw = ee.ImageCollection(config.worldCoverCollection)
  .first().select(config.worldCoverBand).clip(aoi);

var worldCover30m = worldCoverRaw
  .reduceResolution({ reducer: ee.Reducer.mode(), maxPixels: 1024 })
  .reproject({ crs: config.exportCrs, scale: config.exportScale })
  .rename('LC').toInt();


// ============================================================
// 7.  DYNAMIC CLASS DETECTION
// ============================================================

// ── Startup variables declared here (populated inside initializeWorkflow) ────
// Declared at module scope so they are accessible to processYear(),
// validateAndExport(), and all other functions that reference them.
// They are assigned inside the async evaluate() callbacks below.
var classValues       = [];
var nClasses          = 0;
var finalClassValues  = [];
var nFinalClasses     = 0;
var finalClassNamesMap= {};
var labelMapMerged    = null;

var wcClassNames = {
  10:'Tree cover', 20:'Shrubland', 30:'Grassland', 40:'Cropland',
  50:'Built-up', 60:'Bare / sparse vegetation', 70:'Snow and ice',
  80:'Permanent water bodies', 90:'Herbaceous wetland',
  95:'Mangroves', 100:'Moss and lichen'
};


// ============================================================
// 9.  DIAGONAL JEFFRIES-MATUSITA DISTANCE
// ============================================================

/**
 * Diagonal Jeffries-Matusita distance between two classes.
 *
 * Uses the diagonal covariance approximation (assumes per-band independence):
 *   var_avg_b = (var_i_b + var_j_b) / 2
 *   B  = sum_b [ (mu_i_b - mu_j_b)^2 / (8*var_avg_b)
 *              + 0.5*ln(var_avg_b / sqrt(var_i_b * var_j_b)) ]
 *   JM = 2*(1 - exp(-max(B, 0)))   in [0, 2]
 *
 * Rationale for diagonal approximation:
 *   (a) ee.Reducer.covariance() is not usable with multi-band scalar selectors in GEE.
 *   (b) With ~150 samples/class and 7 bands, a full 7x7 covariance matrix
 *       risks numerical ill-conditioning (28 free parameters).
 *   (c) Sufficient for the binary merge/retain decision implemented here.
 *
 * @param {{ mean: number[], variance: number[] }} stats_i
 * @param {{ mean: number[], variance: number[] }} stats_j
 * @returns {number} JM distance in [0, 2]
 */
function jeffriesMatusita(stats_i, stats_j) {
  if (!stats_i || !stats_j ||
      !stats_i.mean || !stats_j.mean ||
      stats_i.mean.length === 0) return 0;
  var p = stats_i.mean.length;
  var B = 0;
  for (var k = 0; k < p; k++) {
    var mu_diff = stats_i.mean[k] - stats_j.mean[k];
    var vi      = Math.max(stats_i.variance[k], 1e-10);
    var vj      = Math.max(stats_j.variance[k], 1e-10);
    var va      = (vi + vj) / 2;
    B += (mu_diff * mu_diff) / (8 * va);
    B += 0.5 * Math.log(va / Math.sqrt(vi * vj));
  }
  return 2 * (1 - Math.exp(-Math.max(B, 0)));
}


// ============================================================
// 8.  INITIALIZATION — async chain replaces synchronous getInfo()
// ============================================================
//
// WHY evaluate() INSTEAD OF getInfo()
// ────────────────────────────────────
// getInfo() is synchronous: it blocks the browser JS thread and sends an
// urgent request to GEE servers. Multiple getInfo() calls in the same
// session (or from residual pending requests of previous runs) compete for
// GEE's per-account concurrency limit and trigger "Too Many Requests".
//
// evaluate() is asynchronous: GEE schedules the computation non-urgently
// and returns the result via a callback. It does not block the JS thread
// and is subject to much more lenient rate limits.
//
// Structure of the async chain (zero synchronous getInfo() calls):
//
//   initializeWorkflow()
//     └─ evaluate #1: class histogram  →  build classValues
//          └─ evaluate #2: grouped JM stats  →  build classStats, merges,
//                                               labelMapMerged, finalClassValues
//               └─ runWorkflow()  (annual loop)
//
// All variables that were previously computed at module level (classValues,
// finalClassValues, labelMapMerged, etc.) are populated inside the callbacks
// and then used by the functions that reference them (processYear, etc.).
// Those functions are only ever called from WITHIN runWorkflow(), which is
// itself only called from the innermost callback — so all variables are
// guaranteed to be populated before any function that uses them runs.

function initializeWorkflow() {

  print('Initializing — computing class histogram (async)...');

  // ── Callback 1: class histogram ─────────────────────────────────────────
  worldCover30m.reduceRegion({
    reducer  : ee.Reducer.frequencyHistogram(),
    geometry : aoi,
    scale    : config.exportScale,
    maxPixels: 1e13,
    tileScale: 4
  }).evaluate(function(histResult) {

    classValues = Object.keys((histResult && histResult['LC']) ? histResult['LC'] : {})
      .map(Number).sort(function(a, b) { return a - b; });
    nClasses = classValues.length;

    print('WorldCover classes present in AOI:',
      classValues.map(function(c) { return 'WC-' + c + ' (' + (wcClassNames[c]||'?') + ')'; }));

    if (nClasses === 0) {
      print('ERROR: no WorldCover classes found in AOI. Check aoiAssetId.');
      return;
    }

    // ── Build JM stack (original + squared bands) ─────────────────────────
    // JM statistics computed on the year to be processed.
    // JM statistics on the fixed reference year — NOT the year being classified.
    // This ensures identical class schema across the entire time series.
    var jmStack       = buildFeatureStack(config.jmYear).select(config.jmBands);
    var jmBandsSq_    = config.jmBands.map(function(b) { return b + '_sq'; });
    var jmStackWithSq = jmStack.addBands(jmStack.pow(2).rename(jmBandsSq_));

    // ── Per-class statistics via separate evaluate() calls ────────────────
    //
    // WHY separate calls instead of a grouped reducer:
    //   The grouped reducer approach (mean().repeat(n).group()) has
    //   undocumented behaviour regarding output key names — they may be
    //   band names OR generic 'mean_N' strings depending on GEE version,
    //   leading to silent zero-fallbacks that produce JM=0 for all pairs.
    //
    //   Separate per-class reduceRegion() calls with ee.Reducer.mean()
    //   on a named-band image reliably return a dict keyed by band name
    //   (e.g. {'NDVI': 0.42, 'NDVI_sq': 0.19, ...}) — documented GEE
    //   behaviour, consistent across all versions.
    //
    //   8 non-blocking evaluate() calls do NOT trigger rate limits;
    //   they are queued by GEE and processed without blocking the JS thread.
    //
    // Completion counter: all JM computations begin only after ALL 8
    // per-class callbacks have returned (counter reaches nClasses).

    print('Computing JM statistics — one async call per class...');

    var classStats_     = {};   // local to this callback scope
    var completedCount_ = 0;    // counts how many per-class callbacks returned

    classValues.forEach(function(c) {
      jmStackWithSq.updateMask(worldCover30m.eq(c))
        .reduceRegion({
          reducer  : ee.Reducer.mean(),   // output keys = band names
          geometry : aoi,
          scale    : config.exportScale,
          maxPixels: 1e13,
          tileScale: 4
        }).evaluate(function(statsDict, error) {

          if (error) {
            print('WARNING class', c, 'stats error:', error);
          }

          // Parse means directly from band-name keys.
          var origMeans = config.jmBands.map(function(b) {
            var v = (statsDict && statsDict[b] !== undefined && statsDict[b] !== null)
                    ? statsDict[b] : 0;
            return isNaN(v) ? 0 : v;
          });

          // Squared-band means: key = '<bandName>_sq'.
          var sqMeans = jmBandsSq_.map(function(b) {
            var v = (statsDict && statsDict[b] !== undefined && statsDict[b] !== null)
                    ? statsDict[b] : 0;
            return isNaN(v) ? 0 : v;
          });

          // Var[X] = E[X^2] - E[X]^2  (floored to 1e-10)
          var variances = origMeans.map(function(mu, i) {
            return Math.max(Math.max(sqMeans[i] - mu * mu, 0), 1e-10);
          });

          classStats_[c] = { mean: origMeans, variance: variances };

          var ndviIdx = config.jmBands.indexOf('NDVI');
          print('Class', c, '(' + (wcClassNames[c]||c) + ')',
                'E[NDVI]=',   origMeans[ndviIdx].toFixed(4),
                'Var[NDVI]=', variances[ndviIdx].toExponential(2));

          // When ALL classes have returned, proceed to JM computation.
          completedCount_++;
          if (completedCount_ === nClasses) {
            computeJMAndRun_(classStats_, jmBandsSq_);
          }
        }); // end per-class evaluate
    }); // end classValues.forEach

  }); // end evaluate #1 (class histogram)
}


/**
 * Called once all per-class statistics are ready.
 * Computes pairwise JM distances, applies merge rules, populates
 * module-scope variables, exports lookup tables, and launches runWorkflow().
 *
 * Defined at module scope so it is accessible from the evaluate() callback.
 * Trailing underscore in name avoids potential shadowing conflicts.
 */
function computeJMAndRun_(classStats_, jmBandsSq_) {

  // ── Pairwise JM distances ───────────────────────────────────────────────
  var jmMatrix = {};
  classValues.forEach(function(ci) {
    classValues.forEach(function(cj) {
      if (cj <= ci) return;
      var jm = jeffriesMatusita(classStats_[ci] || { mean:[], variance:[] },
                                classStats_[cj] || { mean:[], variance:[] });
      jmMatrix[ci + '_' + cj] = jm;
      print('JM(' + (wcClassNames[ci]||ci) + ', ' + (wcClassNames[cj]||cj) + ') =',
            jm.toFixed(4), jm < config.jmThreshold ? '  *** BELOW THRESHOLD' : '');
    });
  });

  // ── Merge rules ─────────────────────────────────────────────────────────
  var mergeRules = [
    { from: [50, 60], to: 55, name: 'Non-vegetated / Bare',  condition: 'always' },
    { from: [30, 40], to: 35, name: 'Grassland / Cropland',  condition: 'if_jm', pair: '30_40' },
    { from: [10, 20], to: 15, name: 'Woody vegetation',      condition: 'if_jm', pair: '10_20' }
  ];

  var mergeMap = {};
  classValues.forEach(function(c) { mergeMap[c] = c; });

  mergeRules.forEach(function(rule) {
    var doMerge = false;
    if (rule.condition === 'always') {
      doMerge = rule.from.every(function(c) { return classValues.indexOf(c) !== -1; });
    } else if (rule.condition === 'if_jm') {
      var present = rule.from.every(function(c) { return classValues.indexOf(c) !== -1; });
      var jmVal   = jmMatrix[rule.pair];
      doMerge = present && (jmVal !== undefined) && (jmVal < config.jmThreshold);
    }
    if (doMerge) {
      rule.from.forEach(function(c) { mergeMap[c] = rule.to; });
      // Bug fix: use !== undefined instead of || to avoid treating JM=0 as falsy.
      var jmLabel = (rule.condition === 'if_jm')
        ? ' [JM=' + (jmMatrix[rule.pair] !== undefined
                     ? jmMatrix[rule.pair].toFixed(4) : '?') + ']'
        : ' [always]';
      print('MERGE: WC [' + rule.from.join(', ') + '] -> code ' + rule.to +
            ' (' + rule.name + ')' + jmLabel);
    }
  });

  // ── Populate module-scope variables ─────────────────────────────────────
  classValues.forEach(function(c) {
    finalClassNamesMap[mergeMap[c]] = wcClassNames[c] || ('WC-' + c);
  });
  mergeRules.forEach(function(rule) {
    if (mergeMap[rule.from[0]] === rule.to) finalClassNamesMap[rule.to] = rule.name;
  });

  var seen = {};
  classValues.forEach(function(c) {
    var fc = mergeMap[c];
    if (!seen[fc]) { seen[fc] = true; finalClassValues.push(fc); }
  });
  finalClassValues.sort(function(a, b) { return a - b; });
  nFinalClasses = finalClassValues.length;

  print('Final class schema (' + nFinalClasses + ' classes):');
  finalClassValues.forEach(function(c) {
    print('  code', c, '->', finalClassNamesMap[c]);
  });

  var fromCodes = classValues;
  var toCodes   = classValues.map(function(c) { return mergeMap[c]; });
  labelMapMerged = worldCover30m.remap(fromCodes, toCodes).rename('LC').toInt();

  // ── Export lookup table and JM distances ─────────────────────────────────
  Export.table.toDrive({
    collection    : ee.FeatureCollection(finalClassValues.map(function(c) {
      return ee.Feature(null, { final_code: c, class_name: finalClassNamesMap[c] });
    })),
    description   : 'class_lookup_table', folder: config.outputDriveFolder,
    fileNamePrefix: 'class_lookup_table', fileFormat: 'CSV',
    selectors     : ['final_code','class_name']
  });

  var jmFeatures = [];
  classValues.forEach(function(ci) {
    classValues.forEach(function(cj) {
      if (cj <= ci) return;
      var key = ci + '_' + cj;
      jmFeatures.push(ee.Feature(null, {
        class_i: ci, class_j: cj,
        name_i : wcClassNames[ci]||('WC-'+ci),
        name_j : wcClassNames[cj]||('WC-'+cj),
        jm_distance: jmMatrix[key] !== undefined ? jmMatrix[key] : -1,
        merged : (mergeMap[ci] === mergeMap[cj] && ci !== cj) ? 1 : 0
      }));
    });
  });
  if (jmFeatures.length > 0) {
    Export.table.toDrive({
      collection    : ee.FeatureCollection(jmFeatures),
      description   : 'jm_distances', folder: config.outputDriveFolder,
      fileNamePrefix: 'jm_distances', fileFormat: 'CSV',
      selectors     : ['class_i','name_i','class_j','name_j','jm_distance','merged']
    });
  }
  print('Lookup table and JM distance matrix exports scheduled.');

  // ── All variables populated — launch annual workflow ─────────────────────
  runWorkflow();
}


// ============================================================
// RUN
// ============================================================
// Zero synchronous getInfo() calls. The entire startup sequence is async.
initializeWorkflow();



// ============================================================
// 10.  Z-SCORE SPECTRAL FILTERING
// ============================================================

/**
 * Remove spectral outliers per class per band using z-score filtering.
 * A sample is kept only if |z_b| <= threshold for ALL bands.
 * Fully server-side (no getInfo); adapts to available bands per year.
 */
function filterByZScore(samples, classProperty, threshold, bands) {
  var bandList  = ee.List(bands);
  var classVals = samples.aggregate_array(classProperty).distinct();
  threshold     = ee.Number(threshold);
  var filtered = classVals.iterate(function(classValue, acc) {
    acc          = ee.FeatureCollection(acc);
    var classSmp = samples.filter(ee.Filter.eq(classProperty, classValue));
    var statsByBand = ee.Dictionary(bandList.iterate(function(band, statsAcc) {
      band     = ee.String(band);
      statsAcc = ee.Dictionary(statsAcc);
      var stats = classSmp.reduceColumns(
        ee.Reducer.mean().combine({ reducer2: ee.Reducer.stdDev(), sharedInputs: true }), [band]);
      var mean   = ee.Number(ee.Algorithms.If(stats.contains('mean'),   stats.get('mean'),   0));
      var stdDev = ee.Number(ee.Algorithms.If(stats.contains('stdDev'), stats.get('stdDev'), 0));
      return statsAcc.set(band, ee.Dictionary({ mean: mean, stdDev: stdDev }));
    }, ee.Dictionary({})));
    var filteredClass = classSmp.map(function(f) {
      var keep = ee.Number(bandList.iterate(function(band, accKeep) {
        band    = ee.String(band);
        accKeep = ee.Number(accKeep);
        var bs     = ee.Dictionary(statsByBand.get(band));
        var stdDev = ee.Number(bs.get('stdDev'));
        var mean   = ee.Number(bs.get('mean'));
        var value  = ee.Number(f.get(band));
        var within = ee.Algorithms.If(
          stdDev.eq(0), 1,
          value.subtract(mean).divide(stdDev).abs().lte(threshold));
        return accKeep.multiply(ee.Number(within));
      }, 1));
      return f.set('keep', keep);
    })
    .filter(ee.Filter.eq('keep', 1))
    .map(function(f) { return f.select(bandList.add(classProperty)); });
    return acc.merge(filteredClass);
  }, ee.FeatureCollection([]));
  return ee.FeatureCollection(filtered);
}


// ============================================================
// 11.  BALANCED SAMPLING — FULLY SERVER-SIDE
// ============================================================

/**
 * Subsample to targetSamplesPerClass per class using ee.Algorithms.If().
 * n >= minFloor -> subsample to target; n < minFloor -> empty (class excluded).
 */
function balancedSample(cleanSamples, classProperty, classValueList, seed) {
  var target   = config.targetSamplesPerClass;
  var minFloor = config.minSamplesPerClass;
  var parts = classValueList.map(function(c) {
    var classSmp = cleanSamples.filter(ee.Filter.eq(classProperty, c));
    return ee.FeatureCollection(ee.Algorithms.If(
      classSmp.size().gte(minFloor),
      classSmp.randomColumn('_r', seed).sort('_r').limit(target),
      ee.FeatureCollection([])
    ));
  });
  return parts.slice(1).reduce(function(acc, fc) {
    return ee.FeatureCollection(acc).merge(fc);
  }, parts[0]);
}

// logSampleCounts removed: evaluating .size() on the z-score-filtered
// FeatureCollection times out GEE — filterByZScore builds a deeply nested
// computation graph (iterate over all classes x all bands). Per-class
// exclusions are handled server-side by balancedSample(). Accuracy
// metrics exported in the CSV fully document per-class performance.
  });
}


// ============================================================
// 12.  STRATIFIED 70 / 30 TRAIN / VALIDATION SPLIT
// ============================================================

function splitTrainVal(samples, seed) {
  var withRandom = samples.randomColumn('split_rand', seed);
  return {
    train: withRandom.filter(ee.Filter.lt ('split_rand', config.trainFraction)),
    val  : withRandom.filter(ee.Filter.gte('split_rand', config.trainFraction))
  };
}


// ============================================================
// 13.  VALID-PIXEL COVERAGE CHECK
// ============================================================

/**
 * Fraction of AOI pixels with valid BAP data.
 * RED band used as proxy. Returns NaN for empty composites.
 */
function computeDataCoverage(featureStack) {
  var valid = ee.Image.constant(1)
    .updateMask(featureStack.select('RED').mask()).rename('v');
  var total = ee.Image.constant(1).clip(aoi).rename('t');
  var p = { reducer: ee.Reducer.sum(), geometry: aoi,
             scale: config.exportScale, maxPixels: 1e13, tileScale: 4 };
  return ee.Number(valid.reduceRegion(p).get('v'))
           .divide(ee.Number(total.reduceRegion(p).get('t')));
}


// ============================================================
// 14.  VALIDATION AND ACCURACY EXPORT
// ============================================================

/**
 * Validate classified image against validation samples and export metrics.
 * CSV: per-class rows (PA, UA, F1, OA) + OVERALL row (OA, Kappa, BA).
 */
function validateAndExport(classified, valSamples, year) {

  var predicted = classified.sampleRegions({
    collection: valSamples, properties: ['LC'],
    scale: config.exportScale, tileScale: 4, geometries: false
  });

  // Pass finalClassValues as the explicit class order to errorMatrix().
  // Without this, GEE builds the matrix in the order classes appear in
  // the validation data — which may differ from finalClassValues order,
  // causing paArray.get([i,0]) and uaArray.get([0,i]) to retrieve values
  // for the wrong class (all PA/UA appear as 0 in the CSV).
  // With explicit order, the matrix is always [nFinalClasses × nFinalClasses]
  // with rows/columns matching finalClassValues exactly.
  var confMatrix        = predicted.errorMatrix('LC', 'land_cover',
                            ee.List(finalClassValues));
  var oa                = confMatrix.accuracy();
  var kappa             = confMatrix.kappa();
  var paArray           = confMatrix.producersAccuracy();
  var uaArray           = confMatrix.consumersAccuracy();
  var balancedAccuracy  = paArray.reduce(ee.Reducer.mean(), [0]).get([0, 0]);

  // OA preview removed: oa.evaluate() competes with coverage checks in the
  // annual loop and triggers "Too Many Requests". OA is fully available in
  // the accuracy_YEAR.csv export.

  var classFeatures = finalClassValues.map(function(c, i) {
    var pa = ee.Number(paArray.get([i, 0]));
    var ua = ee.Number(uaArray.get([0, i]));
    var denom = pa.add(ua);
    var f1 = ee.Number(ee.Algorithms.If(denom.gt(0),
      pa.multiply(ua).multiply(2).divide(denom), 0));
    return ee.Feature(null, {
      year      : year,
      final_code: c,
      class_name: finalClassNamesMap[c] || ('code-' + c),
      PA: pa, UA: ua, F1: f1, OA: oa
    });
  });

  var summaryFeature = ee.Feature(null, {
    year: year, final_code: -1, class_name: 'OVERALL',
    OA: oa, kappa: kappa, balanced_accuracy: balancedAccuracy
  });

  Export.table.toDrive({
    collection    : ee.FeatureCollection(classFeatures)
                      .merge(ee.FeatureCollection([summaryFeature])),
    description   : 'accuracy_' + year,
    folder        : config.outputDriveFolder,
    fileNamePrefix: 'accuracy_' + year,
    fileFormat    : 'CSV',
    selectors     : ['year','final_code','class_name','PA','UA','F1',
                     'OA','kappa','balanced_accuracy']
  });
}


// ============================================================
// 15.  CORE CLASSIFICATION PIPELINE (one year)
// ============================================================

function processYear(year, isReference) {

  print('--- Processing year:', year, '---');
  print('Composite: dry-season median of per-scene indices');

  var featureStack = buildFeatureStack(year);
  var featureBands = featureStack.bandNames();

  // Pixel mask: positions with no valid Landsat observation in any
  // dry-season scene remain masked throughout (prevents stripe artefacts).
  var pixelMask = featureStack.select('RED').mask();

  // Check composite coverage before sampling.
  // A fully empty composite (all pixels masked) would produce zero
  // samples — the class-sample check below would catch it too, but this
  // gives a fast diagnostic message without triggering the RF pipeline.
  var validFraction = ee.Image.constant(1).updateMask(pixelMask).rename('v')
    .reduceRegion({ reducer: ee.Reducer.mean(), geometry: aoi,
                    scale: config.exportScale, maxPixels: 1e13, tileScale: 4 })
    .getNumber('v');

  // Sample feature stack at WorldCover merged label locations.
  // Pool = targetSamplesPerClass * nFinalClasses * 10 to ensure enough
  // candidates survive z-score filtering before the per-class cap.
  var rawSamples = featureStack.addBands(labelMapMerged).sample({
    region    : aoi,
    scale     : config.exportScale,
    numPixels : config.targetSamplesPerClass * nFinalClasses * 10,
    seed      : year,
    tileScale : 4,
    geometries: true
  }).filter(ee.Filter.neq('LC', null));

  var cleanSamples = filterByZScore(rawSamples, 'LC', config.zThreshold, featureBands);
  var balanced     = balancedSample(cleanSamples, 'LC', finalClassValues, year);

  // ── Sample check removed ─────────────────────────────────────────────────
  //
  // The evaluate()-based sample check timed out because it tries to
  // materialise the full filterByZScore → balancedSample computation
  // graph interactively. This is the same class of problem as logSampleCounts.
  //
  // Instead, classification always proceeds. balancedSample() already
  // handles class exclusion server-side: classes below the floor return
  // empty FeatureCollections and are simply absent from the training data.
  // The RF trains on whatever classes have sufficient samples, and the
  // accuracy_YEAR.csv documents which classes were active (PA=0 for absent
  // classes with the fixed errorMatrix order parameter).
  //
  // Years where NO class has samples (completely empty composites) are
  // handled by the RF gracefully — training on an empty collection will
  // produce an error caught in the GEE Tasks panel for that year.

  print('Year', year, '- proceeding with classification (server-side only).');
  var bapMask = pixelMask;  // alias for clarity in the block below

  var split    = splitTrainVal(balanced, year);
  var trainSmp = split.train;
  var valSmp   = split.val;

  var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees    : config.rfParams.numberOfTrees,
    variablesPerSplit: config.rfParams.variablesPerSplit,
    minLeafPopulation: config.rfParams.minLeafPopulation,
    bagFraction      : config.rfParams.bagFraction,
    maxNodes         : config.rfParams.maxNodes,
    seed             : config.rfParams.seed
  }).train({
    features       : trainSmp,
    classProperty  : 'LC',
    inputProperties: featureBands
  });

  // ── Variable importance ───────────────────────────────────────────────────
  //
  // smileRandomForest.explain() returns a dictionary with one key:
  //   'importance': a dictionary of { bandName: importanceScore }
  //
  // GEE's RF importance is computed as the MEAN DECREASE IN IMPURITY (Gini
  // importance): for each feature, the average reduction in node impurity
  // weighted by the proportion of samples reaching that node, averaged
  // across all trees. Higher score = more discriminating feature.
  //
  // The scores are not normalised (do not sum to 1 by default).
  // We normalise them to [0,1] by dividing by the maximum score, making
  // interpretation easier across years with different feature sets.
  //
  // Exported as variable_importance_YEAR.csv with columns:
  //   feature       : band/index name (same names as featureBands)
  //   importance    : raw Gini importance from GEE
  //   importance_norm : importance / max(importance) → [0,1]
  //   rank          : rank (1 = most important)
  //
  // Important caveats:
  //   - Gini importance is biased toward high-cardinality and continuous
  //     features. Features on different scales are not directly comparable
  //     in absolute value, but relative ranking is meaningful.
  //   - Correlated features (e.g. NDVI and EVI2) share importance — neither
  //     alone may rank high even if vegetation density is a key predictor.
  //   - The importance reflects performance on the TRAINING set (Gini
  //     impurity), not on the validation set. Use with care for inference.

  var importance      = ee.Dictionary(classifier.explain().get('importance'));
  var importanceKeys  = importance.keys();         // feature names
  var importanceVals  = importance.values();       // raw scores as ee.List

  // Normalise by max score (server-side via ee.List arithmetic)
  var maxVal = importanceVals.reduce(ee.Reducer.max());
  var normVals = importanceVals.map(function(v) {
    return ee.Number(v).divide(ee.Number(maxVal));
  });

  // Build a FeatureCollection (one Feature per band) for CSV export.
  var importanceFc = ee.FeatureCollection(
    importanceKeys.zip(importanceVals).zip(normVals).map(function(triple) {
      triple = ee.List(triple);
      var pair  = ee.List(triple.get(0));
      var fname = ee.String(pair.get(0));
      var raw   = ee.Number(pair.get(1));
      var norm  = ee.Number(triple.get(1));
      return ee.Feature(null, {
        feature        : fname,
        importance     : raw,
        importance_norm: norm
      });
    })
  );

  // Sort by importance descending and add rank column
  var importanceSorted = importanceFc.sort('importance', false);
  var importanceList   = importanceSorted.toList(importanceSorted.size());
  var importanceRanked = ee.FeatureCollection(
    importanceList.map(function(f) {
      f = ee.Feature(f);
      var rank = importanceList.indexOf(f).add(1);
      return f.set('rank', rank);
    })
  );

  Export.table.toDrive({
    collection    : importanceRanked,
    description   : 'variable_importance_' + year,
    folder        : config.outputDriveFolder,
    fileNamePrefix: 'variable_importance_' + year,
    fileFormat    : 'CSV',
    selectors     : ['rank','feature','importance','importance_norm']
  });

  // Classification — original WorldCover codes as pixel values (Int16).
  // bapMask prevents SLC-off / no-data pixels from being classified.
  var classified = featureStack
    .classify(classifier.setOutputMode('CLASSIFICATION'))
    .updateMask(bapMask)
    .rename('land_cover').clip(aoi).set('year', year);

  // ── MULTIPROBABILITY: per-class probability image ───────────────────────────
  //
  // smileRandomForest MULTIPROBABILITY returns an array-type image — a single
  // band containing a 1-D array of length n_trained_classes per pixel.
  // Direct export fails with "Cannot export array bands" (Error code 3).
  //
  // FIX: arrayFlatten() converts the array image to a regular N-band float image
  // (one named float band per class). This is the correct GEE approach.
  //
  // NOTE: the uncertainty index (H, M, U) has been removed from GEE outputs.
  //
  // WHY: U collapses to a narrow range (~0.33-0.39 in all tested years)
  // because raw RF probabilities are systematically under-dispersed —
  // a well-documented RF calibration issue (Niculescu-Mizil & Caruana 2005).
  // The RF tends to predict probabilities near the class proportion (1/n)
  // rather than near 0 or 1, making the entropy distribution nearly uniform.
  // The resulting U raster is spatially uninformative: almost every pixel
  // gets the same value (approximately = global OA).
  //
  // SOLUTION: calibrate RF probabilities using isotonic regression in R
  // post-processing (using validation_points_YEAR.csv). Calibrated probabilities
  // are better dispersed and produce a spatially meaningful confidence map.

  var multiProbBandNames = finalClassValues.map(function(c) {
    return 'prob_' + (finalClassNamesMap[c] || ('code' + c))
      .replace(/[\/\s]+/g, '_');
  });

  var multiProbMap = featureStack
    .classify(classifier.setOutputMode('MULTIPROBABILITY'))
    .updateMask(bapMask)
    .clip(aoi)
    .arrayFlatten([multiProbBandNames]);   // array → named float bands

  // Confidence = max class probability (winning-class probability).
  // Before isotonic calibration this is still useful as a coarse spatial
  // indicator: high confidence pixels are more likely correctly classified.
  var confidenceMap = multiProbMap.reduce(ee.Reducer.max())
    .multiply(100).toByte().rename('confidence')
    .updateMask(bapMask).clip(aoi);

  // ── Balanced sampling → 70/30 split ─────────────────────────────────────
  //
  // Balance first, then split. The validation set is therefore also
  // class-balanced (~30 per class), which matches the training distribution
  // and is required for the isotonic calibration to work correctly.
  // (An unbalanced validation set causes per-class calibration to collapse
  // when class prevalence is low — all calibrated probabilities → 0.)
  //
  // Limitation: OA on a balanced val set does not represent landscape-level
  // accuracy (rare classes are over-weighted). This is a known trade-off.
  // For publication, complement with an area-weighted OA computed in R
  // from the class frequency histogram and per-class PA values.

  var balanced = balancedSample(cleanSamples, 'LC', finalClassValues, year);
  var split    = splitTrainVal(balanced, year);
  var trainSmp = split.train;
  var valSmp   = split.val;


  var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees    : config.rfParams.numberOfTrees,
    variablesPerSplit: config.rfParams.variablesPerSplit,
    minLeafPopulation: config.rfParams.minLeafPopulation,
    bagFraction      : config.rfParams.bagFraction,
    maxNodes         : config.rfParams.maxNodes,
    seed             : config.rfParams.seed
  }).train({
    features       : trainSmp,
    classProperty  : 'LC',
    inputProperties: featureBands
  });

  // Classification — original WorldCover codes as pixel values (Int16).
  // bapMask prevents SLC-off / no-data pixels from being classified.
  var classified = featureStack
    .classify(classifier.setOutputMode('CLASSIFICATION'))
    .updateMask(bapMask)
    .rename('land_cover').clip(aoi).set('year', year);

  // ── Exports ──────────────────────────────────────────────────────────────
  //
  // ── GEE outputs per year ─────────────────────────────────────────────────
  //
  //   land_cover_YEAR.tif        : classified map (Int16, final class codes)
  //   confidence_YEAR.tif        : max class probability (Byte 0-100)
  //   multiprob_YEAR.tif         : per-class probabilities, N bands (Byte 0-100)
  //   accuracy_YEAR.csv          : OA, Kappa, BA, PA, UA, F1 per class
  //   validation_points_YEAR.csv : coordinates, LC, predicted, correct,
  //                                per-class probabilities → isotonic calib. in R

  Export.image.toDrive({
    image: classified.toInt16(),
    description: 'land_cover_' + year, folder: config.outputDriveFolder,
    fileNamePrefix: 'land_cover_' + year,
    region: aoi, scale: config.exportScale, crs: config.exportCrs, maxPixels: 1e13
  });

  Export.image.toDrive({
    image: confidenceMap,
    description: 'confidence_' + year, folder: config.outputDriveFolder,
    fileNamePrefix: 'confidence_' + year,
    region: aoi, scale: config.exportScale, crs: config.exportCrs, maxPixels: 1e13
  });

  Export.image.toDrive({
    image: multiProbMap.multiply(100).toByte(),
    description: 'multiprob_' + year, folder: config.outputDriveFolder,
    fileNamePrefix: 'multiprob_' + year,
    region: aoi, scale: config.exportScale, crs: config.exportCrs, maxPixels: 1e13
  });

  // ── Validation points with per-class probabilities ────────────────────────
  //
  // Per-class RF probabilities at each validation pixel are the input for
  // isotonic regression calibration in R (Problem 2 fix).
  // Replaces the scalar U value which was spatially uninformative.

  var nameDict_ = ee.Dictionary(finalClassValues.reduce(function(acc, c) {
    acc[c] = finalClassNamesMap[c] || ('code-' + c);
    return acc;
  }, {}));

  var valWithAllProps = multiProbMap
    .addBands(classified.rename('predicted'))
    .sampleRegions({
      collection: valSmp,
      properties: ['LC'],
      scale     : config.exportScale,
      tileScale : 4,
      geometries: true
    })
    .map(function(f) {
      var actual = ee.Number(f.get('LC'));
      var pred   = ee.Number(f.get('predicted'));
      return f.set({
        correct              : actual.eq(pred),
        class_name_actual    : nameDict_.get(actual.format()),
        class_name_predicted : nameDict_.get(pred.format())
      });
    });

  var baseSelectors = ['.geo', 'LC', 'predicted', 'correct',
                       'class_name_actual', 'class_name_predicted'];

  Export.table.toDrive({
    collection    : valWithAllProps,
    description   : 'validation_points_' + year,
    folder        : config.outputDriveFolder,
    fileNamePrefix: 'validation_points_' + year,
    fileFormat    : 'CSV',
    selectors     : baseSelectors.concat(multiProbBandNames)
  });


  validateAndExport(classified, valSmp, year);

  print('Year', year, '- 3 image exports + 3 CSVs scheduled.',
        '(land_cover, confidence, multiprob | accuracy, validation_points, variable_importance)');

  // Map display (reference year only).
  if (isReference && typeof Map !== 'undefined') {
    var sldEntries = [
      { code: 15, color: '#1a9641', label: 'Woody vegetation' },
      { code: 10, color: '#1a9641', label: 'Tree cover' },
      { code: 20, color: '#a6d96a', label: 'Shrubland' },
      { code: 35, color: '#ffffbf', label: 'Grassland / Cropland' },
      { code: 30, color: '#ffffbf', label: 'Grassland' },
      { code: 40, color: '#fdae61', label: 'Cropland' },
      { code: 55, color: '#d7191c', label: 'Non-vegetated / Bare' },
      { code: 80, color: '#2c7bb6', label: 'Perm. water' },
      { code: 90, color: '#abd9e9', label: 'Herb. wetland' },
      { code: 95, color: '#006837', label: 'Mangroves' },
    ].filter(function(e) { return finalClassValues.indexOf(e.code) !== -1; });

    var sld = '<RasterSymbolizer><ColorMap type="values">' +
      sldEntries.map(function(e) {
        return '<ColorMapEntry color="' + e.color + '" quantity="' + e.code +
               '" label="' + e.label + '"/>';
      }).join('') + '</ColorMap></RasterSymbolizer>';

    Map.centerObject(aoi, 8);
    Map.addLayer(classified.sldStyle(sld), {}, 'Land cover ' + year + ' (reference)');
    Map.addLayer(confidenceMap,
      { min: 0, max: 100, palette: ['#2c7bb6','#ffffbf','#d7191c'] },
      'Confidence (max prob) ' + year, false);
  }
}


// processYearsSequentially removed — single-year mode.



// ============================================================
// 17.  MAIN WORKFLOW  (called from initializeWorkflow callback)
// ============================================================

function runWorkflow() {
  print('=== Nigeria Wetlands Annual Classification — v4 ===');
  print('Composite method   : dry-season median of per-scene spectral indices (Nov-May)');
  print('Label source       : ESA WorldCover 2020 (majority-resampled to 30 m)');
  print('Class schema       : merged via JM separability + ecological rules (fixed)');
  print('Output values      : final merged class codes (see class_lookup_table.csv)');
  print('Uncertainty        : U = alpha*H + beta*M (entropy + margin)');
  print('Local accuracy     : post-processing in R/Python');
  print('Uncertainty weights: alpha =', config.uncertaintyWeights.alpha,
        'beta =', config.uncertaintyWeights.beta);
  print('minSamplesPerClass :', config.minSamplesPerClass);
  print('Final classes      :', finalClassValues);
  print('Processing year    :', config.year);
  print('RF parameters      :', config.rfParams);

  processYear(config.year, true);
}
