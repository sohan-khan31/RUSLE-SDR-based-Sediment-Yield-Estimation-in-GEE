// ============================================================================
// LOWER TEESTA: SEASONAL SOIL LOSS & SEDIMENT YIELD, 2018-2025
// Pre-monsoon : March 1-April 30
// Post-monsoon: October 1-November 30
// 16 soil-loss maps + 16 sediment-yield maps
// Two separate charts (Pre vs Post for each variable)
// One combined aggregation for mean / min / max statistics
//
// Original K, LS, R, C, P, and logistic SDR equations retained.
// ============================================================================

// 1. SETTINGS
var START_YEAR = 2018;
var END_YEAR = 2025;
var ANALYSIS_SCALE = 30;

// Optional exports
var EXPORT_MAPS = false;
var EXPORT_CSV = false;


// ============================================================================
// 2. DATA SOURCES
// ============================================================================

var soil = ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02');
var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED');
var dwCollection = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1');
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY');
var basins7 = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_7');
var divisions = ee.FeatureCollection('projects/ee-sohan31du/assets/SHP/BD_Division');


// ============================================================================
// 3. AREA OF INTEREST
// ============================================================================

// IMPORTANT INSRUCTIONS 
// If you want to define AOI by intersecting two geographic boundaries (i.e. clipping a large basin into smaller area), modify the following codes accordingly.
// If you want to use your own AOI asset without intersection, add "//" before the following three variable definitions (basin, rangpur, aoi) from Line 42–Line 49.

var basin = basins7.filter(ee.Filter.eq('HYBAS_ID', 4070912240));

var rangpur = divisions.filter(ee.Filter.eq('ADM1_EN', 'Rangpur'));

var aoi = basin.geometry().intersection(
  rangpur.geometry(),
  ee.ErrorMargin(1)
);

Map.centerObject(aoi, 9);

Map.addLayer(
  aoi,
  {color: '333333'},
  'Teesta Basin in BD',
  false
);


// ============================================================================
// 4. SOIL ERODIBILITY FACTOR (K)
// ============================================================================

var texture = soil
  .select('b0')
  .rename('soil')
  .clip(aoi);

var K = texture.expression(

  "(b('soil') > 11) ? 0.0053" +
  ": (b('soil') > 10) ? 0.0170" +
  ": (b('soil') > 9) ? 0.045" +
  ": (b('soil') > 8) ? 0.050" +
  ": (b('soil') > 7) ? 0.0499" +
  ": (b('soil') > 6) ? 0.0394" +
  ": (b('soil') > 5) ? 0.0264" +
  ": (b('soil') > 4) ? 0.0423" +
  ": (b('soil') > 3) ? 0.0394" +
  ": (b('soil') > 2) ? 0.036" +
  ": (b('soil') > 1) ? 0.0341" +
  ": (b('soil') > 0) ? 0.0288" +
  ": 0"

).rename('K').clip(aoi);


// ============================================================================
// 5. TOPOGRAPHIC FACTOR (LS)
// ============================================================================

// NASADEM elevation
var dem = ee.Image('NASA/NASADEM_HGT/001')
  .select('elevation')
  .clip(aoi);

// Slope in degrees
var slopeDeg = ee.Terrain.slope(dem)
  .rename('slope');

// Convert slope to radians
var slopeRad = slopeDeg.multiply(Math.PI / 180);

// HydroSHEDS flow accumulation
var flowAccum = ee.Image(
  'WWF/HydroSHEDS/15ACC'
)
  .select('b1')
  .clip(aoi)
  .rename('flow_acc');

// Resample flow accumulation
var fa30 = flowAccum
  .resample('bilinear')
  .reproject(dem.projection());

// Cell size
var cellSize = ee.Image.pixelArea().sqrt();

// LS exponents
var m = 0.4;
var n = 1.3;

// Slope-length component
var areaTerm = fa30
  .multiply(cellSize)
  .divide(22.13)
  .max(0.001);

// Slope-steepness component
var slopeTerm = slopeRad
  .sin()
  .divide(0.0896)
  .max(0.001);

// Combined LS factor
var LS = areaTerm.pow(m)
  .multiply(slopeTerm.pow(n))
  .rename('LS')
  .clip(aoi);


// ============================================================================
// 6. SEDIMENT DELIVERY RATIO (SDR)
// ============================================================================

// Original spatial logistic SDR formulation

var Dup = flowAccum
  .add(1)
  .log()
  .multiply(K.add(0.001));

var Ddn = slopeRad
  .sin()
  .add(0.001);

var IC = Dup
  .subtract(Ddn)
  .rename('IC');

var SDR = IC.expression(

  '1 / (1 + exp(-a * (IC - b)))',

  {
    IC: IC,
    a: 1.5,
    b: 0.5
  }

).rename('SDR').clip(aoi);


// ============================================================================
// 7. STATIC FACTOR MAPS
// ============================================================================

Map.addLayer(
  K,
  {
    min: 0,
    max: 0.05,
    palette: ['fff7bc', 'fec44f', 'd95f0e']
  },
  'K Factor',
  false
);

Map.addLayer(
  LS,
  {
    min: 0,
    max: 100,
    palette: ['ffffff', '00ffff', 'ffff00', 'ff0000']
  },
  'LS Factor',
  false
);

Map.addLayer(
  SDR,
  {
    min: 0,
    max: 1,
    palette: ['002a7f', 'ffff99', 'cb181d']
  },
  'SDR Factor',
  false
);


// ============================================================================
// 8. SEASONAL SOIL LOSS AND SEDIMENT YIELD FUNCTION
// ============================================================================

function computeSeason(year, season) {

  var isPre = (season === 'Pre');

  // March-April / October-November
  var startDate = ee.Date.fromYMD(
    year,
    isPre ? 3 : 10,
    1
  );

  // Exclusive end date
  var endDate = ee.Date.fromYMD(
    year,
    isPre ? 5 : 12,
    1
  );


  // ------------------------------------------------------------------------
  // 8.1 RAINFALL EROSIVITY FACTOR (R)
  // ------------------------------------------------------------------------

  var rainfall = chirps
    .filterDate(startDate, endDate)
    .select('precipitation')
    .sum()
    .clip(aoi);

  var R = rainfall
    .multiply(0.35)
    .add(38.5)
    .rename('R');


  // ------------------------------------------------------------------------
  // 8.2 COVER MANAGEMENT FACTOR (C)
  // ------------------------------------------------------------------------

  var s2Composite = s2
    .filterBounds(aoi)
    .filterDate(startDate, endDate)
    .mean();

  var ndvi = s2Composite
    .normalizedDifference(['B8', 'B4']);

  var C = ee.Image(1)
    .subtract(ndvi)
    .divide(2)
    .pow(ndvi.add(1))
    .rename('C');


  // ------------------------------------------------------------------------
  // 8.3 SUPPORT PRACTICE FACTOR (P)
  // ------------------------------------------------------------------------

  var dw = dwCollection
    .filterBounds(aoi)
    .filterDate(startDate, endDate)
    .select('label')
    .mode()
    .rename('lulc')
    .clip(aoi);

  var P = dw.addBands(slopeDeg).expression(

    '(lulc < 4) ? 0.8' +

    ' : (lulc == 4) ? ' +

    '((slope < 2) ? 0.6' +
    ' : (slope < 5) ? 0.5' +
    ' : (slope < 8) ? 0.5' +
    ' : (slope < 12) ? 0.6' +
    ' : (slope < 16) ? 0.7' +
    ' : (slope < 20) ? 0.8' +
    ' : 0.9)' +

    ' : ((lulc >= 5 && lulc <= 8) ? 1 : 1)',

    {
      lulc: dw,
      slope: slopeDeg
    }

  ).rename('P').clip(aoi);


  // ------------------------------------------------------------------------
  // 8.4 RUSLE SOIL LOSS
  // ------------------------------------------------------------------------

  var soilLoss = R
    .multiply(K)
    .multiply(LS)
    .multiply(C)
    .multiply(P)
    .rename('soil_loss')
    .clip(aoi);


  // ------------------------------------------------------------------------
  // 8.5 SEDIMENT YIELD
  // ------------------------------------------------------------------------

  var sedimentYield = soilLoss
    .multiply(SDR)
    .rename('sediment_yield')
    .clip(aoi);


  // Return both outputs
  return {

    soilLoss: soilLoss,

    sedimentYield: sedimentYield

  };

}


// ============================================================================
// 9. VISUALIZATION PARAMETERS
// ============================================================================

// Same color range for each year and season

var soilVis = {

  min: 0,

  max: 20,

  palette: [
    '0066cc',
    '00c5cd',
    'a2d729',
    'ffff00',
    'ff8c00',
    'cc0000'
  ]

};

var yieldVis = {

  min: 0,

  max: 20,

  palette: [
    '008000',
    'ffff00',
    'ff9900',
    'ff0000'
  ]

};


// ============================================================================
// 10. OPTIONAL GEOTIFF EXPORT FUNCTION
// ============================================================================

function exportMap(image, variable, year, season) {

  var name = variable + '_' + year + '_' + season;

  Export.image.toDrive({

    image: image.toFloat(),

    description: name,

    fileNamePrefix: name,

    region: aoi,

    scale: 30,

    crs: 'EPSG:32645',

    maxPixels: 1e13,

    fileFormat: 'GeoTIFF'

  });

}


// ============================================================================
// 11. CREATE ALL 32 SEASONAL MAP LAYERS
// ============================================================================

// No reduceRegion() calls inside the loop.
//
// 16 soil-loss maps
// 16 sediment-yield maps
//
// Images are stacked into 32 bands for statistical analysis.

var combinedStack = null;

for (
  var year = START_YEAR;
  year <= END_YEAR;
  year++
) {

  // Seasonal calculations
  var pre = computeSeason(year, 'Pre');

  var post = computeSeason(year, 'Post');


  // ------------------------------------------------------------------------
  // 11.1 SOIL LOSS MAPS
  // ------------------------------------------------------------------------

  Map.addLayer(

    pre.soilLoss,

    soilVis,

    'Soil Loss ' + year + ' (Pre)',

    false

  );

  Map.addLayer(

    post.soilLoss,

    soilVis,

    'Soil Loss ' + year + ' (Post)',

    false

  );


  // ------------------------------------------------------------------------
  // 11.2 SEDIMENT YIELD MAPS
  // ------------------------------------------------------------------------

  Map.addLayer(

    pre.sedimentYield,

    yieldVis,

    'Sediment Yield ' + year + ' (Pre)',

    false

  );

  Map.addLayer(

    post.sedimentYield,

    yieldVis,

    'Sediment Yield ' + year + ' (Post)',

    year === END_YEAR

  );


  // ------------------------------------------------------------------------
  // 11.3 CREATE UNIQUE BANDS FOR EACH OUTPUT
  // ------------------------------------------------------------------------

  var yearBands = pre.soilLoss
    .rename('SL_' + year + '_Pre')

    .addBands(
      post.soilLoss.rename('SL_' + year + '_Post')
    )

    .addBands(
      pre.sedimentYield.rename('SY_' + year + '_Pre')
    )

    .addBands(
      post.sedimentYield.rename('SY_' + year + '_Post')
    );


  // ------------------------------------------------------------------------
  // 11.4 COMBINE ALL IMAGES
  // ------------------------------------------------------------------------

  if (combinedStack === null) {

    combinedStack = yearBands;

  } else {

    combinedStack = combinedStack.addBands(yearBands);

  }


  // ------------------------------------------------------------------------
  // 11.5 OPTIONAL GEOTIFF EXPORT
  // ------------------------------------------------------------------------

  if (EXPORT_MAPS) {

    exportMap(
      pre.soilLoss,
      'Soil_Loss',
      year,
      'Pre'
    );

    exportMap(
      post.soilLoss,
      'Soil_Loss',
      year,
      'Post'
    );

    exportMap(
      pre.sedimentYield,
      'Sediment_Yield',
      year,
      'Pre'
    );

    exportMap(
      post.sedimentYield,
      'Sediment_Yield',
      year,
      'Post'
    );

  }

}


// ============================================================================
// 12. SINGLE AGGREGATION FOR ALL STATISTICS
// ============================================================================

// Mean, minimum and maximum for all 32 bands.
// One reduceRegion() operation avoids separate
// concurrent regional aggregations.

// Combined reducer
var statsReducer = ee.Reducer.mean().combine({

  reducer2: ee.Reducer.minMax(),

  sharedInputs: true

});

// Compute all statistics together
var stats = combinedStack.reduceRegion({

  reducer: statsReducer,

  geometry: aoi,

  scale: ANALYSIS_SCALE,

  maxPixels: 1e13,

  tileScale: 8

});


// ============================================================================
// 13. SOIL LOSS AND SEDIMENT YIELD STATISTICS TABLES
// ============================================================================

// Each table contains:
// Year
// Pre-monsoon mean
// Post-monsoon mean
// Pre-monsoon minimum
// Pre-monsoon maximum
// Post-monsoon minimum
// Post-monsoon maximum

var soilRows = [];

var yieldRows = [];


for (
  var y = START_YEAR;
  y <= END_YEAR;
  y++
) {

  var slPre = 'SL_' + y + '_Pre_';

  var slPost = 'SL_' + y + '_Post_';

  var syPre = 'SY_' + y + '_Pre_';

  var syPost = 'SY_' + y + '_Post_';


  // ------------------------------------------------------------------------
  // 13.1 SOIL LOSS STATISTICS
  // ------------------------------------------------------------------------

  soilRows.push(

    ee.Feature(null, {

      year: y,

      'Pre-monsoon': stats.get(slPre + 'mean'),

      'Post-monsoon': stats.get(slPost + 'mean'),

      Pre_Min: stats.get(slPre + 'min'),

      Pre_Max: stats.get(slPre + 'max'),

      Post_Min: stats.get(slPost + 'min'),

      Post_Max: stats.get(slPost + 'max')

    })

  );


  // ------------------------------------------------------------------------
  // 13.2 SEDIMENT YIELD STATISTICS
  // ------------------------------------------------------------------------

  yieldRows.push(

    ee.Feature(null, {

      year: y,

      'Pre-monsoon': stats.get(syPre + 'mean'),

      'Post-monsoon': stats.get(syPost + 'mean'),

      Pre_Min: stats.get(syPre + 'min'),

      Pre_Max: stats.get(syPre + 'max'),

      Post_Min: stats.get(syPost + 'min'),

      Post_Max: stats.get(syPost + 'max')

    })

  );

}


// FeatureCollections
var soilLossStats = ee.FeatureCollection(
  soilRows
).sort('year');

var sedimentYieldStats = ee.FeatureCollection(
  yieldRows
).sort('year');


// Display statistics
print(
  'Soil Loss Statistics (t/ha/season)',
  soilLossStats
);

print(
  'Sediment Yield Statistics (t/ha/season)',
  sedimentYieldStats
);


// ============================================================================
// 14. SOIL LOSS COMPARISON CHART (2018-2025)
// ============================================================================

var yearTicks = [
  2018, 2019, 2020, 2021,
  2022, 2023, 2024, 2025
];

var soilChart = ui.Chart.feature.byFeature({

  features: soilLossStats,

  xProperty: 'year',

  yProperties: [
    'Pre-monsoon',
    'Post-monsoon'
  ]

})

.setChartType('LineChart')

.setOptions({

  title: 'Mean Seasonal Soil Loss (2018-2025)',

  hAxis: {

    title: 'Year',

    ticks: yearTicks,

    format: '0'

  },

  vAxis: {

    title: 'Mean Soil Loss (t/ha/season)',

    viewWindow: {
      min: 0
    }

  },

  legend: {
    position: 'bottom'
  },

  lineWidth: 3,

  pointSize: 6,

  series: {

    0: {
      color: '#2471a3'
    },

    1: {
      color: '#cb4335'
    }

  }

});

print(
  'Soil Loss: Pre- vs Post-Monsoon',
  soilChart
);


// ============================================================================
// 15. SEDIMENT YIELD COMPARISON CHART (2018-2025)
// ============================================================================

var yieldChart = ui.Chart.feature.byFeature({

  features: sedimentYieldStats,

  xProperty: 'year',

  yProperties: [
    'Pre-monsoon',
    'Post-monsoon'
  ]

})

.setChartType('LineChart')

.setOptions({

  title: 'Mean Seasonal Sediment Yield (2018-2025)',

  hAxis: {

    title: 'Year',

    ticks: yearTicks,

    format: '0'

  },

  vAxis: {

    title: 'Mean Sediment Yield (t/ha/season)',

    viewWindow: {
      min: 0
    }

  },

  legend: {
    position: 'bottom'
  },

  lineWidth: 3,

  pointSize: 6,

  series: {

    0: {
      color: '#2471a3'
    },

    1: {
      color: '#cb4335'
    }

  }

});

print(
  'Sediment Yield: Pre- vs Post-Monsoon',
  yieldChart
);


// ============================================================================
// 16. OPTIONAL CSV EXPORTS
// ============================================================================

if (EXPORT_CSV) {

  // Soil loss statistics
  Export.table.toDrive({

    collection: soilLossStats,

    description:
      'Teesta_Soil_Loss_Stats_2018_2025',

    fileFormat: 'CSV'

  });


  // Sediment yield statistics
  Export.table.toDrive({

    collection: sedimentYieldStats,

    description:
      'Teesta_Sediment_Yield_Stats_2018_2025',

    fileFormat: 'CSV'

  });

}


// ============================================================================
// END OF SCRIPT
// ============================================================================
