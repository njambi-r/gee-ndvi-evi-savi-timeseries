// Add aoi layer to the map
Map.addLayer(aoi.style({'fillColor': '00000000'}));
Map.centerObject(aoi, 12);
Map.setOptions("ROADMAP");

// Define time of interest
var start_date = ee.Date("2024-01-01");
var end_date = ee.Date("2024-12-31");

// Rainfall from  CHIRPS (Climate Hazards Group InfraRed Precipitation with Station) Data
// Paper reference: https://www.mdpi.com/2072-4292/12/4/709
var rain_bioclim = ee.Image('WORLDCLIM/V1/BIO').select(['bio12']).rename(['rainfall_bioclim']).clip(aoi);

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD');
var filtered = chirps
  .filter(ee.Filter.date(start_date, end_date));
var rain = filtered.reduce(ee.Reducer.sum()); // Calculate total precipitation.
var rain = ee.Image(rain).rename(['rainfall']).clip(aoi);

//filtering chirps data for wet and dry seasons
//dry season
var rain_dry = chirps
  .filter(ee.Filter.date('2024-06-01', '2024-10-31'));
var rain_dry = rain_dry.reduce(ee.Reducer.sum());
var rain_dry = ee.Image(rain_dry).rename(['rainfall_dry']).clip(aoi);

//wet season
var rain_wet = chirps
  .filter(ee.Filter.date('2024-03-01', '2024-05-31'));
var rain_wet = rain_wet.reduce(ee.Reducer.sum());
var rain_wet = ee.Image(rain_wet).rename(['rainfall_wet']).clip(aoi);


// Add all the new layers to the imageCollection
var annual_rain_composite = rain.addBands(rain_dry).addBands(rain_wet)
.addBands(rain_bioclim)
.clip(aoi);

// Add topography and rainfall layers to the map

Map.addLayer(annual_rain_composite, {min: 0, max: 2000, bands: ['rainfall_dry'], palette: ['White', 'Blue', 'Red']},'rainfall_dry');
Map.addLayer(annual_rain_composite, {min: 0, max: 2000, bands: ['rainfall_wet'], palette: ['White', 'Blue', 'Red']}, 'rainfall_wet');

print("annual_rain_composite:", annual_rain_composite);

// --- Monthly Rainfall Chart ---

// Function to sum monthly rainfall
function monthlyRainfall(year, aoi) {
  var months = ee.List.sequence(1, 12);
  return ee.ImageCollection.fromImages(
    months.map(function(m) {
      var start = ee.Date.fromYMD(year, m, 1);
      var end   = start.advance(1, 'month');
      var monthly = chirps.filterDate(start, end)
                          .sum()
                          .clip(aoi)
                          .set('month', m)
                          .set('system:time_start', start);
      return monthly;
    })
  );
}

// Apply for 2024
var monthlyImages = monthlyRainfall(2024, aoi);

// Month names lookup
var monthNames = ee.List([
  'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'
]);

// Reduce over AOI and add month names
var monthlySeries = monthlyImages.map(function(img) {
  var stats = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 5000,
    maxPixels: 1e13
  });
  var m = ee.Number(img.get('month')).subtract(1); // zero-based index
  var monthName = monthNames.get(m);
  return ee.Feature(null, stats
    .set('month', img.get('month'))
    .set('monthName', monthName));
});

// Create chart with month names on x-axis
var chart = ui.Chart.feature.byFeature({
  features: monthlySeries.sort('month'),
  xProperty: 'monthName',
  yProperties: ['precipitation']
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Monthly Total Rainfall (CHIRPS, 2024)',
  hAxis: {title: 'Month'},
  vAxis: {title: 'Rainfall (mm)'},
  legend: {position: 'none'},
  colors: ['#1f77b4']
});

print(chart);

