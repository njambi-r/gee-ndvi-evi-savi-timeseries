/*
Script loosely based on script by Source of Flavor: https://github.com/Source-of-Flavor/google-earth-engine/blob/main/NDVI_Time_Series.j

Goal:Compare NDVI, EVI, and SAVI 

Context: Made for the Geohub blog
By: Rose Njambi
Date: September 2025

Methodology:
- Define the Area of Interest (AOI)
- Load Sentinel-2 imagery: Level-2A Harmonized Surface Reflectance dataset
- Cloud mask: Sentinel-2 Scene Classification Layer (SCL) with Sentinel-2 Cloud Probability data to mask out clouds, cirrus, and shadows.
- Generate monthly composites: Median values of all cloud-free pixels.
- Calculate Vegetation Indices: NDVI, and EVI, for comparison.
- Normalize: Scaling indices between 0–1 for comparability across months.
- Visualization and analysis:
        - Plot image counts per month to check data availability.
        - Generate line charts showing how NDVI and EVI fluctuated through time.
        - Create animated GIFs of NDVI vs. EVI and the true-color RGB composites.
        - Generate histograms of vegetation indices and image bands.
*/

//=====================================================================================================
//                          MAP VIEW & STUDY AREA DEFINITION
//*****************************************************************************************************

// Create a title label to display on the map
var title = ui.Label({
  value: 'NDVI / EVI / SAVI Time Series',
  style: { fontWeight: 'bold', fontSize: '18px' }
});

// Set the position of the title to top-center
title.style().set('position', 'top-center');

// Add the title label to the map
Map.add(title);

// Define Area of Interest (AOI)

//Karura - Forest area
var aoi = Karura_outline;

//Naivasha - Agricultural area - Irrigate
// var aoi = Morendat;

// Center the map view over the AOI with zoom level 10
Map.centerObject(aoi, 10);

//=====================================================================================================
//                          CLOUD MASKING FUNCTION
//*****************************************************************************************************

// Function to safely add cloud probability band
function addCloudProbability(img) {
  var cloudProb = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
    .filterBounds(img.geometry())
    .filterDate(img.date(), img.date().advance(1, 'day'))
    .first();

  // If cloudProb is missing, add a constant 0 band (no clouds)
  cloudProb = ee.Algorithms.If(
    cloudProb,
    ee.Image(cloudProb).select('probability'),
    ee.Image(0).rename('probability').clip(img.geometry())
  );

  return img.addBands(ee.Image(cloudProb));
}

// Function to mask both clouds and shadows
function maskS2clouds(img) {
  var cloudProb = img.select('probability');
  var scl = img.select('SCL'); // Scene Classification Layer
  
  // Cloud mask (>40% probability)
  // Mask out clouds with >40% probability (tweak threshold if needed)
  var cloudMask = cloudProb.lt(40);
  
  // SCL-based masks
  var shadowMask = scl.neq(3);        // Cloud shadows
  var cirrusMask = scl.neq(10);       // Cirrus clouds
  var cloudMask2 = scl.neq(9);        // Dense clouds
 
  // Combined mask
  var mask = cloudMask
    .and(shadowMask)
    .and(cirrusMask)
    .and(cloudMask2);
              
  // Total pixels: count on a reliable band (e.g. B4)
  var totalPixels = img.select('B4').reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: aoi.geometry(),
    scale: 10,
    maxPixels: 1e13
  });
  
  // Clean pixels: count only where mask == 1
  var cleanPixels = img.select('B4').updateMask(mask).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: aoi.geometry(),
    scale: 10,
    maxPixels: 1e13
  });

  // Get the actual count values (assuming 'probability' band for counting)
  var totalCount = ee.Number(totalPixels.get('B4')).max(1); //Avoid division by zero
  var cleanCount = ee.Number(cleanPixels.get('B4'));
  
  // Calculate contamination percentage
  var contamination = totalCount.subtract(cleanCount)
                                .divide(totalCount)
                                .multiply(100)
                                .max(0)  // Ensure non-negative
                                .min(100); // Cap at 100%

  return img.updateMask(mask) // Apply mask
              .divide(10000) // Scale reflectance to 0–1
              .select("B.*") // Keep only spectral bands
              .copyProperties(img, ["system:time_start"]) // Retain timestamp
              .set('contamination', contamination);   // store % contaminated
}

//=====================================================================================================
//                          MULTI-YEAR MONTHLY INDEX FUNCTION (Safe Normalization)
//*****************************************************************************************************
function monthlyIndicesRange(startYear, endYear, aoi) {
  var years = ee.List.sequence(startYear, endYear);

  var allImages = years.map(function(year) {
    year = ee.Number(year);

    var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(aoi)
      .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year, 12, 31))
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
      .map(addCloudProbability)  // join cloud prob
      .map(maskS2clouds);        // apply probability-based mask
            
    var monthly = ee.List.sequence(1, 12).map(function(month) {
      var start = ee.Date.fromYMD(year, month, 1);
      var end   = start.advance(1, 'month');
      var collection = s2.filterDate(start, end);
      var count      = collection.size();
      var contaminationMean = ee.Algorithms.If(
        collection.size().gt(0),
        collection.aggregate_mean('contamination'),
        0
      );
      var composite  = collection.median().clip(aoi).toFloat();
          
      return ee.Algorithms.If(
        count.gt(0),
        (function() {
          // --- Raw indices ---
          var ndvi = composite.normalizedDifference(['B8','B4']).rename('NDVI').toFloat();
          var evi = composite.expression(
            '2.5 * ((NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1))',
            {NIR: composite.select('B8'), RED: composite.select('B4'), BLUE: composite.select('B2')}
          ).rename('EVI').toFloat(); 
          var savi = composite.expression(
            '(1.5 * (NIR - RED)) / (NIR + RED + 0.5)',
            {NIR: composite.select('B8'), RED: composite.select('B4')}
          ).rename('SAVI').toFloat();

          var stacked = ndvi.addBands([evi, savi]);

          // --- Normalization ---
          //[A good rule of thumb is to set min and max to values 
          //that represent the 2nd and 98th percentile of the data 
          //within your area of interest.
          //Source: https://developers.google.com/earth-engine/guides/ic_visualization]
          
          function normalize(img, bandName) {
            // Compute percentiles, but guard against missing results by using default values.
            var stats = ee.Dictionary(img.reduceRegion({
              reducer: ee.Reducer.percentile([2, 98]),
              geometry: aoi,
              scale: 10,
              maxPixels: 1e13
            }));

            /*
            // Get p2 and p98 safely (use default values if null)
            var min = ee.Number(stats.get(bandName + '_p2', 0));
            var max = ee.Number(stats.get(bandName + '_p98', 1));
            */
            var rawMin = stats.get(bandName + '_p2');
            var rawMax = stats.get(bandName + '_p98');
            
            var min = ee.Number(ee.Algorithms.If(rawMin, rawMin, 0));
            var max = ee.Number(ee.Algorithms.If(rawMax, rawMax, 1));
            
            // Ensure non-zero denominator
            var denom = max.subtract(min).max(0.000001);

            // Operate on the image band (ensure image typed properly)
            var bandImg = ee.Image(img.select([bandName]).toFloat());

            return bandImg.subtract(min).divide(denom)
                      .clamp(0, 1) // Ensure values stay 0-1
                      .rename(bandName + '_Normalized')
                      .toFloat();
          }

          var ndviN = normalize(ndvi, 'NDVI');
          var eviN  = normalize(evi, 'EVI');
          var saviN = normalize(savi, 'SAVI');
          var rgbBands = composite.select(['B4','B3','B2']);
          var NIR = composite.select(['B8']);
          
          return stacked.addBands([ndviN, eviN, saviN, rgbBands, NIR]) //B8(NIR), added to study saturation
            .set('year', year)
            .set('month', month)
            .set('count', count)
            .set('contamination', contaminationMean)   // 👈 add this
            .set('noData', 0)
            .set('system:time_start', ee.Date.fromYMD(year, month, 1).millis());
        })(), //end true branch

        // --- No-data branch: return constant 0 bands to keep schema consistent ---
        ee.Image.constant([0,0,0,0,0,0,0,0,0,0])
          .rename(['NDVI','EVI','SAVI',
                   'NDVI_Normalized','EVI_Normalized','SAVI_Normalized',
                   'B4','B3','B2','B8' ])
          .updateMask(ee.Image(0))   // mask out everything
          .toFloat() // Ensure Float type to enable date ranges before 2021 to work
          .set('year', year)
          .set('month', month)
          .set('count', 0)
          .set('contamination', 0)   // 👈 add this
          .set('noData', 1)
          .set('system:time_start', ee.Date.fromYMD(year, month, 1).millis())
      ); //end If
    }); //end monthly map

    return ee.ImageCollection.fromImages(monthly);
  }); // end years.map

  return ee.ImageCollection(allImages.iterate(function(ic, prev) {
    return ee.ImageCollection(prev).merge(ic);
  }, ee.ImageCollection([])));
} // end monthly indices range

//=====================================================================================================
//                          EXPORT MONTHLY RAW IMAGES TO DRIVE (NDVI/EVI/SAVI + RGB)
//*****************************************************************************************************

// Function to export monthly NDVI/EVI/SAVI images
function exportToDrive(ic, folderName) {
  ic.evaluate(function(images) {
    images.features.forEach(function(f) {
      var year = f.properties.year;
      var month = f.properties.month;
      var noData = f.properties.noData;
      if (noData !== 1) {
        var image = ee.Image(f.id);
        Export.image.toDrive({
          image: image,
          description: 'Indices_RGB_' + year + '_' + month,
          folder: folderName,
          fileNamePrefix: 'Indices_RGB_' + year + '_' + month,
          region: aoi,
          scale: 10,
          maxPixels: 1e13
        });
        print('✅ Exporting NDVI for', year + '-' + month);
      } else {
        print('⚠️ No data for ' + year + '-' + month + ', skipping export.');
      }
    });
  });
}

//=====================================================================================================
//                          ADD TIMESTAMP LABEL ON EACH GIF FRAME
//*****************************************************************************************************
// Load Gena's text module
var text = require('users/gena/packages:text');

// Helper: safely place text inside the GIF region (bottom-left here)
function timestampLabel(image, region4326) {
  // Location in % of region bounds (left 3%, down 95%)
  var pt = text.getLocation(region4326, 'left', '1%', '5%');

  // Options object is REQUIRED by text.draw()
  var opts = {
    fontSize: 18,
    textColor: '000000',   // ffffff white
    outlineColor: 'ffffff',// black
    outlineWidth: 3,
    outlineOpacity: 0.6
  };

  var dateStr = ee.Date(image.get('system:time_start')).format('YYYY-MMM');
  return text.draw(dateStr, pt, 20, opts); // 1000 = text scale (m/px). Tweak if you want larger/smaller text.
}

//=====================================================================================================
//                          GENERATE INDICES GIF FUNCTION
//*****************************************************************************************************

// Function to create an animated GIF from monthly NDVI/EVI/SAVI images
function createGif(ic, aoi, band, title) {
  var region4326 = aoi.geometry().transform('EPSG:4326', 1);
  
  //Option 1: Normalised values 0-1 
  var vis = {min: 0, max: 1, palette: ['#d9a679', '#ffffb2', '#78c679', '#238443']};
  
  //Option 2: Normalized for Comparability based on indices max/min values
  //var vis = {min: 0.05, max: 0.25, palette: ['#d9a679', '#ffffb2', '#78c679', '#238443']};
  
  var gifCollection = ic
    .filter(ee.Filter.neq('noData', 1))
    .sort('system:time_start')
    .map(function(img) {
      var bandImg = img.select(band);
      
      // Count unmasked pixels in AOI
      var count = bandImg.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: aoi,
        scale: 10,
        maxPixels: 1e13
      }).getNumber(band);
      
      // Only keep frames with >0 pixels
      return ee.Algorithms.If(
        count.gt(0),
        (function() {
          var visImg = bandImg.visualize(vis);
          var label = timestampLabel(img, region4326);
          return visImg.blend(label).set('system:time_start', img.get('system:time_start'));
        })(),
        null
      );
    })
    .filter(ee.Filter.notNull(['system:time_start'])); // remove dropped frames

  var gifParams = {
    region: region4326,
    dimensions: 600,
    crs: 'EPSG:4326',
    framesPerSecond: 1,
    format: 'gif'
  };

  gifCollection.size().evaluate(function(n) {
    if (n > 0) {
      print(ui.Thumbnail(gifCollection, gifParams, title + ' Animation (with date)'));
      print('📥 Download ' + title + ' GIF:', gifCollection.getVideoThumbURL(gifParams));
    } else {
      print('⚠️ No ' + title + ' data available for animation.');
    }
  });
}

//====================================================================================
//             GENERATE RGB GIF ALIGNED TO INDICES
//====================================================================================
function createRgbGif(ic, aoi, title) {
  var region4326 = aoi.geometry().transform('EPSG:4326', 1);

  // Sentinel-2 true color visualization
  var visRgb = {bands: ['B4','B3','B2'], min: 0, max: 0.3, gamma:1.2};

  var rgbCollection = ic
    .filter(ee.Filter.neq('noData', 1))   // skip months flagged as dummy
    .sort('system:time_start')
    .map(function(img) {
      var rgbImg = img.select(['B4','B3','B2']);

      // Visualize directly (no need for reduceRegion check anymore)
      var visImg = rgbImg.visualize(visRgb);
      var label  = timestampLabel(img, region4326);
      return visImg.blend(label)
                   .set('system:time_start', img.get('system:time_start'));
    });

  var gifParams = {
    region: region4326,
    dimensions: 600,
    crs: 'EPSG:4326',
    framesPerSecond: 1,
    format: 'gif'
  };

  rgbCollection.size().evaluate(function(n) {
    if (n > 0) {
      print(ui.Thumbnail(rgbCollection, gifParams, title + ' RGB Animation (with date)'));
      print('📥 Download ' + title + ' RGB GIF:', rgbCollection.getVideoThumbURL(gifParams));
    } else {
      print('⚠️ No ' + title + ' RGB data available for animation.');
    }
  });
}

//====================================================================================
//   AUTO-TILED MULTI-PANEL COMPARISON GIF
//====================================================================================

// ✅ Visualization styles

/*
//Option 1: Normalised
var visSat  = {bands:['B4','B3','B2'], min:0, max:0.3, gamma:1.3};
var visNDVI = {bands:['NDVI_Normalized'], min:0, max:1, palette:['#d9a679', '#ffffb2', '#78c679', '#238443']};
var visEVI  = {bands:['EVI_Normalized'],  min:0, max:1, palette:['#d9a679', '#ffffb2', '#78c679', '#238443']};
var visSAVI = {bands:['SAVI_Normalized'], min:0, max:1, palette:['#d9a679', '#ffffb2', '#78c679', '#238443']};
*/

//Option 2: Raw
var visSat  = {bands:['B4','B3','B2'], min:0, max:0.3, gamma:1.2};
var visNDVI = {bands:['NDVI'], min:0, max:1, palette:['#d9a679', '#ffffb2', '#78c679', '#238443']};
var visEVI  = {bands:['EVI'],  min:0, max:1, palette:['#d9a679', '#ffffb2', '#78c679', '#238443']};
var visSAVI = {bands:['SAVI'], min:0, max:1, palette:['#d9a679', '#ffffb2', '#78c679', '#238443']};

function createHorizontalGif(ic, title, panelTypes) {
  // panelTypes should be an array like ['sat', 'evi'] or ['sat', 'savi'] or ['sat', 'ndvi']
  var nPanels = panelTypes.length;
  var region4326 = aoi.geometry().transform('EPSG:4326', 1);

  // AOI bounds
  var bounds = region4326.bounds();
  var coords = ee.List(bounds.coordinates().get(0));
  var xmin = ee.Number(ee.List(coords.get(0)).get(0));
  var ymin = ee.Number(ee.List(coords.get(0)).get(1));
  var xmax = ee.Number(ee.List(coords.get(2)).get(0));
  var ymax = ee.Number(ee.List(coords.get(2)).get(1));
  var dx = xmax.subtract(xmin);

  // Expanded region for output
  var outRegion = ee.Geometry.Rectangle([
    xmin, ymin,
    xmin.add(dx.multiply(nPanels)), ymax
  ]);

  var gifParams = {
    region: outRegion,
    dimensions: 500 * nPanels,
    framesPerSecond: 1,
    crs: 'EPSG:4326'
  };

  var combo = ic.filter(ee.Filter.neq('noData', 1))
                .sort('system:time_start')
                .map(function(img) {
    var label = timestampLabel(img, region4326);

    // Create visualization lookup
    var visualizations = {
      'sat':  img.visualize(visSat).clip(region4326),
      'ndvi': img.visualize(visNDVI).clip(region4326),
      'evi':  img.visualize(visEVI ).clip(region4326),
      'savi': img.visualize(visSAVI).clip(region4326)
    };

    // Create label lookup
    var labelNames = {
      'sat': 'RGB',
      'ndvi': 'NDVI', 
      'evi': 'EVI',
      'savi': 'SAVI'
    };

    // Build panels dynamically
    var panelsToMosaic = [];
    
    for (var i = 0; i < nPanels; i++) {
      var panelType = panelTypes[i];
      
      // Get the visualization for this panel
      var panel = visualizations[panelType];
      
      // Shift panel to correct position
      var shiftedPanel = panel.changeProj(panel.projection(), panel.projection().translate(dx.multiply(i), 0));
      panelsToMosaic.push(shiftedPanel);
      
      // Add frame for this panel
      var panelBounds = ee.Geometry.Rectangle([
        xmin.add(dx.multiply(i)), ymin,
        xmin.add(dx.multiply(i+1)), ymax
      ]);
      var panelOutline = panelBounds.bounds();
      var panelFrame = ee.Image.constant([255,255,255]) 
        .updateMask(ee.Image.constant(0).paint(panelOutline, 1, 3))
        .rename(['vis-red','vis-green','vis-blue']);
      panelsToMosaic.push(panelFrame);
      
      // Create and position label for this panel
      // White text for SAT/RGB, black text for vegetation indices
      var labelColor = panelType === 'sat' ? [255,255,255] : [255,255,255]; //[0,0,0] for black
      var panelLabel = ee.Image.constant(labelColor).updateMask(
        text.draw(labelNames[panelType], 
                 ee.Geometry.Point([xmin.add(dx.multiply(i + 0.5)), ymin.add((ymax.subtract(ymin)).multiply(0.1))]), 24)
      ).rename(['vis-red','vis-green','vis-blue']);
      panelsToMosaic.push(panelLabel);
    }
    
    // Add the timestamp label LAST so it stays on top
    panelsToMosaic.push(label);
    
    /* -- displaying but wierd position -- fix how it's displaying later
    // Add legends for vegetation indices
    for (var j = 0; j < nPanels; j++) {
      var panelTypeLg = panelTypes[j];
      
      if (panelTypeLg !== 'sat') {
        // Create color bar legend for vegetation indices
        var legendX = xmin.add(dx.multiply(j)).add(dx.multiply(0.05)); // 5% from left edge of panel
        var legendY = ymax.subtract((ymax.subtract(ymin)).multiply(0.2)); // 20% from bottom
        var legendWidth = dx.multiply(0.2); // 20% of panel width
        var legendHeight = (ymax.subtract(ymin)).multiply(0.15); // 15% of panel height
        
        // Get min/max for this index
        var visParams = panelTypeLg === 'ndvi' ? visNDVI : 
                       panelTypeLg === 'evi' ? visEVI : visSAVI;
        var legendMin = visParams.min;
        var legendMax = visParams.max;
        
        // Create color gradient for legend
        var legendGradient = ee.Image.pixelLonLat()
          .select('latitude')
          .unitScale(legendY, legendY.add(legendHeight))
          .visualize({
            min: 0, 
            max: 1, 
            palette: visParams.palette
          })
          .clip(ee.Geometry.Rectangle([
            legendX, legendY,
            legendX.add(legendWidth), legendY.add(legendHeight)
          ]));
        
        panelsToMosaic.push(legendGradient);
        
        // Add legend labels (min/max values)
        var minLabel = ee.Image.constant([255,255,255]).updateMask(
          text.draw(ee.Number(legendMin).format('%.2f'), 
                   ee.Geometry.Point([legendX.add(legendWidth.multiply(1.1)), legendY]), 12)
        ).rename(['vis-red','vis-green','vis-blue']);
        
        var maxLabel = ee.Image.constant([255,255,255]).updateMask(
          text.draw(ee.Number(legendMax).format('%.2f'), 
                   ee.Geometry.Point([legendX.add(legendWidth.multiply(1.1)), legendY.add(legendHeight)]), 12)
        ).rename(['vis-red','vis-green','vis-blue']);
        
        panelsToMosaic.push(minLabel);
        panelsToMosaic.push(maxLabel);
      }
    }
    */

    return ee.ImageCollection(panelsToMosaic).mosaic()
            .set('system:time_start', img.get('system:time_start'));
  });

  print(ui.Thumbnail(combo, gifParams, title + ' Comparison GIF'));
  print('📥 Download ' + title + ' GIF:', combo.getVideoThumbURL(gifParams));
}

//=====================================================================================================
//                          TIME SERIES CHART FUNCTION
//*****************************************************************************************************

// Function to plot values as a time series chart

// Function to plot values as a time series chart
function createTimeSeriesChart(ic, band, title, color) {
  // Mask out dummy "noData" months before charting
  var maskedIC = ic.map(function(img) {
    var noData = ee.Number(img.get('noData'));
    return ee.Algorithms.If(
      noData.eq(1),
      img.updateMask(ee.Image(0)), // fully masked => chart sees null
      img
    );
  });

  var chart = ui.Chart.image.series({
    imageCollection: ee.ImageCollection(maskedIC).select(band),
    region: aoi,
    reducer: ee.Reducer.mean(),
    scale: 10,
    xProperty: 'system:time_start'
  }).setOptions({
    lineWidth: 2,
    pointSize: 6,
    title: title,
    interpolateNulls: true,
    vAxis: {
      title: band,
      viewWindow: {min: 0, max: 1},
      gridlines: {count: 11}
    },
    hAxis: {
      title: 'Date',
      format: 'YYYY-MMM',
      slantedText: true,
      slantedTextAngle: 90
    },
    series: {0: {color: color, lineWidth: 3}}
  });

  print(chart);
}

// Function to plot all 3 indices in one time series chart
function createCombinedChart(ic, aoi) {
  // Mask out dummy "noData" months
  var maskedIC = ic.map(function(img) {
    var noData = ee.Number(img.get('noData'));
    return ee.Algorithms.If(
      noData.eq(1),
      img.updateMask(ee.Image(0)), // fully masked => null in chart
      img
    );
  });

  // Select all three bands
  var chart = ui.Chart.image.series({
    imageCollection: ee.ImageCollection(maskedIC).select(['NDVI', 'EVI']),//, 'SAVI_Normalized']),
    region: aoi,
    reducer: ee.Reducer.mean(),
    scale: 10,
    xProperty: 'system:time_start'
  }).setChartType('LineChart')
    .setOptions({
      title: 'NDVI and EVI Time Series',//, and SAVI Time Series',
      lineWidth: 2,
      pointSize: 5,
      interpolateNulls: false,
      vAxis: {
        title: 'Vegetation Index Value',
        viewWindow: {min: 0, max: 1},
        gridlines: {count: 11}
      },
      hAxis: {
        title: 'Date',
        format: 'YYYY-MMM',
        slantedText: true,
        slantedTextAngle: 90
      },
      series: {
        0: {color: '#32CD32', pointSize: 5, lineWidth: 3}, // EVI (lime green)
        1: {color: '#006400', pointSize: 5, lineWidth: 3}, // NDVI (forest green)
        //2: {color: '#4575b4', lineWidth: 3}  // SAVI (blue)
      }
    });

  print(chart);
}

// Function to plot VIs and Reflectances in a combined chart
function createVIReflectances(ic, aoi) {
  // Mask out dummy "noData" months
  var maskedIC = ic.map(function(img) {
    var noData = ee.Number(img.get('noData'));
    return ee.Algorithms.If(
      noData.eq(1),
      img.updateMask(ee.Image(0)), // fully masked => null in chart
      img
    );
  });
  // Select all three bands
  var chart = ui.Chart.image.series({
    imageCollection: ee.ImageCollection(maskedIC).select(['NDVI', 'EVI','B4', 'B8']),//, 'SAVI_Normalized']),
    region: aoi,
    reducer: ee.Reducer.mean(),
    scale: 10,
    xProperty: 'system:time_start'
  }).setChartType('LineChart')
    .setOptions({
      title: 'NDVI, EVI, B4-Red, and B8-NIR Time Series',//, and SAVI Time Series'
      interpolateNulls: false,
      vAxis: {
        title: 'VI and Reflectances',
        viewWindow: {min: 0, max: 1},
        gridlines: {count: 11}
      },
      hAxis: {
        title: 'Date',
        format: 'YYYY-MMM',
        slantedText: true,
        slantedTextAngle: 90
      },
      series: {
        0: {color: '#d73027', lineDashStyle: [4, 2], pointSize: 5, lineWidth:2},  // B4 (Red)
        1: {color: '#800080', lineDashStyle: [4, 2], pointSize: 5, lineWidth: 2},  // B8 (Purple)
        2: {color: '#32CD32', pointSize: 5, lineWidth: 3}, // EVI (lime green)
        3: {color: '#006400', pointSize: 5, lineWidth: 3}, // NDVI (forest green)
      }
    });
  print(chart);
}

//=====================================================================================================
//                          IMAGE COUNT CHART FUNCTION
//*****************************************************************************************************

// Function to plot a bar chart showing number of images used per month
function createImageCountChart(ic) {
  var chart = ui.Chart.feature.byFeature(ic, 'system:time_start', ['count'])
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Monthly Image Count',
      hAxis: {
        title: 'Date',
        format: 'YYYY-MM',
        slantedText: true,
        slantedTextAngle: 90
      },
      vAxis: { title: 'Image Count' },
      legend: { position: 'none' },
      colors: ['#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd']
    });
  print(chart);
}

//=====================================================================================================
//                          IMAGE CONTAMINATION WITH CLOUD/SHADOW CHART FUNCTION
//*****************************************************************************************************
function createContaminationChart(ic) {
  var chart = ui.Chart.feature.byFeature(ic, 'system:time_start', ['contamination'])
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Monthly Cloud/Shadow Contamination (%)',
      hAxis: {
        title: 'Date',
        format: 'YYYY-MM',
        slantedText: true,
        slantedTextAngle: 90
      },
      vAxis: { 
        title: 'Contamination (%)',
        //minValue: 0,
        //maxValue: 100
      },
      legend: { position: 'none' },
      colors: ['#d62728']
    });
  print(chart);
  return chart;
}

// ==================================================================
// 📊 SEPARATE HISTOGRAMS PER MONTH (skip empty months, auto bucket size)
// ==================================================================
function monthlyHistograms(ic, band, title, xMin, xMax) {
  var monthNames = [
    'Jan','Feb','Mar','Apr','May','Jun',
    'Jul','Aug','Sep','Oct','Nov','Dec'
  ];
  
  // Loop client-side over 12 months
  for (var m = 1; m <= 12; m++) {
    var monthCol = ic.filter(ee.Filter.eq('month', m))
                    .filter(ee.Filter.neq('noData', 1)); // Filter out no-data images
    var monthName = monthNames[m - 1];
    
    // Check if we have valid data for this month
    var collectionSize = monthCol.size().getInfo();
    
    if (collectionSize > 0) {
      try {
        var img = monthCol.first();
        
        // Additional check: count unmasked pixels to ensure data exists
        var pixelCount = img.select(band).reduceRegion({
          reducer: ee.Reducer.count(),
          geometry: aoi,
          scale: 10,
          maxPixels: 1e13
        }).getNumber(band).getInfo();
        
        if (pixelCount > 0) {
          var chart = ui.Chart.image.histogram({
            image: img.select(band),
            region: aoi,
            scale: 10,
            maxPixels: 1e13
          })
          .setSeriesNames([band])
          .setOptions({
            title: title + ' - ' + monthName,
            hAxis: {title: band, viewWindow: {min: xMin, max: xMax}},
            vAxis: {title: 'Pixel Count'},
            legend: {position: 'none'},
            histogram: {bucketSize: 0.001}   // 👈 force fixed bin size
          });
          print(chart);
        } else {
          print('⚠️ No valid pixels for ' + title + ' (' + monthName + ').');
        }
        
      } catch (error) {
        print('⚠️ Error:' + title + ' (' + monthName + '): ' + error.message);
      }
    } else {
      print('⚠️ No ' + title + ' data available(' + monthName + ').');
    }
  }
}

//=====================================================================================================
//                          EXECUTION BLOCK: MULTI-YEAR RANGE
//*****************************************************************************************************
// Define start and end year of analysis
var startYear = 2024;
var endYear   = 2024;

//See the image collection
var indices = monthlyIndicesRange(startYear, endYear, aoi);

print('Indices collection:', indices);
print('Size:', indices.size());
print('First image:', indices.first());

// EXPORT all valid images to Google Drive
exportToDrive(indices, 'Indices_RGB_Exports');

// CHARTS
createTimeSeriesChart(indices, 'NDVI_Normalized', 'Normalized NDVI Time Series', '#2E8B57');
createTimeSeriesChart(indices, 'EVI_Normalized',  'Normalized EVI Time Series', '#1f77b4');
createTimeSeriesChart(indices, 'SAVI_Normalized', 'Normalized SAVI Time Series', '#ff7f0e');
createCombinedChart(indices, aoi); //combine all VIs time series
createVIReflectances(indices, aoi); //Combine VIs + Band reflectances
createImageCountChart(indices), //number of images used per month
createContaminationChart(indices); //Percentage of image scene contaminated with clouds or shadow

// GIFs
// Individual GIFS
createGif(indices, aoi, 'NDVI', 'NDVI');
createGif(indices, aoi, 'EVI',  'EVI');
createGif(indices, aoi, 'SAVI', 'SAVI');
createRgbGif(indices, aoi, 'Sentinel-2 RGB');// aligned RGB

// Tiled GIFS
// Pass in visualizations you want to compare
// Usage examples for any number of panels:
// 2 panels:
// createHorizontalGif(indices, 'Sat vs EVI', ['sat', 'evi']);
// createHorizontalGif(indices, 'Sat vs SAVI', ['sat', 'savi']);
createHorizontalGif(indices, 'NDVI vs EVI', ['ndvi', 'evi']);
 createHorizontalGif(indices, 'Sat vs NDVI', ['sat', 'ndvi']);

// 3 panels:
 createHorizontalGif(indices, 'Sat NDVI EVI', ['sat', 'ndvi', 'evi']);
// createHorizontalGif(indices, 'Indices Only', ['ndvi', 'evi', 'savi']);

// 4 panels:
 createHorizontalGif(indices, 'All Indices', ['sat', 'ndvi', 'evi', 'savi']);

// 5 panels (if you have more bands):
// createHorizontalGif(indices, 'Extended', ['sat', 'ndvi', 'evi', 'savi', 'another_band']);


//STATS -- To inform indices visualization min/max
// Extract percentiles for visualization

// 1) median image across the whole time series
var medianImg = indices.select(['NDVI','EVI','SAVI']).median();

// 2) compute stats (min/max + percentiles) on that median image
var stats = ee.Dictionary(medianImg.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.percentile([2,5,10,20,30,50,70,95,98]),
    sharedInputs: true
  }),
  geometry: aoi,
  scale: 10,
  maxPixels: 1e13
}));

var ndviMin = ee.Number(stats.get('NDVI_p2'));
var ndviMax = ee.Number(stats.get('NDVI_p98'));
var eviMin = ee.Number(stats.get('EVI_p2'));
var eviMax = ee.Number(stats.get('EVI_p98'));
var saviMin = ee.Number(stats.get('SAVI_p2'));
var saviMax = ee.Number(stats.get('SAVI_p98'));

print('NDVI percentile range:', ndviMin, ndviMax);
print('EVI percentile range:', eviMin, eviMax);
print('SAVI percentile range:', saviMin, saviMax);

// Print results
print('Stats from median composite (NDVI, EVI, SAVI):', stats);

//Histograms
var chart = ui.Chart.image.histogram({
  image: indices.select('NDVI').median(),
  region: aoi,
  scale: 10,
  maxPixels: 1e13
})
.setOptions({
  title: 'NDVI Histogram',
  vAxis: {title: 'Frequency'},
  hAxis: {
    title: 'NDVI',
    viewWindow: {min: 0, max: 1}   // 👈 clamp y-axis between 0 and 1
  },
  legend: {position: 'none'}
});
print(chart);

var chart = ui.Chart.image.histogram({
  image: indices.select('EVI').median(),
  region: aoi,
  scale: 10,
  maxPixels: 1e13
}).setOptions({
  title: 'EVI Histogram',
  vAxis: {title: 'Frequency'},
  hAxis: {
    title: 'EVI',
    viewWindow: {min: 0, max: 1}   // 👈 clamp y-axis between 0 and 1
  },
  legend: {position: 'none'}
});
print(chart);

var chart = ui.Chart.image.histogram({
  image: indices.select('SAVI').median(),
  region: aoi,
  scale: 10,
  maxPixels: 1e13
}).setOptions({
  title: 'SAVI Histogram',
  vAxis: {title: 'Frequency'},
  hAxis: {
    title: 'SAVI',
    viewWindow: {min: 0, max: 1}   // 👈 clamp y-axis between 0 and 1
  },
  legend: {position: 'none'}
});
print(chart);

// Call histograms
monthlyHistograms(indices, 'B8',   'Monthly Histogram - NIR (Band 8)', 0, 0.6);
monthlyHistograms(indices, 'B4',   'Monthly Histogram - Red (Band 4)', 0, 0.6);
monthlyHistograms(indices, 'NDVI', 'Monthly Histogram - NDVI', 0, 1.0);
monthlyHistograms(indices, 'EVI',  'Monthly Histogram - EVI', 0, 1.0);

//Mean
//EVI Stats
var meanImgEvi = indices.select('EVI').mean();

print('EVI stats (mean image)', meanImgEvi.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.percentile([2, 5, 10,20,30, 50, 70, 95, 98]),
    sharedInputs: true
  }),
  geometry: aoi,
  scale: 10,
  maxPixels: 1e13
}));

//NDVI Stats
var meanImgNdvi = indices.select('NDVI').mean();

print('NDVI stats (mean image)', meanImgNdvi.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.percentile([2, 5, 10,20, 30, 50, 70, 95, 98]),
    sharedInputs: true
  }),
  geometry: aoi,
  scale: 10,
  maxPixels: 1e13
}));

//SAVI Stats
var meanImgSavi = indices.select('SAVI').mean();

print('SAVI stats (mean image)', meanImgSavi.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.percentile([2, 5,10,20, 30, 50, 70, 95, 98]),
    sharedInputs: true
  }),
  geometry: aoi,
  scale: 10,
  maxPixels: 1e13
}));


//TEST
// Test individual band visualization to confirm they work:
function testVisualization(ic) {
  var testImg = indices.first();
  var region4326 = aoi.geometry().transform('EPSG:4326', 1);
  
  print('Testing EVI visualization:');
  Map.addLayer(testImg.select('EVI').clip(region4326), visEVI, 'EVI Test');
  
  print('Testing SAVI visualization:');
  Map.addLayer(testImg.select('SAVI').clip(region4326), visSAVI, 'SAVI Test');
}

testVisualization(indices);

//SATURATION IN B4(RED) AND B8(NIR) BANDS 
// Generating histograms like in Huete (2002) paper
var chart = ui.Chart.image.histogram({
  image: indices.select('B4').median(),
  region: aoi,
  scale: 10,
  maxPixels: 1e13
}).setOptions({
  title: 'B4 Histogram',
  vAxis: {title: 'Frequency'},
  hAxis: {
    title: 'B4 - Red',
    //viewWindow: {min: 0, max: 1}   // 👈 clamp y-axis between 0 and 1
  },
  legend: {position: 'none'}
});
print(chart);

var chart = ui.Chart.image.histogram({
  image: indices.select('B8').median(),
  region: aoi,
  scale: 10,
  maxPixels: 1e13
}).setOptions({
  title: 'B8 Histogram',
  vAxis: {title: 'Frequency'},
  hAxis: {
    title: 'B8 - NIR',
    //viewWindow: {min: 0, max: 1}   // 👈 clamp y-axis between 0 and 1
  },
  legend: {position: 'none'}
});
print(chart);