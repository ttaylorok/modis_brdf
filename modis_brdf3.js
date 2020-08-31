var brdf = require('users/ttaylorok/amazon:modis_kernels');

var START_DATE = '01-01'// MM-DD
var END_DATE = '02-01' // MM-DD
var YEARS = ['2018']
var COMBINED = true // combine terra dna aqua
var MIN_SAMPLES = 10
//var PERCENTILE_FILT = [5,95]

var getData = function(start, end, year, combined){
  var c = ee.ImageCollection('MODIS/006/MOD09GA');
  if(COMBINED){
    c = c.merge('MODIS/006/MYD09GA')
  }
  return c.filterDate(year+'-'+start, year+'-'+end);
}

var collection = getData(START_DATE, END_DATE, YEARS[0], COMBINED);
for(var i = 1; i < YEARS.length; i++){
  collection = collection.merge(getData(START_DATE, END_DATE, YEARS[i], COMBINED))
}
var proj = collection.first().select('sur_refl_b01').projection()
print(proj)

var qcBits = function(image){
  var cloudBitMask = 3;
  var shadowBitMask = 4;
  var adjacentBitMask = 1 << 13;
  var cloudInternalBitMask = 1 << 10;
  var snowBitMask = 1 << 15;

  var qa = image.select('state_1km');

  var cloudMask = qa.bitwiseAnd(cloudBitMask).eq(0).rename('cloudMask');
  var shadowMask = qa.bitwiseAnd(shadowBitMask).eq(0).rename('shadowMask');
  var adjacentMask = qa.bitwiseAnd(adjacentBitMask).eq(0).rename('adjacentMask');
  var cloudInternalMask = qa.bitwiseAnd(cloudInternalBitMask).eq(0).rename('cloudInternalMask');
  var snowMask = qa.bitwiseAnd(snowBitMask).eq(0).rename('snowMask');
  //var mask = cloudMask.and(shadowMask).and(adjacentMask).and(cloudInternalMask);
  
  var qc500 = image.select('QC_500m');
  
  var qaBitMask = 3;
  var atmBitMask = 1 << 30;
  var qab1BitMask = 60;
  
  var qaMask = qc500.bitwiseAnd(qaBitMask).eq(0).rename('qaMask');
  var atmMask = qc500.bitwiseAnd(atmBitMask).eq(0).rename('atmMask');
  var qab1Mask = qc500.bitwiseAnd(qab1BitMask).eq(0).rename('qab1Mask');
  
  var totalMask = shadowMask.and(adjacentMask).and(cloudInternalMask).and(snowMask).rename('totalMask');

  return image.addBands([cloudMask,shadowMask,adjacentMask,cloudInternalMask,
    qaMask, atmMask, qab1Mask, totalMask, snowMask]);
}

var collectionQC = collection.map(qcBits);

var visParams = {bands : ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'],
                min : 0, max : 2000};

// 2018-01-27 has two images
Map.addLayer(collection,visParams,'collection');
//.filterDate('2018-01-27')
// point -94.95904983402565,34.621335231913854
// months 01-04 2018
//.first()
//  .filterMetadata('system:index','equals',0)
//.reduce(ee.Reducer.last())
Map.addLayer(collectionQC,{bands:['qab1Mask'],min:0,max:1},'collectionQC');
//.filterDate('2018-01-27')
//.first()
//.reduce(ee.Reducer.last())
//var collection = collection2018//.merge(collection2019).sort('system:time_start');

var chart_qc = ui.Chart.image.doySeries(
  collectionQC.select(['qaMask']),
  mypoint, ee.Reducer.first(),10,ee.Reducer.first())
  .setChartType('ScatterChart')
//print(chart_qc)

var corrected = ee.ImageCollection('MODIS/006/MCD43A4')
  .filterDate('2018-01-01','2018-04-01');
  
Map.addLayer(corrected.median(),
  {bands:['Nadir_Reflectance_Band1','Nadir_Reflectance_Band4','Nadir_Reflectance_Band3'],
    min:0,max:2000
  },'corrected')


var lat = ee.Image.pixelLonLat().select('latitude').multiply(Math.PI/180);

// mask clouds
function maskMODIS(image) {
  var shadowBitMask = 4;
  var adjacentBitMask = 1 << 13;
  var cloudInternalBitMask = 1 << 10;
  var snowBitMask = 1 << 15;

  var qa = image.select('state_1km');
  
  //var valMask_red = image.select('sur_refl_b01').lt(1200)
  //var valMask_green = image.select('sur_refl_b04').lt(1200)
  //var valMask_blue = image.select('sur_refl_b03').lt(1200)

  var shadowMask = qa.bitwiseAnd(shadowBitMask).eq(0);
  var adjacentMask = qa.bitwiseAnd(adjacentBitMask).eq(0);
  var cloudInternalMask = qa.bitwiseAnd(cloudInternalBitMask).eq(0);
  var snowMask = qa.bitwiseAnd(snowBitMask).eq(0);
  var mask = shadowMask
    .and(adjacentMask)
    .and(cloudInternalMask)
    .and(snowMask)
    //.and(valMask_red)
    //.and(valMask_green)
    //.and(valMask_blue);

  return image.updateMask(mask)
}



var cloudsRemoved = collection.map(maskMODIS);
  //.map(function(image){return image.clip(roi3)})
  
var intMean = cloudsRemoved
  .select(['sur_refl_b01','sur_refl_b04','sur_refl_b03'])
  .reduce(ee.Reducer.intervalMean(20,80))
  
print(intMean)

var calcDistance = function(image){
  var d1 = image.select('sur_refl_b01').subtract(intMean.select('sur_refl_b01_mean'))
  var d4 = image.select('sur_refl_b04').subtract(intMean.select('sur_refl_b04_mean'))
  var d3 = image.select('sur_refl_b03').subtract(intMean.select('sur_refl_b03_mean'))
  var dist = ((d1.pow(2)).add(d4.pow(2)).add(d3.pow(2))).sqrt().rename('dist')
  return image.addBands(dist)
}

var withDist = cloudsRemoved.map(calcDistance)
  
var percentiles = withDist
  .select('dist')
  .reduce(ee.Reducer.percentile([90]))

Map.addLayer(percentiles,{},'percentiles')

var filtPercentiles = function(image){
  var distMask = image.select('dist')
    .lte(percentiles.select('dist_p90'));

  return image.updateMask(distMask)
}

var pFilt = withDist.map(filtPercentiles)

//Map.addLayer(pFilt.median(),visParams,'pFilt')

var numPixels = pFilt.count();
var newMask = numPixels.select('sur_refl_b01').gt(MIN_SAMPLES);
var maskFew = function(image){
  return image.updateMask(newMask);
}

var masked = pFilt.map(maskFew);
var masked_median = masked.median()


Map.addLayer(masked_median,visParams,'masked')
Map.addLayer(cloudsRemoved,visParams,'cloudsRemoved')

Export.image.toDrive({
  image: masked_median,
  description: 'masked_median_full_epsg4326',
  scale: 500,
  region: roi,
  crs:'EPSG:4326'
});

Export.image.toDrive({
  image: refl9,
  description: 'reflectance_region_jan_18_pfilt_10-90_full_epsg4326',
  scale: 500,
  region: roi,
  crs:'EPSG:4326'
});


//Map.addLayer(masked.last(),visParams,'masked_last')

var visParams2 = {bands:['red_nadir','green_nadir','blue_nadir'], min:0, max:0.2}
Map.addLayer(refl9, visParams2, 'ref9')
Map.addLayer(refl8, visParams2, 'ref8')

var addKernels = function(image){
  var theta_i = image.select('SolarZenith').multiply(0.01);
  var theta_v = image.select('SensorZenith').multiply(0.01);
  var phi_i = image.select('SolarAzimuth').multiply(0.01);
  var phi_v = image.select('SensorAzimuth').multiply(0.01);
  var phi_r = phi_i.subtract(phi_v); 
  
  // calculate kernels at acquired angles
  var kernels = brdf.calcKernel(theta_i,theta_v, phi_i, phi_v)
    .select(['Kvol','Ksparse']);
  
  // calculate zenith at solar noon
  var date = ee.Date(image.get('system:time_start'));
  var year = date.get('year');
  var yearStart = ee.Date.fromYMD(year,1,1);
  var days = date.difference(yearStart,'day').toDouble();
  
  var d = ((((ee.Image(days).subtract(173)).multiply(0.98563*Math.PI/180))
    .cos()).multiply(0.39795)).asin().rename('d'); 
  var sn = (lat.subtract(d)).multiply(180/(Math.PI)).rename('sn');
  
  // calculate kernels for nadir and solar zenith
  var kernels_nadir = brdf.calcKernel(sn,
    ee.Image(0), ee.Image(0), ee.Image(0))
    .select(['Kvol','Ksparse']).rename(['Kvol_nadir','Ksparse_nadir']);
  
  // constant for fiso
  var ones = ee.Image(1).rename('ones');
  
  return image.addBands(kernels).addBands(kernels_nadir).addBands(ones)
    .addBands([d.multiply(180/Math.PI),sn]);
};


var withKernels = masked.map(addKernels);
//Map.addLayer(withKernels.select(['Ksparse','Kvol']), {}, 'withKernels');

//var countKernels = withKernels.count();
//Map.addLayer(countKernels, {}, 'countKernels');

// Convert to an array. Give the axes names for more readable code.
var imageAxis = 0;
var bandAxis = 1;

// perform inversion
var red = withKernels.select(['sur_refl_b01','ones','Kvol','Ksparse']).toArray();
var x_red = red.arraySlice(bandAxis, 1, 4);
var y_red = red.arraySlice(bandAxis, 0, 1).divide(10000);
var fit_red = x_red.matrixSolve(y_red);

var green = withKernels.select(['sur_refl_b04','ones','Kvol','Ksparse']).toArray();
var x_green = green.arraySlice(bandAxis, 1, 4);
var y_green = green.arraySlice(bandAxis, 0, 1).divide(10000);
var fit_green = x_green.matrixSolve(y_green);

var blue = withKernels.select(['sur_refl_b03','ones','Kvol','Ksparse']).toArray();
var x_blue = blue.arraySlice(bandAxis, 1, 4);
var y_blue = blue.arraySlice(bandAxis, 0, 1).divide(10000);
var fit_blue = x_blue.matrixSolve(y_blue);
//print(fit_red);
//Map.addLayer(fit_red,{},'fit_red');
//Map.addLayer(x_red,{},'x_red');
//Map.addLayer(y,{},'y');

//var mm = x_red.matrixMultiply(fit_red);
//print('mm',mm);
//Map.addLayer(mm,{},'mm');

//var fiso = fit.arrayGet([0, 0]);
//var fvol = fit.arrayGet([1, 0]);
//var fgeo = fit.arrayGet([2, 0]);

  
//Map.addLayer(params.median(), {}, 'params');

// function
// calc solar noon zenith
// calc kernels
// addBands(ones, kvol, ksparse)
// map function to collection


var correctNadir = function(image){
  var kernels_orig = ee.ImageCollection(image.select(['ones','Kvol','Ksparse'])).toArray()
    
  var kernels_nadir = ee.ImageCollection(image.select(['ones','Kvol_nadir','Ksparse_nadir'])).toArray()
  
  var reflectance_fit = kernels_orig.matrixMultiply(fit_red).arrayGet([0,0]).rename(['red_orig'])
    .addBands(kernels_orig.matrixMultiply(fit_green).arrayGet([0,0]).rename(['green_orig']))
    .addBands(kernels_orig.matrixMultiply(fit_blue).arrayGet([0,0]).rename(['blue_orig']));
    
  var reflectance_nadir = kernels_nadir.matrixMultiply(fit_red).arrayGet([0,0]).rename(['red_nadir'])
    .addBands(kernels_nadir.matrixMultiply(fit_green).arrayGet([0,0]).rename(['green_nadir']))
    .addBands(kernels_nadir.matrixMultiply(fit_blue).arrayGet([0,0]).rename(['blue_nadir']));

  return image.addBands(reflectance_fit).addBands(reflectance_nadir);
}

//.arrayTranspose(1,0)//.arrayReshape(ee.Image(3).toArray(),2)
//var fit2 = fit.arrayReshape(ee.Image(3).toArray(),1)
//var as = nadir_matrix.arraySlice(1,0,2)

//Map.addLayer(nadir_matrix,{},'nadir_matrix')
//Map.addLayer(fit2,{},'fit2')
//Map.addLayer(as,{},'as')

var reflectance = withKernels.map(correctNadir)
//Map.addLayer(reflectance.select('red_nadir','green_nadir','blue_nadir'),{min:0,max:0.3},'reflectance')

var chart1 = ui.Chart.image.doySeries(
  reflectance.select('red_orig','red_nadir'),
  mypoint, ee.Reducer.first(),10,ee.Reducer.first())
  .setChartType('ScatterChart')
//print(chart1)

var scaleAngles = function(image){
  return image.addBands([
    image.select('SolarZenith').divide(100).rename('SolarZenith_scaled'),
    image.select('SensorZenith').divide(100).rename('SensorZenith_scaled'),
    image.select('SolarAzimuth').divide(100).rename('SolarAzimuth_scaled'),
    image.select('SensorAzimuth').divide(100).rename('SensorAzimuth_scaled'),
    image.select('sur_refl_b01').divide(10000).rename('sur_refl_b01_scaled'),
    image.select('sur_refl_b04').divide(10000).rename('sur_refl_b04_scaled'),
    image.select('sur_refl_b03').divide(10000).rename('sur_refl_b03_scaled')
  ]);
}
var scaled = reflectance.map(scaleAngles)

var chart1a = ui.Chart.image.doySeries(
  scaled.select('SolarZenith_scaled',
  'SensorZenith_scaled','SolarAzimuth_scaled','SensorAzimuth_scaled'),
  mypoint, ee.Reducer.first(),10,ee.Reducer.first(),0,365)
  .setChartType('ScatterChart').setOptions({pointSize: 3});
//print(chart1a)

//var pix = reflectance.select('red_orig').getRegion(mypoint, 10);
//print('pix',pix.slice(1,-1).toArray())

var xp = reflectance.select('red_orig').toArray()
  .reduceRegion(ee.Reducer.first(),mypoint,10)
var yp = reflectance.select('sur_refl_b01').toArray()//.divide(10000)
  .reduceRegion(ee.Reducer.first(),mypoint,10)

//var chart_1c = ui.Chart.feature.byFeature(pixel, 'green_orig','green_nadir')
//print(chart_1c)

var chart1b = ui.Chart.array.values(yp.toArray(),0,xp.toArray()).setOptions({
      title: 'Fitted Values',
      hAxis: {'title': 'red_orig'},
      vAxis: {'title': 'sur_refl_b01'},
      pointSize: 5,
});
//print(chart1b);

var chart1c = ui.Chart.image.doySeries(
  reflectance.select('sn','d'),
  mypoint, ee.Reducer.first(),10,ee.Reducer.first())
  .setChartType('LineChart')
//print(chart1c)

var chart1d = ui.Chart.image.doySeries(
  scaled.select('sur_refl_b01_scaled','red_nadir'),
  mypoint, ee.Reducer.first(),10,ee.Reducer.first())
  .setChartType('ScatterChart').setOptions({pointSize:3})
print(chart1d)


var chart2 = ui.Chart.image.doySeries(
  scaled.select('sur_refl_b04_scaled','green_nadir'),
  mypoint, ee.Reducer.first(),10,ee.Reducer.first())
  .setChartType('ScatterChart').setOptions({pointSize:3})
print(chart2)

var chart3 = ui.Chart.image.doySeries(
  scaled.select('sur_refl_b03_scaled','blue_nadir'),
  mypoint, ee.Reducer.first(),10,ee.Reducer.first())
  .setChartType('ScatterChart').setOptions({pointSize:3})
print(chart3)

chart1.onClick(function(xValue, yValue, seriesName) {
  if (!xValue) return;  // Selection was cleared.

  var equalDate = ee.Date.fromYMD(2018,1,1).advance(xValue-1,'day')
  var image1 = ee.Image(collection.filterDate(equalDate).first());
  var image2 = ee.Image(masked.filterDate(equalDate).first());
  var image3 = ee.Image(reflectance.filterDate(equalDate).first());
  var image4 = ee.Image(collectionQC.filterDate(equalDate).first());
  var layer1 = ui.Map.Layer(image1, visParams,'collection');
  var layer2 = ui.Map.Layer(image2, visParams,'masked');
  var layer3 = ui.Map.Layer(image3,
    {bands:['red_orig','green_orig','blue_orig'],min:0,max:0.2},'reflectance');
  var layer4 = ui.Map.Layer(image4,imageVisParam2,'collectionQC');
  Map.layers().reset([layer1,layer2,layer3,layer4]);
});

Export.image.toAsset(reflectance.median(),
  "reflectance_region_jan_18_pfilt_10-90_reproj","reflectance_region_jan_18_pfilt_10-90_reproj",null,null,roi,500)

//Map.addLayer(fiso, {}, 'fiso');
//Map.addLayer(fvol, {}, 'fvol');
//Map.addLayer(fgeo, {}, 'fgeo');
