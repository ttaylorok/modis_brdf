var brdf = require('users/ttaylorok/amazon:modis_kernels');

var collection = ee.ImageCollection('MODIS/006/MYD09GA')
  .merge(ee.ImageCollection('MODIS/006/MOD09GA'))
  .filterDate('2018-09-01','2018-10-01')
  .sort('system:time_start');

var visParams = {bands : ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'],
                min : 0, max : 3000};
             
Map.addLayer(collection,visParams,'collection');

var lat = ee.Image.pixelLonLat().select('latitude').multiply(Math.PI/180);

// mask clouds
function maskMODIS(image) {
  var cloudBitMask = 3;
  var shadowBitMask = 4;
  var adjacentBitMask = 1 << 13;
  var cloudInternalBitMask = 1 << 10;

  var qa = image.select('state_1km');

  var cloudMask = qa.bitwiseAnd(cloudBitMask).eq(0);
  var shadowMask = qa.bitwiseAnd(shadowBitMask).eq(0);
  var adjacentMask = qa.bitwiseAnd(adjacentBitMask).eq(0);
  var cloudInternalMask = qa.bitwiseAnd(cloudInternalBitMask).eq(0);
  var mask = cloudMask.and(shadowMask).and(adjacentMask).and(cloudInternalMask);

  return image.updateMask(mask).clip(roi)
      //.copyProperties(image, ["system:time_start"]);
}


var cloudsRemoved = collection.map(maskMODIS);
  //.map(function(image){return image.clip(roi3)})
  
  
Map.addLayer(cloudsRemoved.median(),visParams,'cloudsRemoved')

var numPixels = cloudsRemoved.count();
var newMask = numPixels.select('sur_refl_b01').gt(5);
var maskFew = function(image){
  return image.updateMask(newMask);
}

var masked = cloudsRemoved.map(maskFew);
Map.addLayer(masked,visParams,'masked')

var addKernels = function(image){
  // calculate kernals
  var theta_i = image.select('SolarZenith').multiply(0.01);
  var theta_v = image.select('SensorZenith').multiply(0.01);
  var phi_i = image.select('SolarAzimuth').multiply(0.01);
  var phi_v = image.select('SensorAzimuth').multiply(0.01);
  var phi_r = phi_i.subtract(phi_v); 
  
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
    
  var ones = ee.Image(1).rename('ones');
  
  return image.addBands(kernels).addBands(kernels_nadir).addBands(ones);
};


var withKernels = masked.map(addKernels);
//Map.addLayer(withKernels.select(['Ksparse','Kvol']), {}, 'withKernels');

//var countKernels = withKernels.count();
//Map.addLayer(countKernels, {}, 'countKernels');

// Convert to an array. Give the axes names for more readable code.
var imageAxis = 0;
var bandAxis = 1;

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
Map.addLayer(reflectance.select('red_nadir','green_nadir','blue_nadir'),{min:0,max:0.3},'reflectance')

var chart1 = ui.Chart.image.doySeries(
  reflectance.select('red_orig','red_nadir'),
  mypoint, ee.Reducer.mean(),10,ee.Reducer.mean())
  .setChartType('ScatterChart')
print(chart1)

var chart2 = ui.Chart.image.doySeries(
  reflectance.select('green_orig','green_nadir'),
  mypoint, ee.Reducer.mean(),10,ee.Reducer.mean())
  .setChartType('ScatterChart')
print(chart2)

var chart3 = ui.Chart.image.doySeries(
  reflectance.select('blue_orig','blue_nadir'),
  mypoint, ee.Reducer.mean(),10,ee.Reducer.mean())
  .setChartType('ScatterChart')
print(chart3)



//Map.addLayer(fiso, {}, 'fiso');
//Map.addLayer(fvol, {}, 'fvol');
//Map.addLayer(fgeo, {}, 'fgeo');
