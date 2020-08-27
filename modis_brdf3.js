var brdf = require('users/ttaylorok/amazon:modis_kernels');

var collection = ee.ImageCollection('MODIS/006/MYD09GA')
  .merge(ee.ImageCollection('MODIS/006/MOD09GA'))
  .filterDate('2018-08-01', '2018-08-31');

//var collection = ee.ImageCollection('MODIS/006/MYD09GA')
//  .filterDate('2018-01-01', '2018-12-31');
                
// mask clouds
function maskMODIS(image) {
  // Bits 0-1 and 2 are cloud could and shadow, respectively.
  var cloudBitMask = 3;
  var shadowBitMask = 4;

  // Get the pixel QA band.
  var qa = image.select('state_1km');

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(shadowBitMask).eq(0));

  // Return the masked image, scaled to reflectance, without the QA bands.
  return image.updateMask(mask).clip(roi3)
      .copyProperties(image, ["system:time_start"]);
}


//Map.addLayer(collection.first(),visParams,'modis');


var addKernels = function(image){
  var theta_i = image.select('SolarZenith')//.multiply(0.01*Math.PI/180); // solar zenith angle
  var theta_v = image.select('SensorZenith')//.multiply(0.01*Math.PI/180); // view zenith angle
  var phi_i = image.select('SolarAzimuth')//.multiply(0.01*Math.PI/180);
  var phi_v = image.select('SensorAzimuth')//.multiply(0.01*Math.PI/180);
  var phi_r = phi_i.subtract(phi_v); 
  
  var kernels = brdf.calcKernel(theta_i.multiply(0.01),
    theta_v.multiply(0.01), phi_i.multiply(0.01), phi_v.multiply(0.01));
  
  // NOTE: IF ANY BANDS ARE MASKED, TOARRAY WILL MASK THE WHOLE PIXEL
  var out1 = image.addBands([ee.Image(1).rename('ones'),
    kernels.select(['Kvol', 'Ksparse'])])
  
  return out1
}

var visParams = {bands : ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'],
                min : 0, max : 3000};
                
var cloudsRemoved = collection
  .map(function(image){return image.clip(roi3)})
  .map(maskMODIS);
//var collection2 = cloudsRemoved.map(function(image){
//  return image.clip(roi3)
//});
Map.addLayer(cloudsRemoved.median(),visParams,'cloudsRemoved')

var numPixels = cloudsRemoved.count();
var newMask = numPixels.select('sur_refl_b01').gt(20);
//Map.addLayer(numPixels,{},'numPixels');

var maskFew = function(image){
  return image.updateMask(newMask);
  //return image.unmask(ee.Image(999));
}
var masked = cloudsRemoved.map(maskFew);
//Map.addLayer(masked, {}, 'masked');

var withKernels = masked.map(addKernels); //CHANGE BACK TO MASKED or MASKEDMODIS
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

var mm = x_red.matrixMultiply(fit_red);
print('mm',mm);
//Map.addLayer(mm,{},'mm');

//var fiso = fit.arrayGet([0, 0]);
//var fvol = fit.arrayGet([1, 0]);
//var fgeo = fit.arrayGet([2, 0]);

var params = ee.ImageCollection('MODIS/006/MCD43A1')
  .filterDate('2018-08-01', '2018-12-01')
  .filterBounds(roi3);
  
//Map.addLayer(params.median(), {}, 'params');

// function
// calc solar noon zenith
// calc kernels
// addBands(ones, kvol, ksparse)
// map function to collection

var lat = ee.Image.pixelLonLat().select('latitude').multiply(Math.PI/180);

var calcNadir = function(image){
  var date = ee.Date(image.get('system:time_start'));
  var year = date.get('year');
  var yearStart = ee.Date.fromYMD(year,1,1);
  var days = date.difference(yearStart,'day').toDouble();
  
  var d = ((((ee.Image(days).subtract(173)).multiply(0.98563*Math.PI/180))
    .cos()).multiply(0.39795)).asin().rename('d'); // solar declination
  var sn = (lat.subtract(d)).multiply(180/(Math.PI)).rename('sn'); // zenith at solar noon
  
  return image.addBands([d.multiply(180/Math.PI).rename('d'),sn,ee.Image(days).toDouble().rename('days')])
    .addBands([image.select('SolarZenith').multiply(0.01).rename('szc'),image.select('SensorZenith').multiply(0.01).rename('vzc')])
}

var withNadir = masked.map(calcNadir);
Map.addLayer(withNadir.select(['szc','vzc','d','sn']),{},'withNadir')

var calcNadir = function(image){
  var kernels_nadir = brdf.calcKernel(image.select('sn'),
    ee.Image(0), ee.Image(0), ee.Image(0));
    
  var nadir_matrix = ee.ImageCollection(ee.Image(1)
    .addBands(kernels_nadir.select(['Kvol', 'Ksparse'])))
    .toArray()
    
  var reflectance = nadir_matrix.matrixMultiply(fit_red).arrayGet([0,0])
    .addBands(nadir_matrix.matrixMultiply(fit_green).arrayGet([0,0]))
    .addBands(nadir_matrix.matrixMultiply(fit_blue).arrayGet([0,0]))

  return reflectance;
}

//.arrayTranspose(1,0)//.arrayReshape(ee.Image(3).toArray(),2)
//var fit2 = fit.arrayReshape(ee.Image(3).toArray(),1)
//var as = nadir_matrix.arraySlice(1,0,2)

//Map.addLayer(nadir_matrix,{},'nadir_matrix')
//Map.addLayer(fit2,{},'fit2')
//Map.addLayer(as,{},'as')

var reflectance = withNadir.map(calcNadir)
Map.addLayer(reflectance.mean(),{min:0,max:0.3},'reflectance')

//Map.addLayer(fiso, {}, 'fiso');
//Map.addLayer(fvol, {}, 'fvol');
//Map.addLayer(fgeo, {}, 'fgeo');
