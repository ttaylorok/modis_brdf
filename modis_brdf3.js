var brdf = require('users/ttaylorok/amazon:modis_kernels');

var collection = ee.ImageCollection('MODIS/006/MYD09GA')
  .filterDate('2018-08-01', '2018-12-01');
                
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
//Map.addLayer(cloudsRemoved.first(),visParams,'modis_mask');
 
var modis = collection.first();

var addKernels = function(image){
  var theta_i = image.select('SolarZenith')//.multiply(0.01*Math.PI/180); // solar zenith angle
  var theta_v = image.select('SensorZenith')//.multiply(0.01*Math.PI/180); // view zenith angle
  var phi_i = image.select('SolarAzimuth')//.multiply(0.01*Math.PI/180);
  var phi_v = image.select('SensorAzimuth')//.multiply(0.01*Math.PI/180);
  var phi_r = phi_i.subtract(phi_v); 
  
  var kernels = brdf.calcKernel(theta_i.multiply(0.01),
    theta_v.multiply(0.01), phi_i.multiply(0.01), phi_v.multiply(0.01));
  
  // NOTE: IF ANY BANDS ARE MASKED, TOARRAY WILL MASK THE WHOLE PIXEL
  var out1 = image.select(['sur_refl_b01'])
    .addBands([ee.Image(1),kernels.select(['Kvol', 'Ksparse'])])
  
  return out1.rename(['sur_refl_b01','ones','Kvol','Ksparse'])
}

var visParams = {bands : ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'],
                min : 0, max : 6000};
                
var cloudsRemoved = collection
  .map(function(image){return image.clip(roi3)})
  .map(maskMODIS);
//var collection2 = cloudsRemoved.map(function(image){
//  return image.clip(roi3)
//});
Map.addLayer(cloudsRemoved,{},'cloudsRemoved')

var numPixels = cloudsRemoved.count();
var newMask = numPixels.select('sur_refl_b01').gt(10);
Map.addLayer(numPixels,{},'numPixels');

var maskFew = function(image){
  return image.updateMask(newMask);
  //return image.unmask(ee.Image(999));
}
var masked = cloudsRemoved.map(maskFew);
//Map.addLayer(masked, {}, 'masked');

var withKernels = masked.map(addKernels); //CHANGE BACK TO MASKED or MASKEDMODIS
Map.addLayer(withKernels.select(['Ksparse','Kvol']), {}, 'withKernels');

var countKernels = withKernels.count();
Map.addLayer(countKernels, {}, 'countKernels');

// Convert to an array. Give the axes names for more readable code.
var array = withKernels.toArray();
var imageAxis = 0;
var bandAxis = 1;
 
// Slice off the year and ndvi, and solve for the coefficients.
var x = array.arraySlice(bandAxis, 1, 4);
var y = array.arraySlice(bandAxis, 0, 1);
var fit = x.matrixSolve(y);
print(fit);
Map.addLayer(fit,{},'fit');
Map.addLayer(x,{},'x');
Map.addLayer(y,{},'y');

var mm = x.matrixMultiply(fit);
print('mm',mm);
//Map.addLayer(mm,{},'mm');

var fiso = fit.arrayGet([0, 0]);
var fvol = fit.arrayGet([1, 0]);
var fgeo = fit.arrayGet([2, 0]);

//Map.addLayer(fiso, {}, 'fiso');
//Map.addLayer(fvol, {}, 'fvol');
//Map.addLayer(fgeo, {}, 'fgeo');
