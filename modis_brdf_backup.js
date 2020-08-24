var collection = ee.ImageCollection('MODIS/006/MYD09GA')
  .filterBounds(roi)
  .filterDate('2018-01-01', '2018-4-01');

var visParams = {bands : ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'],
                min : 0, max : 6000};
                

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
  return image.updateMask(mask)
      .copyProperties(image, ["system:time_start"]);
}

var cloudsRemoved = collection.map(maskMODIS);

Map.addLayer(collection.first(),visParams,'modis');
Map.addLayer(cloudsRemoved.first(),visParams,'modis_mask');
 
var modis = collection.first();

var addKernels = function(image){
  var theta_i = image.select('SolarZenith').multiply(0.01*Math.PI/180); // solar zenith angle
  var theta_v = image.select('SensorZenith').multiply(0.01*Math.PI/180); // view zenith angle
  var phi_i = image.select('SolarAzimuth').multiply(0.01*Math.PI/180);
  var phi_v = image.select('SensorAzimuth').multiply(0.01*Math.PI/180);
  var phi_r = phi_i.subtract(phi_v); // WHICH WAY IS POSITIVE???
  
  // phase angle scattering
  var cos_xi = theta_i.cos()
    .multiply(theta_v.cos())
    .add(theta_i.sin()
    .multiply(theta_v.sin())
    .multiply(phi_r.cos()));
    
  var xi = cos_xi.acos();
  
  var Kvol = (((xi.multiply(-1).add(Math.PI/2))
    .multiply(cos_xi)
    .add(xi.sin())).divide(theta_i.cos().add(theta_v.cos())))
    .subtract(Math.PI/4);
    
  // LiSparse kernel
  var h_b = 2;
  var b_r = 1;
  
  var theta_ip = (theta_i.tan().multiply(b_r)).atan();
  var theta_vp = (theta_v.tan().multiply(b_r)).atan();
  
  var D = (theta_ip.tan().multiply(theta_ip.tan())
    .add(theta_vp.tan().multiply(theta_vp.tan()))
    .subtract(theta_ip.tan()
    .multiply(theta_vp.tan())
    .multiply(phi_r.cos())
    .multiply(2))).sqrt();
    
  var cos_t = (D.multiply(D).add(theta_ip.tan()
    .multiply(theta_vp.tan())
    .multiply(phi_r.sin()))
    .pow(2)).sqrt()
    .multiply(h_b)
    .divide(theta_ip.cos().pow(-1).add(theta_vp.cos().pow(-1)))
    
  var t = cos_t.acos()
  
  var O = (t.subtract(t.sin().multiply(cos_t)))
    .multiply(theta_ip.cos().pow(-1).add(theta_vp.cos().pow(-1)))
    .multiply(1/Math.PI)
  
  var cos_xip = theta_ip.cos()
    .multiply(theta_vp.cos())
    .add(theta_ip.sin()
    .multiply(theta_vp.sin())
    .multiply(phi_r.cos()));
  
  var Ksparse = O.subtract(theta_ip.cos().pow(-1))
    .subtract(theta_vp.cos().pow(-1))
    .add((cos_xip.add(1)).multiply(theta_v.cos().pow(-1)).multiply(0.5));
  
  var out1 = image.select(['sur_refl_b01','sur_refl_b02','sur_refl_b03',
  'sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07'])
  .addBands([Kvol,Ksparse])
  
  return out1.rename(['sur_refl_b01','sur_refl_b02','sur_refl_b03',
  'sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07',
  'Kvol','Ksparse'])
}

var withKernels = cloudsRemoved.map(addKernels);
Map.addLayer(withKernels, {}, 'withKernels');

