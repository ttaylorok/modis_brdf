
var calcKernel = function(theta_i1, theta_v1, phi_i1, phi_v1){
  var theta_i = ee.Image(theta_i1).multiply(Math.PI/180); // solar zenith angle
  var theta_v = ee.Image(theta_v1).multiply(Math.PI/180); // view zenith angle
  var phi_i = ee.Image(phi_i1).multiply(Math.PI/180);
  var phi_v = ee.Image(phi_v1).multiply(Math.PI/180);
  var phi_r = phi_v.subtract(phi_i); // WHICH WAY IS POSITIVE???
  
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
  
  //var theta_ip = (theta_i.tan().multiply(b_r)).atan();
  //var theta_vp = (theta_v.tan().multiply(b_r)).atan();

  var theta_ip = theta_i.toDouble();
  var theta_vp = theta_v.toDouble();
  
  var D = (theta_ip.tan().multiply(theta_ip.tan())
    .add(theta_vp.tan().multiply(theta_vp.tan()))
    .subtract(theta_ip.tan()
    .multiply(theta_vp.tan())
    .multiply(phi_r.cos())
    .multiply(2))).sqrt();
    
  var cos_t = (D.multiply(D)
    .add(theta_ip.tan()
    .multiply(theta_vp.tan())
    .multiply(phi_r.sin()).pow(2)))
    .sqrt()
    .multiply(h_b)
    .divide(theta_ip.cos().pow(-1).add(theta_vp.cos().pow(-1)))
    
  var t = cos_t.acos()
  
  var O = (t.subtract(t.sin().multiply(cos_t)))
    .multiply(theta_ip.cos().pow(-1).add(theta_vp.cos().pow(-1)))
    .divide(Math.PI)
  
  var cos_xip = theta_ip.cos()
    .multiply(theta_vp.cos())
    .add(theta_ip.sin()
    .multiply(theta_vp.sin())
    .multiply(phi_r.cos()));
  
  var Ksparse = O.subtract(theta_ip.cos().pow(-1))
    .subtract(theta_vp.cos().pow(-1))
    .add((cos_xip.add(1)).multiply(theta_v.cos().pow(-1)).multiply(0.5)); 
    
  var out = Kvol.addBands(Ksparse)
    .addBands(theta_ip)
    .addBands(theta_vp)
    .addBands(D)
    .addBands(cos_t)
    .addBands(t)
    .addBands(O)
    .addBands(cos_xip)
  return out.rename(['Kvol','Ksparse','theta_ip','theta_vp',
  'D','cost_t','t','O','cos_xip'])
}

var createPlot = function(sza,phi_i1,phi_v1){
  var c1 = calcKernel(sza,-90,phi_i1,phi_v1).set({'angle':-90});
  var c2 = calcKernel(sza,-60,phi_i1,phi_v1).set({'angle':-60});
  var c3 = calcKernel(sza,-30,phi_i1,phi_v1).set({'angle':-30});
  var c4 = calcKernel(sza,0,phi_i1,phi_v1).set({'angle':0});
  var c5 = calcKernel(sza,30,phi_i1,phi_v1).set({'angle':30});
  var c6 = calcKernel(sza,60,phi_i1,phi_v1).set({'angle':60});
  var c7 = calcKernel(sza,90,phi_i1,phi_v1).set({'angle':90});
  //Map.addLayer(test, {}, 'test')
  
  var ic = ee.ImageCollection([c1,c2,c3,c4,c5,c6,c7]);
  var chart2 = ui.Chart.image.series(ic,mypoint,ee.Reducer.first(),10,'angle')
    .setOptions({title: 'Hehehe wheoeoeoe',
      hAxis: {'title': 'orig'},
      vAxis: {'title': 'orig', viewWindow : {'min' : -3, 'max' : 3}},//
      pointSize: 5,
  });
  print(chart2);
  return ic
}

//principal plane, pi = 0 or 180
//cross principal plane, pi = 90 or 270
var ic = createPlot(33,0,180);

Map.addLayer(ee.Image(Math.PI).cos(),{},'power')
Map.addLayer(ic,{},'ic')


//var array = ic..toArray();
//var bc = c1.addBands(c2).addBands(c3).addBands(c4).addBands(c5)
//  .addBands(c6).addBands(c7)


//var vals = bc.reduceRegion(ee.Reducer.first(),mypoint,10).toArray();
//var ksparse = ic.select('Ksparse').reduceRegion(ee.Reducer.first(),mypoint,10).toArray();

//var chart = ui.Chart.array.values(vals,0).setOptions({
//      title: 'Hehehe wheoeoeoe',
//      hAxis: {'title': 'orig'},
//      vAxis: {'title': 'pred'},
//      pointSize: 5,
//});
//print(chart);

//var vals = bc.reduceRegion(ee.Reducer.first(),mypoint,10).toArray();


// NOTE: IF ANY BANDS ARE MASKED, TOARRAY WILL MASK THE WHOLE PIXEL
//var out1 = image.select(['sur_refl_b01']).addBands([ee.Image(1),Kvol,Ksparse])

//return out1.rename(['sur_refl_b01','ones','Kvol','Ksparse'])