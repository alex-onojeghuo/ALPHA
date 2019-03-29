//-----------------------------------------------------------------------------------------
//Random Forest Supervised Classification - GEE
//Compiled by ALEX O. ONOJEGHUO (Geospatial Data Scientist)
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------
//List of variables
//-----------------------------------------------------------------------------------------
var PCA1 = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/S2_PCA1"),
    PCA2 = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/S2_PCA2"),
    PCA3 = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/S2_PCA3"),
    PCA4 = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/S2_PCA4"),
    NDWIdiff = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/NDWI_diff"),
    SWI = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/SWI_LIDAR"),
    SA = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/grassland_input/SA"),
    Clip1 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/grassland_input/Clip1"),
    Clip2 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/grassland_input/Clip2"),
    Clip3 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/grassland_input/Clip3"),
    VVdiff = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/VVdiff_summer_fall"),
    Grid1 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid1"),
    Grid2 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid2"),
    Grid3 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid3"),
    Grid4 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid4"),
    Grid6 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid6"),
    Grid7 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid7"),
    Grid8 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid8"),
    Grid9 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid9"),
    Grid10 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid10"),
    Grid11 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid11"),
    Grid12 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid12"),
    Grid13 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid13"),
    Grid14 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid14"),
    Grid15 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid15"),
    Grid16 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid16"),
    Grid17 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid17"),
    Grid18 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid18"),
    Grid19 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid19"),
    Grid20 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid20"),
    Grid22 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid22"),
    Grid23 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid23"),
    Grid24 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid24"),
    Grid26 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid26"),
    Grid27 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid27"),
    Grid28 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid28"),
    Grid29 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid29"),
    Grid30 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid30"),
    Grid31 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid31"),
    Grid32 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid32"),
    Grid25 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid25"),
    Grid21 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid21_new"),
    Grid5 = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/PGNR_INPUT/Grid5_new"),
    PrjFile = ee.Image("users/Alberta_Sentinel_analysis_v20/Parkland_grassland"),
    AB = ee.FeatureCollection("users/Alberta_Sentinel_analysis_v20/grassland_input/SA"),
    Topo_PC1 = ee.Image("users/Alberta_Sentinel_analysis_v20/grassland_input/Topo_PC1");
	
//-----------------------------------------------------------------------------------------
//Calculate Sentinel-2 optical band differences
//-----------------------------------------------------------------------------------------

//Cloud masking functions
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//parameters
var cloudThresh = 15;
var cloudHeights = ee.List.sequence(200,10000,250);//Height of clouds to use to project cloud shadows
var irSumThresh =0.35;//Sum of IR bands to include as shadows within TDOM and the shadow shift method (lower number masks out less)
var dilatePixels = 2; //Pixels to dilate around clouds
var contractPixels = 1;//Pixels to reduce cloud mask and dark shadows by to reduce inclusion of single-pixel comission errors
//Functions
// var vizParams = {bands: ['red', 'green', 'blue'], min: 0, max: 0.3};
//////////////////////////////////////////////////////////////////////////
var rescale = function(img, exp, thresholds) {
    return img.expression(exp, {img: img})
        .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
  };
  
  var getNotWaterClusterID = function(clusterizedImage){
  var ID = clusterizedImage.reduceRegion({
    reducer:ee.Reducer.mean(),
    geometry:NotWaterPoint,
    scale:30
  });
  ID = ID.get('cluster');
  return ee.Number.parse(ID);
}
var getWaterClusterID = function(clusterizedImage){
  var ID = clusterizedImage.reduceRegion({
    reducer:ee.Reducer.mean(),
    geometry:WaterPoint,
    scale:30
  });
  ID = ID.get('cluster');
  return ee.Number.parse(ID);
}
////////////////////////////////////////
////////////////////////////////////////
// Cloud masking algorithm for Sentinel2
//Built on ideas from Landsat cloudScore algorithm
//Currently in beta and may need tweaking for individual study areas
function sentinelCloudScore(img) {
  

  // Compute several indicators of cloudyness and take the minimum of them.
  var score = ee.Image(1);
  
  // Clouds are reasonably bright in the blue and cirrus bands.
  score = score.min(rescale(img, 'img.blue', [0.1, 0.5]));
  score = score.min(rescale(img, 'img.cb', [0.1, 0.3]));
  score = score.min(rescale(img, 'img.cb + img.cirrus', [0.15, 0.2]));
  
  // Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]));

  
  //Clouds are moist
  var ndmi = img.normalizedDifference(['nir','swir1']);
  score=score.min(rescale(ndmi, 'img', [-0.1, 0.1]));
  // However, clouds are not snow.
  var ndsi = img.normalizedDifference(['green', 'swir1']);
  score=score.min(rescale(ndsi, 'img', [0.8, 0.6]));
  
  score = score.multiply(100).byte();
 
  return img.addBands(score.rename('cloudScore'));
}
//////////////////////////////////////////////////////////////////////////
// Function to mask clouds using the Sentinel-2 QA band.
function maskS2clouds(image) {
  var qa = image.select('QA60').int16();
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);
  
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));

  // Return the masked and scaled data.
  return image.updateMask(mask);
} 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function for finding dark outliers in time series
//Masks pixels that are dark, and dark outliers
function simpleTDOM2(c){
  var shadowSumBands = ['nir','swir1'];
  var irSumThresh = 0.4;
  var zShadowThresh = -1.2;
  //Get some pixel-wise stats for the time series
  var irStdDev = c.select(shadowSumBands).reduce(ee.Reducer.stdDev());
  var irMean = c.select(shadowSumBands).mean();
  var bandNames = ee.Image(c.first()).bandNames();
  
  //Mask out dark dark outliers
  c = c.map(function(img){
    var z = img.select(shadowSumBands).subtract(irMean).divide(irStdDev);
    var irSum = img.select(shadowSumBands).reduce(ee.Reducer.sum());
    var m = z.lt(zShadowThresh).reduce(ee.Reducer.sum()).eq(2).and(irSum.lt(irSumThresh)).not();
    
    return img.updateMask(img.mask().and(m));
  });
  
  return c.select(bandNames);
}
////////////////////////////////////////////////////////
/////////////////////////////////////////////
/***
 * Implementation of Basic cloud shadow shift
 * 
 * Author: Gennadii Donchyts
 * License: Apache 2.0
 */
function projectShadows(cloudMask,image,cloudHeights){
  var meanAzimuth = image.get('MEAN_SOLAR_AZIMUTH_ANGLE');
  var meanZenih = image.get('MEAN_SOLAR_ZENITH_ANGLE');
  ///////////////////////////////////////////////////////
  // print('a',meanAzimuth);
  // print('z',meanZenith)
  
  //Find dark pixels
  var darkPixels = image.select(['nir','swir1','swir2']).reduce(ee.Reducer.sum()).lt(irSumThresh)
    .focal_min(contractPixels).focal_max(dilatePixels)
  ;//.gte(1);
  
  
  //Get scale of image
  var nominalScale = cloudMask.projection().nominalScale();
  //Find where cloud shadows should be based on solar geometry
  //Convert to radians
  var azR =ee.Number(meanAzimuth).add(180).multiply(Math.PI).divide(180.0);
  var zenR  =ee.Number(meanZenith).multiply(Math.PI).divide(180.0);
  
  
 
  //Find the shadows
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
  
    var shadowCastedDistance = zenR.tan().multiply(cloudHeight);//Distance shadow is cast
    var x = azR.sin().multiply(shadowCastedDistance).divide(nominalScale);//X distance of shadow
    var y = azR.cos().multiply(shadowCastedDistance).divide(nominalScale);//Y distance of shadow
    // print(x,y)
   
    return cloudMask.changeProj(cloudMask.projection(), cloudMask.projection().translate(x, y));
    
    
  });
  
  
  var shadowMask = ee.ImageCollection.fromImages(shadows).max();
  // Map.addLayer(cloudMask.updateMask(cloudMask),{'min':1,'max':1,'palette':'88F'},'Cloud mask');
  // Map.addLayer(shadowMask.updateMask(shadowMask),{'min':1,'max':1,'palette':'880'},'Shadow mask');
  
  //Create shadow mask
  shadowMask = shadowMask.and(cloudMask.not());
  shadowMask = shadowMask.and(darkPixels).focal_min(contractPixels).focal_max(dilatePixels);
  
  var cloudShadowMask = shadowMask.or(cloudMask);
  
  image = image.updateMask(cloudShadowMask.not()).addBands(shadowMask.rename(['cloudShadowMask']));
  return image;
}
//////////////////////////////////////////////////////
//Function to bust clouds from S2 image
function bustClouds(img){
  img = sentinelCloudScore(img);
  img = img.updateMask(img.select(['cloudScore']).gt(cloudThresh).focal_min(contractPixels).focal_max(dilatePixels).not());
  return img;
}
//////////////////////////////////////////////////////
//Function for wrapping the entire process to be applied across collection
function wrapIt(img){
  img = sentinelCloudScore(img);
  var cloudMask = img.select(['cloudScore']).gt(cloudThresh)
    .focal_min(contractPixels).focal_max(dilatePixels)

  img = projectShadows(cloudMask,img,cloudHeights);

  return img.clip(geometry3);
}
////////////////////////////////////////////////////// 
//Function to find unique values of a field in a collection
function uniqueValues(collection,field){
    var values  = ee.Dictionary(collection.reduceColumns(ee.Reducer.frequencyHistogram(),[field]).get('histogram')).keys();
    
    return values;
  }
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------
//Define Alberta shape and projection file
//-----------------------------------------------------------------------------------------
var prj = PrjFile.projection();
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//Get static Sentinel-2 image collections
//-----------------------------------------------------------------------------------------
var S2_1 = ee.ImageCollection('COPERNICUS/S2')
  .filterDate('2016-04-01', '2016-08-31')
  .filterBounds(AB)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 50)
  .map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
    return out;
    })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);

var S2_2 = ee.ImageCollection('COPERNICUS/S2') 
  .filterDate('2016-09-01', '2017-03-31')
  .filterBounds(AB)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 50)
  .map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
    return out;
    })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);

var S2_3 = ee.ImageCollection('COPERNICUS/S2') 
  .filterDate('2017-04-01', '2017-08-31')
  .filterBounds(AB)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 50)
  .map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
    return out;
    })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);

var S2_4 = ee.ImageCollection('COPERNICUS/S2') 
  .filterDate('2017-09-01', '2018-03-31')
  .filterBounds(AB)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 50)
  .map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
    return out;
    })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);

var S2 = ee.ImageCollection(S2_1.merge(S2_2));
var S2 = ee.ImageCollection(S2.merge(S2_3));
var S2 = ee.ImageCollection(S2.merge(S2_4));

//print(S2_1);
//print(S2_2);
//print(S2_3);
//print(S2_4);
//print(S2);

//-----------------------------------------------------------------------------------------
//Apply cloud masking algorithms
//-----------------------------------------------------------------------------------------
//Bust clouds using BQA method
var S2QA_1 = S2_1.map(maskS2clouds);
var S2QA_2 = S2_2.map(maskS2clouds);
var S2QA_3 = S2_3.map(maskS2clouds);
var S2QA_4 = S2_4.map(maskS2clouds);

//Bust clouds using cloudScore method
var S2_1 = S2_1.map(bustClouds);
var S2_2 = S2_2.map(bustClouds);
var S2_3 = S2_3.map(bustClouds);
var S2_4 = S2_4.map(bustClouds);

//Bust clouds using cloudScore and shadows using TDOM
var S2_1 = simpleTDOM2(S2_1);
var S2_2 = simpleTDOM2(S2_2);
var S2_3 = simpleTDOM2(S2_3);
var S2_4 = simpleTDOM2(S2_4);

//get median of all bands
var S2_1 = S2_1.median();
var S2_2 = S2_2.median();
var S2_3 = S2_3.median();
var S2_4 = S2_4.median();

//-----------------------------------------------------------------------------------------
//get bands as varibles
//----------------------------------------------------------------------------------------- 
var B2_1 = S2_1.select(['blue']);
var B3_1 = S2_1.select(['green']);
var B4_1 = S2_1.select(['red']);
var B5_1 = S2_1.select(['re1']);
var B6_1 = S2_1.select(['re2']);
var B7_1 = S2_1.select(['re3']);
var B8_1 = S2_1.select(['nir']);
var B8A_1 = S2_1.select(['nir2']);
var B11_1 = S2_1.select(['swir1']);
var B12_1 = S2_1.select(['swir2']);

var B2_2 = S2_2.select(['blue']);
var B3_2 = S2_2.select(['green']);
var B4_2 = S2_2.select(['red']);
var B5_2 = S2_2.select(['re1']);
var B6_2 = S2_2.select(['re2']);
var B7_2 = S2_2.select(['re3']);
var B8_2 = S2_2.select(['nir']);
var B8A_2 = S2_2.select(['nir2']);
var B11_2 = S2_2.select(['swir1']);
var B12_2 = S2_2.select(['swir2']);

var B2_3 = S2_3.select(['blue']);
var B3_3 = S2_3.select(['green']);
var B4_3 = S2_3.select(['red']);
var B5_3 = S2_3.select(['re1']);
var B6_3 = S2_3.select(['re2']);
var B7_3 = S2_3.select(['re3']);
var B8_3 = S2_3.select(['nir']);
var B8A_3 = S2_3.select(['nir2']);
var B11_3 = S2_3.select(['swir1']);
var B12_3 = S2_3.select(['swir2']);

var B2_4 = S2_4.select(['blue']);
var B3_4 = S2_4.select(['green']);
var B4_4 = S2_4.select(['red']);
var B5_4 = S2_4.select(['re1']);
var B6_4 = S2_4.select(['re2']);
var B7_4 = S2_4.select(['re3']);
var B8_4 = S2_4.select(['nir']);
var B8A_4 = S2_4.select(['nir2']);
var B11_4 = S2_4.select(['swir1']);
var B12_4 = S2_4.select(['swir2']);
//-----------------------------------------------------------------------------------------

var b2_1 = S2_1.select(['blue']);
var b3_1 = S2_1.select(['green']);
var b4_1 = S2_1.select(['red']);
var b8_1 = S2_1.select(['nir']);

var b2_2 = S2_2.select(['blue']);
var b3_2 = S2_2.select(['green']);
var b4_2 = S2_2.select(['red']);
var b8_2 = S2_2.select(['nir']);

var b2_3 = S2_3.select(['blue']);
var b3_3 = S2_3.select(['green']);
var b4_3 = S2_3.select(['red']);
var b8_3 = S2_3.select(['nir']);

var b2_4 = S2_4.select(['blue']);
var b3_4 = S2_4.select(['green']);
var b4_4 = S2_4.select(['red']);
var b8_4 = S2_4.select(['nir']);

//-----------------------------------------------------------------------------------------
//Calculate the image band difference
var b2_diff = b2_3.subtract(b2_4);
var b3_diff = b3_3.subtract(b3_4);
var b4_diff = b4_3.subtract(b4_4);
var b8_diff = b8_3.subtract(b8_4);
 
//print(b2_diff);
//print(b3_diff);
//print(b4_diff);
//print(b8_diff);

//1. Function to define training layer - binary training data
//var train = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/AMWI_HR_TrainData_PGNR_v8'); 
var train = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/AMWI_HR_TrainData_5Class'); 

var train = train.select(['b1']);
var Apal = ['#08306b', '#006d2c', '#8c510a', '#b8e186', '#fdae61'];
//var Apal =['#0057e7', '#ffa700', '#a1a0a0', '#db3236'];
 
//2. Specify input training area
var Plots = ee.FeatureCollection('users/Alberta_Sentinel_analysis_v20/grassland_input/amwi_boundary_final');

//3. Create image stack of radar, optical, and topographic variables for Random Forest Classification
var ImgStack1 = ee.Image([VVdiff,SWI]);


//print(ImgStack);
var sen2_composite = ee.Image([b4_diff,b8_diff]).rename(['b1_2','b1_3']);
//print (sen2_composite);
var ImgStack = ImgStack1.addBands(sen2_composite);
print(ImgStack);


//4. Function to create stratified sample points for image classification training / validation 
var points = ImgStack.addBands(train).stratifiedSample({
  numPoints: 2000, 
  classBand: "b1_4", 
  region: Plots, 
  scale: 10,
  tileScale: 4
}).randomColumn(); 

var training = points.filter(ee.Filter.lt('random', 0.7));
var validation = points.filter(ee.Filter.gte('random', 0.9)); 

//5. Apply Random Forest (RF) algorithim to train classifier
var classifier = ee.Classifier.randomForest({numberOfTrees: 500}).train(training, "b1_4", ImgStack.bandNames());

//6. Clip classified output to area of interest 
var pred = ImgStack.classify(classifier);
var pred_clip = pred.clip(SA);
var pred_clip1 = pred.clip(Grid1);
var pred_clip2 = pred.clip(Grid2);
var pred_clip3 = pred.clip(Grid3);
var pred_clip4 = pred.clip(Grid4);
var pred_clip5 = pred.clip(Grid5);
var pred_clip6 = pred.clip(Grid6);
var pred_clip7 = pred.clip(Grid7);
var pred_clip8 = pred.clip(Grid8);
var pred_clip9 = pred.clip(Grid9);
var pred_clip10 = pred.clip(Grid10);
var pred_clip11 = pred.clip(Grid11);
var pred_clip12 = pred.clip(Grid12);
var pred_clip13 = pred.clip(Grid13);
var pred_clip14 = pred.clip(Grid14);
var pred_clip15 = pred.clip(Grid15);
var pred_clip16 = pred.clip(Grid16);
var pred_clip17 = pred.clip(Grid17);
var pred_clip18 = pred.clip(Grid18);
var pred_clip19 = pred.clip(Grid19);
var pred_clip20 = pred.clip(Grid20);
var pred_clip21 = pred.clip(Grid21);
var pred_clip22 = pred.clip(Grid22);
var pred_clip23 = pred.clip(Grid23);
var pred_clip24 = pred.clip(Grid24);
var pred_clip25 = pred.clip(Grid25);
var pred_clip26 = pred.clip(Grid26);
var pred_clip27 = pred.clip(Grid27);
var pred_clip28 = pred.clip(Grid28);
var pred_clip29 = pred.clip(Grid29);
var pred_clip30 = pred.clip(Grid30);
var pred_clip31 = pred.clip(Grid31);
var pred_clip32 = pred.clip(Grid32); 

//7. Add results to map
//Map.addLayer(pred_clip, {min:1, max:5, palette:Apal}, 'Prediction');
Map.addLayer(train, {min:1, max:5, palette:Apal}, 'training data');
Map.addLayer(Plots);


//---------------------------------------------------------------------------------------------//  
//Export classified output
// -------------------------------------------------------------------------------------------//

Export.image.toDrive({
  image: pred_clip1,
  description: 'S2_VV_B4_B8_SWI_clip1',
  scale: 10,
  region: Grid1,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip2,
  description: 'S2_VV_B4_B8_SWI_clip2',
  scale: 10,
  region: Grid2,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip3,
  description: 'S2_VV_B4_B8_SWI_clip3',
  scale: 10,
  region: Grid3,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip4,
  description: 'S2_VV_B4_B8_SWI_clip4',
  scale: 10,
  region: Grid4,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip5,
  description: 'S2_VV_B4_B8_SWI_clip5',
  scale: 10,
  region: Grid5,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip6,
  description: 'S2_VV_B4_B8_SWI_clip6',
  scale: 10,
  region: Grid6,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip7,
  description: 'S2_VV_B4_B8_SWI_clip7',
  scale: 10,
  region: Grid7,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip8,
  description: 'S2_VV_B4_B8_SWI_clip8',
  scale: 10,
  region: Grid8,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip9,
  description: 'S2_VV_B4_B8_SWI_clip9',
  scale: 10,
  region: Grid9,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip10,
  description: 'S2_VV_B4_B8_SWI_clip10',
  scale: 10,
  region: Grid10,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip11,
  description: 'S2_VV_B4_B8_SWI_clip11',
  scale: 10,
  region: Grid11,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip12,
  description: 'S2_VV_B4_B8_SWI_clip12',
  scale: 10,
  region: Grid12,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip13,
  description: 'S2_VV_B4_B8_SWI_clip13',
  scale: 10,
  region: Grid13,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip14,
  description: 'S2_VV_B4_B8_SWI_clip14',
  scale: 10,
  region: Grid14,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip15,
  description: 'S2_VV_B4_B8_SWI_clip15',
  scale: 10,
  region: Grid15,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip16,
  description: 'S2_VV_B4_B8_SWI_clip16',
  scale: 10,
  region: Grid16,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip17,
  description: 'S2_VV_B4_B8_SWI_clip17',
  scale: 10,
  region: Grid17,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip18,
  description: 'S2_VV_B4_B8_SWI_clip18',
  scale: 10,
  region: Grid18,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip19,
  description: 'S2_VV_B4_B8_SWI_clip19',
  scale: 10,
  region: Grid19,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip20,
  description: 'S2_VV_B4_B8_SWI_clip20',
  scale: 10,
  region: Grid20,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip21,
  description: 'S2_VV_B4_B8_SWI_clip21',
  scale: 10,
  region: Grid21,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip22,
  description: 'S2_VV_B4_B8_SWI_clip22',
  scale: 10,
  region: Grid22,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip23,
  description: 'S2_VV_B4_B8_SWI_clip23',
  scale: 10,
  region: Grid23,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip24,
  description: 'S2_VV_B4_B8_SWI_clip24',
  scale: 10,
  region: Grid24,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip25,
  description: 'S2_VV_B4_B8_SWI_clip25',
  scale: 10,
  region: Grid25,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip26,
  description: 'S2_VV_B4_B8_SWI_clip26',
  scale: 10,
  region: Grid26,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip27,
  description: 'S2_VV_B4_B8_SWI_clip27',
  scale: 10,
  region: Grid27,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip28,
  description: 'S2_VV_B4_B8_SWI_clip28',
  scale: 10,
  region: Grid28,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip29,
  description: 'S2_VV_B4_B8_SWI_clip29',
  scale: 10,
  region: Grid29,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip30,
  description: 'S2_VV_B4_B8_SWI_clip30',
  scale: 10,
  region: Grid30,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip31,
  description: 'S2_VV_B4_B8_SWI_clip31',
  scale: 10,
  region: Grid31,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

Export.image.toDrive({
  image: pred_clip32,
  description: 'S2_VV_B4_B8_SWI_clip32',
  scale: 10,
  region: Grid32,
  folder: 'RF_PGNR_V16', 
  maxPixels: 10E10
});

// -------------------------------------------------------------------------------------------//
