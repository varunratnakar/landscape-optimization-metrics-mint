import os
from osgeo import gdal
from osgeo import ogr
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rs
from rasterio.features import shapes
from shapely.geometry import Polygon, MultiPolygon
import warnings
import sys

warnings.filterwarnings("ignore")

intersection_file = sys.argv[1]
flamelen_dir = sys.argv[2]
bldgs_tiles_json_file = sys.argv[3]
footprint_polygons_file = sys.argv[4]

output_damage_response_dir = sys.argv[5]
output_building_damage_dir = sys.argv[6]
output_tmp_dir = sys.argv[7]

def make_damage_response_rasters(filename,in_path,out_path, prefix):

    """
    This script takes fire footprint raster files for flame length and maps the values to 
    building damage response values (0, 25, 40, 55, 70, 85, 100) and writes a new raster

    inputs: 
    filename (str)
    in_path (str)
    out_path (str)
    prefix (str)

    outputs:
    damage response rasters
    """
   
    #open tif with gdal and convert to array
    # print(f'input raster: {os.path.join(in_path, filename)}')
    ds = gdal.Open(os.path.join(in_path, filename))
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    band = ds.GetRasterBand(1)
    array = band.ReadAsArray()
    
    # map flame length to building damage response values
    # approach as per the Wildfire Risk to Communities paper, https://www.fs.usda.gov/rds/archive/Catalog/RDS-2020-0016
    
    #0' flame --> 0 building response function value
    #0' < flame < 2' --> 25 
    #2-4' flame --> 40
    #4-6' flame --> 55
    #6-8' flame --> 70
    #8-12' flame --> 85
    #>12' flame --> 100
    
    array = np.where((array >= 12),100,array)
    array = np.where(((array >= 8) & (array < 12)),85,array)
    array = np.where(((array >= 6) & (array < 8)),70,array)
    array = np.where(((array >= 4) & (array < 6)),55,array)
    array = np.where(((array >= 2) & (array < 4)),40,array)
    array = np.where(((array > 0) & (array < 2)),25,array)
    array = np.where((array == 0),0,array)
    
    #binmask = np.where((array > 0),1,0)  # keep all the values that are greater than 0

    # export
    driver = gdal.GetDriverByName("GTiff")
    driver.Register()
    
    new_filename = f"{prefix}-{filename[-22:]}"
    # bin_filename = f"{prefix}-{filename[9:]}"
    # print(f'output raster: {new_filename}')
    # return
    #out_path='/landscape-optimization/test/'+new_filename
    outds = driver.Create(os.path.join(out_path, new_filename), xsize = array.shape[1],
                      ysize = array.shape[0], bands = 1, 
                      eType = gdal.GDT_Int16)
    outds.SetGeoTransform(gt)
    outds.SetProjection(proj)
    outband = outds.GetRasterBand(1)
    outband.WriteArray(array)
    outband.SetNoDataValue(np.nan)
    outband.FlushCache()

    # close your datasets and bands!
    outband = None
    outds = None


def groupby_multipoly(df, by, aggfunc="first"):
    
    """
    This function make multipolygons from polygons for fire footprints

    inputs:
    geopandas dataframe with polygons

    outputs:
    geopandas dataframe with multipolygons
    """
    
    data = df.drop(labels=df.geometry.name, axis=1)
    aggregated_data = data.groupby(by=by).agg(aggfunc)

    # Process spatial component
    def merge_geometries(block):
        return MultiPolygon(block.values)

    g = df.groupby(by=by, group_keys=False)[df.geometry.name].agg(
        merge_geometries
    )

    # Aggregate
    aggregated_geometry = gpd.GeoDataFrame(g, geometry=df.geometry.name, crs=df.crs)
    # Recombine
    aggregated = aggregated_geometry.join(aggregated_data)
    return aggregated

def make_raster(filename,shp_filename,in_path,out_path, prefix,out_path_tmp):

    """
    This function takes in damage response rasters and creates binary building 
    damage (tmp) rasters and building damage rasters with varying intensity.

    inputs:
    
    filename (str)
    shp_filename (str)
    in_path (str)
    out_path (str)
    prefix (str)
    out_path_tmp (str)

    outputs:

    raster files
    """
    
    fn_ras = os.path.join(in_path, 'damage_response-'+filename)
    # print('input raster: ',fn_ras) 
    ras_ds = gdal.Open(fn_ras)
    intensity_array = ras_ds.GetRasterBand(1).ReadAsArray()
    driver = ogr.GetDriverByName('ESRI Shapefile')
    vec_ds = driver.Open(shp_filename)
    lyr = vec_ds.GetLayer() 
    geot = ras_ds.GetGeoTransform()
    geo_proj = ras_ds.GetProjection() 
    
    # Setup the New Raster
    drv_tiff = gdal.GetDriverByName("GTiff") 
    out_net=os.path.join(out_path_tmp, f'{prefix}-binary-{filename}')
    # print(out_net)
    chn_ras_ds = drv_tiff.Create(out_net, ras_ds.RasterXSize, ras_ds.RasterYSize, 1, gdal.GDT_Float32)
    chn_ras_ds.SetGeoTransform(geot)
    chn_ras_ds.SetProjection(geo_proj) 
    gdal.RasterizeLayer(chn_ras_ds, [1], lyr, burn_values=[1])
    chn_ras_ds.GetRasterBand(1).SetNoDataValue(np.nan) 
    chn_ras_ds = None
    
    # open binary bldg damage raster back up
    ds = gdal.Open(out_net)
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    band = ds.GetRasterBand(1)
    bin_array = band.ReadAsArray()
    
    # get array and use as mask to get values from damage_response raster array
    damage_array = intensity_array*bin_array.astype(int)
    
    #print(damage_array.shape)
    
    # write new raster using chn_ras_ds.GetRasterBand(1).WriteArray(binmask)
    # export
    driver_ = gdal.GetDriverByName("GTiff")
    driver_.Register()
    out=os.path.join(out_path, f'{prefix}-{filename}')
    # print(f'output raster: {out}')
    outds = driver_.Create(out, xsize = damage_array.shape[1],
                      ysize = damage_array.shape[0], bands = 1, 
                      eType = gdal.GDT_Int16)
    outds.SetGeoTransform(gt)
    outds.SetProjection(proj)
    outband = outds.GetRasterBand(1)
    outband.WriteArray(damage_array)
    outband.SetNoDataValue(np.nan)
    outband.FlushCache()

    # close datasets and bands
    outband = None
    outds = None
    

if __name__ == "__main__":
    os.makedirs(output_damage_response_dir, exist_ok=True)
    os.makedirs(output_building_damage_dir, exist_ok=True)
    os.makedirs(output_tmp_dir, exist_ok=True)

    data = pd.read_csv(intersection_file)
    filenames = data['filename'].values

    # Create Damage response rasters

    in_path = flamelen_dir
    out_path = output_damage_response_dir

    prefix='damage_response'

    gdf = gpd.GeoDataFrame(columns=['feature'],geometry='feature',crs='EPSG:32610')
    count = 0
    for fn in filenames:
        # lp = fn.split('/')[-1]
        count += 1
        # print("current iter:", count)
        f = fn+'.tif'
        make_damage_response_rasters('flamelen-'+f,in_path,out_path,prefix)


    footprint_polygons = gpd.read_file(footprint_polygons_file)
    footprint_polygons = footprint_polygons.set_crs('EPSG:32610',allow_override=True)

    # make multipolygons from polygons for fire footprints
    grouped = groupby_multipoly(footprint_polygons, by='filename').reset_index()

    # load bldgs
    bldgs = gpd.read_file(bldgs_tiles_json_file).to_crs('EPSG:32610') 

    # intersect fire footprint MPs with bldgs polygons
    intersection = gpd.overlay(grouped, bldgs, how='intersection')

    # for given fire footprint, merge (essentially a groupby) bldg polygons into multipolygon
    gdf = intersection.dissolve(by='filename').reset_index()

    in_path = output_damage_response_dir
    out_path = output_building_damage_dir

    prefix='building_damage'
    
    out_path_tmp = output_tmp_dir

    for i in range(len(gdf)):
        #create shp files
        shp_filename=os.path.join(out_path_tmp, f"{prefix}-{gdf.iloc[i].filename[4:22]}.shp")
        gdf.iloc[[i]].to_file(driver = 'ESRI Shapefile', filename=shp_filename)
        
        root_filename=f"{gdf.iloc[i].filename[4:22]}.tif"
        #print(root_filename)
        make_raster(root_filename,shp_filename,in_path,out_path, prefix,out_path_tmp)                                                                                                                                                                                                              


    # Quality Control: test building damage rasters for all zero values (edge case)
    in_path = output_building_damage_dir
    out_path = output_tmp_dir

    files = [file for file in os.listdir(in_path) if file.endswith(".tif")]

    zero_rasters=[]

    for file in files:
        ds = gdal.Open(os.path.join(in_path, file))
        gt = ds.GetGeoTransform()
        proj = ds.GetProjection()
        band = ds.GetRasterBand(1)
        array = band.ReadAsArray()
        
        if array.sum()==0:
            zero_rasters.append(file)
            os.rename(os.path.join(in_path, file), os.path.join(out_path, file))
            
    print('# of zero rasters: ',len(zero_rasters))
    print('# of files checked: ',len(files))     
