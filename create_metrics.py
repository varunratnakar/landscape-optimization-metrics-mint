# Imports
import sys
import os
import shutil
from glob import glob
from zipfile import ZipFile

from osgeo import gdal
from osgeo import ogr
import geopandas as gpd
import pandas as pd
import numpy as np
import rioxarray as rxr

from shapely.geometry import MultiPolygon
from geocube.vector import vectorize


def unzip(zipfile, todir):
    if not os.path.exists(todir):
        with ZipFile(zipfile) as zip_file:
            filenames = zip_file.namelist()
            zip_maindir = None
            if filenames[0].endswith("/"):
                zip_maindir = filenames[0]
            
            for filename in zip_file.namelist():
                if filename.startswith("__"):
                    continue
                if zip_maindir and not filename.startswith(zip_maindir):
                    continue

                tofilename = filename
                if zip_maindir:
                    tofilename = filename.replace(zip_maindir, todir + "/", 1)

                if tofilename.endswith("/"):
                    os.makedirs(tofilename)
                else:
                    source = zip_file.open(filename)
                    target = open(tofilename, "wb")
                    with source, target:
                        shutil.copyfileobj(source, target)


def _make_polygons(filename, in_path, out_path, prefix):
    """
    Function takes fire footprint raster, converts the intensity values to binary values (burn/not burn) and writes new raster.
    It also creates a geopandas dataframe of footprint polygons

    Parameters
    ----------
    filename (str)
    in_path (str)
    out_path (str)
    prefix (str) ex. burned_area

    Returns
    -------
    geopandas dataframe with footprint polygon
    """

    #open tif with gdal and convert to array
    print(f'input raster: {filename}')
    ds = gdal.Open(in_path+ '/' + filename)
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    band = ds.GetRasterBand(1)
    array = band.ReadAsArray()
    
    # create burn/not burn binary mask

    binmask = np.where((array > 0),1,0)  # keep all the values that are greater than 0

    # export
    driver = gdal.GetDriverByName("GTiff")
    driver.Register()
    
    bin_filename = f"{prefix}-{filename[4:]}"
    print(f'output raster: {bin_filename}')
    outds = driver.Create(out_path+ '/' + bin_filename, xsize = binmask.shape[1],
                      ysize = binmask.shape[0], bands = 1, 
                      eType = gdal.GDT_Int16)
    outds.SetGeoTransform(gt)
    outds.SetProjection(proj)
    outband = outds.GetRasterBand(1)
    outband.WriteArray(binmask)
    outband.SetNoDataValue(np.nan)
    outband.FlushCache()

    # close your datasets and bands!
    outband = None
    outds = None
    
    #open bin_mask and polygonize it
    bin_ = rxr.open_rasterio(out_path+ '/' + bin_filename).squeeze('band', drop=True)

    polygons = vectorize(bin_)
    polygons.rename(columns={None: "value"},inplace=True)
    
    # polygons = polygonize(bin_, gt)
    
    perimeter = polygons[polygons['value']==1.0]  # select outside polygon
    
    # returns polygon in geodataframe
    perimeter['filename'] = filename
    return perimeter


def _make_damage_response_rasters(filename,in_path,out_path, prefix):
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
    ds = gdal.Open(in_path+ '/' + filename)
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    band = ds.GetRasterBand(1)
    array = band.ReadAsArray()
    
    # map flame length to building damage response values
    # approach as per the Wildfire Risk to Communities paper, https://www.fs.usda.gov/rds/archive/Catalog/RDS-2020-0016
    
    #<10' flame --> 0 building response function value
    #>=10' flame --> 100
    
    # New modified binary threshold for calculating building damage rasters
    array = np. where((array < 10),0, array)
    array = np.where((array >= 10), 100, array)
    
    #binmask = np.where((array > 0),1,0)  # keep all the values that are greater than 0

    # export
    driver = gdal.GetDriverByName("GTiff")
    driver.Register()
    
    new_filename = f"{prefix}-{filename[-22:]}"
    #bin_filename = f"damage_response-{filename[9:]}"
    print(new_filename)
    
    #out_path='../test/'+new_filename
    outds = driver.Create(out_path+ '/' + new_filename, xsize = array.shape[1],
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

def _groupby_multipoly(df, by, aggfunc="first"):
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

def _make_bldg_damage_raster(filename,shp_filename,in_path,out_path, prefix,out_path_tmp):

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
    
    fn_ras = in_path+ '/' + 'damage_response-'+filename
    print('input raster: ',fn_ras) 
    ras_ds = gdal.Open(fn_ras)
    intensity_array = ras_ds.GetRasterBand(1).ReadAsArray()
    driver = ogr.GetDriverByName('ESRI Shapefile')
    vec_ds = driver.Open(shp_filename)
    lyr = vec_ds.GetLayer() 
    geot = ras_ds.GetGeoTransform()
    geo_proj = ras_ds.GetProjection() 
    
    # Setup the New Raster
    drv_tiff = gdal.GetDriverByName("GTiff") 
    out_net=f'{out_path_tmp+ "/" + prefix}-binary-{filename}'
    print(out_net)
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
    damage_array = np.where((damage_array > 0), 1, 0) # keep all the values that are greater than 0
    
    #print(damage_array.shape)
    
    # write new raster using chn_ras_ds.GetRasterBand(1).WriteArray(binmask)
    # export
    driver_ = gdal.GetDriverByName("GTiff")
    driver_.Register()
    out=f'{out_path+ "/" + prefix}-{filename}'
    print(f'output raster: {out}')
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

def _make_raster(filename, shp_filename, in_path, out_path, prefix):
    """
    This function writes raster files

    parameters:
    
        filename (str)
        shp_filename (str)
        in_path (str)
        out_path (str)
        prefix (str)

    outputs:

        raster files
    """
    
    fn_ras = in_path+'/' + 'toa-'+filename
    ras_ds = gdal.Open(fn_ras)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    vec_ds = driver.Open(shp_filename) 
    lyr = vec_ds.GetLayer() 
    geot = ras_ds.GetGeoTransform()
    geo_proj = ras_ds.GetProjection() 
    
    # Setup the New Raster
    drv_tiff = gdal.GetDriverByName("GTiff") 
    out_net=out_path+ '/' + prefix+filename
    print('writing output raster: ',out_net)
    chn_ras_ds = drv_tiff.Create(out_net, ras_ds.RasterXSize, ras_ds.RasterYSize, 1, gdal.GDT_Float32)
    chn_ras_ds.SetGeoTransform(geot)
    chn_ras_ds.SetProjection(geo_proj) 
    gdal.RasterizeLayer(chn_ras_ds, [1], lyr, burn_values=[1])
    chn_ras_ds.GetRasterBand(1).SetNoDataValue(np.nan) 
    chn_ras_ds = None

def make_burned_area_rasters(toa_dir, out_burned_area_dir, footprints_polygons_path):
    prefix='burned_area'

    files = [file for file in os.listdir(toa_dir) if file.endswith(".tif")]
    gdf = gpd.GeoDataFrame(columns=['feature'],geometry='feature',crs='EPSG:32610')
    for file in files:
        new_gdf = _make_polygons(file, toa_dir, out_burned_area_dir, prefix)
        gdf = pd.concat([gdf, new_gdf], ignore_index=True)
    gdf = gdf.drop(['feature'], axis=1)
    gdf = gdf.set_geometry("geometry")
    gdf.to_file(footprints_polygons_path, driver="GeoJSON")   

def make_damage_response_rasters(flamelen_dir, out_damage_response_dir):
    prefix='damage_response'
    files = [file for file in os.listdir(flamelen_dir) if file.endswith(".tif")]
    for file in files:
        _make_damage_response_rasters(file, flamelen_dir, out_damage_response_dir, prefix)


def make_building_damage_rasters(damage_response_dir, ms_buildings_tiles_path, footprints_polygons_path, tmp_dir, out_building_damage_dir):
    # footprint_polygons

    footprint_polygons = gpd.read_file(footprints_polygons_path)
    footprint_polygons = footprint_polygons.set_crs('EPSG:32610',allow_override=True)

    # make multipolygons from polygons for fire footprints
    grouped = _groupby_multipoly(footprint_polygons, by='filename').reset_index()

    # load bldgs
    bldgs = gpd.read_file(ms_buildings_tiles_path).to_crs('EPSG:32610') 

    # intersect fire footprint MPs with bldgs polygons
    intersection = gpd.overlay(grouped, bldgs, how='intersection')

    # for given fire footprint, merge (essentially a groupby) bldg polygons into multipolygon
    gdf = intersection.dissolve(by='filename').reset_index()

    prefix='building_damage'

    for i in range(len(gdf)):
        #create shp files
        shp_filename=f"{tmp_dir + '/' + prefix}-{gdf.iloc[i].filename[4:22]}.shp"
        gdf.iloc[[i]].to_file(driver = 'ESRI Shapefile', filename=shp_filename)
        
        root_filename=f"{gdf.iloc[i].filename[4:22]}.tif"
        #print(root_filename)
        _make_bldg_damage_raster(root_filename,shp_filename,damage_response_dir,out_building_damage_dir, prefix, tmp_dir)  

def make_habitat_damage_raster(toa_dir, critical_habitats_dir, footprints_polygons_path, tmp_dir, out_habitat_damage_dir):
    # footprint_polygons
    footprint_polygons = gpd.read_file(footprints_polygons_path)
    footprint_polygons = footprint_polygons.set_crs('EPSG:32610',allow_override=True)

    # make multipolygons from polygons for fire footprints
    grouped = _groupby_multipoly(footprint_polygons, by='filename').reset_index()
  
    # load critical habitat
    habitat_shape_file = glob(os.path.join(critical_habitats_dir, '*.shp'))[0]
    habitat = gpd.read_file(habitat_shape_file).to_crs('EPSG:32610') 

    # intersect fire footprint MPs with habitat polygons
    intersection = gpd.overlay(grouped, habitat, how='intersection')

    # for given fire footprint, merge (essentially a groupby) bldg polygons into multipolygon
    gdf = intersection.dissolve(by='filename').reset_index()

    prefix='habitat_damage-'

    for i in range(len(gdf)):
        #create shp files
        shp_filename=f"{tmp_dir + '/' + prefix}-{gdf.iloc[i].filename[4:22]}.shp"
        gdf.iloc[[i]].to_file(driver = 'ESRI Shapefile', filename=shp_filename)
        
        root_filename=f"{gdf.iloc[i].filename[4:22]}.tif"
        #print(root_filename)
        _make_raster(root_filename,shp_filename,toa_dir,out_habitat_damage_dir, prefix)




if __name__ == "__main__":
    # check if there is a config file in the arguments
    # fetch data
    # - ms_buildings_tiles_file
    # - critical_habitats_file
    # - flamelens_file
    # - toas_file

    pass