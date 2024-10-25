from osgeo import gdal
from osgeo import ogr
import geopandas as gpd
import numpy as np
import rasterio as rs
from rasterio.features import shapes
from shapely.geometry import Polygon, MultiPolygon
import shutil
import sys
import os
import pandas as pd

toa_dir = sys.argv[1]
critical_habitats_shape_file = sys.argv[2]
footprint_polygons_file = sys.argv[3]

output_habitat_damage_dir = sys.argv[4]
output_tmp_dir = sys.argv[5]

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

def make_raster(filename,shp_filename,in_path,out_path, prefix):
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
    
    fn_ras = os.path.join(in_path, 'toa-'+filename)
    
    ras_ds = gdal.Open(fn_ras)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    vec_ds = driver.Open(shp_filename) 
    lyr = vec_ds.GetLayer() 
    geot = ras_ds.GetGeoTransform()
    geo_proj = ras_ds.GetProjection() 
    
    # Setup the New Raster
    drv_tiff = gdal.GetDriverByName("GTiff") 
    out_net=os.path.join(out_path, f'{prefix}-{filename}')

    # print('writing output raster: ',out_net)
    chn_ras_ds = drv_tiff.Create(out_net, ras_ds.RasterXSize, ras_ds.RasterYSize, 1, gdal.GDT_Float32)
    chn_ras_ds.SetGeoTransform(geot)
    chn_ras_ds.SetProjection(geo_proj) 
    gdal.RasterizeLayer(chn_ras_ds, [1], lyr, burn_values=[1])
    chn_ras_ds.GetRasterBand(1).SetNoDataValue(np.nan) 
    chn_ras_ds = None


if __name__ == "__main__":
    os.makedirs(output_habitat_damage_dir, exist_ok=True)
    os.makedirs(output_tmp_dir, exist_ok=True)

    # footprint_polygons
    footprint_polygons = gpd.read_file(footprint_polygons_file)
    footprint_polygons = footprint_polygons.set_crs('EPSG:32610',allow_override=True)

    # make multipolygons from polygons for fire footprints
    grouped = groupby_multipoly(footprint_polygons, by='filename').reset_index()
  
    # load critical habitat
    habitat = gpd.read_file(critical_habitats_shape_file).to_crs('EPSG:32610') 

    # intersect fire footprint MPs with habitat polygons
    intersection = gpd.overlay(grouped, habitat, how='intersection')

    # for given fire footprint, merge (essentially a groupby) bldg polygons into multipolygon
    gdf = intersection.dissolve(by='filename').reset_index()

    gdf_tmp_file = os.path.join(output_tmp_dir, 'habitat_gdf_old.geojson')
    gdf.to_file(gdf_tmp_file, driver='GeoJSON')
    gdf_read = gpd.read_file(gdf_tmp_file)

    in_path = toa_dir
    out_path = output_habitat_damage_dir

    prefix='habitat_damage'
    
    count = 0
    for i in range(len(gdf_read)):
        count += 1
        #create shp files
        shp_filename=os.path.join(output_tmp_dir, f"{prefix}-{gdf_read.iloc[i].filename[4:22]}.shp")
        gdf.iloc[[i]].to_file(driver = 'ESRI Shapefile', filename=shp_filename)
        root_filename=f"{gdf_read.iloc[i].filename[4:22]}.tif"
        # print("filenames", shp_filename, root_filename)
        # print(count, "f:", root_filename)
        make_raster(root_filename,shp_filename,in_path,out_path, prefix)