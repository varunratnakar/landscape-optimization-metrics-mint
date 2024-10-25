import os
from osgeo import gdal
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rs
import rasterio.features
import xarray as xr
import affine
import shapely.geometry as sg
import warnings
import sys
warnings.filterwarnings("ignore")

intersection_file = sys.argv[1]
toa_dir = sys.argv[2]

output_burned_area_dir = sys.argv[3]
output_footprints_polygon_file = sys.argv[4]


def polygonize(da: xr.DataArray) -> gpd.GeoDataFrame:
    """
    Polygonize a 2D-DataArray into a GeoDataFrame of polygons.

    Parameters
    ----------
    da : xr.DataArray

    Returns
    -------
    polygonized : geopandas.GeoDataFrame
    """
    if da.dims != ("y", "x"):
        raise ValueError('Dimensions must be ("y", "x")')

    values = da.values
    transform = da.attrs.get("transform", None)
    if transform is None:
        raise ValueError("transform is required in da.attrs")
    transform = affine.Affine(*transform)
    shapes = rasterio.features.shapes(values, transform=transform)

    geometries = []
    colvalues = []
    for (geom, colval) in shapes:
        geometries.append(sg.Polygon(geom["coordinates"][0]))
        colvalues.append(colval)

    gdf = gpd.GeoDataFrame({"value": colvalues, "geometry": geometries})
    gdf.crs = da.attrs.get("crs")
    return gdf

def make_polygons(filename, in_path, out_path, prefix):
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
    # print(f'input raster: {os.path.join(in_path, filename)}')

    ds = gdal.Open(os.path.join(in_path, filename))
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
    # print(f'output raster: {bin_filename}')
    outds = driver.Create(os.path.join(out_path, bin_filename), xsize = binmask.shape[1],
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
    bin_ = xr.open_rasterio(os.path.join(out_path, bin_filename)).squeeze('band', drop=True)

    polygons = polygonize(bin_)
    perimeter = polygons[polygons['value']==1.0]  # select outside polygon
    
    # returns polygon in geodataframe
    perimeter['filename'] = filename
    return perimeter


if __name__ == "__main__":
    os.makedirs(output_burned_area_dir, exist_ok=True)
    
    data = pd.read_csv(intersection_file)
    filenames = data['filename'].values

    in_path = toa_dir
    out_path = output_burned_area_dir

    prefix='burned_area'
    gdf = gpd.GeoDataFrame(columns=['feature'],geometry='feature',crs='EPSG:32610')
    count = 0
    for fn in filenames:
        # lp = fn.split('/')[-1]
        count += 1
        # print("current iter:", count)
        f = fn+'.tif'
        new_gdf = make_polygons('toa-'+f, in_path, out_path, prefix)
        gdf = pd.concat([gdf, new_gdf])

    gdf = gdf.drop(['feature'], axis=1)
    gdf = gdf.set_geometry("geometry")
    gdf.to_file(output_footprints_polygon_file, driver="GeoJSON")