from time import time_ns
import yaml
import geopandas as gpd
import pandas as pd
import shapely.speedups

shapely.speedups.enable()

import sys
args = sys.argv
# print(args)
print('Running intersection script')

# Required arguments
# Argument 1: path to burn unit polygons
# Argument 2: path to ignition points csv
# Argument 3: path to output file

rx_burn_units_path = args[1]
full_ignitions_file_path = args[2]


full_ignitions_df = pd.read_csv(full_ignitions_file_path)
ign_gdf = gpd.GeoDataFrame(full_ignitions_df, geometry=gpd.points_from_xy(full_ignitions_df['x_ignition'], full_ignitions_df['y_ignition']))
# print(ign_gdf)
ign_gdf = ign_gdf.set_crs('EPSG:32610',allow_override=True)
# ign_gdf = ign_gdf.to_crs('EPSG:32610')
# print(ign_gdf)
rx_burn_units = gpd.read_file(rx_burn_units_path)
rx_burn_units = rx_burn_units.to_crs('epsg:32610')
# rx_burn_units.set_crs('EPSG:32610',allow_override=True)


# intersection = gpd.overlay(rx_burn_units, ign_gdf, how='intersection')
intersection = gpd.overlay(ign_gdf, rx_burn_units, how='intersection')
# print(intersection)
new_gdf = intersection.dissolve(by='filename').reset_index()
# print(new_gdf)
new_gdf.drop('geometry', axis=1).to_csv(args[3], index=False)
# new_gdf.to_json()
filenames = new_gdf['filename'].values
print('Total intersections= ', len(filenames))

# tiles = [fn.split('-')[0] for fn in filenames]
# print(set(tiles))


