import geopandas
import pandas as pd
import shapely
import numpy as np
import yaml
import datetime
import os
import pickle

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
flaml_dir = config["StorageDir"] + "/pkl/flaml/"
regridded_dir = config["StorageDir"] + "/geojson/"

def regrid(filename,start_datetime,end_datetime,lon_lower,lon_upper,lon_res,lat_lower,lat_upper,lat_res):
    
    '''
    Use geopandas to regrid tropomi/gosat pairs (using tropomi pixel centers)

    Arguments
        filename [str] : name of regridded geojson file to save to regridded_dir
        start_datetime [datetime] : beginning of period to regrid
        end_datetime [datetime] : end of period to regrid
        lon_lower [float] : longitude of left bound of regridded area
        lon_upper [float] : longitude of right bound of regridded area
        lon_res [float] : longitude resolution of final regridded
        lat_lower [float] : latitude of lower bound of regridded area
        lat_upper [float] : latitude of upper bound of regridded area
        lat_res [float] : latitude resolution of final regridded

    Returns
        none [none] : nothing explicity returned but a geojson file is saved to regridded_dir
    '''

    # Do not run if file already exists
    if os.path.exists(regridded_dir+filename):
        return
    
    # Read in TROPOMI/GOSAT pairs
    pairs = pd.read_pickle(flaml_dir+"filtered_tropomi_gosat_pairs.pkl")

    # Add column for predicted delta_GOSAT_TROPOMI
    X = pairs[["tropomi_sza","tropomi_raa","tropomi_surface_altitude","tropomi_surface_altitude_stdv","tropomi_u10","tropomi_v10",\
               "tropomi_cirrus_reflectance","tropomi_xch4_precision","tropomi_xch4_apriori","tropomi_fluorescence","tropomi_co_column",\
               "tropomi_co_column_precision","tropomi_h2o_column","tropomi_h2o_column_precision","tropomi_aerosol_size","tropomi_aerosol_size_precision","tropomi_aerosol_column",\
               "tropomi_aerosol_column_precision","tropomi_aerosol_altitude","tropomi_aerosol_altitude_precision","tropomi_nir_surface_albedo","tropomi_swir_surface_albedo",\
               "tropomi_nir_surface_albedo_precision","tropomi_swir_surface_albedo_precision","tropomi_nir_aerosol_optical_thickness","tropomi_swir_aerosol_optical_thickness",\
               "tropomi_nir_chi_squared_band","tropomi_swir_chi_squared_band","tropomi_ground_pixel","tropomi_landflag"]]
    with open(config["StorageDir"] + f"/pkl/flaml/model_{config['Model']}.pkl", "rb") as handle:
        model = pickle.load(handle)
    pairs["predicted_delta_tropomi_gosat"] = (model.predict(X)*config["a"] + config["b"])

    if "Land" in filename:
        pairs = pairs[(pairs["tropomi_time"] >= start_datetime) & (pairs["tropomi_time"] <= end_datetime) &\
                      ((pairs["tropomi_landflag"] == 0) | (pairs["tropomi_landflag"] == 2))].reset_index(drop=True).copy()
    elif "Water" in filename:
        pairs = pairs[(pairs["tropomi_time"] >= start_datetime) & (pairs["tropomi_time"] <= end_datetime) &\
                      ((pairs["tropomi_landflag"] == 1) | (pairs["tropomi_landflag"] == 3))].reset_index(drop=True).copy()
    else:
        pairs = pairs[(pairs["tropomi_time"] >= start_datetime) & (pairs["tropomi_time"] <= end_datetime)].reset_index(drop=True).copy()
    
    # Put paired data into a geodataframe
    pairs_gdf = geopandas.GeoDataFrame(pairs, geometry=geopandas.points_from_xy(pairs.tropomi_longitude, pairs.tropomi_latitude), crs="EPSG:4326")

    # Build a grid to average onto
    lon_min,lon_max= lon_lower+lon_res/2, lon_upper-lon_res/2
    lat_min,lat_max= lat_lower+lat_res/2, lat_upper-lat_res/2
    grid_points = []

    for lon in np.arange(lon_min, lon_max+lon_res, lon_res):
        for lat in np.arange(lat_min, lat_max+lat_res, lat_res):
            xmin = lon - lon_res/2
            xmax = lon + lon_res/2
            ymin = lat - lat_res/2
            ymax = lat + lat_res/2
            grid_points.append(shapely.geometry.box(xmin,ymin,xmax,ymax))

    grid = geopandas.GeoDataFrame(grid_points,columns=['geometry'],crs="EPSG:4326")

    # Associate a grid cell with each paired obs.
    merged_pairs = geopandas.sjoin(pairs_gdf, grid, how="left", predicate="within")

    # Average the paired measurements within each grid cell
    for idx in grid.index:
        # Loop through each grid cell in the standard grid 
        # For each grid cell, take every paired measurement that is within it and average/count them
        pairs_subset = merged_pairs[merged_pairs["index_right"] == idx]
        grid.loc[idx,"true_delta_tropomi_gosat"] = pairs_subset["delta_tropomi_gosat"].mean()
        grid.loc[idx,"predicted_delta_tropomi_gosat"] = pairs_subset["predicted_delta_tropomi_gosat"].mean()
        grid.loc[idx,"residual_delta_tropomi_gosat"] = (pairs_subset["delta_tropomi_gosat"] - pairs_subset["predicted_delta_tropomi_gosat"]).mean()
        grid.loc[idx, "paired_n"] = len(pairs_subset)
        
    # Save geojson file
    grid.to_file(regridded_dir+filename, driver='GeoJSON')
    
if __name__ == "__main__":
    
    # 2021 (Land, Ocean, All) [2 x 2.5, 0.25 x 0.3125]
    regrid("Paired_Global_Land_2021_2_25deg.geojson",datetime.datetime(2021,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,2.5,-90,90,2)
    regrid("Paired_Global_Water_2021_2_25deg.geojson",datetime.datetime(2021,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,2.5,-90,90,2)
    regrid("Paired_Global_All_2021_2_25deg.geojson",datetime.datetime(2021,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,2.5,-90,90,2)
    
    regrid("Paired_Global_Land_2021_025_03125deg.geojson",datetime.datetime(2021,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,0.3125,-90,90,0.25)
    regrid("Paired_Global_Water_2021_025_03125deg.geojson",datetime.datetime(2021,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,0.3125,-90,90,0.25)
    regrid("Paired_Global_All_2021_025_03125deg.geojson",datetime.datetime(2021,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,0.3125,-90,90,0.25)
    
    # 2018-2021 (Land, Ocean, All) [2 x 2.5, 0.25 x 0.3125]
    regrid("Paired_Global_Land_2018_2021_2_25deg.geojson",datetime.datetime(2018,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,2.5,-90,90,2)
    regrid("Paired_Global_Water_2018_2021_2_25deg.geojson",datetime.datetime(2018,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,2.5,-90,90,2)
    regrid("Paired_Global_All_2018_2021_2_25deg.geojson",datetime.datetime(2018,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,2.5,-90,90,2)
    
    regrid("Paired_Global_Land_2018_2021_025_03125deg.geojson",datetime.datetime(2018,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,0.3125,-90,90,0.25)
    regrid("Paired_Global_Water_2018_2021_025_03125deg.geojson",datetime.datetime(2018,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,0.3125,-90,90,0.25)
    regrid("Paired_Global_All_2018_2021_025_03125deg.geojson",datetime.datetime(2018,1,1,0,0,0),datetime.datetime(2021,12,31,23,59,59),-180,180,0.3125,-90,90,0.25)