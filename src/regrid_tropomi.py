import geopandas
import pandas as pd
import shapely
import os
import numpy as np
import yaml

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
tropomi_processed_dir = config["StorageDir"] + "/pkl/tropomi/"
blended_tropomi_gosat_dir = config["StorageDir"] + "/pkl/blended_tropomi_gosat/"
regridded_dir = config["StorageDir"] + "/geojson/"

def regrid(filename,start_orbit,end_orbit,lon_lower,lon_upper,lon_res,lat_lower,lat_upper,lat_res):

    '''
    Use geopandas to regrid tropomi onto a regular grid (using tropomi pixel centers)

    Arguments
        filename [str] : name of regridded geojson file to save to regridded_dir
        start_orbit [int] : beginning orbit number to include in regrid
        end_orbit [int] : end orbit number to include in regrid
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

    print(filename, flush=True)

    tropomi = None
    filelist = os.listdir(tropomi_processed_dir)
    filelist.sort()

    # Build a dataframe of TROPOMI data (including the blended product)
    if type(start_orbit) == int: # if consecutive orbit numbers (e.g., to do a full year)
        for file in filelist:
            orbit_number = int(file[16:21])
            if (orbit_number >= start_orbit) & (orbit_number <= end_orbit):
                df1 = pd.read_pickle(tropomi_processed_dir+file)
                if len(df1) != 0:
                    df2 = pd.read_pickle(blended_tropomi_gosat_dir+file) # corrected with ML
                    assert len(df1) == len(df2)
                    df = pd.concat([df1[["latitude","longitude","xch4_corrected","swir_surface_albedo","nir_surface_albedo","aerosol_size","landflag"]],df2], axis=1) # add the corrected data
                    df = df[(df["latitude"] >= lat_lower) & (df["latitude"] <= lat_upper) & (df["longitude"] >= lon_lower) & (df["longitude"] <= lon_upper)].reset_index(drop=True)

                    if "Water" in filename:
                        df = df[(df["landflag"] == 1) | (df["landflag"] == 3)].reset_index(drop=True)
                    elif "Land" in filename:
                        df = df[(df["landflag"] == 0) | (df["landflag"] == 2)].reset_index(drop=True)

                    # This is to avoid pd.concat chaning the dtype to object
                    if len(df) != 0:
                        if tropomi is None:
                            tropomi = df
                        else:
                            tropomi = pd.concat([tropomi,df], ignore_index=True)

    if type(start_orbit) == list: # if there are gaps in the orbit numbers (e.g., to do seasons over a year)
        for file in filelist:
            orbit_number = int(file[16:21])
            if any([(orbit_number >= start_orbit[i]) and (orbit_number <= end_orbit[i]) for i in range(len(start_orbit))]):
                df1 = pd.read_pickle(tropomi_processed_dir+file)
                if len(df1) != 0:
                    df2 = pd.read_pickle(blended_tropomi_gosat_dir+file) # corrected with ML
                    assert len(df1) == len(df2)
                    df = pd.concat([df1[["latitude","longitude","xch4_corrected","swir_surface_albedo","nir_surface_albedo","aerosol_size","landflag"]],df2], axis=1) # add the corrected data
                    df = df[(df["latitude"] >= lat_lower) & (df["latitude"] <= lat_upper) & (df["longitude"] >= lon_lower) & (df["longitude"] <= lon_upper)].reset_index(drop=True)

                    if "Water" in filename:
                        df = df[(df["landflag"] == 1) | (df["landflag"] == 3)].reset_index(drop=True)
                    elif "Land" in filename:
                        df = df[(df["landflag"] == 0) | (df["landflag"] == 2)].reset_index(drop=True)

                    # This is to avoid pd.concat chaning the dtype to object
                    if len(df) != 0:
                        if tropomi is None:
                            tropomi = df
                        else:
                            tropomi = pd.concat([tropomi,df], ignore_index=True)

                        
    # Put TROPOMI data into a geodataframe
    tropomi_gdf = geopandas.GeoDataFrame(tropomi, geometry=geopandas.points_from_xy(tropomi.longitude, tropomi.latitude), crs="EPSG:4326")

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

    # Associate a grid cell with each GOSAT and each TROPOMI measurement
    merged_tropomi = geopandas.sjoin(tropomi_gdf, grid, how="left", predicate="within")

    # Average the TROPOMI measurements within each grid cell
    for idx in grid.index:
        # Loop through each grid cell in the standard grid 
        # For each grid cell, take every tropomi measurement that is within it and average/count them
        tropomi_subset = merged_tropomi[merged_tropomi["index_right"] == idx]
        for var in ["xch4_corrected","xch4_blended_tropomi_gosat","swir_surface_albedo","nir_surface_albedo","aerosol_size"]:
            grid.loc[idx, "tropomi_"+var] = tropomi_subset[var].mean()
        grid.loc[idx, "tropomi_n"] = len(tropomi_subset)
        
    grid.to_file(regridded_dir+filename, driver='GeoJSON')
    
if __name__ == "__main__":
    regrid("TROPOMI_Global_2021_1deg.geojson",16679,21856,-180,180,1,-90,90,1)
    regrid("TROPOMI_Land_2021_1deg.geojson",16679,21856,-180,180,1,-90,90,1)
    regrid("TROPOMI_Water_2021_1deg.geojson",16679,21856,-180,180,1,-90,90,1)
    regrid("TROPOMI_Global_Jan_Mar_2021_1deg.geojson",16679,17956,-180,180,1,-90,90,1)
    regrid("TROPOMI_Global_Apr_Jun_2021_1deg.geojson",17957,19247,-180,180,1,-90,90,1)
    regrid("TROPOMI_Global_Jul_Sep_2021_1deg.geojson",19248,20551,-180,180,1,-90,90,1)
    regrid("TROPOMI_Global_Oct_Dec_2021_1deg.geojson",20552,21856,-180,180,1,-90,90,1)