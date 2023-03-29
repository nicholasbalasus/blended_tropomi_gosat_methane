import time
import os
import glob
import yaml
from geopy import distance
import pandas as pd
import numpy as np
from scipy import interpolate
import multiprocessing
from netCDF4 import Dataset
import pickle

from utilities import get_gosat_df, get_tccon_df

with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# TCCON station altitudes
tccon_station_altitudes = {"easttroutlake01": 0.50, "reunion01": 0.09, "hefei01": 0.03, "pasadena01": 0.23, "saga01": 0.01,
                           "edwards01": 0.70, "xianghe01": 0.04, "parkfalls01": 0.44, "rikubetsu01": 0.38, "jpl02": 0.39,
                           "nicosia01": 0.19, "burgos01": 0.04, "sodankyla01": 0.19, "paris01": 0.06, "eureka01": 0.61,
                           "orleans01": 0.13, "nyalesund01": 0.02, "lauder02": 0.37,"tsukuba02": 0.03, "izana01": 2.37,
                           "lamont01": 0.32, "lauder03": 0.37, "bremen01": 0.03, "garmisch01": 0.74, "karlsruhe01": 0.12} # [km a.s.l.]

# Co-location criteria
max_distance_diff = {x:500 for x in tccon_station_altitudes.keys() if "edwards01" not in x} # [km]
max_distance_diff["edwards01"] = 50 # [km]
max_elevation_diff = {x:250 for x in tccon_station_altitudes.keys()} # [m]
max_time_diff = {x:pd.Timedelta(hours=2) for x in tccon_station_altitudes.keys()} # [h]

# GOSAT data
gosat_files = glob.glob(os.path.join(config['StorageDir'], 'gosat', "*.nc"))
gosat_df = get_gosat_df(gosat_files)

# Function that return a dataframe of GOSAT+TCCON pairs corresponding to a input TCCON file
def get_gosat_tccon_pairs(tccon_file):
    
    # Get station name so we can access the dicts for altitude, max_distance_diff, max_elevation_diff, and max_time_diff
    with Dataset(tccon_file) as ds:
        tccon_station_name = ds.long_name
        
    tccon_df = get_tccon_df(tccon_file)
    
    gosat_tccon_pairs_list = []
    for idx in gosat_df.index:
        
        gosat_instance = gosat_df.loc[idx:idx]
        
        # If GOSAT observation is not within the allowable altitude difference, return None
        if np.abs(1000*tccon_station_altitudes[tccon_station_name] - gosat_instance["surface_altitude"].iloc[0]) > max_elevation_diff[tccon_station_name]:
            continue
        
        # Subset the TCCON to less than allowable time difference from GOSAT observation
        tccon_subset_temporal = tccon_df[(tccon_df["time"] >= gosat_instance["time"].iloc[0] - max_time_diff[tccon_station_name]) & 
                                         (tccon_df["time"] <= gosat_instance["time"].iloc[0] + max_time_diff[tccon_station_name])].copy()
        
        if len(tccon_subset_temporal) == 0:
            continue
        
        # Subset TCCON to less than max_distance_diff from this GOSAT observation
        tccon_subset_temporal["gosat_distance"] = tccon_subset_temporal.apply(lambda x: distance.distance((gosat_instance['latitude'].iloc[0], gosat_instance['longitude'].iloc[0]), (x['latitude'], x['longitude'])).km, axis=1)
        tccon_subset_temporal_spatial = tccon_subset_temporal[tccon_subset_temporal["gosat_distance"] < max_distance_diff[tccon_station_name]]
        
        if len(tccon_subset_temporal_spatial) == 0:
            continue
        
        # Identify the closest matching TCCON observation and concatenate it with the GOSAT row
        closest_tccon = (tccon_subset_temporal_spatial["time"] - gosat_instance["time"].iloc[0]).abs().idxmin()
        gosat_tccon_pair_df = pd.concat([gosat_instance.add_prefix("gosat_").reset_index(drop=True),
                                         tccon_subset_temporal_spatial.loc[closest_tccon:closest_tccon].add_prefix("tccon_").reset_index(drop=True)], axis=1)
        gosat_tccon_pair = gosat_tccon_pair_df.iloc[0]
        
        # Calculate ∆(TROPOMI-GOSAT) using equations (A1-A4)
        gosat_mask = gosat_tccon_pair["gosat_ch4_profile_apriori"][::-1] != -9999.99 # when there are 19 pressure levels
        gosat_pressure_levels = 100*gosat_tccon_pair["gosat_pressure_levels"][::-1][gosat_mask] # [Pa]
        tccon_pressure_levels = gosat_tccon_pair["tccon_prior_pressure"][::-1] # [Pa]
        tccon_prior = (1e9*gosat_tccon_pair["tccon_prior_wet_ch4"] * (1 + gosat_tccon_pair["tccon_prior_wet_h2o"]/(1-gosat_tccon_pair["tccon_prior_wet_h2o"])))[::-1] # [ppb], have to "dry" the profile
        f_interp_tccon_prior_to_gosat_pressure_grid = interpolate.interp1d(tccon_pressure_levels, tccon_prior, bounds_error=False, fill_value="extrapolate")

        c_Gr = gosat_tccon_pair["gosat_xch4"] # [ppb]
        h_G = gosat_tccon_pair["gosat_pressure_weight"][::-1][gosat_mask] # [unitless]
        A_G = gosat_tccon_pair["gosat_xch4_averaging_kernel"][::-1][gosat_mask] # [unitless]
        x_Fa = f_interp_tccon_prior_to_gosat_pressure_grid(np.array(gosat_pressure_levels)) # [ppb]
        x_Ga = gosat_tccon_pair["gosat_ch4_profile_apriori"][::-1][gosat_mask] # [ppb]
        c_Fr = gosat_tccon_pair["tccon_xch4"] # [ppb]
        c_Fa = np.sum(h_G*x_Fa) # [ppb]
        x_Fr = x_Fa * (c_Fr/c_Fa) # [ppb]

        delta_gosat_tccon = (c_Gr+np.sum(h_G*(1-A_G)*(x_Fa-x_Ga))) - (np.sum(h_G*(x_Fa+(A_G*(x_Fr-x_Fa)))))
        gosat_tccon_pair_df.loc[0,"delta_gosat_tccon"] = delta_gosat_tccon # [ppb]
        
        gosat_tccon_pairs_list.append(gosat_tccon_pair_df)
        
    if len(gosat_tccon_pairs_list) == 0:
        return None, None
    else:
        gosat_tccon_pairs_df = pd.concat(gosat_tccon_pairs_list, ignore_index=True)
        return tccon_station_name, gosat_tccon_pairs_df

if __name__ == "__main__":

    # Get list of all TCCON files that intersect the time period of 30 Apr 2018 to 31 Dec 2021
    tccon_dir = os.path.join(config['StorageDir'], 'tccon')
    tccon_files = glob.glob(os.path.join(tccon_dir, "*.nc"))
    tccon_files = [tccon_file for tccon_file in tccon_files
               if pd.Interval(pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").min(),
                               pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").max())
                   .overlaps(pd.Interval(pd.to_datetime('2018-04-30 00:00'),
                                          pd.to_datetime('2021-12-31 23:59')))]

    # Get dataframe for each station in parallel
    num_processes = multiprocessing.cpu_count()  # Use all available CPU cores
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(get_gosat_tccon_pairs, tccon_files)
        pool.close()
        pool.join()
        
    # Form dictionary of form {station name: pandas DataFrame of ∆(GOSAT-TCCON)}
    dict_delta_gosat_tccon = {x:y for x,y in results if x is not None}

    # Save dictionary
    with open(os.path.join(config['StorageDir'], 'processed', 'delta_gosat_tccon.pkl'), "wb") as f:
        pickle.dump(dict_delta_gosat_tccon, f)