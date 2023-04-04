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

from utilities import get_blended_df, get_tccon_df

with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# TCCON station altitudes
tccon_station_altitudes = {"easttroutlake01": 0.50, "reunion01": 0.09, "hefei01": 0.03, "pasadena01": 0.23, "saga01": 0.01,
                           "edwards01": 0.70, "xianghe01": 0.04, "parkfalls01": 0.44, "rikubetsu01": 0.38, "jpl02": 0.39,
                           "nicosia01": 0.19, "burgos01": 0.04, "sodankyla01": 0.19, "paris01": 0.06, "eureka01": 0.61,
                           "orleans01": 0.13, "nyalesund01": 0.02, "lauder02": 0.37,"tsukuba02": 0.03, "izana01": 2.37,
                           "lamont01": 0.32, "lauder03": 0.37, "bremen01": 0.03, "garmisch01": 0.74, "karlsruhe01": 0.12} # [km a.s.l.]

# Co-location criteria
max_distance_diff = {x:100 for x in tccon_station_altitudes.keys() if "edwards01" not in x} # [km]
max_distance_diff["edwards01"] = 50 # [km]
max_elevation_diff = {x:250 for x in tccon_station_altitudes.keys()} # [m]
max_time_diff = {x:pd.Timedelta(hours=1) for x in tccon_station_altitudes.keys()} # [h]

def get_tropomi_tccon_pairs(tropomi_file, tccon_file):

    print(tropomi_file, flush=True)

    # Get station name so we can access the dicts for altitude, max_distance_diff, max_elevation_diff, and max_time_diff
    with Dataset(tccon_file) as ds:
        tccon_station_name = ds.long_name
    
    # Dataframe of TCCON data for this station
    tccon_df = get_tccon_df(tccon_file)

    # Dataframe of this blended file
    tropomi_df = get_blended_df(tropomi_file)
    
    # List to append single row dataframes of TROPOMI+TCCON pairs to
    tropomi_tccon_pairs_list = []
    
    for idx in tropomi_df.index:
            
        tropomi_instance = tropomi_df.loc[idx:idx]
        
        # If TROPOMI observation not within the allowable elevation difference, continue to next TROPOMI observation
        if np.abs(1000*tccon_station_altitudes[tccon_station_name] - tropomi_instance["surface_altitude"].iloc[0]) > max_elevation_diff[tccon_station_name]:
            continue
            
        # Subset the TCCON to less than the allowable time difference from the TROPOMI observation
        tccon_subset_temporal = tccon_df[(tccon_df["time"] >= tropomi_instance["time"].iloc[0] - max_time_diff[tccon_station_name]) &
                                            (tccon_df["time"] <= tropomi_instance["time"].iloc[0] + max_time_diff[tccon_station_name])].copy()
        
        if len(tccon_subset_temporal) == 0:
            continue
            
        # Subset TCCON to less than the max_distance_diff from this TROPOMI observation
        tccon_subset_temporal["tropomi_distance"] = tccon_subset_temporal.apply(lambda x: distance.distance((tropomi_instance['latitude'].iloc[0], tropomi_instance['longitude'].iloc[0]), (x['latitude'], x['longitude'])).km, axis=1)
        tccon_subset_temporal_spatial = tccon_subset_temporal[tccon_subset_temporal["tropomi_distance"] < max_distance_diff[tccon_station_name]]
        
        if len(tccon_subset_temporal_spatial) == 0:
            continue
            
        # Identify the closest matching TCCON observation and concatenate it with the TROPOMI observation
        closest_tccon = (tccon_subset_temporal_spatial["time"] - tropomi_instance["time"].iloc[0]).abs().idxmin()
        tropomi_tccon_pair_df = pd.concat([tropomi_instance.add_prefix("tropomi_").reset_index(drop=True),
                                            tccon_subset_temporal_spatial.loc[closest_tccon:closest_tccon].add_prefix("tccon_").reset_index(drop=True)], axis=1)
        tropomi_tccon_pair = tropomi_tccon_pair_df.iloc[0]
        
        # Calculate ∆(TROPOMI-TCCON) using equation (A6)
        tropomi_pressure_levels = [tropomi_tccon_pair["tropomi_surface_pressure"] - i*tropomi_tccon_pair["tropomi_pressure_interval"] for i in np.arange(0,13)][::-1] # [Pa]
        tropomi_pressure_layers = np.array([(tropomi_pressure_levels[i]+tropomi_pressure_levels[i+1])/2 for i in np.arange(0,12)])
        tccon_pressure_levels = tropomi_tccon_pair["tccon_prior_pressure"][::-1] # [Pa]
        tccon_prior = (1e9*tropomi_tccon_pair["tccon_prior_wet_ch4"] * (1 + tropomi_tccon_pair["tccon_prior_wet_h2o"]/(1-tropomi_tccon_pair["tccon_prior_wet_h2o"])))[::-1] # [ppb], have to "dry" the profile
        f_interp_tccon_prior_to_tropomi_pressure_grid = interpolate.interp1d(tccon_pressure_levels, tccon_prior, bounds_error=False, fill_value="extrapolate")
        
        c_Tr = tropomi_tccon_pair["tropomi_xch4_corrected"] # [ppb]
        h_T = tropomi_tccon_pair["tropomi_dry_air_subcolumns"]/np.sum(tropomi_tccon_pair["tropomi_dry_air_subcolumns"]) # [unitless]
        A_T = tropomi_tccon_pair["tropomi_column_averaging_kernel"] # [unitless]
        x_Fa = f_interp_tccon_prior_to_tropomi_pressure_grid(tropomi_pressure_layers) # [ppb]
        x_Ta = 1e9*tropomi_tccon_pair["tropomi_methane_profile_apriori"]/tropomi_tccon_pair["tropomi_dry_air_subcolumns"] # [ppb]
        c_Fr = tropomi_tccon_pair["tccon_xch4"] # [ppb]
        c_Fa = np.sum(h_T*x_Fa) # [ppb]
        x_Fr = x_Fa * (c_Fr/c_Fa) # [ppb]
        
        delta_tropomi_tccon = (c_Tr + np.sum(h_T*(1-A_T)*(x_Fa-x_Ta))) - (np.sum(h_T*(x_Fa+(A_T*(x_Fr-x_Fa)))))
        tropomi_tccon_pair_df.loc[0,"delta_tropomi_tccon"] = delta_tropomi_tccon # [ppb]
        tropomi_tccon_pair_df.loc[0,"delta_blended_tccon"] = delta_tropomi_tccon - tropomi_tccon_pair["tropomi_xch4_corrected"] + tropomi_tccon_pair["tropomi_xch4_blended"]
        
        tropomi_tccon_pairs_list.append(tropomi_tccon_pair_df)
        
    if len (tropomi_tccon_pairs_list) == 0:
        return None
    else:
        tropomi_tccon_pairs_df = pd.concat(tropomi_tccon_pairs_list, ignore_index=True)
        return tropomi_tccon_pairs_df

if __name__ == "__main__":

    # Get list of all TCCON files that intersect the time period of 30 Apr 2018 to 31 Dec 2021
    tccon_dir = os.path.join(config['StorageDir'], 'tccon')
    tccon_files = glob.glob(os.path.join(tccon_dir, "*.nc"))
    tccon_files = [tccon_file for tccon_file in tccon_files
               if pd.Interval(pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").min(),
                               pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").max())
                   .overlaps(pd.Interval(pd.to_datetime('2018-04-30 00:00'),
                                          pd.to_datetime('2021-12-31 23:59')))]

    # Dict to hold TROPOMI+TCCON pairs for each station
    dict_delta_tropomi_tccon = {}

    # Find pairs for each TCCON station
    for tccon_file in tccon_files:
        
        with Dataset(tccon_file) as ds:
            tccon_station_name = ds.long_name
            print(f"Starting {tccon_station_name}", flush=True)

        # Get list of blended files that intersect the TCCON file time range
        start_tccon_dt = pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").min()
        end_tccon_dt = pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").max()
        all_blended_files = glob.glob(os.path.join(config["StorageDir"], "blended", "*.nc"))
        tropomi_files = [file for file in all_blended_files 
                            if pd.Interval(pd.to_datetime(file.split("_")[12]), pd.to_datetime(file.split("_")[13]))
                            .overlaps(pd.Interval(start_tccon_dt,end_tccon_dt))]
        tropomi_files.sort()

        # Get list of inputs for this TCCON station then run get_tropomi_tccon_pairs in parallel
        inputs = [(tropomi_file, tccon_file) for tropomi_file in tropomi_files]
        num_processes = multiprocessing.cpu_count()  # Use all available CPU cores
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.starmap(get_tropomi_tccon_pairs, inputs)
            pool.close()
            pool.join()
        results = [r for r in results if r is not None]
        if len(results) != 0:
            dict_delta_tropomi_tccon[tccon_station_name] = pd.concat(results, ignore_index=True)

        print(f"Finished {tccon_station_name}", flush=True)

    # Save dictionary
    with open(os.path.join(config['StorageDir'], 'processed', 'delta_tropomi_tccon.pkl'), "wb") as f:
        pickle.dump(dict_delta_tropomi_tccon, f)