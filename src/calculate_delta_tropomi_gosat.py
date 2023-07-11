import time
import os
import glob
import yaml
from geopy import distance
import pandas as pd
import numpy as np
from scipy import interpolate
import multiprocessing

from utilities import get_gosat_df, get_tropomi_df

with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

tropomi_files = glob.glob(os.path.join(config["StorageDir"], "tropomi", "*.nc"))
gosat_files = glob.glob(os.path.join(config['StorageDir'], 'gosat', "*.nc"))
gosat_df = get_gosat_df(gosat_files)

# Function that return a dataframe of TROPOMI+GOSAT pairs corresponding to a input TROPOMI file
def get_tropomi_gosat_pairs(tropomi_file):
    
    print(tropomi_file, flush=True)

    tropomi_df = get_tropomi_df(tropomi_file)
    
    tropomi_gosat_pairs_list = []
    for idx in tropomi_df.index:
        tropomi_instance = tropomi_df.loc[idx:idx]
        
        # Subset GOSAT to less than 1 hour apart from this TROPOMI observation
        gosat_subset_temporal = gosat_df[(gosat_df["time"] >= (tropomi_instance["time"].iloc[0] - pd.Timedelta(hours=1))) & 
                                         (gosat_df["time"] <= (tropomi_instance["time"].iloc[0] + pd.Timedelta(hours=1)))].copy()
        if len(gosat_subset_temporal) == 0:
            continue
        
        # Subset GOSAT to less than 5 km from this TROPOMI observation
        gosat_subset_temporal["tropomi_distance"] = gosat_subset_temporal.apply(lambda x: distance.distance((tropomi_instance['latitude'].iloc[0], tropomi_instance['longitude'].iloc[0]), (x['latitude'], x['longitude'])).km, axis=1)
        gosat_subset_temporal_spatial = gosat_subset_temporal[gosat_subset_temporal["tropomi_distance"] < 5]
        
        if len(gosat_subset_temporal_spatial) == 0:
            continue
        
        # Identify the closest matching GOSAT observation and concatenate it with the TROPOMI row
        closest_gosat = gosat_subset_temporal_spatial["tropomi_distance"].idxmin()
        tropomi_gosat_pair_df =  pd.concat([tropomi_instance.add_prefix("tropomi_").reset_index(drop=True),
                                            gosat_subset_temporal_spatial.loc[closest_gosat:closest_gosat].add_prefix("gosat_").reset_index(drop=True)], axis=1)
        tropomi_gosat_pair = tropomi_gosat_pair_df.iloc[0]
        
        # Calculate ∆(TROPOMI-GOSAT) using equations (A1-A4)
        gosat_mask = tropomi_gosat_pair["gosat_ch4_profile_apriori"][::-1] != -9999.99 # when there are 19 pressure levels
        gosat_prior = tropomi_gosat_pair["gosat_ch4_profile_apriori"][::-1][gosat_mask] # [ppb]
        gosat_pressure_levels = 100*tropomi_gosat_pair["gosat_pressure_levels"][::-1][gosat_mask] # [Pa]
        f_interp_gosat_prior_to_tropomi_pressure_grid = interpolate.interp1d(gosat_pressure_levels, gosat_prior, bounds_error=False, fill_value="extrapolate")
        tropomi_pressure_levels = [tropomi_gosat_pair["tropomi_surface_pressure"] - i*tropomi_gosat_pair["tropomi_pressure_interval"] for i in np.arange(0,13)][::-1] # [Pa]
        tropomi_pressure_layers = np.array([(tropomi_pressure_levels[i]+tropomi_pressure_levels[i+1])/2 for i in np.arange(0,12)])

        c_Tr = tropomi_gosat_pair["tropomi_xch4_corrected"] # [ppb]
        h_T = tropomi_gosat_pair["tropomi_dry_air_subcolumns"]/np.sum(tropomi_gosat_pair["tropomi_dry_air_subcolumns"]) # [unitless]
        A_T = tropomi_gosat_pair["tropomi_column_averaging_kernel"] # [unitless]
        x_Ga = f_interp_gosat_prior_to_tropomi_pressure_grid(tropomi_pressure_layers) # [ppb]
        x_Ta = 1e9*tropomi_gosat_pair["tropomi_methane_profile_apriori"]/tropomi_gosat_pair["tropomi_dry_air_subcolumns"] # [ppb]
        c_Gr = tropomi_gosat_pair["gosat_xch4"] - config["GlobalOffsetGOSAT"] # [ppb], adjust to GOSAT having a global mean bias of 0 relative to GGG2020
        c_Ga = np.sum(h_T*x_Ga) # [ppb]
        x_Gr = x_Ga * (c_Gr/c_Ga) # [ppb]

        c_T_star = c_Tr + np.sum(h_T*(1-A_T)*(x_Ga-x_Ta)) # [ppb]
        c_G_star = np.sum(h_T*(x_Ga+(A_T*(x_Gr-x_Ga)))) # [ppb]
        delta_tropomi_gosat = c_T_star - c_G_star # [ppb]
        tropomi_gosat_pair_df.loc[0,"delta_tropomi_gosat"] = delta_tropomi_gosat
        
        tropomi_gosat_pairs_list.append(tropomi_gosat_pair_df)
        
    if len(tropomi_gosat_pairs_list) == 0:
        return None
    else:
        tropomi_gosat_pairs_df = pd.concat(tropomi_gosat_pairs_list, ignore_index=True)
        return tropomi_gosat_pairs_df

if __name__ == "__main__":
    num_processes = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(get_tropomi_gosat_pairs, tropomi_files)
        pool.close()
        pool.join()
    results = [r for r in results if r is not None]
    tropomi_gosat_pairs = pd.concat(results, ignore_index=True)
    tropomi_gosat_pairs.to_pickle(os.path.join(config['StorageDir'], 'processed', 'tropomi_gosat_pairs.pkl'))