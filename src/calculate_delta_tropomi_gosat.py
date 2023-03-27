import os
import pandas as pd
import numpy as np
import glob
import yaml
from netCDF4 import Dataset
from multiprocessing import Pool
from geopy import distance
from scipy import interpolate
import sys
import pickle

# Load configuration file
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Get number of available cores
n_jobs = len(os.sched_getaffinity(0))

def find_closest_gosat(tropomi_row):

    # Set the maximum time difference and distance threshold for matching TROPOMI and GOSAT observations.
    max_time_diff = pd.Timedelta(hours=1)
    max_distance = 5  # km
    
    # Subset the GOSAT DataFrame to include only observations within the time window.
    gosat_subset_temporal = gosat_df[(gosat_df["time"] >= tropomi_row["time"].iloc[0] - max_time_diff) & 
                                     (gosat_df["time"] <= tropomi_row["time"].iloc[0] + max_time_diff)].copy()
    
    # If no GOSAT observations are found within the time window, return None.
    if len(gosat_subset_temporal) == 0:
        return None
        
    # Calculate the distances between the TROPOMI row and each GOSAT row in the temporal subset.
    gosat_subset_temporal["tropomi_distance"] = gosat_subset_temporal.apply(lambda x: distance.distance((tropomi_row['latitude'].iloc[0], tropomi_row['longitude'].iloc[0]), (x['latitude'], x['longitude'])).km, axis=1)
    
    # Subset the GOSAT DataFrame to include only observations within the distance threshold.
    gosat_subset_temporal_spatial = gosat_subset_temporal[gosat_subset_temporal["tropomi_distance"] < max_distance]
    
    # If no GOSAT observations are found within the distance threshold, return None.
    if len(gosat_subset_temporal_spatial) == 0:
        return None
    
    # Identify the closest matching GOSAT observation and concatenate it with the TROPOMI row.
    else:
        closest_gosat = gosat_subset_temporal_spatial["tropomi_distance"].idxmin()
        tropomi_gosat_pair_df =  pd.concat([tropomi_row.add_prefix("tropomi_").reset_index(drop=True),
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
        
        return tropomi_gosat_pair_df
    
def find_tropomi_gosat_pairs_for_batch_of_tropomi_instances(tropomi_df_index_batch):
    
    # Create an empty list to hold dataframes for each TROPOMI+GOSAT pair
    list_batch_dfs = []
    
    # For each TROPOMI instance, append a dataframe of itself and its GOSAT pair to list_batch_dfs (if a pair exists)
    for i in tropomi_df_index_batch:
        df = find_closest_gosat(tropomi_chunk_this_task.loc[i:i]) # weird notation is to give a single row dataframe instead of a series
        if df is not None:
            list_batch_dfs.append(df)
            
    # If this batch of TROPOMI instances had any pairs, put them in a single dataframe and return it
    if len(list_batch_dfs) > 0:
        batch_df = pd.concat(list_batch_dfs, ignore_index=True)
        return batch_df
    else:
        return None

if __name__ == "__main__":

    # This processing is running across an array of SLURM jobs (1 job per 1 chunk file generated in process_tropomi_gosat.py)
    SLURM_ARRAY_TASK_ID = int(sys.argv[1])
    SLURM_ARRAY_TASK_COUNT = int(sys.argv[2])

    # Only run if file hasn't already been made
    if not os.path.isfile(os.path.join(config['StorageDir'], 'tmp', f'tropomi_gosat_pairs_chunk{SLURM_ARRAY_TASK_ID}of{SLURM_ARRAY_TASK_COUNT-1}.pkl')):

        # Get GOSAT dataframe from previous module
        gosat_df = pd.read_pickle(os.path.join(config['StorageDir'], 'processed', 'gosat.pkl'))

        # Get list of all temporoary chunks of processed TROPOMI data
        tropomi_chunk_dir = os.path.join(config['StorageDir'], 'tmp')
        tropomi_chunk_files = glob.glob(os.path.join(tropomi_chunk_dir, 'chunk*.pkl'))
        tropomi_chunk_files.sort()
        assert(len(tropomi_chunk_files) == SLURM_ARRAY_TASK_COUNT)
        print(tropomi_chunk_files[SLURM_ARRAY_TASK_ID], flush=True)
        tropomi_chunk_this_task = pd.read_pickle(tropomi_chunk_files[SLURM_ARRAY_TASK_ID])
        tropomi_df_index_batches = np.array_split(tropomi_chunk_this_task.index, n_jobs)

        # Get a dataframe of all of the TROPOMI/GOSAT pairs from this chunk of TROPOMI data
        with Pool(n_jobs) as p:
            all_dfs = p.map(find_tropomi_gosat_pairs_for_batch_of_tropomi_instances, tropomi_df_index_batches)
        all_dfs = [x for x in all_dfs if x is not None]
        if len(all_dfs) > 0:
            tropomi_gosat_pairs_this_task = pd.concat(all_dfs, ignore_index=True)
            tropomi_gosat_pairs_this_task.to_pickle(os.path.join(config['StorageDir'], 'tmp', f'tropomi_gosat_pairs_chunk{SLURM_ARRAY_TASK_ID}of{SLURM_ARRAY_TASK_COUNT-1}.pkl'), protocol=pickle.HIGHEST_PROTOCOL)
        else:
            # Need to create an empty dataframe with the correct column names so as to not mess up the bottom portion
            empty_df = pd.DataFrame({}, columns=list(tropomi_chunk_this_task.add_prefix("tropomi_").columns) + list(gosat_df.add_prefix("gosat_").columns) + ["delta_tropomi_gosat"])
            empty_df.to_pickle(os.path.join(config['StorageDir'], 'tmp', f'tropomi_gosat_pairs_chunk{SLURM_ARRAY_TASK_ID}of{SLURM_ARRAY_TASK_COUNT-1}.pkl'), protocol=pickle.HIGHEST_PROTOCOL)

        # We have created 128 files of 'tropomi_gosat_pairs_chunkxof127' (where x starts at 0). We need to merge them to one dataframe.
        # Only can do this after all 128 files are written. This will only be true for one of the 128 tasks.
        all_files_created = [f'tropomi_gosat_pairs_chunk{x}of{SLURM_ARRAY_TASK_COUNT-1}.pkl' for x in range(SLURM_ARRAY_TASK_COUNT-1)]
        if all([os.path.isfile(os.path.join(config["StorageDir"], 'tmp', f)) for f in all_files_created]):
            tropomi_gosat_pairs_all_tasks = pd.concat([pd.read_pickle(os.path.join(config["StorageDir"], 'tmp', f)) for f in all_files_created], ignore_index=True)
            tropomi_gosat_pairs_all_tasks = tropomi_gosat_pairs_all_tasks.sort_values("tropomi_time").reset_index(drop=True)
            tropomi_gosat_pairs_all_tasks.to_pickle(os.path.join(config['StorageDir'], 'processed', 'tropomi_gosat_pairs.pkl'))