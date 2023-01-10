import pandas as pd
from geopy import distance
import datetime
import multiprocessing
import numpy as np
import sys
import yaml
import os

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
gosat_processed_dir = config["StorageDir"] + "/pkl/gosat/"
tropomi_processed_dir = config["StorageDir"] + "/pkl/tropomi/"
paired_dir = config["StorageDir"] + "/pkl/paired/"

# Read in dataframe of all gosat data
gosat_data = pd.read_pickle(gosat_processed_dir+"gosat_data.pkl")
gosat_data = gosat_data.add_prefix("gosat_")

def f_pair_tropomi_gosat(file):

    '''
    Function that looks at every dataframe of TROPOMI data in tropomi_process_dir
    For each TROPOMI observation in that file, search for a GOSAT observation pair in time & space (of the "same" atmosphere)

    Arguments
        file [str] : name of the TROPOMI processed dataframe

    Returns
        none [none] : saves the paired TROPOMI and GOSAT data to paired_dir but the function doesn't return anything
    '''
    
    # Read in dataframe of the TROPOMI observations
    tropomi_data = pd.read_pickle(tropomi_processed_dir+file)
    tropomi_data = tropomi_data.add_prefix("tropomi_")
    
    # Dataframe that will contain paired TROPOMI and GOSAT observations
    paired_data = pd.DataFrame()

    # If there is any TROPOMI data available
    if len(tropomi_data) != 0:

        for idx in tropomi_data.index:
            tropomi_time = tropomi_data.loc[idx,"tropomi_time"]

            # GOSAT observations within TimeThreshold minutes of tropomi observation idx
            gosat_time_subset = gosat_data[(gosat_data["gosat_time"] >= (tropomi_time - datetime.timedelta(minutes=config["TimeThreshold"]))) &
                                            (gosat_data["gosat_time"] <= (tropomi_time + datetime.timedelta(minutes=config["TimeThreshold"])))].copy()

            # If there are any GOSAT observations that meet the time threshold, search the distance threshold
            if len(gosat_time_subset) != 0:

                # GOSAT observations within TimeThreshold minutes & DistanceThreshold km of tropomi observation idx
                gosat_time_subset["tropomi_gosat_distance_km"] = gosat_time_subset.apply(lambda row: distance.distance((row["gosat_latitude"],row["gosat_longitude"]),(tropomi_data.loc[idx,"tropomi_latitude"],tropomi_data.loc[idx,"tropomi_longitude"])).km, axis=1)
                gosat_time_distance_subset = gosat_time_subset[gosat_time_subset["tropomi_gosat_distance_km"] <= config["DistanceThreshold"]].copy()

                # If any GOSAT observations remain, pair the closest one in distance with the TROPOMI observation; add to paired_data
                if len(gosat_time_distance_subset) != 0:
                    min_distance_idx = gosat_time_distance_subset.index[gosat_time_distance_subset["tropomi_gosat_distance_km"] == gosat_time_distance_subset["tropomi_gosat_distance_km"].min()]
                    assert (len(min_distance_idx) == 1)
                    min_distance_idx = min_distance_idx[0]
                    temp_paired_data = pd.concat([gosat_time_distance_subset.loc[[min_distance_idx]].reset_index(drop=True),tropomi_data.loc[[idx]].reset_index(drop=True)], axis=1)
                    paired_data = pd.concat([paired_data,temp_paired_data], ignore_index=True)

    # Save paired data with one .pkl file per one .nc file from TROPOMI
    paired_data.to_pickle(paired_dir+"paired_"+file)
        
if __name__ == "__main__":

    # open list of files to be processed
    with open("tropomi.txt") as f:
        filelist = [line.strip() for line in f]

    # split list of files into number of slurm array tasks
    # sys.argv[1] = SLURM_ARRAY_TASK_ID
    # sys.argv[2] = SLURM_ARRAY_TASK_COUNT
    filelist = np.array_split(filelist,int(sys.argv[2]))[int(sys.argv[1])]
    
    # set up so that jobs can be requeued and progress is not lost
    filelist_not_done_yet = []
    for file in filelist:
        if ("paired_"+file) not in os.listdir(paired_dir):
            filelist_not_done_yet.append(file)

    pool = multiprocessing.Pool(config["Cores"])
    pool.map(f_pair_tropomi_gosat, filelist_not_done_yet)
    pool.close()