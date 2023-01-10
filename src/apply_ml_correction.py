import pickle
import pandas as pd
import numpy as np
import multiprocessing
import yaml
import sys
import pickle

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
blended_tropomi_gosat_dir = config["StorageDir"] + "/pkl/blended_tropomi_gosat/"
tropomi_processed_dir = config["StorageDir"] + "/pkl/tropomi/"

# Read in the machine learning model
with open(config["StorageDir"] + f"/pkl/flaml/model_{config['Model']}.pkl", "rb") as handle:
    model = pickle.load(handle)

def f_apply_ml_correction(file):

    # Read in the TROPOMI file
    df = pd.read_pickle(tropomi_processed_dir+file)

    # If there are TROPOMI observations, write an equivalent dataframe with just the blended TROPOMI+GOSAT xch4 variables
    if len(df) != 0:
        
        X = df[["sza","raa","surface_altitude","surface_altitude_stdv","u10","v10",\
               "cirrus_reflectance","xch4_precision","xch4_apriori","fluorescence","co_column",\
               "co_column_precision","h2o_column","h2o_column_precision","aerosol_size","aerosol_size_precision","aerosol_column",\
               "aerosol_column_precision","aerosol_altitude","aerosol_altitude_precision","nir_surface_albedo","swir_surface_albedo",\
               "nir_surface_albedo_precision","swir_surface_albedo_precision","nir_aerosol_optical_thickness","swir_aerosol_optical_thickness",\
               "nir_chi_squared_band","swir_chi_squared_band","ground_pixel","landflag"]].add_prefix("tropomi_")

        # Predict the delta(TROPOMI-GOSAT) then remove from TROPOMI
        y = df["xch4_corrected"] - (model.predict(X)*config["a"] + config["b"])

        final_df = pd.DataFrame(index=df.index, data={"xch4_blended_tropomi_gosat":y})

    # If there are not, just have an empty dataframe
    else:
        final_df = pd.DataFrame(columns=["xch4_blended_tropomi_gosat"])

    final_df.to_pickle(blended_tropomi_gosat_dir+file)
        
if __name__ == "__main__":
    
    # open list of files to be processed
    with open("tropomi.txt") as f:
        filelist = [line.strip() for line in f]
    
    # split list of files into number of slurm array tasks
    # sys.argv[1] = SLURM_ARRAY_TASK_ID
    # sys.argv[2] = SLURM_ARRAY_TASK_COUNT
    pool = multiprocessing.Pool(config["Cores"])
    pool.map(f_apply_ml_correction, np.array_split(filelist,int(sys.argv[2]))[int(sys.argv[1])])
    pool.close()
