import pandas as pd
import pickle
from flaml import AutoML
import datetime
import yaml

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
flaml_dir = config["StorageDir"] + "/pkl/flaml/"

def run_flaml():
    
    '''
    Function that trains a machine learning model to predict delta(TROPOMI-GOSAT)
    Uses FLAML for tuning and gives 1 hour to each of lightgbm, random forest, and xgboost

    Arguments
        none [none] : no arguments are explicity fed to the function; the filtered tropomi gosat pairs are used though

    Returns
        none [none] : nothing is explicity returned from the function, though the best model for lgbm, rf, and xgboost are saved to flaml_dir
    '''

    # Read in processed TROPOMI/GOSAT pairs
    pairs = pd.read_pickle(flaml_dir+"filtered_tropomi_gosat_pairs.pkl").sort_values("tropomi_time").reset_index(drop=True)

    # Define predictor variables
    X = pairs[["tropomi_sza","tropomi_raa","tropomi_surface_altitude","tropomi_surface_altitude_stdv","tropomi_u10","tropomi_v10",\
               "tropomi_cirrus_reflectance","tropomi_xch4_precision","tropomi_xch4_apriori","tropomi_fluorescence","tropomi_co_column",\
               "tropomi_co_column_precision","tropomi_h2o_column","tropomi_h2o_column_precision","tropomi_aerosol_size","tropomi_aerosol_size_precision","tropomi_aerosol_column",\
               "tropomi_aerosol_column_precision","tropomi_aerosol_altitude","tropomi_aerosol_altitude_precision","tropomi_nir_surface_albedo","tropomi_swir_surface_albedo",\
               "tropomi_nir_surface_albedo_precision","tropomi_swir_surface_albedo_precision","tropomi_nir_aerosol_optical_thickness","tropomi_swir_aerosol_optical_thickness",\
               "tropomi_nir_chi_squared_band","tropomi_swir_chi_squared_band","tropomi_ground_pixel","tropomi_landflag"]]

    # Define variable to be predicted
    y = pairs["delta_tropomi_gosat"]

    # Training set is 2018-2020
    train_index = pairs[pairs["tropomi_time"] < datetime.datetime(2021,1,1,0,0,0)].index
    X_train = X.loc[train_index]
    y_train = y.loc[train_index]

    # For each algorithm give it the time specific in config.yml to train using 8 cores
    for estimator in ["lgbm", "rf", "xgboost"]:
        automl = AutoML()
        automl.fit(X_train, y_train, task="regression", metric="mse", time_budget=config["TimeFLAML"], n_jobs=8, estimator_list = [estimator], eval_method="cv", split_type="time", n_splits=10)
        # Save the model
        with open(flaml_dir + f"model_{estimator}.pkl", "wb") as handle:
            pickle.dump(automl, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
if __name__ == "__main__":
    run_flaml()
