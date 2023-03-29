from flaml import AutoML
import pandas as pd
import os
import yaml
import pickle

# Load configuration file
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

def f_run_flaml():

    # Read in TROPOMI/GOSAT pairs
    tropomi_gosat_pairs = pd.read_pickle(os.path.join(config["StorageDir"], "processed", "tropomi_gosat_pairs.pkl"))

    # Define predictor variables
    X = tropomi_gosat_pairs[["tropomi_solar_zenith_angle","tropomi_relative_azimuth_angle","tropomi_across_track_pixel_index",
                         "tropomi_surface_classification","tropomi_surface_altitude","tropomi_surface_altitude_precision",
                         "tropomi_eastward_wind","tropomi_northward_wind","tropomi_xch4_apriori","tropomi_reflectance_cirrus_VIIRS_SWIR",
                         "tropomi_xch4_precision","tropomi_fluorescence","tropomi_co_column","tropomi_co_column_precision",
                         "tropomi_h2o_column","tropomi_h2o_column_precision","tropomi_aerosol_size","tropomi_aerosol_size_precision",
                         "tropomi_aerosol_height","tropomi_aerosol_height_precision","tropomi_aerosol_column","tropomi_aerosol_column_precision",
                         "tropomi_surface_albedo_SWIR","tropomi_surface_albedo_SWIR_precision","tropomi_surface_albedo_NIR",
                         "tropomi_surface_albedo_NIR_precision","tropomi_aerosol_optical_thickness_SWIR","tropomi_aerosol_optical_thickness_NIR",
                         "tropomi_chi_square_SWIR","tropomi_chi_square_NIR"]]

    # Defin variable to be predicted
    y = tropomi_gosat_pairs["delta_tropomi_gosat"]

    # Training set is 2018-2020
    train_index = tropomi_gosat_pairs[tropomi_gosat_pairs["tropomi_time"] < pd.to_datetime("2021-01-01")].index
    X_train = X.loc[train_index]
    y_train = y.loc[train_index]

    # For each algorithm give it the time specific in config.yml to train using 8 cores
    for estimator in ["lgbm", "rf", "xgboost"]:
        automl = AutoML()
        automl.fit(X_train, y_train, task="regression", metric="mse", time_budget=config["TimeFLAML"], n_jobs=8, estimator_list = [estimator], eval_method="cv", split_type="time", n_splits=10)
        # Save the model
        with open(os.path.join(config["StorageDir"], "processed", f"model_{estimator}.pkl"), "wb") as handle:
            pickle.dump(automl, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
if __name__ == "__main__":
    f_run_flaml()
