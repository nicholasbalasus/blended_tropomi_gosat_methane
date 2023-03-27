import shap
import pickle
import yaml
import pandas as pd
import os

# Load configuration file
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

def f_shap_explainer():
    
    # Read in the predictive model
    with open(os.path.join(config["StorageDir"], "processed", f"model_{config['Model']}.pkl"), "rb") as handle:
        model = pickle.load(handle)

    # Read in the pairs and subset to the training data (2018-2020)
    tropomi_gosat_pairs = pd.read_pickle(os.path.join(config["StorageDir"], "processed", "tropomi_gosat_pairs.pkl"))
    X = tropomi_gosat_pairs[["tropomi_solar_zenith_angle","tropomi_relative_azimuth_angle","tropomi_across_track_pixel_index",
                         "tropomi_surface_classification","tropomi_surface_altitude","tropomi_surface_altitude_precision",
                         "tropomi_eastward_wind","tropomi_northward_wind","tropomi_xch4_apriori","tropomi_reflectance_cirrus_VIIRS_SWIR",
                         "tropomi_xch4_precision","tropomi_fluorescence","tropomi_co_column","tropomi_co_column_precision",
                         "tropomi_h2o_column","tropomi_h2o_column_precision","tropomi_aerosol_size","tropomi_aerosol_size_precision",
                         "tropomi_aerosol_height","tropomi_aerosol_height_precision","tropomi_aerosol_column","tropomi_aerosol_column_precision",
                         "tropomi_surface_albedo_SWIR","tropomi_surface_albedo_SWIR_precision","tropomi_surface_albedo_NIR",
                         "tropomi_surface_albedo_NIR_precision","tropomi_aerosol_optical_thickness_SWIR","tropomi_aerosol_optical_thickness_NIR",
                         "tropomi_chi_square_SWIR","tropomi_chi_square_NIR"]]
    train_index = tropomi_gosat_pairs[tropomi_gosat_pairs["tropomi_time"] < pd.to_datetime("2021-01-01")].index
    X_train = X.loc[train_index]

    # Form and save explainer
    explainer = shap.TreeExplainer(model.model.estimator)
    with open(os.path.join(config["StorageDir"], "shap_explainer.pkl"), "wb") as handle:
        pickle.dump(explainer, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Calculate and save SHAP values for train data
    shap_values_train = explainer(X_train)
    with open(os.path.join(config["StorageDir"], "shap_values_train.pkl"), "wb") as handle:
        pickle.dump(shap_values_train, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
if __name__ == "__main__":
    f_shap_explainer()