import shap
import pickle
import yaml
import pandas as pd
import datetime

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
flaml_dir = config["StorageDir"] + "/pkl/flaml/"

def f_shap_explainer():

    '''
    Using the training data, determines a SHAP explainer and SHAP values

    Arugments
        none [none] : no explicit arguments are taken, but the trained model and filtered TROPOMI+GOSAT pairs are used as inputs

    Returns
        none [none] : nothing is explicity returned by the SHAP explainer and values are saved to flaml_dir
    '''    
    
    # Read in the predictive model
    with open(config["StorageDir"] + f"/pkl/flaml/model_{config['Model']}.pkl", "rb") as handle:
        model = pickle.load(handle)

    # Read in the pairs and subset to the training data (2018-2020)
    pairs = pd.read_pickle(flaml_dir+"filtered_tropomi_gosat_pairs.pkl").sort_values("tropomi_time").reset_index(drop=True)
    pairs["tropomi_landflag"] = pairs["tropomi_landflag"].astype("category")
    train_index = pairs[pairs["tropomi_time"] < datetime.datetime(2021,1,1,0,0,0)].index
    X = pairs[["tropomi_sza","tropomi_raa","tropomi_surface_altitude","tropomi_surface_altitude_stdv","tropomi_u10","tropomi_v10",\
               "tropomi_cirrus_reflectance","tropomi_xch4_precision","tropomi_xch4_apriori","tropomi_fluorescence","tropomi_co_column",\
               "tropomi_co_column_precision","tropomi_h2o_column","tropomi_h2o_column_precision","tropomi_aerosol_size","tropomi_aerosol_size_precision","tropomi_aerosol_column",\
               "tropomi_aerosol_column_precision","tropomi_aerosol_altitude","tropomi_aerosol_altitude_precision","tropomi_nir_surface_albedo","tropomi_swir_surface_albedo",\
               "tropomi_nir_surface_albedo_precision","tropomi_swir_surface_albedo_precision","tropomi_nir_aerosol_optical_thickness","tropomi_swir_aerosol_optical_thickness",\
               "tropomi_nir_chi_squared_band","tropomi_swir_chi_squared_band","tropomi_ground_pixel","tropomi_landflag"]]
    X_train = X.loc[train_index]

    # Form and save explainer
    explainer = shap.TreeExplainer(model.model.estimator)
    with open(flaml_dir+"shap_explainer.pkl", "wb") as handle:
        pickle.dump(explainer, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Calculate and save SHAP values for train data
    shap_values_train = explainer(X_train)
    with open(flaml_dir+"shap_values_train.pkl", "wb") as handle:
        pickle.dump(shap_values_train, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
if __name__ == "__main__":
    f_shap_explainer()