from netCDF4 import Dataset
import pandas as pd
import os
import numpy as np
import glob
import yaml
import pickle

from utilities import get_tropomi_df

with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

def predict_blended_xch4(tropomi_file, model):

    df = get_tropomi_df(tropomi_file)
    # Get rid of the non-predictor variables
    df = df.drop(["latitude","longitude","time","latitude_bounds","xch4","xch4_corrected","pressure_interval","surface_pressure","dry_air_subcolumns","methane_profile_apriori","column_averaging_kernel"], axis=1) 
    df = df.add_prefix("tropomi_")
    with Dataset(tropomi_file) as ds:
        mask = ds["PRODUCT/qa_value"][:] == 1.0
        assert len(df) == np.sum(mask)
        blended = ds["PRODUCT/methane_mixing_ratio_bias_corrected"][:][mask] - (config["a"]*model.predict(df) + config["b"])
    
    return blended

def f_write_blended_files(src_file):
    
    # new file will have the same name but with BLND as the acronym and the creation time changed
    dst_file = src_file.split("/")[-1].replace("RPRO","BLND")
    dst_file = os.path.join(config["StorageDir"], "blended", dst_file[:dst_file.rfind("_")+1], pd.Timestamp.now().strftime("%Y%m%dT%H%M%S"), ".nc")

    with Dataset(src_file) as src, Dataset(dst_file, "w") as dst:
    
        # Set global attributes
        dst.setncatts({
            "Title": "Blended TROPOMI+GOSAT Methane Product",
            "Contact": "Nicholas Balasus (nicholasbalasus@g.harvard.edu)"
        })

        # Create a mask to select data
        mask = src["PRODUCT/qa_value"][:] == 1.0
        
        # Create nobs dimension (dimension across which qa_value == 1.0) and layer dimension
        dst.createDimension("nobs", np.sum(mask))
        dst.createDimension("layer", len(src["PRODUCT/layer"]))
        dst.createDimension("corner", len(src["PRODUCT/corner"]))
        
        # Copy over variables and their attributes from the PRODUCT group
        vars_to_keep_in_PRODUCT = ["qa_value","latitude","longitude","methane_mixing_ratio","methane_mixing_ratio_precision",
                                "methane_mixing_ratio_bias_corrected"]
        for var in vars_to_keep_in_PRODUCT:
            dst.createVariable(var, src["PRODUCT/"+var].datatype, ('nobs'))
            dst[var].setncatts(src["PRODUCT/"+var].__dict__)
            dst[var][:] = src["PRODUCT/"+var][:][mask]
            
        # Time originally has only scanline dimensions, so need to expand it
        dst.createVariable("time_utc", np.str_, ('nobs'), fill_value='')
        dst["time_utc"].setncatts({'long_name': 'Time of observation as ISO 8601 date-time string'})
        dst["time_utc"][:] = np.array(np.expand_dims(np.tile(src["PRODUCT/time_utc"][:][0,:], (mask.shape[2],1)).T, axis=0), dtype=np.str_)[mask]
        
        # Copy over variables and their attributes from the PRODUCT/SUPPORT_DATA/GEOLOCATIONS group
        vars_to_keep_in_PRODUCT_SUPPORT_DATA_GEOLOCATIONS = ["latitude_bounds","longitude_bounds"]
        for var in vars_to_keep_in_PRODUCT_SUPPORT_DATA_GEOLOCATIONS:
            dst.createVariable(var, src["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/"+var].datatype, ('nobs', 'corner'))
            dst[var].setncatts(src["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/"+var].__dict__)
            dst[var][:] = src["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/"+var][:][mask]
            
        # Copy over variables and their attributes from the PRODUCT/SUPPORT_DATA/DETAILED_RESULTS group
        vars_to_keep_in_PRODUCT_SUPPORT_DATA_DETAILED_RESULTS = ["chi_square_SWIR","surface_albedo_SWIR","surface_albedo_NIR"]
        for var in vars_to_keep_in_PRODUCT_SUPPORT_DATA_DETAILED_RESULTS:
            dst.createVariable(var, src["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/"+var].datatype, ('nobs'))
            dst[var].setncatts(src["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/"+var].__dict__)
            dst[var][:] = src["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/"+var][:][mask]
            
        vars_to_keep_in_PRODUCT_SUPPORT_DATA_DETAILED_RESULTS = ["column_averaging_kernel"]
        for var in vars_to_keep_in_PRODUCT_SUPPORT_DATA_DETAILED_RESULTS:
            dst.createVariable(var, src["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/"+var].datatype, ('nobs','layer'))
            dst[var].setncatts(src["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/"+var].__dict__)
            dst[var][:] = src["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/"+var][:][mask]
        
        # Copy over variables and their attributes from the PRODUCT/SUPPORT_DATA/INPUT_DATA group
        vars_to_keep_in_PRODUCT_SUPPORT_DATA_INPUT_DATA = ["surface_altitude","surface_classification","surface_pressure",
                                                        "pressure_interval","reflectance_cirrus_VIIRS_SWIR"]
        for var in vars_to_keep_in_PRODUCT_SUPPORT_DATA_INPUT_DATA:
            dst.createVariable(var, src["PRODUCT/SUPPORT_DATA/INPUT_DATA/"+var].datatype, ('nobs'))
            dst[var].setncatts(src["PRODUCT/SUPPORT_DATA/INPUT_DATA/"+var].__dict__)
            dst[var][:] = src["PRODUCT/SUPPORT_DATA/INPUT_DATA/"+var][:][mask]
            
        vars_to_keep_in_PRODUCT_SUPPORT_DATA_INPUT_DATA = ["methane_profile_apriori","dry_air_subcolumns"]
        for var in vars_to_keep_in_PRODUCT_SUPPORT_DATA_INPUT_DATA:
            dst.createVariable(var, src["PRODUCT/SUPPORT_DATA/INPUT_DATA/"+var].datatype, ('nobs', 'layer'))
            dst[var].setncatts(src["PRODUCT/SUPPORT_DATA/INPUT_DATA/"+var].__dict__)
            dst[var][:] = src["PRODUCT/SUPPORT_DATA/INPUT_DATA/"+var][:][mask]
            
        # Edit some of the attributes to be consistent with the new file format
        for var in dst.variables:
            if var == "time_utc":
                dst[var].setncattr("coordinates", "longitude latitude")
                
            for attribute in dst[var].ncattrs():
                if attribute == "coordinates":
                    dst[var].setncattr(attribute, "longitude latitude")
                if attribute == "bounds" and var == "latitude":
                    dst[var].setncattr(attribute, "latitude_bounds")
                if attribute == "bounds" and var == "longitude":
                    dst[var].setncattr(attribute, "longitude_bounds")
                if attribute == "comment" and var == "methane_mixing_ratio_bias_corrected":
                    del dst[var].comment
                    del dst[var].ancillary_variables
                if attribute == "ancillary_variables" and var == "methane_mixing_ratio":
                    del dst[var].ancillary_variables
                    
        # Add blended xch4
        dst.createVariable("methane_mixing_ratio_blended", src["PRODUCT/methane_mixing_ratio"].dtype, ('nobs'))
        dst["methane_mixing_ratio_blended"].setncatts(src["PRODUCT/methane_mixing_ratio"].__dict__)
        dst["methane_mixing_ratio_blended"].setncattr("comment", "produced as described in Balasus et al. (2023)")
        with open(os.path.join(config["StorageDir"], "processed", f"model_{config['Model']}.pkl"), "rb") as handle:
            model = pickle.load(handle)
        dst["methane_mixing_ratio_blended"][:] = predict_blended_xch4(src, model)

if __name__ == "__main__":
    
    src_files = glob.glob(os.path.join(config["StorageDir"], "tropomi", "*.nc"))
    num_processes = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(f_write_blended_files, src_files)
        pool.close()
        pool.join()