from netCDF4 import Dataset
import numpy as np
import pandas as pd
import os
import multiprocessing
import yaml
import sys

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
flaml_dir = config["StorageDir"] + "/pkl/flaml/"
tropomi_raw_dir =  config["StorageDir"] + "/tropomi/"
blended_nc_dir = config["StorageDir"] + "/blended_tropomi_gosat/"
tropomi_processed_dir = config["StorageDir"] + "/pkl/tropomi/"
blended_dir = config["StorageDir"] + "/pkl/blended_tropomi_gosat/"

def f_find_valid_n_obs_idx(ds):

    '''
    Function that returns the same array of valid_idx that was used in process_tropomi.py
    Criteria are that all variables are available and that qa_value == 1

    Arguments
        ds [Dataset] : netCDF4 dataset of TROPOMI data

    Returns
        valid_idx [np.array] : array of True and False for valid and invalid data points in the TROPOMI dataset
    '''
    subset_data = {}
    variables = {"/instrument/latitude_center":"latitude",\
                    "/instrument/longitude_center":"longitude",\
                    "/instrument/solar_zenith_angle":"sza",\
                    "/instrument/viewing_zenith_angle":"vza",\
                    "/instrument/relative_azimuth_angle":"raa",\
                    "/instrument/glintflag":"glintflag",\
                    "/instrument/ground_pixel":"ground_pixel",\
                    "/meteo/surface_altitude":"surface_altitude",\
                    "/meteo/surface_altitude_stdv":"surface_altitude_stdv",\
                    "/meteo/dp":"dp",\
                    "/meteo/u10":"u10",\
                    "/meteo/v10":"v10",\
                    "/meteo/fluorescence_apriori":"fluorescence_apriori",\
                    "/meteo/weak_h2o_column":"weak_h2o_column",\
                    "/meteo/strong_h2o_column":"strong_h2o_column",\
                    "/meteo/weak_ch4_column":"weak_ch4_column",\
                    "/meteo/strong_ch4_column":"strong_ch4_column",\
                    "/meteo/cirrus_reflectance":"cirrus_reflectance",\
                    "/meteo/stdv_h2o_ratio":"stdv_h2o_ratio",\
                    "/meteo/stdv_ch4_ratio":"stdv_ch4_ratio",\
                    "/meteo/surface_pressure":"surface_pressure",\
                    "/target_product/xch4":"xch4",\
                    "/target_product/xch4_precision":"xch4_precision",\
                    "/target_product/xch4_apriori":"xch4_apriori",\
                    "/target_product/xch4_corrected":"xch4_corrected",\
                    "/side_product/fluorescence":"fluorescence",\
                    "/side_product/co_column":"co_column",\
                    "/side_product/co_column_precision":"co_column_precision",\
                    "/side_product/h2o_column":"h2o_column",\
                    "/side_product/h2o_column_precision":"h2o_column_precision",\
                    "/side_product/aerosol_size":"aerosol_size",\
                    "/side_product/aerosol_size_precision":"aerosol_size_precision",\
                    "/side_product/aerosol_column":"aerosol_column",\
                    "/side_product/aerosol_column_precision":"aerosol_column_precision",\
                    "/side_product/aerosol_altitude":"aerosol_altitude",\
                    "/side_product/aerosol_altitude_precision":"aerosol_altitude_precision",\
                    "/diagnostics/chi_squared":"chi_squared",\
                    "/diagnostics/degrees_of_freedom":"degrees_of_freedom",\
                    "/diagnostics/degrees_of_freedom_ch4":"degrees_of_freedom_ch4",\
                    "/diagnostics/degrees_of_freedom_aerosol":"degrees_of_freedom_aerosol",\
                    "/diagnostics/rms":"rms",\
                    "/diagnostics/qa_value":"qa_value"}

    # Get the variables we are interested in (those with only a dimension of nobs)
    for idx,v in enumerate(variables.keys()):

        subset_data[variables[v]] = ds[v][:]
        
        # Start with every pixel being valid, then remove those that are masked
        if idx == 0:
            validmask = np.ones_like(subset_data[variables[v]], dtype="bool")
        validmask &= ~subset_data[variables[v]].mask

    # Get variables with NIR ans SWIR components (nobs,nwin)
    for var in ["spectral_shift","surface_albedo","surface_albedo_precision","aerosol_optical_thickness","reflectance_max"]:
        subset_data["nir_"+var] = ds["/side_product/"+var][:,0]
        subset_data["swir_"+var] = ds["/side_product/"+var][:,1]
        validmask &= (~subset_data["nir_"+var].mask & ~subset_data["swir_"+var].mask)
        
    for var in ["chi_squared_band","number_of_spectral_points_in_retrieval","signal_to_noise_ratio"]:
        subset_data["nir_"+var] = ds["/diagnostics/"+var][:,0]
        subset_data["swir_"+var] = ds["/diagnostics/"+var][:,1]
        validmask &= (~subset_data["nir_"+var].mask & ~subset_data["swir_"+var].mask)
        
    # Get variables with FOV components (nobs,nfov)
    for idx in [0,1,2,3]:
        subset_data["fov_"+str(idx)] = ds["/meteo/cloud_fraction"][:,idx]
        validmask &= (~subset_data["fov_"+str(idx)].mask)

    # Remove pixels without qa_value == 1
    validmask &= (subset_data["qa_value"] == 1)

    return validmask

def f_write_nc(file):

    '''
    Function that creates a netCDF file that mimics the SRON file exactly, but with the blended TROPOMI+GOSAT data added
    Only includes data with a qa_value == 1 since this is the only data the machine learning is applied to

    Arguments
        file [str] : name of netCDF file in tropomi_raw_dir that is to be replicated

    Returns
        none [none] : nothing explicity returned, but a netCDF file is saved to blended_nc_dir
    '''

    file = file.replace(".pkl",".nc")

    # If file exists, remove it
    if os.path.exists(blended_nc_dir+file):
        os.remove(blended_nc_dir+file)

    # Original TROPOMI data is src, blended product is destination
    with Dataset(tropomi_raw_dir+file) as src, Dataset(blended_nc_dir+file, "w") as dst:

        # Filter the tropomi data the same way we did in process_tropomi.py
        valid_nobs_idx = f_find_valid_n_obs_idx(src)

        # Add attributes to distinguish from SRONv19
        dst.setncatts({"Title":"Blended TROPOMI+GOSAT Methane Product",\
                       "Contact":"Nicholas Balasus (nicholasbalasus@g.harvard.edu)"})

        # Replicate dimensions of src
        for name, dimension in src.dimensions.items():
            if name == "nobs":
                if len(valid_nobs_idx) == np.sum(valid_nobs_idx.mask):
                    dst.createDimension(name, 0)
                else:
                    dst.createDimension(name, np.sum(valid_nobs_idx))
            else:
                dst.createDimension(name, len(dimension))

        # Replicate groups and variables of src
        for groupname, _ in src.groups.items():
            dst.createGroup(groupname)
            for varname, variable in src[groupname].variables.items():
                dst[groupname].createVariable(varname, variable.datatype, variable.dimensions)
                dst[groupname][varname].setncatts(src[groupname][varname].__dict__)

                # This is okay because 'nobs' is always the 1st dim if it exists
                if 'nobs' in variable.dimensions:
                    dst[groupname][varname][:] = src[groupname][varname][:][valid_nobs_idx]
                else:
                    dst[groupname][varname][:] = src[groupname][varname][:]
                    
                if varname == "xch4_corrected":
                    dst[groupname].createVariable("xch4_blended", variable.datatype, variable.dimensions)
                    dst[groupname]["xch4_blended"].setncatts(src[groupname][varname].__dict__)
        
        # Add a variable for the blended product
        if len(pd.read_pickle(blended_dir+file.replace(".nc",".pkl"))["xch4_blended_tropomi_gosat"]) != 0:
            dst["target_product/xch4_blended"][:] = pd.read_pickle(blended_dir+file.replace(".nc",".pkl"))["xch4_blended_tropomi_gosat"]

        # Checks
        df = pd.read_pickle(tropomi_processed_dir+file.replace(".nc",".pkl"))
        df2 = pd.read_pickle(blended_dir+file.replace(".nc",".pkl"))
        if len(df) != 0:
            assert len(df) == np.sum(valid_nobs_idx)
            assert np.array_equal(dst["target_product/xch4_corrected"][:],df["xch4_corrected"])
            assert np.array_equal(dst["meteo/dp"][:],df["dp"])
            assert len(df) == len(df2)

if __name__ == "__main__":
    
    # Open list of files to be processed
    with open("tropomi.txt") as f:
        filelist = [line.strip() for line in f]
            
    # split list of files into number of slurm array tasks
    # sys.argv[1] = SLURM_ARRAY_TASK_ID
    # sys.argv[2] = SLURM_ARRAY_TASK_COUNT
    pool = multiprocessing.Pool(config["Cores"])
    filelist = np.array_split(filelist,int(sys.argv[2]))[int(sys.argv[1])]
    pool.map(f_write_nc, filelist)
    pool.close()