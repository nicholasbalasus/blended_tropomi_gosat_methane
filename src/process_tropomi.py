import os
import pandas as pd
import numpy as np
import datetime
from netCDF4 import Dataset
from geopy import distance
import multiprocessing
import yaml
import sys
from shapely.geometry import Polygon
import geopandas

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
tropomi_raw_dir = config["StorageDir"] + "/tropomi/"
tropomi_processed_dir = config["StorageDir"] + "/pkl/tropomi/"

# Location of TCCON stations (latitude, longitude)
tccon = {
    "Bremen": (53.1,8.85),
    "Burgos": (18.5325,120.6496),
    "Pasadena": (34.136,-118.127),
    "Edwards": (34.959917, -117.881069),
    "East_Trout_Lake": (54.353738,-104.986667),
    "Eureka": (80.05,-86.42),
    "Garmisch-Partenkirchen": (47.476,11.063),
    "Hefei": (31.9, 117.17),
    "Izana": (28.3,-16.5),
    "Jet_Propulsion_Laboratory_2": (34.202,-118.175),
    "Saga": (33.24104,130.28818),
    "Karlsruhe": (49.1,8.438),
    "Lauder_2": (-45.038,169.684),
    "Lauder_3": (-45.038,169.684),
    "Nicosia": (35.141,33.381),
    "Ny-Ålesund": (78.923,11.923),
    "Lamont": (36.604,-97.486),
    "Orléans": (47.97,2.113),
    "Park_Falls": (45.945,-90.273),
    "Paris": (48.846,2.356),
    "Réunion_Island": (-20.901,55.485),
    "Rikubetsu": (43.46,143.77),
    "Sodankylä": (67.3668,26.631),
    "Tsukuba_2": (36.0513,140.1215),
    "Xianghe": (39.75,116.96)
}

def f_process_tropomi(file):
    
    '''
    Function that converts the netCDF TROPOMI orbit files into pandas datafamres
    Reduces to only the variables we are interested in and subsets to quality assurance value of 1

    Arguments
        file [str] : name of TROPOMI level 2 file in tropomi_raw_dir

    Returns
        none [none] : subset_data dataframe is saved to tropomi_processed_dir as a pickle file, but nothing is returned from the function
    '''

    # Read in the netCDF file
    ds = Dataset(tropomi_raw_dir+file)
    
    # dict of path to variables --> name of variable for TROPOMI files
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
                 "/meteo/landflag":"landflag",\
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
    subset_data = {}
    for idx,v in enumerate(variables.keys()):

        subset_data[variables[v]] = ds[v][:]
        
        # Start with every pixel being valid, then remove those that are masked
        if idx == 0:
            validmask = np.ones_like(subset_data[variables[v]], dtype="bool")
        validmask &= ~subset_data[variables[v]].mask
    
    # Get variables with NIR ans SWIR components (dimensions of nobs,nwin)
    for var in ["spectral_shift","surface_albedo","surface_albedo_precision","aerosol_optical_thickness","reflectance_max"]:
        subset_data["nir_"+var] = ds["/side_product/"+var][:,0]
        subset_data["swir_"+var] = ds["/side_product/"+var][:,1]
        validmask &= (~subset_data["nir_"+var].mask & ~subset_data["swir_"+var].mask)
        
    for var in ["chi_squared_band","number_of_spectral_points_in_retrieval","signal_to_noise_ratio"]:
        subset_data["nir_"+var] = ds["/diagnostics/"+var][:,0]
        subset_data["swir_"+var] = ds["/diagnostics/"+var][:,1]
        validmask &= (~subset_data["nir_"+var].mask & ~subset_data["swir_"+var].mask)
        
    # Get variables with FOV components (dimensions of nobs,nfov)
    for idx in [0,1,2,3]:
        subset_data["fov_"+str(idx)] = ds["/meteo/cloud_fraction"][:,idx]
        validmask &= (~subset_data["fov_"+str(idx)].mask)
    
    # Remove pixels without qa_value == 1
    validmask &= (subset_data["qa_value"] == 1)
        
    # Use the validmask we have built to actually subset the data
    for key in subset_data.keys():
        assert np.shape(subset_data[key]) == np.shape(validmask)
        # Convert to regular array, there should be no masked values left
        subset_data[key] = subset_data[key][validmask].filled(fill_value=-9990999)
        assert ~np.any(subset_data[key] == -9990999)
        
    # Get the time
    time_as_separate_array = ds["/instrument/time"][:][validmask]
    time_to_fill_subset = np.empty_like(subset_data["latitude"], dtype="object")
    for idx,_ in enumerate(time_to_fill_subset):
        time_to_fill_subset[idx] = datetime.datetime(*time_as_separate_array[idx])
    subset_data["time"] = time_to_fill_subset
    
    # Add variables with vertical dimension (ch4_profile_apriori, altitude_levels, dry_air_subcolumns, xch4_column_averaging_kernel)
    profile_variables = {"/target_product/ch4_profile_apriori":"ch4_profile_apriori",\
                         "/target_product/xch4_column_averaging_kernel":"xch4_column_averaging_kernel",\
                         "/meteo/altitude_levels":"altitude_levels",\
                         "/meteo/dry_air_subcolumns":"dry_air_subcolumns"}
    
    for v in profile_variables.keys():
        subset_data[profile_variables[v]] = ds[v][:][validmask].filled(fill_value=-9990999)
        assert ~np.any(subset_data[profile_variables[v]] == -9990999)
        subset_data[profile_variables[v]] = list(subset_data[profile_variables[v]]) # so it can go in a DataFrame
    
    # Get the distance to TCCON sites
    for site in tccon.keys():
        subset_data[f"{site}_distance_km"] = np.empty_like(subset_data["latitude"])
        for idx,_ in enumerate(subset_data[f"{site}_distance_km"]):
            subset_data[f"{site}_distance_km"][idx] = distance.distance(tccon[site],(subset_data["latitude"][idx],subset_data["longitude"][idx])).km
    
    # All vlaues of qa are 1, so delete them
    assert(all(subset_data["qa_value"] == 1))
    del(subset_data["qa_value"])
    
    # Close the file
    ds.close()
    
    # Convert subset_Data to dataframe, then save it to a pickle file
    pd.DataFrame(subset_data).to_pickle(tropomi_processed_dir+file.replace(".nc",".pkl"))

if __name__ == "__main__":
    
    # Open list of files to be processed
    with open("tropomi.txt") as f:
        filelist = [line.strip() for line in f]
    
    # split list of files into number of slurm array tasks
    # sys.argv[1] = SLURM_ARRAY_TASK_ID
    # sys.argv[2] = SLURM_ARRAY_TASK_COUNT
    pool = multiprocessing.Pool(config["Cores"])
    filelist = np.array_split(filelist,int(sys.argv[2]))[int(sys.argv[1])]
    filelist_not_done_yet = []
    for file in filelist:
        if (file.replace(".nc", ".pkl")) not in os.listdir(tropomi_processed_dir):
            filelist_not_done_yet.append(file)
    pool.map(f_process_tropomi, filelist_not_done_yet)
    pool.close()