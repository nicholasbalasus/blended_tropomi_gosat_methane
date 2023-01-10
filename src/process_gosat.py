import os
import pandas as pd
import numpy as np
import datetime
from netCDF4 import Dataset
from geopy import distance
import yaml

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
gosat_raw_dir = config["StorageDir"] + "/gosat/"
gosat_processed_dir = config["StorageDir"] + "/pkl/gosat/"

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

def f_process_gosat(file):
    
    '''
    Function that converts the daily netCDF GOSAT files into dictionaries of numpy arrays
    Reduces to only the variables we are interested in and subsets to quality flag of 0

    Arguments
        file [str] : name of GOSAT level 2 file in gosat_raw_dir

    Returns
        subset_data [dict] : dictionary of the subset of GOSAT data w/ lat,lon, xch4, model_xco2_range, model_xco2, time, surface_altitude, tccon_distances
    '''

    # Read in the netCDF file
    ds = Dataset(gosat_raw_dir+file)
    
    # Get the variables we are interested in
    subset_data = {}
    variables = ["latitude","longitude","xch4","xch4_quality_flag","model_xco2_range","model_xco2","surface_altitude","time"]

    for idx,v in enumerate(variables):
        subset_data[v] = ds[v][:]
        
        # Start with every pixel being valid, then remove those that are masked
        if idx == 0:
            validmask = np.ones_like(subset_data[v], dtype="bool")
        validmask &= ~subset_data[v].mask
    
    # Removed variables without quality flag == 0
    validmask &= (subset_data["xch4_quality_flag"] == 0)
    
    # Use the validmask we have built to actually subset the data
    for key in subset_data.keys():
        assert np.shape(subset_data[key]) == np.shape(validmask)
        # Convert to regular array, there should be no masked values left
        subset_data[key] = subset_data[key][validmask].filled(fill_value=-9990999)
        assert ~np.any(subset_data[key] == -9990999)
        
    # Convert time from seconds since 1970-01-01 00:00:00 to datetime (UTC)
    reference_time = datetime.datetime(1970,1,1,0,0,0)
    time_to_fill_subset = np.empty(np.shape(subset_data["time"]), dtype="object")
    for idx,_ in enumerate(time_to_fill_subset):
        time_to_fill_subset[idx] = reference_time + datetime.timedelta(seconds=subset_data["time"][idx])
    subset_data["time"] = time_to_fill_subset

    # Test to make sure all times are on the day that the file represents
    begin_day = datetime.datetime.strptime(file.replace("UoL-GHG-L2-CH4-GOSAT-OCPR-","").replace("-fv9.0.nc",""), "%Y%m%d")
    end_day = begin_day + datetime.timedelta(hours=25)
    assert(all((subset_data["time"] >= begin_day) & (subset_data["time"] < end_day)))
    
    # Add variables with a vertical dimension (i.e., profiles)
    for v in ["ch4_profile_apriori", "xch4_averaging_kernel", "pressure_levels", "pressure_weight"]:
        subset_data[v] = ds[v][:][validmask].filled(fill_value=-9990999)
        assert ~np.any(subset_data[v] == -9990999)
        subset_data[v] = list(subset_data[v]) # so it can go in a DataFrame
            
    # Close the file
    ds.close()
    
    # Calculate distance to each TCCON site
    for site in tccon.keys():
        subset_data[f"{site}_distance_km"] = np.empty_like(subset_data["latitude"])
        for idx,_ in enumerate(subset_data[f"{site}_distance_km"]):
            subset_data[f"{site}_distance_km"][idx] = distance.distance(tccon[site],(subset_data["latitude"][idx],subset_data["longitude"][idx])).km
        
    # Remove quality flag since they are all 0
    del(subset_data["xch4_quality_flag"])
    
    return subset_data

if __name__ == "__main__":
    
    # List of all GOSAT files
    filelist = os.listdir(gosat_raw_dir)
    filelist.sort() 

    # Convert every GOSAT netCDF file into a dataframe, then put them all together into a single dataframe
    gosat_data = pd.DataFrame()
    for file in filelist:
        gosat_one_file = pd.DataFrame(f_process_gosat(file))
        gosat_data = pd.concat([gosat_data,gosat_one_file], ignore_index=True)

    # Save a single dataframe that has all of the valid GOSAT data
    gosat_data.to_pickle(gosat_processed_dir+"gosat_data.pkl")