from netCDF4 import Dataset
import numpy as np
import pandas as pd
import datetime
import multiprocessing
import scipy.interpolate
import yaml
import os

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
tropomi_processed_dir = config["StorageDir"] + "/pkl/tropomi/"
tccon_raw_dir = config["StorageDir"] + "/tccon/"
tccon_tropomi_processed_dir = config["StorageDir"] + "/pkl/tccon/tropomi/"
blended_tropomi_gosat_dir = config["StorageDir"] + "/pkl/blended_tropomi_gosat/"

# TCCON station names and their associated files in tccon_raw_dir
tccon_station_filenames = {
    "Bremen": "br20090106_20210624.public.qc.nc",
    "Burgos": "bu20170303_20210430.public.qc.nc",
    "Pasadena": "ci20120920_20220901.public.qc.nc",
    "Edwards": "df20130720_20220901.public.qc.nc",
    "East_Trout_Lake": "et20161003_20220813.public.qc.nc",
    "Eureka": "eu20100724_20200707.public.qc.nc",
    "Garmisch-Partenkirchen": "gm20070716_20211018.public.qc.nc",
    "Hefei": "hf20160108_20201231.public.qc.nc",
    "Izana": "iz20140102_20220730.public.qc.nc",
    "Jet_Propulsion_Laboratory_2": "jf20110519_20180514.public.qc.nc",
    "Saga": "js20110728_20220330.public.qc.nc",
    "Karlsruhe": "ka20140115_20211222.public.qc.nc",
    "Lauder_2": "ll20130102_20180930.public.qc.nc",
    "Lauder_3": "lr20181002_20220630.public.qc.nc",
    "Nicosia": "ni20190903_20210601.public.qc.nc",
    "Ny-Ålesund": "ny20050314_20210919.public.qc.nc",
    "Lamont": "oc20080706_20221031.public.qc.nc",
    "Orléans": "or20090829_20210303.public.qc.nc",
    "Park_Falls": "pa20040526_20220831.public.qc.nc",
    "Paris": "pr20140923_20210616.public.qc.nc",
    "Réunion_Island": "ra20150301_20200718.public.qc.nc",
    "Rikubetsu": "rj20140624_20210630.public.qc.nc",
    "Sodankylä": "so20090516_20220614.public.qc.nc",
    "Tsukuba_2": "tk20140328_20210628.public.qc.nc",
    "Xianghe": "xh20180614_20211130.public.qc.nc"
}

# Altitude of TCCON stations in km (https://tccon-wiki.caltech.edu/Main/TCCONSites)
tccon_station_altitude = {
    "Bremen": 0.027,
    "Burgos": 0.035,
    "Pasadena": 0.23,
    "Edwards": 0.699,
    "East_Trout_Lake": 0.5108,
    "Eureka": 0.61,
    "Garmisch-Partenkirchen": 0.74,
    "Hefei": 0.029,
    "Izana": 2.37,
    "Jet_Propulsion_Laboratory_2": 0.39,
    "Saga": 0.007,
    "Karlsruhe": 0.116,
    "Lauder_2": 0.37,
    "Lauder_3": 0.37,
    "Nicosia": 0.185,
    "Ny-Ålesund": 0.02,
    "Lamont": 0.32,
    "Orléans": 0.13,
    "Park_Falls": 0.44,
    "Paris": 0.06,
    "Réunion_Island": 0.087,
    "Rikubetsu": 0.38,
    "Sodankylä": 0.188,
    "Tsukuba_2": 0.03,
    "Xianghe": 0.036
}

def f_pair_tropomi_tccon(args):

    '''
    For the given TCCON station, determine a set of valid TROPOMI and TCCON pairs for comparison
    Metric used is distance of 100 km spatially, 250 m vertically, and 1 hour temporally
    Calculate delta(TROPOMI-TCCON) accounting for differences in priors & averaging kernels
    Equations from https://climate.esa.int/media/documents/PVIR_GHG-CCIp_v3_reducedsize.pdf (last accessed 5 Dec 2022)

    Arguments
        station [str] : name of TCCON station to pair with TROPOMI observations
        blended [bool] : if True, use blended product; if False, use SRON v19

    Returns
        none [none] : function doesn't return anything but saves a dataframe of TCCON/TROPOMI pairs to tccon_tropomi_processed_dir
    '''
    
    # Using a list as a single input because of multiprocessing
    station,blended = args

    # Read in each TROPOMI dataframe
    # Subset TROPOMI to that within 100 km spatially and 250 m of altitude of the TCCON station
    tropomi_files = os.listdir(tropomi_processed_dir)
    tropomi_files.sort()
    tropomi_subset = pd.DataFrame()
    for file in tropomi_files:
        tropomi = pd.concat([pd.read_pickle(tropomi_processed_dir + file),pd.read_pickle(blended_tropomi_gosat_dir + file)], axis=1) # regular TROPOMI + blended product var
        if station == "Edwards":
            valid_idx = (tropomi[station+"_distance_km"] <= 50) & (np.abs(tropomi["surface_altitude"] - tccon_station_altitude[station]*1000) <= 250)
        else:
            valid_idx = (tropomi[station+"_distance_km"] <= 100) & (np.abs(tropomi["surface_altitude"] - tccon_station_altitude[station]*1000) <= 250)
        tropomi_subset_tmp = tropomi[valid_idx].reset_index(drop=True)
        if len(tropomi_subset) == 0: # i.e., if this is the first file
            tropomi_subset = tropomi_subset_tmp
        else:
            tropomi_subset = pd.concat([tropomi_subset,tropomi_subset_tmp], ignore_index=True)

    # Dataframe that will contain our pairs of TROPOMI/TCCON
    # Include variables for TROPOMI, TCCON, and delta so we can see if our analysis is robust to accounting for differences in priors/avg. kernels
    all_pairs = pd.DataFrame()
    all_pairs.index = tropomi_subset.index
    for var in ["time","tropomi_xch4","tccon_xch4","delta_tropomi_tccon","tropomi_swir_surface_albedo","tropomi_aerosol_size"]:
        all_pairs[var] = np.nan
        
    # Read in netCDF of TCCON data and form a dataframe from it
    tccon_data = Dataset(tccon_raw_dir+tccon_station_filenames[station])
    tccon_time = tccon_data["time"][:] # seconds since 1970-01-01 00:00:00 UTC
    tccon_xch4 = tccon_data["xch4"][:]*1000 # ppm to ppb
    tccon_datetime = np.empty_like(tccon_time,dtype="object")
    for idx,time in enumerate(tccon_time):
        tccon_datetime[idx] = datetime.datetime(1970,1,1,0,0,0) + datetime.timedelta(seconds=time)
    tccon_f_wetH2O = list(tccon_data["prior_h2o"][:])
    tccon_f_wetCH4 = list(tccon_data["prior_ch4"][:]) # these profiles are wet (i.e., molec. CH4/(molec. H2O + molec. air))
    tccon_apriori_pressure = list(tccon_data["prior_pressure"][:]) # [atm]
    tccon = pd.DataFrame({"datetime":tccon_datetime,"xch4":tccon_xch4,"f_wetH2O":tccon_f_wetH2O,\
                            "f_wetCH4":tccon_f_wetCH4,"apriori_pressure":tccon_apriori_pressure})
    
    for idx in tropomi_subset.index:

        # Get all TROPOMI data for this specific observation
        tropomi_time = tropomi_subset.loc[idx,"time"]
        p_level_tropomi = np.array([0.1+(x*tropomi_subset.loc[idx,"dp"]) for x in range(13)]) # pressure levels that define boundaries of 12 pressure layers of TROPOMI [hPa]
        p_layer_tropomi = [(p_level_tropomi[i] + p_level_tropomi[i+1])/2 for i in range(12)] # pressure layer average (defined by pressure level bounds) [hPa]
        tropomi_dry_air_subcolumns = tropomi_subset.loc[idx,"dry_air_subcolumns"] # profile of TROPOMI dry air in each layer [molec. cm-2]
        x_A_tropomi_ppb = 1e9*tropomi_subset.loc[idx,"ch4_profile_apriori"]/tropomi_dry_air_subcolumns # profile of prior CH4 for TROPOMI [ppb]
        avg_k_tropomi = tropomi_subset.loc[idx,"xch4_column_averaging_kernel"] # profile of TROPOMI averaging kernel [unitless]
        h_TROPOMI = tropomi_subset.loc[idx,"dry_air_subcolumns"]/np.sum(tropomi_subset.loc[idx,"dry_air_subcolumns"]) # pressure weight TROPOMI
        if blended == True:
            tropomi_c = tropomi_subset.loc[idx,"xch4_blended_tropomi_gosat"]
        elif blended == False:
            tropomi_c = tropomi_subset.loc[idx,"xch4_corrected"]
        tropomi_swir_surface_albedo = tropomi_subset.loc[idx,"swir_surface_albedo"]
        tropomi_aerosol_size = tropomi_subset.loc[idx,"aerosol_size"]
        
        # Subset TCCON to be within 1 hour of this specific TROPOMI observation
        tccon_subset = tccon[(tccon["datetime"] >= (tropomi_time - datetime.timedelta(hours=1))) & (tccon["datetime"] < (tropomi_time + datetime.timedelta(hours=1)))]
        
        # If there are any TCCON measurements within 1 hour of this specific TROPOMI observation
        if len(tccon_subset) != 0:
            
            # Choose TCCON measurement closest in time
            closest_tccon = (tccon_subset["datetime"] - tropomi_time).map(lambda t: np.abs(t.total_seconds())).argmin()
            tccon_subset = tccon_subset.iloc[closest_tccon]
            
            # Dry the vertical profiles as in https://tccon-wiki.caltech.edu/Main/AuxiliaryDataGGG2020 (last accessed 5 Dec 2022)
            f_wetH2O = tccon_subset["f_wetH2O"][::-1]
            f_dryH2O = (1/f_wetH2O - 1)**-1
            f_wetCH4 = tccon_subset["f_wetCH4"][::-1]
            tccon_apriori_ch4 = (f_wetCH4/(1-f_dryH2O)) # tccon prior ch4 profile [ppb]
            tccon_apriori_pressure = (tccon_subset["apriori_pressure"]*1013.25)[::-1] # atm to hPa

            # Put the TCCON ch4 prior profile on the TROPOMI pressure grid
            f = scipy.interpolate.interp1d(tccon_apriori_pressure,tccon_apriori_ch4,bounds_error=False,fill_value="extrapolate") # maps tccon pressure to tccon ch4 prior
            x_A_tccon_ppb_tropomi_p_layers = f(p_layer_tropomi)

            # TROPOMI XCH4 (column amount) that would have been retrieved using the TCCON CH4 prior profile
            c_tropomi_on_tccon_prior = tropomi_c + np.sum(h_TROPOMI*(1-avg_k_tropomi)*(x_A_tccon_ppb_tropomi_p_layers-x_A_tropomi_ppb))

            # TCCON doesn't report retrieved profile, so we scale the prior profile
            tccon_c_A = np.sum(h_TROPOMI*x_A_tccon_ppb_tropomi_p_layers)
            tccon_x_ppb = x_A_tccon_ppb_tropomi_p_layers*(tccon_subset["xch4"]/tccon_c_A)

            # TCCON XCH4 (column amount) that would have been retrieved if TCCON had the vertical sensitivity of TROPOMI
            c_tccon_on_tropomi_avg_k = np.sum(h_TROPOMI*(x_A_tccon_ppb_tropomi_p_layers+(avg_k_tropomi*(tccon_x_ppb-x_A_tccon_ppb_tropomi_p_layers))))
            
            all_pairs.loc[idx,"time"] = tropomi_time
            all_pairs.loc[idx,"tropomi_xch4"] = tropomi_c
            all_pairs.loc[idx,"tccon_xch4"] = tccon_subset["xch4"]
            all_pairs.loc[idx,"delta_tropomi_tccon"] = c_tropomi_on_tccon_prior - c_tccon_on_tropomi_avg_k
            all_pairs.loc[idx,"tropomi_swir_surface_albedo"] = tropomi_swir_surface_albedo
            all_pairs.loc[idx,"tropomi_aerosol_size"] = tropomi_aerosol_size
        
    # Drop any row that has a NaN
    all_pairs = all_pairs[~all_pairs.isnull().any(axis=1)].reset_index(drop=True)
    
    # Save the DataFrame
    if blended == True:
        all_pairs.to_pickle(tccon_tropomi_processed_dir + station + f"_blended.pkl")
    elif blended == False:
        all_pairs.to_pickle(tccon_tropomi_processed_dir + station + f"_SRONv19.pkl")
    
if __name__ == "__main__":
    
    # Across 50 cores, run f_pair_tropomi_tccon for each of the 25 TCCON stations
    # Run once with blended product (True) and once with original TROPOMI data (False)
    stations = list(tccon_station_filenames.keys())
    args = []
    for station in stations:
        args.append([station,True])
        args.append([station,False])
    pool = multiprocessing.Pool(50)
    pool.map(f_pair_tropomi_tccon, args)
    pool.close()