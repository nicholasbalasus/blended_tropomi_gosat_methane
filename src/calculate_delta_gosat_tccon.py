import os
import pandas as pd
import numpy as np
import glob
import yaml
from netCDF4 import Dataset
from multiprocessing import Pool
from geopy import distance
from scipy import interpolate
import pickle

# Load configuration file
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Get number of available cores
n_jobs = len(os.sched_getaffinity(0))

# TCCON station altitudes
tccon_station_altitudes = {"easttroutlake01": 0.50, "reunion01": 0.09, "hefei01": 0.03, "pasadena01": 0.23, "saga01": 0.01,
                           "edwards01": 0.70, "xianghe01": 0.04, "parkfalls01": 0.44, "rikubetsu01": 0.38, "jpl02": 0.39,
                           "nicosia01": 0.19, "burgos01": 0.04, "sodankyla01": 0.19, "paris01": 0.06, "eureka01": 0.61,
                           "orleans01": 0.13, "nyalesund01": 0.02, "lauder02": 0.37,"tsukuba02": 0.03, "izana01": 2.37,
                           "lamont01": 0.32, "lauder03": 0.37, "bremen01": 0.03, "garmisch01": 0.74, "karlsruhe01": 0.12} # [km a.s.l.]

# Co-location criteria
max_distance_diff = {x:500 for x in tccon_station_altitudes.keys() if "edwards01" not in x} # [km]
max_distance_diff["edwards01"] = 50 # [km]
max_elevation_diff = {x:250 for x in tccon_station_altitudes.keys()} # [m]
max_time_diff = {x:pd.Timedelta(hours=2) for x in tccon_station_altitudes.keys()} # [h]

# Define a function to create find the closest TCCON observation to a given GOSAT observation (if there is one).
def find_closest_tccon(gosat_row, tccon_df, tccon_station):
    
    # If GOSAT observation is not within the allowable altitude difference, return None.
    if np.abs(1000*tccon_station_altitudes[tccon_station] - gosat_row.iloc[0]["surface_altitude"]) > max_elevation_diff[tccon_station]:
        return None
    
    # Subset the TCCON DataFrame to include only observations within the allowable trime difference from the GOSAT observation.
    tccon_subset_temporal = tccon_df[(tccon_df["time"] >= gosat_row["time"].iloc[0] - max_time_diff[tccon_station]) & 
                                     (tccon_df["time"] <= gosat_row["time"].iloc[0] + max_time_diff[tccon_station])].copy()
    
    # If no TCCON observations are found within the time window, return None.
    if len(tccon_subset_temporal) == 0:
        return None
        
    # Calculate the distances between the GOSAT row and each TCCON row in the temporal subset.
    tccon_subset_temporal["gosat_distance"] = tccon_subset_temporal.apply(lambda x: distance.distance((gosat_row['latitude'].iloc[0], gosat_row['longitude'].iloc[0]), (x['latitude'], x['longitude'])).km, axis=1)
    
    # Subset the TCCON DataFrame to include only observations within the distance threshold.
    tccon_subset_temporal_spatial = tccon_subset_temporal[tccon_subset_temporal["gosat_distance"] < max_distance_diff[tccon_station]]
    
    # If no TCCON observations are found within the distance threshold, return None.
    if len(tccon_subset_temporal_spatial) == 0:
        return None
    
    # Identify the closest matching TCCON observation in time and concatenate it with the GOSAT row.
    else:
        closest_tccon = (tccon_subset_temporal_spatial["time"] - gosat_row["time"].iloc[0]).abs().idxmin()
        tccon_gosat_pair =  pd.concat([gosat_row.add_prefix("gosat_").reset_index(drop=True),
                                         tccon_subset_temporal_spatial.loc[closest_tccon:closest_tccon].add_prefix("tccon_").reset_index(drop=True)], axis=1)
        return tccon_gosat_pair
    
# Define a function to create a TCCON-GOSAT dataframe for a given TCCON station
def make_tccon_gosat_dataframe_for_a_tccon_station(tccon_file):

    with Dataset(tccon_file) as ds:

        # Make a dataframe of the TCCON data in the time range we are interested in.
        mask = (pd.to_datetime(ds["time"][:], unit="s") > pd.to_datetime("2018-04-30 00:00")) & (pd.to_datetime(ds["time"][:], unit="s") < pd.to_datetime("2021-12-31 23:59"))
        tccon_df = pd.DataFrame({
                            "time": pd.to_datetime(ds["time"][:][mask], unit="s"),
                            "latitude": ds["lat"][:][mask],
                            "longitude": ds["long"][:][mask],
                            "prior_wet_ch4": list(1e-9*ds["prior_ch4"][:][mask]), # [mol CH4/mol wet air]
                            "xch4": 1000*ds["xch4"][:][mask], # [ppb]
                            "prior_pressure": list(101325*ds["prior_pressure"][:][mask]), # [Pa]
                            "prior_wet_h2o": list(ds["prior_h2o"][:][mask]), # [mol H2O/mol wet air]
                            "ak_xch4": list(ds["ak_xch4"][:][mask])
        })
        
        # Hefei files incorrectly have longitude as 119.17 when they should be 117.17
        if (ds.long_name == "hefei01") and (len(np.unique(tccon_df["longitude"])) == 1) and (tccon_df["longitude"].iloc[0] == np.float32(119.17)):
            tccon_df["longitude"] -= 2.0

        # For each GOSAT observation, see if there is a TCCON observation within our distance, altitude, and time ranges (and add it to list_tccon_gosat_pairs_dfs if there is).
        list_tccon_gosat_pairs_dfs = []
        for i in gosat_df.index:
            list_tccon_gosat_pairs_dfs.append(find_closest_tccon(gosat_df.loc[i:i], tccon_df, ds.long_name))

        # Put all of the TCCON/GOSAT pairs into a single list (for this station).
        list_tccon_gosat_pairs_dfs = [x for x in list_tccon_gosat_pairs_dfs if x is not None]

        # If this list has dataframes, put them together and then add a column for ∆(GOSAT-TCCON)
        if len(list_tccon_gosat_pairs_dfs) > 0:
            station_tccon_gosat_pairs_df = pd.concat(list_tccon_gosat_pairs_dfs, ignore_index=True)

            # Calculate ∆(GOSAT-TCCON) using equation (A5)
            for i in station_tccon_gosat_pairs_df.index:

                tccon_gosat_pair = station_tccon_gosat_pairs_df.loc[i]

                gosat_mask = tccon_gosat_pair["gosat_ch4_profile_apriori"][::-1] != -9999.99 # when there are 19 pressure levels
                gosat_pressure_levels = 100*tccon_gosat_pair["gosat_pressure_levels"][::-1][gosat_mask] # [Pa]
                tccon_pressure_levels = tccon_gosat_pair["tccon_prior_pressure"][::-1] # [Pa]
                tccon_prior = (1e9*tccon_gosat_pair["tccon_prior_wet_ch4"] * (1 + tccon_gosat_pair["tccon_prior_wet_h2o"]/(1-tccon_gosat_pair["tccon_prior_wet_h2o"])))[::-1] # [ppb], have to "dry" the profile
                f_interp_tccon_prior_to_gosat_pressure_grid = interpolate.interp1d(tccon_pressure_levels, tccon_prior, bounds_error=False, fill_value="extrapolate")

                c_Gr = tccon_gosat_pair["gosat_xch4"] # [ppb]
                h_G = tccon_gosat_pair["gosat_pressure_weight"][::-1][gosat_mask] # [unitless]
                A_G = tccon_gosat_pair["gosat_xch4_averaging_kernel"][::-1][gosat_mask] # [unitless]
                x_Fa = f_interp_tccon_prior_to_gosat_pressure_grid(np.array(gosat_pressure_levels)) # [ppb]
                x_Ga = tccon_gosat_pair["gosat_ch4_profile_apriori"][::-1][gosat_mask] # [ppb]
                c_Fr = tccon_gosat_pair["tccon_xch4"] # [ppb]
                c_Fa = np.sum(h_G*x_Fa) # [ppb]
                x_Fr = x_Fa * (c_Fr/c_Fa) # [ppb]

                delta_gosat_tccon = (c_Gr+np.sum(h_G*(1-A_G)*(x_Fa-x_Ga))) - (np.sum(h_G*(x_Fa+(A_G*(x_Fr-x_Fa)))))
                station_tccon_gosat_pairs_df.loc[i,"delta_gosat_tccon"] = delta_gosat_tccon # [ppb]
                
            return ds.long_name, station_tccon_gosat_pairs_df
        
        # If no TCCON/GOSAT pairs were found, returning None
        else:
            return None, None

if __name__ == "__main__":
    
    # Get GOSAT dataframe from previous module
    gosat_df = pd.read_pickle(os.path.join(config['StorageDir'], 'processed', 'gosat.pkl'))

    # Get list of all TCCON files
    tccon_dir = os.path.join(config['StorageDir'], 'tccon')
    tccon_files = glob.glob(os.path.join(tccon_dir, "*.nc"))

    # Reduce the list to only those stations that intersect the time period of 30 Apr 2018 to 31 Dec 2021
    tccon_files = [tccon_file for tccon_file in tccon_files
               if pd.Interval(pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").min(),
                               pd.to_datetime(Dataset(tccon_file)["time"][:], unit="s").max())
                   .overlaps(pd.Interval(pd.to_datetime('2018-04-30 00:00'),
                                          pd.to_datetime('2021-12-31 23:59')))]

    # Get dataframe for each station in parallel
    with Pool(n_jobs) as p:
        l = p.map(make_tccon_gosat_dataframe_for_a_tccon_station, tccon_files)
        
    # Form dictionary of form {station name: pandas DataFrame of ∆(GOSAT-TCCON)}
    dict_delta_gosat_tccon = {x:y for x,y in l if x is not None}

    # Save dictionary
    with open(os.path.join(config['StorageDir'], 'processed', 'delta_gosat_tccon.pkl'), "wb") as f:
        pickle.dump(dict_delta_gosat_tccon, f)