import os
import pandas as pd
import numpy as np
import glob
import yaml
from netCDF4 import Dataset
from multiprocessing import Pool
import sys
import pickle

# Load configuration file
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Get number of available cores
n_jobs = len(os.sched_getaffinity(0))

def convert_batch_of_gosat_nc_files_to_dataframe(gosat_file_batch):
    # Create an empty list to hold dataframes for each file
    list_batch_dfs = []

    # Loop through each file in the batch
    for gosat_file in gosat_file_batch:
        # Read in the file and filter data with quality_flag of 0 into a dataframe
        with Dataset(gosat_file) as ds:
            mask = ds["xch4_quality_flag"][:] == 0
            df = pd.DataFrame({
                "latitude": ds["latitude"][:][mask],
                "longitude": ds["longitude"][:][mask],
                "xch4": ds["xch4"][:][mask],
                "model_xco2_range": ds["model_xco2_range"][:][mask],
                "surface_altitude": ds["surface_altitude"][:][mask],
                "time": ds["time"][:][mask],
                "retr_flag": ds["retr_flag"][:][mask],
                "ch4_profile_apriori": list(ds["ch4_profile_apriori"][:][mask]),
                "xch4_averaging_kernel": list(ds["xch4_averaging_kernel"][:][mask]),
                "pressure_levels": list(ds["pressure_levels"][:][mask]),
                "pressure_weight": list(ds["pressure_weight"][:][mask])
            })

        # Convert column of seconds since 1970-01-01 00:00:00 to datetime
        df["time"] = pd.to_datetime(df["time"], unit="s")

        # Add dataframe for the current file to the list
        list_batch_dfs.append(df)

    # Concatenate all dataframes in the list into a single dataframe
    batch_df = pd.concat(list_batch_dfs, ignore_index=True)

    return batch_df


def convert_batch_of_tropomi_nc_files_to_dataframe(tropomi_file_batch):
    # Create an empty list to hold dataframes for each file
    list_batch_dfs = []
    
    # Loop through each file in the batch
    for tropomi_file in tropomi_file_batch:
        # Read in the file and filter data with quality_flag of 0 into a dataframe
        with Dataset(tropomi_file) as ds:
            mask = ds["PRODUCT/qa_value"][:] == 1.0
            df = pd.DataFrame({
                               # non-predictor variables
                               "latitude": ds["PRODUCT/latitude"][:][mask],
                               "longitude": ds["PRODUCT/longitude"][:][mask],
                               "time": np.expand_dims(np.tile(ds["PRODUCT/time_utc"][:][0,:], (mask.shape[2],1)).T, axis=0)[mask],
                               "latitude_bounds": list(ds["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds"][:][mask]),
                               "longitude_bounds": list(ds["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds"][:][mask]),
                               "xch4": ds["PRODUCT/methane_mixing_ratio"][:][mask],
                               "xch4_corrected": ds["PRODUCT/methane_mixing_ratio_bias_corrected"][:][mask],
                               "pressure_interval": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/pressure_interval"][:][mask],
                               "surface_pressure": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure"][:][mask],
                               "dry_air_subcolumns": list(ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/dry_air_subcolumns"][:][mask]),
                               "methane_profile_apriori": list(ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/methane_profile_apriori"][:][mask]),
                               "column_averaging_kernel": list(ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel"][:][mask]),
                               # predictor variables
                               "solar_zenith_angle": ds["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle"][:][mask],
                               "relative_azimuth_angle": np.abs(180 - np.abs(ds["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_azimuth_angle"][:][mask] - ds["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_azimuth_angle"][:][mask])),
                               "across_track_pixel_index": np.expand_dims(np.tile(ds["PRODUCT/ground_pixel"][:], (mask.shape[1],1)), axis=0)[mask],
                               "surface_classification": (ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification"][:][mask] & 0x03).astype(int),
                               "surface_altitude": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude"][:][mask],
                               "surface_altitude_precision": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude_precision"][:][mask],
                               "eastward_wind": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/eastward_wind"][:][mask],
                               "northward_wind": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/northward_wind"][:][mask],
                               "xch4_apriori": np.sum(ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/methane_profile_apriori"][:][mask]/np.expand_dims(np.sum(ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/dry_air_subcolumns"][:][mask], axis=1),axis=1), axis=1)*1e9,
                               "reflectance_cirrus_VIIRS_SWIR": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/reflectance_cirrus_VIIRS_SWIR"][:][mask],
                               "xch4_precision": ds["PRODUCT/methane_mixing_ratio_precision"][:][mask],
                               "fluorescence": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/fluorescence"][:][mask],
                               "co_column": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/carbonmonoxide_total_column"][:][mask],
                               "co_column_precision": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/carbonmonoxide_total_column_precision"][:][mask],
                               "h2o_column": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_total_column"][:][mask],
                               "h2o_column_precision": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_total_column_precision"][:][mask],
                               "aerosol_size": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_size"][:][mask],
                               "aerosol_size_precision": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_size_precision"][:][mask],
                               "aerosol_height": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_mid_altitude"][:][mask],
                               "aerosol_height_precision": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_mid_altitude_precision"][:][mask],
                               "aerosol_column": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_number_column"][:][mask],
                               "aerosol_column_precision": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_number_column_precision"][:][mask],
                               "surface_albedo_SWIR": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_SWIR"][:][mask],
                               "surface_albedo_SWIR_precision": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_SWIR_precision"][:][mask],
                               "surface_albedo_NIR": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_NIR"][:][mask],
                               "surface_albedo_NIR_precision": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_NIR_precision"][:][mask],
                               "aerosol_optical_thickness_SWIR": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_optical_thickness_SWIR"][:][mask],
                               "aerosol_optical_thickness_NIR": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_optical_thickness_NIR"][:][mask],
                               "chi_square_SWIR": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/chi_square_SWIR"][:][mask],
                               "chi_square_NIR": ds["PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/chi_square_NIR"][:][mask]
                              })

        # Convert column of strings to datetime
        df["time"] = pd.to_datetime(df["time"], format="%Y-%m-%dT%H:%M:%S.%fZ")
        
        # Add dataframe for the current file to the list
        list_batch_dfs.append(df)

    # Concatenate all dataframes in the list into a single dataframe
    batch_df = pd.concat(list_batch_dfs, ignore_index=True)
    
    return batch_df

if __name__ == "__main__":

    # TROPOMI processing is running across an array of SLURM jobs
    SLURM_ARRAY_TASK_ID = int(sys.argv[1])
    SLURM_ARRAY_TASK_COUNT = int(sys.argv[2])
    print(SLURM_ARRAY_TASK_ID, SLURM_ARRAY_TASK_COUNT, flush=True)

    # Process GOSAT Data only in the first task (this only takes a few minutes anyway)
    if SLURM_ARRAY_TASK_ID == 0:
        gosat_dir = os.path.join(config['StorageDir'], 'gosat')
        gosat_files = glob.glob(os.path.join(gosat_dir, "*.nc"))
        gosat = convert_batch_of_gosat_nc_files_to_dataframe(gosat_files)
        gosat.to_pickle(os.path.join(config['StorageDir'], 'processed', 'gosat.pkl'))

    # Process TROPOMI Data for the given task #
    tropomi_dir = os.path.join(config['StorageDir'], 'tropomi')
    tropomi_files = glob.glob(os.path.join(tropomi_dir, "*.nc"))
    tropomi_files.sort()
    tropomi_file_batches_by_task = np.array_split(tropomi_files, SLURM_ARRAY_TASK_COUNT)
    tropomi_file_batch_this_task = tropomi_file_batches_by_task[SLURM_ARRAY_TASK_ID]
    tropomi_file_batches_this_task = np.array_split(tropomi_file_batch_this_task, n_jobs)
    with Pool(n_jobs) as p:
        all_dfs = p.map(convert_batch_of_tropomi_nc_files_to_dataframe, tropomi_file_batches_this_task)
    tropomi_this_task = pd.concat(all_dfs, ignore_index=True)
    tropomi_this_task.to_pickle(os.path.join(config['StorageDir'], 'tmp', f'chunk{str(SLURM_ARRAY_TASK_ID).zfill(3)}of{SLURM_ARRAY_TASK_COUNT-1}.pkl'), protocol=pickle.HIGHEST_PROTOCOL)