import pandas as pd
from netCDF4 import Dataset
import numpy as np
    
# Function to turn a list of netCDF GOSAT files into one pandas dataframe
def get_gosat_df(gosat_files):

    gosat_dfs = []
    for gosat_file in gosat_files:

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

        df["time"] = pd.to_datetime(df["time"], unit="s")
        gosat_dfs.append(df)

    gosat_df = pd.concat(gosat_dfs, ignore_index=True)
    
    return gosat_df

# Function to turn one netCDF TROPOMI file into one pandas dataframe
def get_tropomi_df(tropomi_file):
    
    with Dataset(tropomi_file) as ds:
        mask = ds["PRODUCT/qa_value"][:] == 1.0
        tropomi_df = pd.DataFrame({
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
                           "relative_azimuth_angle": np.abs(180 - np.abs(ds["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_azimuth_angle"][:][mask] -
                                                                         ds["PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_azimuth_angle"][:][mask])),
                           "across_track_pixel_index": np.expand_dims(np.tile(ds["PRODUCT/ground_pixel"][:], (mask.shape[1],1)), axis=0)[mask],
                           "surface_classification": (ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification"][:][mask] & 0x03).astype(int),
                           "surface_altitude": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude"][:][mask],
                           "surface_altitude_precision": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude_precision"][:][mask],
                           "eastward_wind": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/eastward_wind"][:][mask],
                           "northward_wind": ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/northward_wind"][:][mask],
                           "xch4_apriori": np.sum(ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/methane_profile_apriori"][:][mask]/
                                                  np.expand_dims(np.sum(ds["PRODUCT/SUPPORT_DATA/INPUT_DATA/dry_air_subcolumns"][:][mask], axis=1),axis=1), axis=1)*1e9,
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
    tropomi_df["time"] = pd.to_datetime(tropomi_df["time"], format="%Y-%m-%dT%H:%M:%S.%fZ")
    
    return tropomi_df

# Function to turn one TCCON netCDF file into one pandas dataframe
def get_tccon_df(tccon_file):

    with Dataset(tccon_file) as ds:
        tccon_df = pd.DataFrame({
                            "time": pd.to_datetime(ds["time"][:], unit="s"),
                            "latitude": ds["lat"][:],
                            "longitude": ds["long"][:],
                            "prior_wet_ch4": list(1e-9*ds["prior_ch4"][:]), # [mol CH4/mol wet air]
                            "xch4": 1000*ds["xch4"][:], # [ppb]
                            "prior_pressure": list(101325*ds["prior_pressure"][:]), # [Pa]
                            "prior_wet_h2o": list(ds["prior_h2o"][:]), # [mol H2O/mol wet air]
                            "ak_xch4": list(ds["ak_xch4"][:])
        })
        
        # Hefei files incorrectly have longitude as 119.17 when they should be 117.17
        if (ds.long_name == "hefei01") and (len(np.unique(tccon_df["longitude"])) == 1) and (tccon_df["longitude"].iloc[0] == np.float32(119.17)):
            tccon_df["longitude"] -= 2.0

    return tccon_df

# Function to turn one netCDF Blended file (written in write_blended_files.py) into one pandas dataframe
def get_blended_df(blended_file):
    
    with Dataset(blended_file) as ds:
        blended_df = pd.DataFrame({
                           "latitude": ds["latitude"][:],
                           "longitude": ds["longitude"][:],
                           "time": pd.to_datetime(ds["time_utc"][:], format="%Y-%m-%dT%H:%M:%S.%fZ"),
                           "xch4_corrected": ds["methane_mixing_ratio_bias_corrected"][:],
                           "xch4_blended": ds["methane_mixing_ratio_blended"][:],
                           "pressure_interval": ds["pressure_interval"][:],
                           "surface_pressure": ds["surface_pressure"][:],
                           "dry_air_subcolumns": list(ds["dry_air_subcolumns"][:]),
                           "methane_profile_apriori": list(ds["methane_profile_apriori"][:]),
                           "column_averaging_kernel": list(ds["column_averaging_kernel"][:]),
                           "surface_altitude": ds["surface_altitude"][:],
                           "surface_albedo_SWIR": ds["surface_albedo_SWIR"][:],
                           "aerosol_size": ds["aerosol_size"][:]
                          })
    
    return blended_df