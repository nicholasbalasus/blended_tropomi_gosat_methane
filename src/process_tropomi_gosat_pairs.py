import os
import pandas as pd
import numpy as np
import scipy.interpolate
import yaml

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
paired_dir = config["StorageDir"] + "/pkl/paired/"
flaml_dir = config["StorageDir"] + "/pkl/flaml/"

def f_process_tropomi_gosat_pairs():

    '''
    Process all of the TROPOMI and GOSAT pairs
    This involes two things: putting them all into a single dataframe and calculating delta(TROPOMI-GOSAT)
    Using equations from https://climate.esa.int/media/documents/PVIR_GHG-CCIp_v3_reducedsize.pdf (last accessed 5 Dec 2022)

    Arguments
        none [none] : no arguments are explicity taken into the function, though it uses all of the files in paired_dir

    Returns
        none [none] : nothing is returned from the function, though it saves the single dataframe formed to flaml_dir
    '''

    # Concatenate all TROPOMI and GOSAT pairs together into a single dataframe
    all_pairs = pd.concat([pd.read_pickle(paired_dir+f) for f in os.listdir(paired_dir)], ignore_index=True)

    # Add a variable for the difference between GOSAT and TROPOMI accounting for prior & averaging kernel differences
    all_pairs["delta_tropomi_gosat"] = np.nan

    # Loop through each individual TROPOMI/GOSAT pair
    for idx in all_pairs.index:
        
        pair = all_pairs.loc[idx]
        
        # Get GOSAT data (and flip it around for interpolation)
        gosat_valid_idx = (pair.gosat_ch4_profile_apriori[::-1] != -9999.99) # when there are 19 pressure levels
        x_A_gosat_ppb = pair.gosat_ch4_profile_apriori[::-1][gosat_valid_idx] # profile of prior CH4 for GOSAT [ppb]
        p_levels_gosat = pair.gosat_pressure_levels[::-1][gosat_valid_idx] # profile of GOSAT pressures [hPa]

        # Get TROPOMI data
        p_level_tropomi = np.array([0.1+(x*pair.tropomi_dp) for x in range(13)]) # pressure levels that define boundaries of 12 pressure layers of TROPOMI [hPa]
        p_layer_tropomi = [(p_level_tropomi[i] + p_level_tropomi[i+1])/2 for i in range(12)] # pressure layer average (defined by pressure level bounds) [hPa]
        tropomi_dry_air_subcolumns = pair.tropomi_dry_air_subcolumns # profile of TROPOMI dry air in each layer [molec. cm-2]
        x_A_tropomi_ppb = 1e9*pair.tropomi_ch4_profile_apriori/tropomi_dry_air_subcolumns # profile of prior CH4 for TROPOMI [ppb]
        avg_k_tropomi = pair.tropomi_xch4_column_averaging_kernel # profile of TROPOMI averaging kernel [unitless]
        h_TROPOMI = pair.tropomi_dry_air_subcolumns/np.sum(pair.tropomi_dry_air_subcolumns) # pressure weight TROPOMI

        # Convert GOSAT prior CH4 profile to TROPOMI pressure grid
        f = scipy.interpolate.interp1d(p_levels_gosat,x_A_gosat_ppb,bounds_error=False,fill_value="extrapolate") # function that maps gosat pressure --> gosat prior CH4
        x_A_gosat_ppb_tropomi_p_layers = f(p_layer_tropomi) # gosat prior CH4 profile but mapped to TROPOMI pressure levels now [ppb]

        # Calculate TROPOMI XCH4 (column amount) that would have been retrieved if GOSAT's prior CH4 profile was used 
        c_tropomi_on_gosat_prior = pair.tropomi_xch4_corrected + np.sum(h_TROPOMI*(1-avg_k_tropomi)*(x_A_gosat_ppb_tropomi_p_layers-x_A_tropomi_ppb))

        # Calculate GOSAT XCH4 (column amount) that would have been retrived if GOSAT had the vertical sensitivity of TROPOMI
        c_gosat_A_ppb = np.sum(h_TROPOMI*x_A_gosat_ppb_tropomi_p_layers) # GOSAT prior XCH4 (column amount) [ppb]
        c_gosat_ppb = pair.gosat_xch4 - config["GlobalOffsetGOSAT"] # adjust all GOSAT values to have a global mean bias of 0 relative to TCCON GGG2020
        x_gosat_ppb = x_A_gosat_ppb_tropomi_p_layers * (c_gosat_ppb/c_gosat_A_ppb) # because GOSAT doesn't report CH4 at vertical levels, scale the prior CH4 profile
        c_gosat_on_tropomi_avg_k = np.sum(h_TROPOMI*(x_A_gosat_ppb_tropomi_p_layers + (avg_k_tropomi*(x_gosat_ppb-x_A_gosat_ppb_tropomi_p_layers))))
        
        # Calculate delta(GOSAT_XCH4 - TROPOMI_XCH4)
        all_pairs.loc[idx,"delta_tropomi_gosat"] =  c_tropomi_on_gosat_prior - c_gosat_on_tropomi_avg_k
        
    all_pairs.to_pickle(flaml_dir + "filtered_tropomi_gosat_pairs.pkl")

if __name__ == "__main__":
    f_process_tropomi_gosat_pairs()