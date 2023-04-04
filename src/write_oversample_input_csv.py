import geopandas
import pandas as pd
import os
import numpy as np
import yaml
from netCDF4 import Dataset
from shapely.geometry import Polygon
import multiprocessing
import glob

from utilities import get_blended_df

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Write a csv file that can be used as input to tools/oversample for all satellite pixels in this time range (start_dt,end_dt) and this domain (extent)
def csv_oversample(filename,start_dt,end_dt,extent,filter_bool):
    
    print(filename, flush=True)

    # Dict to hold all satellite data that is in the valid time range and intersects with the oversampling region
    satellite = {} 
    for var in ["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","methane_mixing_ratio_bias_corrected","methane_mixing_ratio_blended","methane_mixing_ratio_precision","surface_albedo_SWIR","surface_albedo_SWIR_precision","aerosol_size","aerosol_size_precision"]:
        satellite[var] = np.array([])

    # Polygon for entire oversampling region
    lon_lower,lon_upper,lat_lower,lat_upper = extent
    oversampling_region = Polygon(zip([lon_lower,lon_lower,lon_upper,lon_upper],[lat_lower,lat_upper,lat_upper,lat_lower]))

    # Get list of all blended files that overlap our time period
    blended_files = [f for f in glob.glob(os.path.join(config["StorageDir"], "blended", "*.nc")) if pd.Interval(pd.to_datetime(f.split("_")[12]), pd.to_datetime(f.split("_")[13])).overlaps(pd.Interval(start_dt,end_dt))]

    # Compile all of the valid satellite observations into the satellite{} dict
    for blended_file in blended_files:
        with Dataset(blended_file) as ds:
            latc = ds["latitude_bounds"][:]
            lonc = ds["longitude_bounds"][:]
            polygons = geopandas.GeoSeries([Polygon(zip(lonc[i],latc[i])) for i in range(len(latc))]) # form polygons from the pixel corners
            valid_idx = polygons.intersects(oversampling_region) # satellite pixels that intersect with the oversampling region
            
            # Remove pixels that cross antimeridian because they mess up regridding - this is not a perfect solution
            for i in range(len(valid_idx)):
                valid_idx[i] &= (np.max(lonc[i]) - np.min(lonc[i])) < 180
            
            # Remove coastal pixels if necessary (landflag == 3 AND landflag == 2 with SWIR chi2 > 20)
            if filter_bool:
                tropomi_surface_classification = (ds["surface_classification"][:] & 0x03).astype(int)
                valid_idx &= ~((tropomi_surface_classification == 3) | ((tropomi_surface_classification == 2) & (ds["chi_square_SWIR"][:] > 20)))
            
            # As we loop through each file, append the satellite retrievals that intersect with the oversampling region
            for var in ["methane_mixing_ratio_bias_corrected","methane_mixing_ratio_blended","methane_mixing_ratio_precision",
                        "surface_albedo_SWIR","surface_albedo_SWIR_precision","aerosol_size","aerosol_size_precision"]:
                satellite[var] = np.append(satellite[var],ds[var][:][valid_idx])
            for i in [0,1,2,3]:
                satellite[f"lat{i}"] = np.append(satellite[f"lat{i}"],ds["latitude_bounds"][:,i][valid_idx])
                satellite[f"lon{i}"] = np.append(satellite[f"lon{i}"],ds["longitude_bounds"][:,i][valid_idx])
            satellite["lat"] = np.append(satellite["lat"],ds["latitude"][:][valid_idx])
            satellite["lon"] = np.append(satellite["lon"],ds["longitude"][:][valid_idx])

    # Form the CSV files
    satellite = pd.DataFrame.from_dict(satellite)
    csv_operational = satellite[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","methane_mixing_ratio_bias_corrected","methane_mixing_ratio_precision"]]
    csv_blended = satellite[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","methane_mixing_ratio_blended","methane_mixing_ratio_precision"]]
    csv_sa = satellite[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","surface_albedo_SWIR","surface_albedo_SWIR_precision"]]
    csv_aerosol = satellite[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","aerosol_size","aerosol_size_precision"]]

    # Save the CSV files
    csv_operational.to_csv(os.path.join(config["StorageDir"], "oversample", "input", filename+"_Operational.csv"), sep=',', float_format='%.6f', header=False, index=True)
    csv_blended.to_csv(os.path.join(config["StorageDir"], "oversample", "input", filename+"_Blended.csv"), sep=',', float_format='%.6f', header=False, index=True)
    csv_sa.to_csv(os.path.join(config["StorageDir"], "oversample", "input", filename+"_SWIR_Surface_Albedo.csv"), sep=',', float_format='%.6f', header=False, index=True)
    csv_aerosol.to_csv(oos.path.join(config["StorageDir"], "oversample", "input", filename+"_Aerosol_Size.csv"), sep=',', float_format='%.6f', header=False, index=True)
    
if __name__ == "__main__":

    inputs = [("TROPOMI_NorthAfrica_2021", pd.to_datetime("2021-01-01 00:00:00"), pd.to_datetime("2022-01-01 00:00:00"), (-20.,65.,0.,40.), False),
              ("TROPOMI_NorthAfrica_Filter_2021", pd.to_datetime("2021-01-01 00:00:00"), pd.to_datetime("2022-01-01 00:00:00"), (-20.,65.,0.,40.), True)]
    
    with multiprocessing.Pool() as pool:
        results = pool.starmap(csv_oversample, inputs)
        pool.close()
        pool.join()