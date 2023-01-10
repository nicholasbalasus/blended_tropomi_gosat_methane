import geopandas
import pandas as pd
import os
import numpy as np
import yaml
from netCDF4 import Dataset
from shapely.geometry import Polygon
import multiprocessing

# Config file with settings
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Key directories
blended_tropomi_gosat_nc_dir = config["StorageDir"] + "/blended_tropomi_gosat/"
oversample_input_dir = config["StorageDir"] + "/oversample/input/"

def csv_oversample(args):

    '''
    Function that writes a csv file that can be used in the tools/oversample code for oversampling

    Arguments
        filename [str] : file prefix that describes the area and time that the csv file will contain
        start_orbit [int] : beginning TROPOMI orbit to include
        end_orbit [int] : end TROPOMI orbit to include
        lon_lower [float] : lowest longitude of oversampling region
        lon_upper [float] : highest longitude of oversampling region
        lat_lower [float] : lowest latitude of oversampling region
        lat_upper [float] : highest latitude of oversampling region

    Returns
        none [none] : nothing explicity returned but four csv files are saved to oversample_input_dir
    '''

    # Using a list as a single input because of multiprocessing
    filename,start_orbit,end_orbit,lon_lower,lon_upper,lat_lower,lat_upper = args
    
    # Do not run if files already exists
    if os.path.exists(oversample_input_dir+filename+"_Blended.csv"):
        return

    print(filename, flush=True)

    # Dict to hold all tropomi dat that is between start_orbit and end_orbit and intersects with the oversampling region
    tropomi = {} 
    for var in ["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","xch4_corrected","xch4_blended","xch4_precision","swir_surface_albedo","swir_surface_albedo_precision","surface_altitude","surface_altitude_stdv","aerosol_size","aerosol_size_precision"]:
        tropomi[var] = np.array([])

    # Polygon for entire oversampling region
    oversampling_region = Polygon(zip([lon_lower,lon_lower,lon_upper,lon_upper],[lat_lower,lat_upper,lat_upper,lat_lower]))

    filelist = os.listdir(blended_tropomi_gosat_nc_dir)
    filelist.sort()

    # Compile all of the valid tropomi observations into the tropomi{} dict
    for file in filelist:
        orbit_number = int(file[16:21])
        if (orbit_number >= start_orbit) & (orbit_number <= end_orbit):
            ds = Dataset(blended_tropomi_gosat_nc_dir+file)
            latc = ds["instrument/latitude_corners"][:]
            lonc = ds["instrument/longitude_corners"][:]
            polygons = geopandas.GeoSeries([Polygon(zip(lonc[i],latc[i])) for i in range(len(latc))]) # form polygons from the pixel corners
            valid_idx = polygons.intersects(oversampling_region) # tropomi pixels that intersect with the oversampling region
            
            # Remove pixels that cross antimeridian because they mess up regridding - this is not a perfect solution
            cross_antimeridian = np.ones_like(ds["target_product/xch4_blended"], dtype="bool")
            for i in range(len(cross_antimeridian)):
                cross_antimeridian[i] = not((len(np.where(lonc[i] > 179)[0]) > 0) & (len(np.where(lonc[i] < -179)[0]) > 0))
            valid_idx &= cross_antimeridian
            
            # Remove coastal pixels if necessary
            if "No_Coastal" in filename:
                valid_idx &= (ds["meteo/landflag"][:] != 3)
            elif "No_Mixed" in filename:
                valid_idx &= (ds["meteo/landflag"][:] != 2)
                valid_idx &= (ds["meteo/landflag"][:] != 3)
            
            # As we loop through each file, append the tropomi retrievals that intersect with the oversampling region
            for var in ["xch4_corrected","xch4_blended","xch4_precision"]:
                tropomi[var] = np.append(tropomi[var],ds[f"target_product/{var}"][:][valid_idx])
            for var in ["swir_surface_albedo","swir_surface_albedo_precision"]:
                tropomi[var] = np.append(tropomi[var],ds[f"side_product/{var.replace('swir_','')}"][:,1][valid_idx])
            for var in ["surface_altitude","surface_altitude_stdv"]:
                tropomi[var] = np.append(tropomi[var],ds[f"meteo/{var}"][:][valid_idx])
            for var in ["aerosol_size","aerosol_size_precision"]:
                tropomi[var] = np.append(tropomi[var],ds[f"side_product/{var}"][:][valid_idx])
            for i in [0,1,2,3]:
                tropomi[f"lat{i}"] = np.append(tropomi[f"lat{i}"],ds[f"instrument/latitude_corners"][:,i][valid_idx])
                tropomi[f"lon{i}"] = np.append(tropomi[f"lon{i}"],ds[f"instrument/longitude_corners"][:,i][valid_idx])
            tropomi["lat"] = np.append(tropomi["lat"],ds["instrument/latitude_center"][:][valid_idx])
            tropomi["lon"] = np.append(tropomi["lon"],ds["instrument/longitude_center"][:][valid_idx])


    # Form the CSV files
    tropomi = pd.DataFrame.from_dict(tropomi)
    csv_sron = tropomi[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","xch4_corrected","xch4_precision"]]
    csv_blended = tropomi[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","xch4_blended","xch4_precision"]]
    csv_sa = tropomi[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","swir_surface_albedo","swir_surface_albedo_precision"]]
    csv_h = tropomi[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","surface_altitude","surface_altitude_stdv"]]
    csv_aerosol = tropomi[["lat0","lat1","lat2","lat3","lat","lon0","lon1","lon2","lon3","lon","aerosol_size","aerosol_size_precision"]]

    # Save the CSV files
    csv_sron.to_csv(oversample_input_dir+filename+"_SRON.csv", sep=',', float_format='%.6f', header=False, index=True)
    csv_blended.to_csv(oversample_input_dir+filename+"_Blended.csv", float_format='%.6f', header=False, index=True)
    csv_sa.to_csv(oversample_input_dir+filename+"_SWIR_Surface_albedo.csv", float_format='%.6f', header=False, index=True)
    csv_h.to_csv(oversample_input_dir+filename+"_Surface_altitude.csv", float_format='%.6f', header=False, index=True)
    csv_aerosol.to_csv(oversample_input_dir+filename+"_Aerosol_size.csv", float_format='%.6f', header=False, index=True)
    
if __name__ == "__main__":
    
    args = [["TROPOMI_NorthAfrica_2021",16679,21856,-20,65,0,40],\
            ["TROPOMI_NorthAfrica_No_Coastal_2021",16679,21856,-20,65,0,40],\
            ["TROPOMI_NorthAfrica_No_Mixed_2021",16679,21856,-20,65,0,40],\
            ["TROPOMI_NorthAmerica_2021",16679,21856,-135,-60,15,60],\
            ["TROPOMI_NorthAmerica_No_Coastal_2021",16679,21856,-135,-60,15,60],\
            ["TROPOMI_NorthAmerica_No_Mixed_2021",16679,21856,-135,-60,15,60]]
    
    pool = multiprocessing.Pool(6)
    pool.map(csv_oversample, args)
    pool.close()