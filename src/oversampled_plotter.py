import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import yaml
import pandas as pd
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import sys
import multiprocessing

sys.path.insert(1,"./tools/")
from tol_colors import tol_cmap,tol_cset

with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

def f_plotter(args):

    basename,res,extent = args
    
    blended = pd.read_fwf(f"{config['StorageDir']}/oversample/output/{basename}_Blended_oversampled{res}.csv", header=None, widths=[6,6,12,12,15,6], names=["row","col","lat","lon","xch4_blended","n"], usecols=["lat","lon","xch4_blended","n"])
    sron = pd.read_fwf(f"{config['StorageDir']}/oversample/output/{basename}_SRON_oversampled{res}.csv", header=None, widths=[6,6,12,12,15,6], names=["row","col","lat","lon","xch4_corrected","n"], usecols=["lat","lon","xch4_corrected","n"])
    albedo = pd.read_fwf(f"{config['StorageDir']}/oversample/output/{basename}_SWIR_Surface_albedo_oversampled{res}.csv", header=None, widths=[6,6,12,12,15,6], names=["row","col","lat","lon","swir_surface_albedo","n"], usecols=["lat","lon","swir_surface_albedo","n"])
    altitude = pd.read_fwf(f"{config['StorageDir']}/oversample/output/{basename}_Surface_altitude_oversampled{res}.csv", header=None, widths=[6,6,12,12,15,6], names=["row","col","lat","lon","surface_altitude","n"], usecols=["lat","lon","surface_altitude","n"])
    aerosol = pd.read_fwf(f"{config['StorageDir']}/oversample/output/{basename}_Aerosol_size_oversampled{res}.csv", header=None, widths=[6,6,12,12,15,6], names=["row","col","lat","lon","aerosol_size","n"], usecols=["lat","lon","aerosol_size","n"])

    patches = []

    for i in blended.index:
        lon = blended.loc[i,"lon"]
        lat = blended.loc[i,"lat"]
        xmin = lon - res/2
        xmax = lon + res/2
        ymin = lat - res/2
        ymax = lat + res/2
        patches.append(Polygon(np.array(((xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin)))))

    sron_collection = PatchCollection(patches, array=sron["xch4_corrected"], cmap=tol_cmap("rainbow_PuBr"))
    sron_collection.set_clim(vmin=1850,vmax=1925)

    blended_collection = PatchCollection(patches, array=blended["xch4_blended"], cmap=tol_cmap("rainbow_PuBr"))
    blended_collection.set_clim(vmin=1850,vmax=1925)

    albedo_collection = PatchCollection(patches, array=albedo["swir_surface_albedo"], cmap=tol_cmap("rainbow_PuBr"))
    albedo_collection.set_clim(vmin=0,vmax=0.6)
    
    altitude_collection = PatchCollection(patches, array=altitude["surface_altitude"], cmap=tol_cmap("rainbow_PuBr"))
    altitude_collection.set_clim(vmin=0,vmax=2500)
    
    aerosol_collection = PatchCollection(patches, array=aerosol["aerosol_size"], cmap=tol_cmap("rainbow_PuBr"))
    aerosol_collection.set_clim(vmin=3,vmax=4.5)

    n_collection = PatchCollection(patches, array=albedo["n"], cmap=tol_cmap("rainbow_PuBr"))

    fig,ax = plt.subplots(figsize=(15,15), subplot_kw={"projection":ccrs.PlateCarree()})
    ax.imshow(np.tile(np.array([[[200, 200, 200]]],dtype=np.uint8), [2, 2, 1]),origin='upper',transform=ccrs.PlateCarree(),extent=[-180, 180, -90, 90])
    ax.add_collection(sron_collection)
    ax.set_title("TROPOMI SRON v19",fontsize=13)
    plt.colorbar(sron_collection, ax=ax, label=f"XCH$_4$ [ppb]", extend="both")
    ax.coastlines(linewidth=0.1,linestyle="--")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, zorder=-1, linewidth=0)
    gl.top_labels = gl.right_labels = False
    fig.savefig(f"{config['RunDir']}/notebooks/paper_figures/oversampled_figs/{basename}_TROPOMI.png", dpi=300, bbox_inches="tight")

    fig,ax = plt.subplots(figsize=(15,15), subplot_kw={"projection":ccrs.PlateCarree()})
    ax.imshow(np.tile(np.array([[[200, 200, 200]]],dtype=np.uint8), [2, 2, 1]),origin='upper',transform=ccrs.PlateCarree(),extent=[-180, 180, -90, 90])
    ax.add_collection(blended_collection)
    ax.set_title("Blended TROPOMI+GOSAT",fontsize=13)
    plt.colorbar(blended_collection, ax=ax, label=f"XCH$_4$ [ppb]", extend="both")
    ax.coastlines(linewidth=0.1,linestyle="--")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, zorder=-1, linewidth=0)
    gl.top_labels = gl.right_labels = False
    fig.savefig(f"{config['RunDir']}/notebooks/paper_figures/oversampled_figs/{basename}_Blended.png", dpi=300, bbox_inches="tight")

    fig,ax = plt.subplots(figsize=(15,15), subplot_kw={"projection":ccrs.PlateCarree()})
    ax.imshow(np.tile(np.array([[[200, 200, 200]]],dtype=np.uint8), [2, 2, 1]),origin='upper',transform=ccrs.PlateCarree(),extent=[-180, 180, -90, 90])
    ax.add_collection(albedo_collection)
    ax.set_title("TROPOMI SWIR Surface Albedo",fontsize=13)
    plt.colorbar(albedo_collection, ax=ax, label=f"SWIR Surface Albedo", extend="both")
    ax.coastlines(linewidth=0.1,linestyle="--")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, zorder=-1, linewidth=0)
    gl.top_labels = gl.right_labels = False
    fig.savefig(f"{config['RunDir']}/notebooks/paper_figures/oversampled_figs/{basename}_Albedo.png", dpi=300, bbox_inches="tight")
    
    fig,ax = plt.subplots(figsize=(15,15), subplot_kw={"projection":ccrs.PlateCarree()})
    ax.imshow(np.tile(np.array([[[200, 200, 200]]],dtype=np.uint8), [2, 2, 1]),origin='upper',transform=ccrs.PlateCarree(),extent=[-180, 180, -90, 90])
    ax.add_collection(altitude_collection)
    ax.set_title("TROPOMI Surface Altitude",fontsize=13)
    plt.colorbar(altitude_collection, ax=ax, label=f"Surface Altitude", extend="both")
    ax.coastlines(linewidth=0.1,linestyle="--")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, zorder=-1, linewidth=0)
    gl.top_labels = gl.right_labels = False
    fig.savefig(f"{config['RunDir']}/notebooks/paper_figures/oversampled_figs/{basename}_Altitude.png", dpi=300, bbox_inches="tight")
    
    fig,ax = plt.subplots(figsize=(15,15), subplot_kw={"projection":ccrs.PlateCarree()})
    ax.imshow(np.tile(np.array([[[200, 200, 200]]],dtype=np.uint8), [2, 2, 1]),origin='upper',transform=ccrs.PlateCarree(),extent=[-180, 180, -90, 90])
    ax.add_collection(aerosol_collection)
    ax.set_title("TROPOMI Aerosol Size",fontsize=13)
    plt.colorbar(aerosol_collection, ax=ax, label=f"Aerosol Size", extend="both")
    ax.coastlines(linewidth=0.1,linestyle="--")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, zorder=-1, linewidth=0)
    gl.top_labels = gl.right_labels = False
    fig.savefig(f"{config['RunDir']}/notebooks/paper_figures/oversampled_figs/{basename}_Aerosol.png", dpi=300, bbox_inches="tight")

    fig,ax = plt.subplots(figsize=(15,15), subplot_kw={"projection":ccrs.PlateCarree()})
    ax.imshow(np.tile(np.array([[[200, 200, 200]]],dtype=np.uint8), [2, 2, 1]),origin='upper',transform=ccrs.PlateCarree(),extent=[-180, 180, -90, 90])
    ax.add_collection(n_collection)
    ax.set_title("n",fontsize=13)
    plt.colorbar(n_collection, ax=ax, label=f"n", extend="both")
    ax.coastlines(linewidth=0.1,linestyle="--")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, zorder=-1, linewidth=0)
    gl.top_labels = gl.right_labels = False
    fig.savefig(f"{config['RunDir']}/notebooks/paper_figures/oversampled_figs/{basename}_n.png", dpi=300, bbox_inches="tight")
    
if __name__ == "__main__":
    args = [["TROPOMI_NorthAfrica_2018",0.01,[-20,65,0,40]],\
            ["TROPOMI_NorthAfrica_2019",0.01,[-20,65,0,40]],\
            ["TROPOMI_NorthAfrica_2020",0.01,[-20,65,0,40]],\
            ["TROPOMI_NorthAfrica_2021",0.01,[-20,65,0,40]],\
            ["TROPOMI_NorthAfrica_No_Mixed_2021",0.01,[-20,65,0,40]],\
            ["TROPOMI_NorthAmerica_2021",0.01,[-135,-60,15,60]],\
            ["TROPOMI_NorthAmerica_No_Mixed_2021",0.01,[-135,-60,15,60]]]
    
    pool = multiprocessing.Pool(7)
    pool.map(f_plotter, args)
    pool.close()
