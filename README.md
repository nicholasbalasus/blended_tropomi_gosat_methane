# Blended TROPOMI GOSAT Methane Product
The entire project can be run by running `run.sh` (e.g., `sbatch -J bash -p seas_compute -t 3-00:00 --mem 1000 --wrap "bash run.sh" --output run.out`). The files that run, where they run, where the outputs are saved, and everything else is controlled by `config.yml`. All of the code that I have written is in `src/`, while code from others is in `tools/`.

The project is broken into modules for downloading data, processing data, and writing data. Each of these modules is broken down below with approximations for their run times (run on the `serial_requeue` partion of Harvard's Cannon cluster), number of cores requested, and total amount of memory requested. Because I was using `serial_requeue`, the resources requested are large. These can be reduced in exchanged for longer run times. At the end of the project, the storage directory specificed in `config.yml` will be ~1.3 TB. After all of the modules have been run, `notebooks/paper.ipynb` can be run to make the figures.

* Module 1: Download data
    * **Download_GOSAT**: download GOSAT level 2 data from UoL for 2018-2021 (~5 minutes, 1 core, 4 GB).
    * **Download_TROPOMI**: download TROPOMI level 2 data for 2018-2021 from the SRON ftp (~200 minutes, 8 cores, 32 GB).
    * **Download_TCCON**: download TCCON data from tccondata.org (~1 minute, 1 core, 4 GB).

* Module 2: Process data
    * **Process_GOSAT**: process all daily netCDF GOSAT data to one dataframe (~210 minutes, 1 core, 8 GB).
    * **Process_TROPOMI**: process each netCDF TROPOMI file to a pickled dataframe (~100 minutes, 1024 cores, 3072 GB).
    * **Pair_TROPOMI_GOSAT**: pair TROPOMI and GOSAT measurements with time and distance thresholds specificed in `config.yml` (~1000 minutes, 1024 cores, 3072 GB).
    * **Process_TROPOMI_GOSAT_Pairs**: concatenate all pairs and calculate delta(TROPOMI-GOSAT) (~10 minutes, 1 core, 160 GB).
    * **Pair_GOSAT_TCCON**: make dataframes of GOSAT/TCCON pairs (with and without global GOSAT offset) for each TCCON site (~2 minutes, 25 cores, 160 GB).
    * **Run_FLAML**: train models to predict delta(TROPOMI-GOSAT) (~90 minutes, 8 cores, 64 GB).
    * **Predict_Delta_GOSAT_TROPOMI**: predict and remove delta(TROPOMI-GOSAT) from all TROPOMI data (~10 minutes, 1024 cores, 3072 GB).
    * **Pair_TROPOMI_TCCON**: make dataframes of TROPOMI/TCCON pairs and Blended/TCCON pairs for each TCCON site (~100 minutes, 50 cores, 250 GB).
    * **SHAP_Explainer**: make SHAP explainer and calculate shap values for train data (~320 minutes, 1 core, 64 GB).

* Module 3: Write data
    * **Write_NetCDF**: write netCDF files that mimic the original TROPOMI data but add a variable for the blended product (~5 minutes, 512 cores, 1536 GB).
    * **Paired_Regrid**: regrid the TROPOMI and GOSAT pairs to a standard grid (~180 minute, 1 core, 128 GB).
    * **TROPOMI_Regrid**: regrid the TROPOMI data to a standard grid (~1100 minutes, 1 core, 200 GB).
    * **Oversample_TROPOMI**: for specific regions, oversample TROPOMI data to 0.01 degrees (~600 minutes, 4 cores, 1200 GB). 
