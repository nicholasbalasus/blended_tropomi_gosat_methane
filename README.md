# Blended TROPOMI+GOSAT Methane Product
The entire project can be run by running `run.sh` (e.g., `sbatch -J bash -p huce_intel -t 14-00:00 --mem 32000 --wrap "bash run.sh" --output run.out`). The files that run, where they run, where the outputs are saved, and everything else is controlled by `config.yml`. All of the code that I have written is in `src/`, while code from others is in `tools/`.

The project is broken into modules, each with approximations for their run times (run on the `huce_ice` partion of Harvard's Cannon cluster), number of cores requested, and total amount of memory requested. At the end of the project, the storage directory specificed in `config.yml` will be ~1 TB. After all of the modules have been run, `notebooks/paper.ipynb` can be run to make the figures.

0. **Make_Conda_Env**: make or update the conda environment specified in `environment.yml`.
1. **Download_GOSAT**: download GOSAT data from UoL for 2018-2021 (~5 minutes, 1 core, 4 GB).
2. **Download_TROPOMI**: download operational TROPOMI data for 2018-2021 from the copernicus hub (~1.5 days, 8 cores, 32 GB).
3. **Download_TCCON**: download TCCON data from tccondata.org (~1 minute, 1 core, 4 GB).
4. **Calculate_Delta_GOSAT_TCCON**: for each TCCON station, find GOSAT pairs for TCCON observations and calculate delta(GOSAT-TCCON) (~3.5 hours, 25 cores, 500 GB).
5. **Calculate_Delta_TROPOMI_GOSAT**: pair TROPOMI and GOSAT measurements and calculate delta(TROPOMI-GOSAT) (~4 days, 64 cores, 500 GB).
6. **Run_FLAML_SHAP**: train models to predict delta(TROPOMI-GOSAT) then run SHAP (~4.5 hours, 8 cores, 64 GB).
7. **Write_Blended_Files**: write netCDF files with an added variable for `methane_mixing_ratio_blended` (~1.5 hours, 16 jobs at 64 cores, 500 GB)
8. **Calculate_Delta_TROPOMI_TCCON**: for each TCCON station, find TROPOMI pairs for TCCON observations and calculate delta(TROPOMI-TCCON) (~10 hours, 25 jobs at 64 cores, 500 GB).
9. **Oversample_TROPOMI**: oversample the TROPOMI and Blended data to a 0.01 degree grid for 2021 (~8.5 hours, 2 cores, 500 GB).
