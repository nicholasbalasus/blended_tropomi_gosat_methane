#!/bin/bash

#############################################################################
# Script Name: run.sh
# Description: modules to make/evaluate the blended TROPOMI+GOSAT methane product
# Args: controlled by config.yml
# Author: Nicholas Balasus
# Email: nicholasbalasus@g.harvard.edu
#############################################################################

# To run this: sbatch -J bash -p seas_compute -t 3-00:00 --mem 32000 --wrap "bash run.sh" --output run.out

# Source .bashrc to gain access to miniconda
source ~/.bashrc

# Read in the setings from the configuration file
source tools/parse_yaml.sh
eval $(parse_yaml config.yml)

# Make necessary directories if they don't exist
mkdir -p $StorageDir/gosat $StorageDir/tropomi $StorageDir/blended_tropomi_gosat $StorageDir/tccon $StorageDir/geojson $StorageDir/pkl $StorageDir/oversample
mkdir -p $StorageDir/pkl/gosat $StorageDir/pkl/tropomi $StorageDir/pkl/blended_tropomi_gosat $StorageDir/pkl/paired $StorageDir/pkl/flaml $StorageDir/pkl/tccon
mkdir -p $StorageDir/pkl/tccon/gosat $StorageDir/pkl/tccon/tropomi
mkdir -p $StorageDir/oversample/input $StorageDir/oversample/output

#############################################################################
# Module Number: 0 
# Module Name: Make_Conda_Env
# Description: make a conda env from the packages in environment.yml
# Files Used: environment.yml
#############################################################################

if "$Make_Conda_Env"; then
    # if conda env exists, update it; if it does not exist, create it
    if { conda env list | grep 'ch4_env'; } >/dev/null 2>&1; then conda env update -f environment.yml --prune; else conda env create -f environment.yml; fi
fi

#############################################################################
# Module Number: 1A
# Module Name: Download_GOSAT
# Description: download GOSAT level 2 data from UoL for 2018-2021
# Files used: 
#############################################################################

if "$Download_GOSAT"; then
    printf "\n=== Downloading GOSAT Data ===\n"
    start_time=$(date +%s)

    # if directory is not empty, exit the script
    cd $StorageDir/gosat
    [ "$(ls -A)" ] && printf "ERROR: The GOSAT directory is not empty. Remove everything in ${PWD} and try again.\n" && exit 1

    # download all GOSAT data as a tar file - last access 5 Dec 2022
    sbatch -W -J mod_1A -p $Partition -t 60 --mem 4000 --wrap "wget -q https://dap.ceda.ac.uk/neodc/gosat/data/ch4/nceov1.0/CH4_GOS_OCPR/CH4_GOS_OCPR_v9.0_final_nceo_2009_2021.tar.gz"; wait;
    rm slurm*.out
    # extract the tar file then remove it
    tar xzf CH4_GOS_OCPR_v9.0_final_nceo_2009_2021.tar.gz
    rm CH4_GOS_OCPR_v9.0_final_nceo_2009_2021.tar.gz
    # move all of the individual files we want into $StorageDir/gosat
    cd CH4_GOS_OCPR
    mv 2018 2019 2020 2021 -t ..
    cd ..
    rm -r CH4_GOS_OCPR
    mv 2018/* 2019/* 2020/* 2021/* -t .
    rm -r 2018 2019 2020 2021

    end_time=$(date +%s)
    printf "=== Finished Downloading GOSAT Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 1B
# Module Name: Download_TROPOMI
# Description: download TROPOMI data for 2018-2021 from the SRON ftp
# Files Used: 
#############################################################################

if "$Download_TROPOMI"; then
    printf "\n=== Downloading TROPOMI Data ===\n"
    start_time=$(date +%s)
    cd $StorageDir/tropomi

    # if directory is not empty, exit the script
    [ "$(ls -A)" ] && printf "ERROR: The TROPOMI directory is not empty. Remove everything in ${PWD} and try again.\n" && exit 1

    # write a text file of all of the urls to download (some of these orbit numbers might not exist in the ftp which is okay with wget)
    # orbit numbers 01130-21856 corresponds to all of 2018-2021 - last access 5 Dec 2022
    touch tropomi_urls.txt
    for i in {01130..21856}
    do
        echo "ftp://ftp.sron.nl/open-access-data-2/TROPOMI/tropomi/ch4/19_446/s5p_l2_ch4_0446_$i.nc" >> tropomi_urls.txt
    done; wait;

    # download tropomi data with a maximum of 8 wget running simultaneously
    sbatch -W -J mod_1B -p $Partition -t 300 -c 8 --mem 32000 --wrap "cat tropomi_urls.txt | xargs -n 1 -P 8 wget"; wait;
    rm slurm*.out

    # remove text file of urls
    rm tropomi_urls.txt

    end_time=$(date +%s)
    printf "=== Finished Downloading TROPOMI Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 1C
# Module Name: Download_TCCON
# Description: download TCCON data from tccondata.org
# Files Used: 
#############################################################################

if "$Download_TCCON"; then
    printf "\n=== Downloading TCCON Data ===\n"
    start_time=$(date +%s)
    cd $StorageDir/tccon

    # if directory is not empty, exit the script
    [ "$(ls -A)" ] && printf "ERROR: The TCCON directory is not empty. Remove everything in ${PWD} and try again.\n" && exit 1

    # download all TCCON data as a tgz file
    # data last accessed at this link on 5 Dec 2022
    sbatch -W -J mod_1C -p $Partition -t 30 --mem 4000 --wrap "wget https://renc.osn.xsede.org/ini210004tommorrell/10.14291/TCCON.GGG2020/tccon.latest.public.tgz --max-redirect=2 --trust-server-names --content-disposition -q"; wait;
    rm slurm*.out
    
    # extract the tar file then remove it
    tar xzf tccon.latest.public.tgz
    rm tccon.latest.public.tgz

    end_time=$(date +%s)
    printf "=== Finished Downloading TCCON Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2A
# Module Name: Process_GOSAT
# Description: process GOSAT data to one dataframe using process_gosat.py
# Files Used: process_gosat.py
#############################################################################

if "$Process_GOSAT"; then
    printf "\n=== Processing GOSAT Data ===\n"
    start_time=$(date +%s)

    # if directory is not empty, exit the script
    cd $StorageDir/pkl/gosat
    [ "$(ls -A)" ] && printf "ERROR: The GOSAT directory is not empty. Remove everything in ${PWD} and try again.\n" && exit 1

    cd $RunDir
    sbatch -W --export=NONE -J mod_2A -p $Partition -t 360 --mem 8000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/process_gosat.py"; wait;
    rm slurm*.out

    end_time=$(date +%s)
    printf "=== Finished Processing GOSAT Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2B
# Module Name: Process_TROPOMI
# Description: process each netCDF TROPOMI file to a pickled dataframe
# Files Used: process_tropomi.py
#############################################################################

if "$Process_TROPOMI"; then
    printf "\n=== Processing TROPOMI Data ===\n"
    start_time=$(date +%s)
    
    # write tropomi.txt which has all of the files to process
    cd $StorageDir/tropomi
    printf '%s\n' * > tropomi.txt
    mv tropomi.txt $RunDir
    cd $RunDir

    # run process_tropomi.py 64 times (each running 1/64 of the files listed in tropomi.txt) with 3 GB per core
    # had to escape the slurm id vars (from https://stackoverflow.com/questions/63516186/accessing-task-id-for-array-jobs)
    sbatch -W --export=NONE -J mod_2B -p $Partition -t 1440 --mem $(($Cores*3000)) -c $Cores --array=0-63 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/process_tropomi.py \${SLURM_ARRAY_TASK_ID} \${SLURM_ARRAY_TASK_COUNT}"; wait;
    rm slurm*.out

    # remove text file of file names
    rm tropomi.txt

    end_time=$(date +%s)
    printf "=== Finished Processing TROPOMI Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2C
# Module Name: Pair_TROPOMI_GOSAT
# Description: pair TROPOMI and GOSAT measurements
# Files Used: pair_tropomi_gosat.py
#############################################################################

if "$Pair_TROPOMI_GOSAT"; then
    printf "\n=== Pairing TROPOMI and GOSAT Data ===\n"
    start_time=$(date +%s)

    # write tropomi.txt which has all of the files to process
    cd $StorageDir/pkl/tropomi
    printf '%s\n' * > tropomi.txt
    mv tropomi.txt $RunDir
    cd $RunDir

    # run pair_tropomi_gosat.py 64 times (each running 1/64 of the files listed in tropomi.txt) with 3 GB per core
    # this loops through each tropomi file in $StorageDir/pkl/tropomi and makes another dataframe in $StorageDir/pkl/paired of only the points with a gosat pair
    # had to escape the slurm id vars (from https://stackoverflow.com/questions/63516186/accessing-task-id-for-array-jobs)
    sbatch -W --export=NONE -J mod_2C -p $Partition -t 2880 --mem $(($Cores*3000)) -c $Cores --array=0-63 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/pair_tropomi_gosat.py \${SLURM_ARRAY_TASK_ID} \${SLURM_ARRAY_TASK_COUNT}"; wait;
    rm slurm*.out

    # remove text file of file names
    rm tropomi.txt

    end_time=$(date +%s)
    printf "=== Finished Pairing TROPOMI and GOSAT Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2D
# Module Name: Process_TROPOMI_GOSAT_Pairs
# Description: concatenate all pairs and calculate delta(GOSAT-TROPOMI)
# Files Used: process_tropomi_gosat_pairs.py
#############################################################################

if "$Process_TROPOMI_GOSAT_Pairs"; then
    printf "\n=== Processing TROPOMI and GOSAT Pairs ===\n"
    start_time=$(date +%s)

    sbatch -W --export=NONE -J mod_2D -p $Partition -t 240 --mem 160000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/process_tropomi_gosat_pairs.py"; wait;
    rm slurm*.out

    end_time=$(date +%s)
    printf "=== Finished Processing TROPOMI and GOSAT Pairs in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2E
# Module Name: Pair_GOSAT_TCCON
# Description: make dataframes of GOSAT/TCCON pairs for each TCCON site
# Files Used: pair_gosat_tccon.py
#############################################################################

if "$Pair_GOSAT_TCCON"; then
    printf "\n=== Pairing GOSAT and TCCON Data ===\n"
    start_time=$(date +%s)

    # use 25 cores (# of TCCON stations)
    sbatch -W --export=NONE -J mod_2E -p $Partition -t 240 -c 25 --mem 160000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/pair_gosat_tccon.py"; wait;
    rm slurm*.out

    end_time=$(date +%s)
    printf "=== Finished Pairing GOSAT and TCCON Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2F
# Module Name: Run_FLAML
# Description: Run AutoML to build model to predict delta(GOSAT-TROPOMI)
# Files Used: run_flaml.py
#############################################################################

if "$Run_FLAML"; then
    printf "\n=== Running FLAML ===\n"
    start_time=$(date +%s)

    # use 8 cores to run FLAML
    sbatch -W --export=NONE -J mod_2F -p $Partition -t 300 -c 8 --mem 64000 --output flaml.out --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/run_flaml.py"; wait;

    end_time=$(date +%s)
    printf "=== Finished Running FLAML in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2G
# Module Name: Predict_Delta_GOSAT_TROPOMI
# Description: predict and apply delta(GOSAT-TROPOMI)to all TROPOMI data
# Files Used: apply_ml_correction.py
#############################################################################

if "$Predict_Delta_GOSAT_TROPOMI"; then
    printf "\n=== Predicting delta(GOSAT-TROPOMI) for all TROPOMI Data ===\n"
    start_time=$(date +%s)

    # write tropomi.txt which has all of the files to process
    cd $StorageDir/pkl/tropomi
    printf '%s\n' * > tropomi.txt
    mv tropomi.txt $RunDir
    cd $RunDir

    # run apply_ml_correction.y 64 times (each running 1/64 of the files listed in tropomi.txt) with 3 GB per core
    sbatch -W --export=NONE -J mod_2G -p $Partition -t 120 --mem $(($Cores*3000)) -c $Cores --array=0-63 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/apply_ml_correction.py \${SLURM_ARRAY_TASK_ID} \${SLURM_ARRAY_TASK_COUNT}"; wait;
    rm slurm*.out

    # remove text file of file names
    rm tropomi.txt

    end_time=$(date +%s)
    printf "=== Finished Predicting delta(GOSAT-TROPOMI) for all TROPOMI Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2H
# Module Name: Pair_TROPOMI_TCCON
# Description: make dataframes of TROPOMIT/TCCON pairs for each TCCON site
# Files Used: pair_tropomi_tccon.py
#############################################################################

if "$Pair_TROPOMI_TCCON"; then
    printf "\n=== Pairing TROPOMI and TCCON Data ===\n"
    start_time=$(date +%s)

    # use 50 cores (# of TCCON stations * (blended + operational))
    sbatch -W --export=NONE -J mod_2H -p $Partition -t 1080 -c 50 --mem 250000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/pair_tropomi_tccon.py"; wait;
    rm slurm*.out

    end_time=$(date +%s)
    printf "=== Finished Pairing TROPOMI and TCCON Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2I
# Module Name: SHAP_Explainer
# Description: make SHAP explainer and calculate shap values for train data
# Files Used: shap_explainer.py
#############################################################################

if "$SHAP_Explainer"; then
    printf "\n=== Making SHAP Explainer ===\n"
    start_time=$(date +%s)

    sbatch -W --export=NONE -J mod_2I -p $Partition -t 5760 --mem 64000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/shap_explainer.py"; wait;
    rm slurm*.out

    end_time=$(date +%s)
    printf "=== Finished Making SHAP Explainer in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 3A
# Module Name: Write_NetCDF
# Description: write netCDF files for the blended product
# Files Used: write_nc_files.py
#############################################################################

if "$Write_NetCDF"; then
    printf "\n=== Writing netCDF files for the blended product ===\n"
    start_time=$(date +%s)

    # write tropomi.txt which has all of the files to process
    cd $StorageDir/pkl/blended_tropomi_gosat
    printf '%s\n' * > tropomi.txt
    mv tropomi.txt $RunDir
    cd $RunDir

    # run apply_ml_correction.y 64 times (each running 1/64 of the files listed in tropomi.txt) with 3 GB per core
    sbatch -W --export=NONE -J mod_3A -p $Partition -t 1-00:00 --mem $(($Cores*3000)) -c $Cores --array=0-63 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/write_nc_files.py \${SLURM_ARRAY_TASK_ID} \${SLURM_ARRAY_TASK_COUNT}"; wait;
    rm slurm*.out

    # remove text file of file names
    rm tropomi.txt

    end_time=$(date +%s)
    printf "=== Finished writing netCDF files for the blended product in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 3B
# Module Name: Paired_Regrid
# Description: regrid paired GOSAT/TROPOMI data to standard grid
# Files Used: regrid_pairs.py
#############################################################################

if "$Paired_Regrid"; then
    printf "\n=== Regridding Paired Data ===\n"
    start_time=$(date +%s)

    sbatch -W --export=NONE -J mod_3B -p $Partition -t 1440 --mem 128000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/regrid_pairs.py"; wait;
    rm slurm*.out

    end_time=$(date +%s)
    printf "=== Finished Regridding Paired Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 3C
# Module Name: TROPOMI_Regrid
# Description: regrid TROPOMI data to standard grid
# Files Used: regrid_tropomi.py
#############################################################################

if "$TROPOMI_Regrid"; then
    printf "\n=== Regridding TROPOMI Data ===\n"
    start_time=$(date +%s)

    sbatch -W --export=NONE -J mod_3C -p $Partition -t 1440 --mem 200000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/regrid_tropomi.py"; wait;
    rm slurm*.out

    end_time=$(date +%s)
    printf "=== Finished Regridding TROPOMI Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 3D
# Module Name: Oversample_TROPOMI
# Description: oversample TROPOMI data to 0.01 degrees
# Files Used: csv_tropomi_oversample.py, tools/oversample/*
#############################################################################

if "$Oversample_TROPOMI"; then
    printf "\n=== Oversampling TROPOMI data ===\n"
    start_time=$(date +%s)

    # write csv files of the TROPOMI data for the requested regions, time frames, and variables (all specified in csv_tropomi_oversample.py)
    sbatch -W --export=NONE -J mod_3Da -p $Partition -t 2-00:00 --mem 1800000 -c 6 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/csv_tropomi_oversample.py"; wait;
    rm slurm*.out

    # oversample the data
    cd $RunDir/tools/oversample
    sbatch -W -J mod_3Db -p $Partition -t 1-00:00 --mem 300000 --wrap "./oversampling.sh"
    rm slurm*.out
    cd $RunDir

    end_time=$(date +%s)
    printf "=== Finished oversampling TROPOMI data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
    
    sbatch -W --export=NONE -J mod_3Dc -p $Partition -t 1-00:00 --mem 1800000 -c 6 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/plotter.py"
fi