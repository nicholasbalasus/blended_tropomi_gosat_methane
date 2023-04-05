#!/bin/bash

#############################################################################
# Script Name: run.sh
# Description: modules to make/evaluate the blended TROPOMI+GOSAT methane product
# Args: controlled by config.yml
# Author: Nicholas Balasus
# Email: nicholasbalasus@g.harvard.edu
#############################################################################

# To run this: sbatch -J bash -p huce_intel -t 14-00:00 --mem 32000 --wrap "bash run.sh" --output run.out

# Source .bashrc to gain access to miniconda
source ~/.bashrc

# Read in the setings from the configuration file
source tools/parse_yaml.sh
eval $(parse_yaml config.yml)

# Make necessary directories if they don't exist
mkdir -p $StorageDir/gosat $StorageDir/tropomi $StorageDir/tccon $StorageDir/processed $StorageDir/blended $StorageDir/oversample/input $StorageDir/oversample/output $RunDir/notebooks/paper_figures

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
# Module Number: 1
# Module Name: Download_GOSAT
# Description: download GOSAT level 2 data from UoL for 2018-2021
# Files used: 
#############################################################################

if "$Download_GOSAT"; then
    printf "=== Downloading GOSAT Data ===\n"
    start_time=$(date +%s)

    # if directory is not empty, exit the script
    cd $StorageDir/gosat
    [ "$(ls -A)" ] && printf "ERROR: The GOSAT directory is not empty. Remove everything in ${PWD} and try again.\n" && exit 1

    # download all GOSAT data as a tar file - last access 18 Mar 2023
    sbatch -W -J mod_1 -p $Partition -t 7-00:00 --mem 4000 --wrap "wget -q https://dap.ceda.ac.uk/neodc/gosat/data/ch4/nceov1.0/CH4_GOS_OCPR/CH4_GOS_OCPR_v9.0_final_nceo_2009_2021.tar.gz"; wait;
    if ! $Debug; then rm slurm*.out; fi
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
    # remove the first four months of 2018 where we don't have TROPOMI data
    rm UoL-GHG-L2-CH4-GOSAT-OCPR-201801*-fv9.0.nc
    rm UoL-GHG-L2-CH4-GOSAT-OCPR-201802*-fv9.0.nc
    rm UoL-GHG-L2-CH4-GOSAT-OCPR-201803*-fv9.0.nc
    mv UoL-GHG-L2-CH4-GOSAT-OCPR-20180430-fv9.0.nc -t ..
    rm UoL-GHG-L2-CH4-GOSAT-OCPR-201804*-fv9.0.nc
    mv ../UoL-GHG-L2-CH4-GOSAT-OCPR-20180430-fv9.0.nc -t .

    end_time=$(date +%s)
    printf "=== Finished Downloading GOSAT Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 2
# Module Name: Download_TROPOMI
# Description: download TROPOMI data for 2018-2021 from copernicus
# Files Used: 
#############################################################################

if "$Download_TROPOMI"; then
    printf "=== Downloading TROPOMI Data ===\n"
    start_time=$(date +%s)
    cd $StorageDir/tropomi

    # if directory is not empty, exit the script
    [ "$(ls -A)" ] && printf "ERROR: The TROPOMI directory is not empty. Remove everything in ${PWD} and try again.\n" && exit 1

    # write a text file of all of the urls to download
    # we can only get information about 100 files at a time, so run a while loop until we run out of files (using prev_length variable)
    touch files.txt
    prev_length=-1
    s=0
    while true; do
        wget -q --no-check-certificate --user=s5pguest --password=s5pguest --output-document="output.txt" "https://s5phub.copernicus.eu/dhus/search?q=platformname:Sentinel-5 AND producttype:L2__CH4___ AND processingmode:Reprocessing AND processorversion:020400 AND endposition:[2018-01-01T00:00:00.000Z TO 2021-12-31T23:59:59.999Z]&rows=100&start=$s"
        grep "link href=\"https" output.txt >> files.txt
        length=$(wc -l < files.txt)
        if [ "$length" -eq "$prev_length" ]; then
            break
        fi
        prev_length=$length
        s=$((s+100))
    done
    sed -i -e 's/<link href=//g' files.txt
    sed -i -e 's/\/>//g' files.txt

    # download tropomi data with a maximum of 16 wget running simultaneously
    sbatch -W -J mod_2 -p $Partition -t 7-00:00 -c 16 --mem 32000 --wrap "xargs -n 1 -P 16 wget --content-disposition --continue --no-check-certificate --tries=5 --user=s5pguest --password=s5pguest < files.txt"; wait;
    if ! $Debug; then rm slurm*.out; fi

    # remove text files
    rm files.txt
    rm output.txt

    end_time=$(date +%s)
    printf "=== Finished Downloading TROPOMI Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 3
# Module Name: Download_TCCON
# Description: download TCCON data from tccondata.org
# Files Used: 
#############################################################################

if "$Download_TCCON"; then
    printf "=== Downloading TCCON Data ===\n"
    start_time=$(date +%s)
    cd $StorageDir/tccon

    # if directory is not empty, exit the script
    [ "$(ls -A)" ] && printf "ERROR: The TCCON directory is not empty. Remove everything in ${PWD} and try again.\n" && exit 1

    # download all TCCON data as a tgz file
    # data last accessed at this link on 18 Mar 2023
    sbatch -W -J mod_3 -p $Partition -t 7-00:00 --mem 4000 --wrap "wget https://renc.osn.xsede.org/ini210004tommorrell/10.14291/TCCON.GGG2020/tccon.latest.public.tgz --max-redirect=2 --trust-server-names --content-disposition -q"; wait;
    if ! $Debug; then rm slurm*.out; fi
    
    # extract the tar file then remove it
    tar xzf tccon.latest.public.tgz
    rm tccon.latest.public.tgz

    end_time=$(date +%s)
    printf "=== Finished Downloading TCCON Data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 4
# Module Name: Calculate_Delta_GOSAT_TCCON
# Description: for each TCCON station, find ∆(GOSAT-TCCON)
# Files Used: calculate_delta_gosat_tccon.py, utilities.py
#############################################################################

if "$Calculate_Delta_GOSAT_TCCON"; then
    printf "=== Calculating delta(GOSAT-TCCON) ===\n"
    start_time=$(date +%s)

    cd $RunDir
    sbatch -W --export=NONE -J mod_4 -p $Partition -t 7-00:00 --mem 500000 -c 25 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/calculate_delta_gosat_tccon.py"; wait;
    if ! $Debug; then rm slurm*.out; fi

    end_time=$(date +%s)
    printf "=== Finished Calculating delta(GOSAT-TCCON) in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 5
# Module Name: Calculate_Delta_TROPOMI_GOSAT
# Description: create a dataframe of TROPOMI/GOSAT pairs
# Files Used: calculate_delta_tropomi_gosat.py, utilities.py
#############################################################################

if "$Calculate_Delta_TROPOMI_GOSAT"; then
    printf "=== Calculating delta(TROPOMI-GOSAT) ===\n"
    start_time=$(date +%s)

    cd $RunDir
    sbatch -W --export=NONE -J mod_5 -p $Partition -t 7-00:00 --mem 500000 -c 64 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/calculate_delta_tropomi_gosat.py"; wait;
    if ! $Debug; then rm slurm*.out; fi

    end_time=$(date +%s)
    printf "=== Finished Calculating delta(TROPOMI-GOSAT) in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 6
# Module Name: Run_FLAML_SHAP
# Description: Run AutoML to predict delta(TROPOMI-GOSAT) then SHAP
# Files Used: run_flaml_and_shap.py
#############################################################################

if "$Run_FLAML_SHAP"; then
    printf "\n=== Running FLAML and SHAP ===\n"
    start_time=$(date +%s)

    # use 8 cores to run FLAML
    sbatch -W --export=NONE -J mod_6 -p $Partition -t 7-00:00 -c 8 --mem 64000 --output flaml.out --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/run_flaml_and_shap.py"; wait;

    end_time=$(date +%s)
    printf "=== Finished Running FLAML and SHAP in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 7
# Module Name: Write_Blended_Files
# Description: write blended TROPOMI+GOSAT netCDF files
# Files Used: write_blended_files.py, utilities.py
#############################################################################

if "$Write_Blended_Files"; then
    printf "=== Writing Blended TROPOMI+GOSAT Files ===\n"
    start_time=$(date +%s)

    cd $RunDir
    sbatch -W --export=NONE -J mod_7 -p $Partition -t 7-00:00 --mem 500000 -c 64 --array=0-15 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/write_blended_files.py \${SLURM_ARRAY_TASK_ID} \${SLURM_ARRAY_TASK_COUNT}"; wait;
    if ! $Debug; then rm slurm*.out; fi

    end_time=$(date +%s)
    printf "=== Finished Writing Blended TROPOMI+GOSAT Files in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 8
# Module Name: Calculate_Delta_TROPOMI_TCCON
# Description: for each TCCON station, find ∆(TROPOMI-TCCON)
# Files Used: calculate_delta_tropomi_tccon.py, utilities.py
#############################################################################

if "$Calculate_Delta_TROPOMI_TCCON"; then
    printf "=== Calculating delta(TROPOMI-TCCON) ===\n"
    start_time=$(date +%s)

    cd $RunDir
    sbatch -W --export=NONE -J mod_8 -p $Partition -t 7-00:00 --mem 500000 -c 64 --array=0-24 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/calculate_delta_tropomi_tccon.py \${SLURM_ARRAY_TASK_ID} \${SLURM_ARRAY_TASK_COUNT}"; wait;
    if ! $Debug; then rm slurm*.out; fi

    end_time=$(date +%s)
    printf "=== Finished Calculating delta(TROPOMI-TCCON) in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi

#############################################################################
# Module Number: 9
# Module Name: Oversample_TROPOMI
# Description: oversample TROPOMI and Blended data to 0.01 degrees
# Files Used: write_oversample_input_csv.py, tools/oversample/*
#############################################################################

if "$Oversample_TROPOMI"; then
    printf "\n=== Oversampling TROPOMI data ===\n"
    start_time=$(date +%s)

    # write csv files of the TROPOMI data for the requested regions, time frames, and variables (all specified in csv_tropomi_oversample.py)
    sbatch -W --export=NONE -J mod_9a -p $Partition -t 7-00:00 --mem 500000 -c 2 --wrap "source ~/.bashrc; conda activate $CondaEnv; python ./src/write_oversample_input_csv.py"; wait;
    if ! $Debug; then rm slurm*.out; fi

    # oversample the data
    cd $RunDir/tools/oversample
    sbatch -W -J mod_9b -p $Partition -t 1-00:00 --mem 500000 --wrap "./oversampling.sh $StorageDir"
    if ! $Debug; then rm slurm*.out; fi
    cd $RunDir

    end_time=$(date +%s)
    printf "=== Finished oversampling TROPOMI data in $(( ($end_time- $start_time)/60 )) Minutes ===\n"
fi