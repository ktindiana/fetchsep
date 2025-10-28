#!/bin/bash

#PLEASE READ IN ENTIRETY BEFORE RUNNING

#SCRIPT TO GENERATE THE CLEAR BENCHMARK DATASET for Linux
#The CLEAR Benchmark Dataset was generated with python 3.10.17.

#REPLACE THE "python" CALL TO REFLECT YOUR REQUIRED COMMAND.

#The steps in this procedure will download all necessary data, process it,
#and create individual curated SEP lists for each GOES spacecraft.
#These lists may be combined to a final list by choosing the primary
#spacecraft for each time period.

#SPECIFIC PYTHON CALL
#--------------------
#The steps provided here work on Linux where the command to call python is python.
#The python command in the script may need to be modified for your specific call,
#e.g. python3.10.

#CONFIG FILE
#-----------
#!!!!If you already have a fetchsep.cfg file in the main fetchsep directory
#that you have customized, copy it to a different name (e.g. fetchsep.cfg.bak)
#so that it will not be overwritten. !!!!

#This script assumes the data directory is located in fetchsep/data and all GOES data
#will be downloaded here. If you want your data directory elsewhere, then
#fetchsep/reference/CLEAR/fetchsep_CLEAR.cfg
#must be edited with the location of your desired data path.

#Output files are stored in fetchsep/CLEAR/output and fetchsep/CLEAR/plots and should
#not be modified to ensure the script works properly.

#RUNNING THE SCRIPT
#------------------
#Run the script in the top directory of the repository, fetchsep/
#You may need to change your execution permissions:
#chmod u+x ./fetchsep/reference/CLEAR/deploy_CLEAR_Linux.sh
#Run in the terminal as:
#./fetchsep/reference/CLEAR/deploy_CLEAR_Linux.sh

#COMMAND LINE ARGUMENTS
#----------------------
#In the command line, the argument LISTS may be added to skip the calculation of
#the mean background with IDSEP IF IT WAS ALREADY DONE PREVIOUSLY.
#The background calculation is extremely time consuming and it needn't be repeated.
#If the user only wants to regenerate the SEP lists,
#the argument LISTS will skip to that part of the script.

#When run without an argument, will default to startpoint ALL and the script will:
#    - Set up environment
#    - Download data and calculate mean background and sigma with IDSEP (12+ hours)
#    - Copy CLEAR curated batch files
#    - Generate SEP lists (2+ hours)

#When run with the argument LISTS, the script will:
#    - Set up environment
#    - Copy CLEAR curated batch files
#    - Generate SEP lists (2+ hours)

startpoint=${1:-"ALL"}
if [ -z "$1" ]; then
  echo "No argument provided for startpoint. Using default ALL to generate the CLEAR dataset from scratch."
  startpoint="ALL"
else
  startpoint="$1"
fi

echo "Starting at $startpoint"

######################################################################
############ SET UP ENVIRONMENT AND CREATE DIRECTORIES ###############
######################################################################
# set up environment and create directories
date '+%Y-%m-%d %H:%M:%S'
echo "Setting up environment"
export PYTHONPATH="$PYTHONPATH:$PWD"
cp fetchsep/reference/CLEAR/fetchsep_CLEAR.cfg ./fetchsep.cfg
python fetchsep/utils/config.py


declare -A start_date=(
   [GOES-06]=1986-01-01
   [GOES-07]=1987-03-01
   [GOES-08]=1995-01-01
   [GOES-10]=1998-07-09
   [GOES-11]=2003-06-01
   [GOES-13]=2010-05-01
   [GOES-15]=2011-01-01
   [GOES-RT/primary]=2020-03-08
   [GOES-RT/secondary]="2021-09-23 15:45:00"
)

declare -A end_date=(
   [GOES-06]=1994-12-01
   [GOES-07]=1996-09-01
   [GOES-08]=2003-06-18
   [GOES-10]=2004-07-01
   [GOES-11]=2011-03-01
   [GOES-13]=2018-01-01
   [GOES-15]=2020-03-05
   [GOES-RT/primary]=2025-09-11
   [GOES-RT/secondary]=2025-09-11
)


######################################################################
################# CALCULATE BACKGROUNDS WITH IDSEP ###################
######################################################################
#Calculate the mean background and sigma solutions for each spacecraft
#This step is very time consuming and may take 12 hours or more.
if [[ "${startpoint}" = "ALL" ]]; then
    for n in 06 07 08 10 11 13 15 RT; do
       if [[ $n != RT ]]; then
          date '+%Y-%m-%d %H:%M:%S'
          echo "[GOES-${n}] Calculate background with idsep"
          python bin/idsep \
             --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
             --Experiment GOES-${n} --FluxType integral --RemoveAbove 10 --saveplot \
             >CLEAR/output/GOES-${n}_integral_idsep.log

          echo
       else
          for type in primary secondary; do
             date '+%Y-%m-%d %H:%M:%S'
             echo "[GOES-${n}/${type}] Calculate background with idsep"
             python bin/idsep \
                --StartDate "${start_date[GOES-${n}/${type}]}" --EndDate "${end_date[GOES-${n}/${type}]}" \
                --Experiment GOES_${n} --Spacecraft $type --FluxType integral --RemoveAbove 10 --saveplot \
                >CLEAR/output/GOES_${n}_integral_${type}_idsep.log

             echo
          done
       fi

       # Sandberg and Bruno calibrations available only for GOES-13 and GOES-15
       [[ $n == 13 || $n == 15 ]] && {
          date '+%Y-%m-%d %H:%M:%S'
          echo "[GOES-${n}/uncor_S14_B17] Calculate background with idsep"
          python bin/idsep \
             --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
             --Experiment GOES-${n} --FluxType differential --RemoveAbove 10 --saveplot \
             --options "S14;Bruno2017;uncorrected" \
             >CLEAR/output/GOES-${n}_differential_uncor_S14_B17_idsep.log

          echo
       }
    done
fi

######################################################################
############## BATCH OPSEP USING CURATED CLEAR LISTS #################
######################################################################
#Generating the SEP event lists for each spacecraft (~2 hours).
#If no changes are needed to be made to the background solutions, but there
#are updates desired for the SEP event lists, then the script may be
#run starting at this point to regenerate the SEP event lists.
# MAKE SURE TO INCLUDE --StartPoint BATCH or the batch file will be
#overwritten by a new run of idsep within fetchsep_prepare_obs
if [[ "${startpoint}" = "ALL" ]] || [[ "${startpoint}" = "LISTS" ]]; then
    for n in 06 07 08 10 11 13 15 RT; do
       if [[ $n != RT ]]; then
          date '+%Y-%m-%d %H:%M:%S'
          echo "[GOES-${n}] Copy curated batch files"
          cp fetchsep/reference/CLEAR/batch_event_list_GOES-${n}_integral_enhance_idsep_CLEAR.txt \
             CLEAR/output/idsep/GOES-${n}_integral/

          date '+%Y-%m-%d %H:%M:%S'
          echo "[GOES-${n}] Batch opsep using curated CLEAR lists"
          python bin/fetchsep_prepare_obs \
             --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
             --Experiment GOES-${n} --FluxType integral --Threshold "30,1;50,1" \
             --BatchFile batch_event_list_GOES-${n}_integral_enhance_idsep_CLEAR.txt \
             --IDSEPEnhancement --Associations --StartPoint BATCH \
             >CLEAR/output/GOES-${n}_integral_batch.log
          echo
       else
          for type in primary secondary; do
             date '+%Y-%m-%d %H:%M:%S'
             echo "[GOES-${n}/${type}] Copy curated batch files"
             cp fetchsep/reference/CLEAR/batch_event_list_GOES_RT_integral_${type}_enhance_idsep_CLEAR.txt \
                CLEAR/output/idsep/GOES_RT_integral_${type}/

             date '+%Y-%m-%d %H:%M:%S'
             echo "[GOES-${n}/${type}] Batch opsep using curated CLEAR lists"
             python bin/fetchsep_prepare_obs \
                --StartDate "${start_date[GOES-${n}/${type}]}" --EndDate "${end_date[GOES-${n}/${type}]}" \
                --Experiment GOES_${n} --Spacecraft $type --FluxType integral --Threshold "30,1;50,1" \
                --BatchFile batch_event_list_GOES_${n}_integral_${type}_enhance_idsep_CLEAR.txt \
                --IDSEPEnhancement --Associations --StartPoint BATCH \
                >CLEAR/output/GOES_${n}_integral_${type}_batch.log
             echo
          done
       fi

       # Sandberg and Bruno calibrations available only for GOES-13 and GOES-15
       [[ $n == 13 || $n == 15 ]] && {
          date '+%Y-%m-%d %H:%M:%S'
          echo "[GOES-${n}/uncor_S14_B17 Copy curated batch files"
          cp fetchsep/reference/CLEAR/batch_event_list_GOES-${n}_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt \
             CLEAR/output/idsep/GOES-${n}_differential_uncor_S14_B17/

          date '+%Y-%m-%d %H:%M:%S'
          echo "[GOES-${n}/uncor_S14_B17] Batch opsep using curated CLEAR lists"
          python bin/fetchsep_prepare_obs \
             --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
             --Experiment GOES-${n} --FluxType differential --Threshold "30,1;50,1" \
             --BatchFile batch_event_list_GOES-${n}_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt \
             --IDSEPEnhancement --IDSEPSubtractBG --Associations --StartPoint BATCH --options "S14;Bruno2017;uncorrected" \
             >CLEAR/output/GOES-${n}_differential_uncor_S14_B17_bgsub_batch.log
          echo
       }
    done
fi

#Remove CLEAR config file so that running FetchSEP will not overwrite files in the
#CLEAR/ directory. Return to FetchSEP defaults.
echo "Returning config file to FetchSEP defaults."
rm fetchsep.cfg
python fetchsep/utils/config.py
date '+%Y-%m-%d %H:%M:%S'
echo "Completed"
