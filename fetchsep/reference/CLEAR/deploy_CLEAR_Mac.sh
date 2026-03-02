#!/bin/bash

#PLEASE READ IN ENTIRETY BEFORE RUNNING

#STEP-BY-STEP PROCEDURE TO GENERATE THE CLEAR BENCHMARK DATASET for Mac
#The benchmark dataset is created using FetchSEP run in python3.10.17.

#REPLACE THE "python" CALL FOR YOUR REQUIRED COMMAND, e.g. py or python3.10

#The steps in this procedure will download all necessary data, process it,
#and create individual curated SEP lists for each GOES spacecraft.
#These lists may be combined to a final list by choosing the primary
#spacecraft for each time period.

#SPECIFIC PYTHON CALL
#--------------------
#The steps provided here work on Mac where the command to call python is for python3.10.
#Best to work in a virtual environment for python3.10.

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
#chmod u+x ./fetchsep/reference/CLEAR/deploy_CLEAR_Mac.sh
#Run in the terminal as:
#./fetchsep/reference/CLEAR/deploy_CLEAR_Mac.sh

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
#In the top level fetchsep directory
date '+%Y-%m-%d %H:%M:%S'
echo "Setting up environment"
export PYTHONPATH="$PYTHONPATH:$PWD"

outpath="./CLEAR/output"
plotpath="./CLEAR/plots"
listpath="./CLEAR/lists"

python fetchsep/utils/config.py --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}"

idsep_nsigma=3
init_win=150
sliding_win=5
percent_points=0.4
opsep_nsigma=3

######################################################################
################# CALCULATE BACKGROUNDS WITH IDSEP ###################
######################################################################

if [[ "${startpoint}" = "ALL" ]]; then
    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-06] Calculate background with idsep"
    python bin/idsep --StartDate 1986-01-01 --EndDate 1994-12-01 --Experiment GOES-06 --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-06_integral_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-07] Calculate background with idsep"
    python bin/idsep --StartDate 1987-03-01 --EndDate 1996-09-01 --Experiment GOES-07 --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-07_integral_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-08] Calculate background with idsep"
    python bin/idsep --StartDate 1995-01-01 --EndDate 2003-06-18 --Experiment GOES-08 --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-08_integral_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-10] Calculate background with idsep"
    python bin/idsep --StartDate 1998-07-09 --EndDate 2004-07-01 --Experiment GOES-10 --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-10_integral_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-11] Calculate background with idsep"
    python bin/idsep --StartDate 2003-06-01 --EndDate 2011-03-01 --Experiment GOES-11 --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-11_integral_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-13] Calculate background with idsep"
    python bin/idsep --StartDate 2010-05-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-13_integral_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-15] Calculate background with idsep"
    python bin/idsep --StartDate 2011-01-01 --EndDate 2020-03-05 --Experiment GOES-15 --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-15_integral_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-RT/primary] Calculate background with idsep"
    python bin/idsep --StartDate 2020-03-08 --EndDate 2025-09-11 --Experiment GOES-RT --Spacecraft primary --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-RT_integral_primary_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-RT/secondary] Calculate background with idsep"
    python bin/idsep --StartDate "2021-09-23 15:45:00" --EndDate 2025-09-11 --Experiment GOES-RT --Spacecraft secondary --FluxType integral --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-RT_integral_secondary_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-13 energy bin calibrated] Calculate background with idsep"
    python bin/idsep --StartDate 2010-05-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType differential --options "S14;Bruno2017;uncorrected" --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-13_differential_uncor_S14_B17_idsep.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-15 energy bin calibrated] Calculate background with idsep"
    python bin/idsep --StartDate 2011-01-01 --EndDate 2020-03-05 --Experiment GOES-15 --FluxType differential --options "S14;Bruno2017;uncorrected" --RemoveAbove 10 --saveplot --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" > "${outpath}"/GOES-15_differential_uncor_S14_B17_idsep.log

fi

if [[ "${startpoint}" = "ALL" ]] || [[ "${startpoint}" = "LISTS" ]]; then
    ######################################################################
    ################### MOVE CURATED BATCH FILES #########################
    ######################################################################
    date '+%Y-%m-%d %H:%M:%S'
    echo "Copying CLEAR curated batch event lists to proper directories."
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-06_integral_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-06_integral/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-07_integral_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-07_integral/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-08_integral_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-08_integral/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-10_integral_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-10_integral/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-11_integral_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-11_integral/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-13_integral_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-13_integral/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-15_integral_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-15_integral/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-RT_integral_primary_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-RT_integral_primary/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-RT_integral_secondary_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-RT_integral_secondary/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-13_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-13_differential_uncor_S14_B17/.
    cp fetchsep/reference/CLEAR/batch_event_list_GOES-15_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt "${outpath}"/idsep/GOES-15_differential_uncor_S14_B17/.

fi

if [[ "${startpoint}" = "ALL" ]] || [[ "${startpoint}" = "LISTS" ]] || [[ "${startpoint}" = "BATCH" ]]; then
    ######################################################################
    ############## BATCH OPSEP USING CURATED CLEAR LISTS #################
    ######################################################################
    # MAKE SURE TO INCLUDE --StartPoint BATCH or the batch file will be
    #overwritten by a new run of idsep within fetchsep_automated_pipeline
    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-06] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 1986-01-01 --EndDate 1994-12-01 --Experiment GOES-06 --FluxType integral --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-06_integral_enhance_idsep_CLEAR.txt  --IDSEPEnhancement --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 1 --NoGOESColors > "${outpath}"/GOES-06_integral_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-07] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 1987-03-01 --EndDate 1996-09-01 --Experiment GOES-07 --FluxType integral  --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-07_integral_enhance_idsep_CLEAR.txt  --IDSEPEnhancement --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 2 --NoGOESColors > "${outpath}"/GOES-07_integral_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-08] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 1995-01-01 --EndDate 2003-06-18 --Experiment GOES-08 --FluxType integral  --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-08_integral_enhance_idsep_CLEAR.txt  --IDSEPEnhancement --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 3 --NoGOESColors > "${outpath}"/GOES-08_integral_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-10] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 1998-07-09 --EndDate 2004-07-01 --Experiment GOES-10 --FluxType integral  --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-10_integral_enhance_idsep_CLEAR.txt  --IDSEPEnhancement --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 4 --NoGOESColors > "${outpath}"/GOES-10_integral_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-11] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 2003-06-01 --EndDate 2011-03-01 --Experiment GOES-11 --FluxType integral  --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-11_integral_enhance_idsep_CLEAR.txt  --IDSEPEnhancement --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 5 --NoGOESColors > "${outpath}"/GOES-11_integral_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-13] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 2010-05-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType integral  --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-13_integral_enhance_idsep_CLEAR.txt  --IDSEPEnhancement --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 9 --NoGOESColors > "${outpath}"/GOES-13_integral_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-15] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 2011-01-01 --EndDate 2020-03-05 --Experiment GOES-15 --FluxType integral  --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-15_integral_enhance_idsep_CLEAR.txt  --IDSEPEnhancement --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 10 --NoGOESColors > "${outpath}"/GOES-15_integral_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-RT/primary] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 2020-03-08 --EndDate 2025-09-11 --Experiment GOES-RT --Spacecraft primary --FluxType integral --BatchFile batch_event_list_GOES-RT_integral_primary_enhance_idsep_CLEAR.txt --IDSEPEnhancement --Associations --Threshold "30,1;50,1"  --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}" > "${outpath}"/GOES-RT_integral_primary_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-RT/secondary] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate "2021-09-23 15:45:00" --EndDate 2025-09-11 --Experiment GOES-RT --Spacecraft secondary --FluxType integral  --BatchFile batch_event_list_GOES-RT_integral_secondary_enhance_idsep_CLEAR.txt --IDSEPEnhancement --Threshold "30,1;50,1" --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}" > "${outpath}"/GOES-RT_integral_secondary_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-13 energy bin calibrated] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 2010-05-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType differential --options "S14;Bruno2017;uncorrected" --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-13_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt --IDSEPEnhancement --IDSEPSubtractBG --Associations --StartPoint BATCH  --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 1 > "${outpath}"/GOES-13_differential_uncor_S14_B17_bgsub_batch.log

    date '+%Y-%m-%d %H:%M:%S'
    echo "[GOES-15 energy bin calibrated] Generating SEP events lists with opsep"
    python bin/fetchsep_automated_pipeline --StartDate 2011-01-01 --EndDate 2020-03-05 --Experiment GOES-15 --FluxType differential --options "S14;Bruno2017;uncorrected"  --Threshold "30,1;50,1" --BatchFile batch_event_list_GOES-15_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt --IDSEPEnhancement --IDSEPSubtractBG --Associations --StartPoint BATCH --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}"  --ColorScheme 2 > "${outpath}"/GOES-15_differential_S14_B17_uncor_bgsub_batch.log

fi

######################################################################
################ COMPILE CURATED LISTS INTO PRIMARY ##################
######################################################################
# MAKE SURE TO INCLUDE --StartPoint COMPILE if all batch runs of OPSEP
# are already complete and you want to extract to a single list.

#Create a single SEP list by extracting SEP events for the primary GOES satellite at the time
date '+%Y-%m-%d %H:%M:%S'
echo "[GOES PRIMARY list] Creating single list of SEP events selecting the primary GOES spacecraft at the time"
find "${outpath}"/opsep/* -name "*_integral_*_sep_events.csv" > CLEARlists.txt
python bin/make_primary_goes_list --Prefix GOES_integral --Filename CLEARlists.txt
rm CLEARlists.txt

date '+%Y-%m-%d %H:%M:%S'
echo "[GOES PRIMARY energy bin calibrated list] Creating single list of SEP events selecting the primary GOES spacecraft at the time"
find "${outpath}"/opsep/* -name "*_differential_*_sep_events.csv" > CLEARlists.txt
python bin/make_primary_goes_list --Prefix GOES_differential_energy_bin_calibrated --Filename CLEARlists.txt
rm CLEARlists.txt

date '+%Y-%m-%d %H:%M:%S'
echo "Completed"
