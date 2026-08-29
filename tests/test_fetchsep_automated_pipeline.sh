#!/bin/bash

#RUNNING THE SCRIPT
#------------------
#Run the script in the top directory of the repository, fetchsep/
#You may need to change your execution permissions:
#chmod u+x ./fetchsep/reference/CLEAR/deploy_CLEAR_Mac.sh
#Run in the terminal as:
#./tests/test_fetchsep_automated_pipeline.sh


######################################################################
############ SET UP ENVIRONMENT AND CREATE DIRECTORIES ###############
######################################################################
#In the top level fetchsep directory
date '+%Y-%m-%d %H:%M:%S'
echo "Setting up environment"
export PYTHONPATH="$PYTHONPATH:$PWD"

#DO NOT CHANGE THESE FROM TEST DIRECTORIES
#DIRECTORIES WILL BE COMPLETELY DELETED IN THE NEXT STEP!!!!!
datapath="tests/data"
outpath="tests/output"
plotpath="tests/plots"
listpath="tests/lists"

echo "Removing test directories"
rm -rf "${datapath}"
rm -rf "${outpath}"
rm -rf "${plotpath}"
rm -rf "${listpath}"

echo "Creating test directories"
mkdir "${datapath}"
mkdir "${outpath}"
mkdir "${plotpath}"
mkdir "${listpath}"

echo "Making directories and copying GOES data from tests/files/"
mkdir "${datapath}"/GOES
mkdir "${datapath}"/GOES/EPEAD
mkdir "${datapath}"/GOES/HEPAD
cp tests/files/data/GOES/EPEAD/* "${datapath}"/GOES/EPEAD/.
cp tests/files/data/GOES/HEPAD/* "${datapath}"/GOES/HEPAD/.
cp tests/files/data/GOES/fetchsep_data_manager.csv "${datapath}"/.

idsep_nsigma=3
init_win=150
sliding_win=5
percent_points=0.4
opsep_nsigma=3

echo "Running fetchsep_automated_pipeline. Output will be generated while the pipeline is running."
echo "Print messages are stored in ${outpath}/GOES-13_integral_batch_TEST.log"
echo "This test may take 5 minutes to run."
# "[GOES-13] Generating SEP events lists with opsep"
python bin/fetchsep_automated_pipeline --StartDate "2011-10-01 00:10:00" --EndDate 2012-06-01 --Experiment GOES-13 --FluxType integral  --Threshold "30,1;50,1" --RemoveAbove 10  --idsep_nsigma "${idsep_nsigma}" --init_win "${init_win}" --sliding_win "${sliding_win}" --percent_points "${percent_points}" --ReferenceEnergyBin "10,-1" --IDSEPEnhancement --Associations --datapath "${datapath}" --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --opsep_nsigma "${opsep_nsigma}" > "${outpath}"/GOES-13_integral_batch_TEST.log
echo "============================="
echo "Last messages printed to log:"
tail "${outpath}"/GOES-13_integral_batch_TEST.log
echo "============================="
echo "Test completed."
date '+%Y-%m-%d %H:%M:%S'
