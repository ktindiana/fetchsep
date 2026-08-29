#!/bin/bash

#RUNNING THE SCRIPT
#------------------
#Run the script in the top directory of the repository, fetchsep/
#You may need to change your execution permissions:
#chmod u+x ./tests/test_prepsep.sh
#Run in the terminal as:
#./tests/test_prepsep.sh

step=${1:-"1"}
if [ -z "$1" ]; then
  echo "No argument provided. Starting at step ONE."
  step="1"
else
  step="$1"
fi

echo "Starting at step $1"


######################################################################
############ SET UP ENVIRONMENT AND CREATE DIRECTORIES ###############
######################################################################

#DO NOT CHANGE THESE FROM TEST DIRECTORIES
#DIRECTORIES WILL BE COMPLETELY DELETED IN THE NEXT STEP!!!!!
datapath="tests/data"
outpath="tests/output"
plotpath="tests/plots"
listpath="tests/lists"
targetpath="tests/target"


if [[ "${step}" = "1" ]]; then
    echo "Running test_prepsep.sh for step 1."
    date '+%Y-%m-%d %H:%M:%S'
    echo "Setting up environment"
    export PYTHONPATH="$PYTHONPATH:$PWD"

    echo "Removing test directories"
    rm -rf "${datapath}"
    rm -rf "${outpath}"
    rm -rf "${plotpath}"
    rm -rf "${listpath}"
    rm -rf "${targetpath}"

    echo "Creating test directories"
    mkdir "${datapath}"
    mkdir "${outpath}"
    mkdir "${plotpath}"
    mkdir "${listpath}"
    mkdir "${targetpath}"

    echo "Making directories and copying GOES data from tests/files/"
    mkdir "${datapath}"/GOES
    mkdir "${datapath}"/GOES/EPEAD
    mkdir "${datapath}"/GOES/HEPAD
    cp tests/files/data/GOES/EPEAD/* "${datapath}"/GOES/EPEAD/.
    cp tests/files/data/GOES/HEPAD/* "${datapath}"/GOES/HEPAD/.
    cp tests/files/data/GOES/fetchsep_data_manager.csv "${datapath}"/.

    echo "Running prepsep for step 1. Output will be generated while the pipeline is running."
    python bin/prepsep --StartDate 2012-05-12 --EndDate 2012-05-17 --Experiment GOES-13 --FluxType integral  --Threshold "30,1;50,1" --Associations --datapath "${datapath}" --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --TargetDir "${targetpath}"

    python bin/prepsep --StartDate 2012-05-17 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral  --Threshold "30,1;50,1" --Associations --datapath "${datapath}" --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --TargetDir "${targetpath}"

    echo "User must interact with prepsep output to approve an SEP event and run again with an argument of 2 to fully test."
    echo "Add SEP events to tests/target/approved)SEP.csv and rerun: ./tests/test_prepsep.sh 2"
    date '+%Y-%m-%d %H:%M:%S'

fi

if [[ "${step}" = "2" ]]; then
    echo "Running test_prepsep.sh for step 2."
    python bin/prepsep --StartDate 2012-05-17 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral  --Threshold "30,1;50,1" --Associations --datapath "${datapath}" --outpath "${outpath}" --plotpath "${plotpath}" --listpath "${listpath}" --TargetDir "${targetpath}"


    echo "Test completed."
    date '+%Y-%m-%d %H:%M:%S'

fi
