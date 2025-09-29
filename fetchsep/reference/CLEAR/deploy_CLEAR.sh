#!/bin/bash

# The steps in this procedure will download all necessary data, process it, and
# create individual curated SEP lists for each GOES spacecraft. These lists may
# be combined to a final list by choosing the primary spacecraft for each time
# period.

#This script using the python call "python". You should modify this to
#python3 or python3.10 depending on the needs of your system.

#This script assumes the data directory is located in fetchsep/data. If you want
#your data directory elsewhere, then fetchsep/reference/CLEAR/fetchsep_CLEAR_GOES-07.cfg
#and fetchsep_CLEAR.cfg must be edited with the location of your desired data path
#(and output, plots, etc).

#Execute this script in the top directory of the repository, fetchsep/


# set up environment and create directories
export PYTHONPATH="$PYTHONPATH:$PWD"
python fetchsep/utils/config.py
#mkdir -p output/idsep/GOES-{06,07,08,10,11,13,15}_integral \
   output/idsep/GOES_RT_integral_{primary,secondary} \
   output/idsep/GOES-{13,15}_differential_uncor_S14_B17


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

for n in 06 07 08 10 11 13 15 RT; do
   # copy customized fetchsep settings
   [[ $n == 07 ]] && {
      cp fetchsep/reference/CLEAR/fetchsep_CLEAR_GOES-${n}.cfg ./fetchsep.cfg
   } || {
      cp fetchsep/reference/CLEAR/fetchsep_CLEAR.cfg ./fetchsep.cfg
   }

   if [[ $n != RT ]]; then
      echo "[GOES-${n}] Calculate background with idsep"
      python bin/idsep \
         --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
         --Experiment GOES-${n} --FluxType integral --RemoveAbove 10 --saveplot \
         >output/GOES-${n}_integral_idsep.log

      echo
      echo "[GOES-${n}] Copy curated batch files"
      cp fetchsep/reference/CLEAR/batch_event_list_GOES-${n}_integral_enhance_idsep_CLEAR.txt \
         output/idsep/GOES-${n}_integral/

      echo "[GOES-${n}] Batch opsep using curated CLEAR lists"
      python bin/fetchsep_prepare_obs \
         --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
         --Experiment GOES-${n} --FluxType integral --Threshold "30,1;50,1" \
         --BatchFile batch_event_list_GOES-${n}_integral_enhance_idsep_CLEAR.txt \
         --IDSEPEnhancement --StartPoint BATCH \
         >output/GOES-${n}_integral_batch.log
      echo
   else
      for type in primary secondary; do
         echo "[GOES-${n}/${type}] Calculate background with idsep"
         python bin/idsep \
            --StartDate "${start_date[GOES-${n}/${type}]}" --EndDate "${end_date[GOES-${n}/${type}]}" \
            --Experiment GOES_${n} --Spacecraft $type --FluxType integral --RemoveAbove 10 --saveplot \
            >output/GOES_${n}_integral_${type}_idsep.log

         echo
         echo "[GOES-${n}/${type}] Copy curated batch files"
         cp fetchsep/reference/CLEAR/batch_event_list_GOES_RT_integral_${type}_enhance_idsep_CLEAR.txt \
            output/idsep/GOES_RT_integral_${type}/

         echo "[GOES-${n}/${type}] Batch opsep using curated CLEAR lists"
         python bin/fetchsep_prepare_obs \
            --StartDate "${start_date[GOES-${n}/${type}]}" --EndDate "${end_date[GOES-${n}/${type}]}" \
            --Experiment GOES_${n} --Spacecraft $type --FluxType integral --Threshold "30,1;50,1" \
            --BatchFile batch_event_list_GOES_${n}_integral_${type}_enhance_idsep_CLEAR.txt \
            --IDSEPEnhancement --StartPoint BATCH \
            >output/GOES_${n}_integral_${type}_batch.log
         echo
      done
   fi

   # Sandberg and Bruno calibrations available only for GOES-13 and GOES-15
   [[ $n == 13 || $n == 15 ]] && {
      echo "[GOES-${n}/uncor_S14_B17] Calculate background with idsep"
      python bin/idsep \
         --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
         --Experiment GOES-${n} --FluxType differential --RemoveAbove 10 --saveplot \
         --options "S14;Bruno2017;uncorrected" \
         >output/GOES-${n}_differential_uncor_S14_B17_idsep.log

      echo
      echo "[GOES-${n}/uncor_S14_B17 Copy curated batch files"
      cp fetchsep/reference/CLEAR/batch_event_list_GOES-${n}_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt \
         output/idsep/GOES-${n}_differential_uncor_S14_B17/

      echo "[GOES-${n}/uncor_S14_B17] Batch opsep using curated CLEAR lists"
      python bin/fetchsep_prepare_obs \
         --StartDate "${start_date[GOES-${n}]}" --EndDate "${end_date[GOES-${n}]}" \
         --Experiment GOES-${n} --FluxType differential --Threshold "30,1;50,1" \
         --BatchFile batch_event_list_GOES-${n}_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt \
         --IDSEPEnhancement --IDSEPSubtractBG --StartPoint BATCH --options "S14;Bruno2017;uncorrected" \
         >output/GOES-${n}_differential_uncor_S14_B17_bgsub_batch.log
      echo
   }
done
