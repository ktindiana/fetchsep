@ECHO OFF
setlocal EnableDelayedExpansion

REM The CLEAR Benchmark Dataset was generated with python 3.10.17.

REM SCRIPT TO GENERATE THE BENCHMARK DATASET FROM SCRATCH
REM REPLACE THE "python" CALL TO REFLECT YOUR REQUIRED COMMAND.

REM The steps in this procedure will download all necessary data, process it, and
REM create individual curated SEP lists for each GOES spacecraft. These lists may
REM be combined to a final list by choosing the primary spacecraft for each time
REM period.

REM This script using the python call "python". You should modify this to
REM python3 or python3.10 or py depending on the needs of your system.

REM This script assumes the data directory is located in fetchsep/data. If you want
REM your data directory elsewhere, then fetchsep/reference/CLEAR/fetchsep_CLEAR.cfg
REM must be edited with the location of your desired data path (and outpath, plotpath).
REM Output files are stored in fetchsep/CLEAR/output and fetchsep/CLEAR/plots.

REM Execute this script in the top directory of the repository, fetchsep/


REM set up environment and create directories
$env:PYTHONPATH = "$env:PYTHONPATH;$PWD"
copy .\fetchsep\reference\CLEAR\fetchsep_CLEAR.cfg .\fetchsep.cfg
python .\fetchsep\utils\config.py


set start_date.GOES-06=1986-01-01
set start_date.GOES-07=1987-03-01
set start_date.GOES-08=1995-01-01
set start_date.GOES-10=1998-07-09
set start_date.GOES-11=2003-06-01
set start_date.GOES-13=2010-05-01
set start_date.GOES-15=2011-01-01
set start_date.GOES-RT.primary=2020-03-08
set "start_date.GOES-RT.secondary=2021-09-23 15:45:00"

set end_date.GOES-06=1994-12-01
set end_date.GOES-07=1996-09-01
set end_date.GOES-08=2003-06-18
set end_date.GOES-10=2004-07-01
set end_date.GOES-11=2011-03-01
set end_date.GOES-13=2018-01-01
set end_date.GOES-15=2020-03-05
set end_date.GOES-RT.primary=2025-09-11
set end_date.GOES-RT.secondary=2025-09-11


for %%G in (06 07 08 10 11 13 15 RT) do (
   if not "%%G" == "RT" (
      echo "[GOES-%%G] Calculate background with idsep"
      python .\bin\idsep ^
         --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
         --Experiment GOES-%%G --FluxType integral --RemoveAbove 10 --saveplot ^
         > .\CLEAR\output\GOES-%%G_integral_idsep.log

      echo
      echo "[GOES-%%G] Copy curated batch files"
      copy .\fetchsep\reference\CLEAR\batch_event_list_GOES-%%G_integral_enhance_idsep_CLEAR.txt ^
         .\CLEAR\output\idsep\GOES-%%G_integral\

      echo "[GOES-%%G] Batch opsep using curated CLEAR lists"
      python .\bin\fetchsep_prepare_obs ^
         --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
         --Experiment GOES-%%G --FluxType integral --Threshold "30,1;50,1" ^
         --BatchFile batch_event_list_GOES-%%G_integral_enhance_idsep_CLEAR.txt ^
         --IDSEPEnhancement --StartPoint BATCH ^
         > .\CLEAR\output\GOES-%%G_integral_batch.log
      echo
   ) else (
      for %%T in (primary secondary) do (
         echo "[GOES-%%G/%%T] Calculate background with idsep"
         python .\bin\idsep ^
            --StartDate "!start_date.GOES-%%G.%%T!" --EndDate "!end_date.GOES-%%G.%%T!" ^
            --Experiment GOES_%%G --Spacecraft %%T --FluxType integral --RemoveAbove 10 --saveplot ^
            > .\CLEAR\output\GOES_%%G_integral_%%T_idsep.log

         echo
         echo "[GOES-%%G/%%T] Copy curated batch files"
         copy .\fetchsep\reference\CLEAR\batch_event_list_GOES_RT_integral_%%T_enhance_idsep_CLEAR.txt ^
            .\CLEAR\output\idsep\GOES_RT_integral_%%T\

         echo "[GOES-%%G/%%T] Batch opsep using curated CLEAR lists"
         python .\bin\fetchsep_prepare_obs ^
            --StartDate "!start_date.GOES-%%G.%%T!" --EndDate "!end_date.GOES-%%G.%%T!" ^
            --Experiment GOES_%%G --Spacecraft %%T --FluxType integral --Threshold "30,1;50,1" ^
            --BatchFile batch_event_list_GOES_%%G_integral_%%T_enhance_idsep_CLEAR.txt ^
            --IDSEPEnhancement --StartPoint BATCH ^
            > .\CLEAR\output\GOES_%%G_integral_%%T_batch.log
         echo
      )
   )
)

REM Sandberg and Bruno calibrations available only for GOES-13 and GOES-15
for %%G in (13 15) do (
  echo "[GOES-%%G/uncor_S14_B17] Calculate background with idsep"
  python .\bin\idsep ^
	 --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
	 --Experiment GOES-%%G --FluxType differential --RemoveAbove 10 --saveplot ^
	 --options "S14;Bruno2017;uncorrected" ^
	 > .\CLEAR\output\GOES-%%G_differential_uncor_S14_B17_idsep.log

  echo
  echo "[GOES-%%G/uncor_S14_B17 Copy curated batch files"
  copy .\fetchsep\reference\CLEAR\batch_event_list_GOES-%%G_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt ^
	 .\CLEAR\output\idsep\GOES-%%G_differential_uncor_S14_B17\

  echo "[GOES-%%G/uncor_S14_B17] Batch opsep using curated CLEAR lists"
  python .\bin\fetchsep_prepare_obs ^
	 --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
	 --Experiment GOES-%%G --FluxType differential --Threshold "30,1;50,1" ^
	 --BatchFile batch_event_list_GOES-%%G_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt ^
	 --IDSEPEnhancement --IDSEPSubtractBG --StartPoint BATCH --options "S14;Bruno2017;uncorrected" ^
	 > .\CLEAR\output\GOES-%%G_differential_uncor_S14_B17_bgsub_batch.log
  echo

)
