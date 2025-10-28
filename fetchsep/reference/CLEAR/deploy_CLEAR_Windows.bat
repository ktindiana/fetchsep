@ECHO OFF
setlocal EnableDelayedExpansion

REM PLEASE READ IN ENTIRETY BEFORE RUNNING

REM STEP-BY-STEP PROCEDURE TO GENERATE THE CLEAR BENCHMARK DATASET for Windows
REM The benchmark dataset is created using FetchSEP run in python3.10.17.

REM REPLACE THE "python" CALL FOR YOUR REQUIRED COMMAND, if necessary.

REM The steps in this procedure will download all necessary data, process it,
REM and create individual curated SEP lists for each GOES spacecraft.
REM These lists may be combined to a final list by choosing the primary
REM spacecraft for each time period.

REM SPECIFIC PYTHON CALL
REM --------------------
REM The steps provided here work on Windows where the command to call python is python.
REM The python command in the script may need to be modified for your specific call,
REM e.g. python3.10 or py

REM CONFIG FILE
REM -----------
REM !!!!If you already have a fetchsep.cfg file in the main fetchsep directory
REM that you have customized, copy it to a different name (e.g. fetchsep.cfg.bak)
REM so that it will not be overwritten. !!!!

REM This script assumes the data directory is located in fetchsep\data and all GOES data
REM will be downloaded here. If you want your data directory elsewhere, then edit
REM fetchsep\reference\CLEAR\fetchsep_CLEAR.cfg
REM with the location of your desired data path. Do not modify outpath and plotpath.

REM Output files are stored in fetchsep\CLEAR\output and fetchsep\CLEAR\plots and should
REM not be modified to ensure the script works properly.

REM RUNNING THE SCRIPT
REM ------------------
REM Run the script in the top directory of the repository, fetchsep\
REM You may need to change your execution permissions for:
REM .\fetchsep\reference\CLEAR\deploy_CLEAR_Windows.bat
REM Run in the terminal as:
REM .\fetchsep\reference\CLEAR\deploy_CLEAR_Windows.bat

REM COMMAND LINE ARGUMENTS
REM ----------------------
REM In the command line, the argument LISTS may be added to skip the calculation of
REM the mean background with IDSEP IF IT WAS ALREADY DONE PREVIOUSLY.
REM The background calculation is extremely time consuming and it needn't be repeated.
REM If the user only wants to regenerate the SEP lists,
REM the argument LISTS will skip to that part of the script.

REM When run without an argument, will default to startpoint ALL and the script will:
REM     - Set up environment
REM     - Download data and calculate mean background and sigma with IDSEP (12+ hours)
REM     - Copy CLEAR curated batch files
REM     - Generate SEP lists (2+ hours)

REM When run with the argument LISTS, the script will:
REM     - Set up environment
REM     - Copy CLEAR curated batch files
REM     - Generate SEP lists (2+ hours)


SET "startpoint=ALL"
IF NOT "%1"=="" (
	SET "startpoint=%1"
)
echo "Starting at %startpoint%"



REM #####################################################################
REM ########### SET UP ENVIRONMENT AND CREATE DIRECTORIES ###############
REM #####################################################################
Get-Date -Format "yyyy-MM-dd HH:mm:ss"
echo "Setting up environment"
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


REM #####################################################################
REM ################ CALCULATE BACKGROUNDS WITH IDSEP ###################
REM #####################################################################
if "%startpoint%" == "ALL" (
	for %%G in (06 07 08 10 11 13 15 RT) do (
	   if not "%%G" == "RT" (
		  Get-Date -Format "yyyy-MM-dd HH:mm:ss"
		  echo "[GOES-%%G] Calculate background with idsep"
		  python .\bin\idsep ^
			 --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
			 --Experiment GOES-%%G --FluxType integral --RemoveAbove 10 --saveplot ^
			 > .\CLEAR\output\GOES-%%G_integral_idsep.log

		  echo
	   ) else (
		  for %%T in (primary secondary) do (
		     Get-Date -Format "yyyy-MM-dd HH:mm:ss"
			 echo "[GOES-%%G/%%T] Calculate background with idsep"
			 python .\bin\idsep ^
				--StartDate "!start_date.GOES-%%G.%%T!" --EndDate "!end_date.GOES-%%G.%%T!" ^
				--Experiment GOES_%%G --Spacecraft %%T --FluxType integral --RemoveAbove 10 --saveplot ^
				> .\CLEAR\output\GOES_%%G_integral_%%T_idsep.log

			 echo
		  )
	   )
	)

	REM Sandberg and Bruno calibrations available only for GOES-13 and GOES-15
	for %%G in (13 15) do (
	  Get-Date -Format "yyyy-MM-dd HH:mm:ss"
	  echo "[GOES-%%G/uncor_S14_B17] Calculate background with idsep"
	  python .\bin\idsep ^
		 --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
		 --Experiment GOES-%%G --FluxType differential --RemoveAbove 10 --saveplot ^
		 --options "S14;Bruno2017;uncorrected" ^
		 > .\CLEAR\output\GOES-%%G_differential_uncor_S14_B17_idsep.log

	  echo

	)
)

REM #####################################################################
REM ############# BATCH OPSEP USING CURATED CLEAR LISTS #################
REM #####################################################################
REM MAKE SURE TO INCLUDE --StartPoint BATCH or the batch file will be
REM overwritten by a new run of idsep within fetchsep_prepare_obs

for %%G in (06 07 08 10 11 13 15 RT) do (
   if not "%%G" == "RT" (
      Get-Date -Format "yyyy-MM-dd HH:mm:ss"
      echo "[GOES-%%G] Copy curated batch files"
      copy .\fetchsep\reference\CLEAR\batch_event_list_GOES-%%G_integral_enhance_idsep_CLEAR.txt ^
         .\CLEAR\output\idsep\GOES-%%G_integral\

      Get-Date -Format "yyyy-MM-dd HH:mm:ss"
      echo "[GOES-%%G] Batch opsep using curated CLEAR lists"
      python .\bin\fetchsep_prepare_obs ^
         --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
         --Experiment GOES-%%G --FluxType integral --Threshold "30,1;50,1" ^
         --BatchFile batch_event_list_GOES-%%G_integral_enhance_idsep_CLEAR.txt ^
         --IDSEPEnhancement --Associations --StartPoint BATCH ^
         > .\CLEAR\output\GOES-%%G_integral_batch.log
      echo
   ) else (
      for %%T in (primary secondary) do (
		 Get-Date -Format "yyyy-MM-dd HH:mm:ss"
         echo "[GOES-%%G/%%T] Copy curated batch files"
         copy .\fetchsep\reference\CLEAR\batch_event_list_GOES_RT_integral_%%T_enhance_idsep_CLEAR.txt ^
            .\CLEAR\output\idsep\GOES_RT_integral_%%T\

	     Get-Date -Format "yyyy-MM-dd HH:mm:ss"
         echo "[GOES-%%G/%%T] Batch opsep using curated CLEAR lists"
         python .\bin\fetchsep_prepare_obs ^
            --StartDate "!start_date.GOES-%%G.%%T!" --EndDate "!end_date.GOES-%%G.%%T!" ^
            --Experiment GOES_%%G --Spacecraft %%T --FluxType integral --Threshold "30,1;50,1" ^
            --BatchFile batch_event_list_GOES_%%G_integral_%%T_enhance_idsep_CLEAR.txt ^
            --IDSEPEnhancement --Associations --StartPoint BATCH ^
            > .\CLEAR\output\GOES_%%G_integral_%%T_batch.log
         echo
      )
   )
)

REM Sandberg and Bruno calibrations available only for GOES-13 and GOES-15
for %%G in (13 15) do (
  Get-Date -Format "yyyy-MM-dd HH:mm:ss"
  echo "[GOES-%%G/uncor_S14_B17 Copy curated batch files"
  copy .\fetchsep\reference\CLEAR\batch_event_list_GOES-%%G_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt ^
	 .\CLEAR\output\idsep\GOES-%%G_differential_uncor_S14_B17\

  Get-Date -Format "yyyy-MM-dd HH:mm:ss"
  echo "[GOES-%%G/uncor_S14_B17] Batch opsep using curated CLEAR lists"
  python .\bin\fetchsep_prepare_obs ^
	 --StartDate "!start_date.GOES-%%G!" --EndDate "!end_date.GOES-%%G!" ^
	 --Experiment GOES-%%G --FluxType differential --Threshold "30,1;50,1" ^
	 --BatchFile batch_event_list_GOES-%%G_differential_uncor_S14_B17_bgsub_enhance_idsep_CLEAR.txt ^
	 --IDSEPEnhancement --IDSEPSubtractBG --Associations --StartPoint BATCH --options "S14;Bruno2017;uncorrected" ^
	 > .\CLEAR\output\GOES-%%G_differential_uncor_S14_B17_bgsub_batch.log
  echo

)

REM Remove CLEAR config file so that running FetchSEP will not overwrite files in the
REM CLEAR/ directory. Return to FetchSEP defaults.

echo "Returning config file to FetchSEP defaults."
del fetchsep.cfg
python .\fetchsep\utils\config.py
Get-Date -Format "yyyy-MM-dd HH:mm:ss"
echo "Completed"
