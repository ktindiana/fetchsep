========
FetchSEP
========


.. image:: https://img.shields.io/pypi/v/fetchsep.svg
        :target: https://pypi.python.org/pypi/fetchsep

.. image:: https://img.shields.io/travis/ktindiana/fetchsep.svg
        :target: https://travis-ci.com/ktindiana/fetchsep

.. image:: https://readthedocs.org/projects/fetchsep/badge/?version=latest
        :target: https://fetchsep.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Identify SEP elevations above background in a time series (`idsep`) and analyze events individually (`opsep`)


* Free software: MIT license
* This documentation is DEPRECATED. WILL UPDATE. Documentation: https://fetchsep.readthedocs.io.
* Best documentation is listed below in this README.

Description
===========

FetchSEP analyzes particle fluxes (tested mainly with protons) for solar energetic particle (SEP) events. FetchSEP is composed of two main packages:

* `idsep`: Calculates the mean background and expected level of variation (sigma) and produces rough SEP event lists identifying enhanced time periods.
* `opsep`: Analyzes individual SEP events and calculates timing, flux, and fluence values.

The mean background and sigma calculated for each day in `idsep` may be used by `opsep` for background subtraction and identification of enhancements above background when processing individual SEP events.

FetchSEP can run `idsep` followed by `opsep` to identify and analyze all SEP events in a time series. The results of this automated process will not be perfect as some identified enhanced periods may consist of multiple SEP events in succession. A small amount of human intervention to appropriately split those time periods can lead to an automatically generated high-quality list of SEP events and their characteristics. More guidance provided in the fetchsep_prepare_obs section.

FetchSEP may be used to download particle data for supported spacecraft and output the fluxes into a simple, easy-to-read csv file.

FetchSEP supported data:

* GOES-05 to GOES-15 (excluding GOES-09) - fluxes are taken from the west-facing detector
* GOES-R real-time integral (CCMC iSWA archive) and differential (NOAA NCEI) fluxes
* SOHO/EPHIN 4 energy bins and SOHO/EPHIN REleASE high resolution energy bins
* SOHO/ERNE
* STEREO-A and STEREO-B
* ACE/SIS
* ACE/EPAM electrons (most energetic)
* IMP-8/CPME
* SEPEM v2 and v3 (if the user has the file)
* CalGOES (NASA JSC SRAG dataset, if the user has the file)
* Users may input their own time series files in the appropriate format

Set Up
======

FetchSEP and requirements are defined for python version 3.10. Recommend setting up a virtual environment then install libraries in requirements.txt with pip.

Add the fetchsep path to your python path. In the top level fetchsep
path, fetchsep/, run:

Windows:

    | $env:PYTHONPATH = "$env:PYTHONPATH;$PWD"
    
Mac:

    | export PYTHONPATH="$PYTHONPATH:$PWD"

Export the fetchsep.cfg file from fetchsep > utils > fetchsep.cfg to
the fetchsep path:

    | python fetchsep/utils/config.py
or

    | python bin/opsep --ExportConfig

In fetchsep.cfg, edit the location of the data, output, plots, and lists directories as needed. By default they will be in fetchsep directory. The data/ directory will hold all of the spacecraft data downloaded by fetchsep.

By default, `idsep` and `opsep` will create the necessary output directories in the current working directory where the command is executed.  The directories that will be created are `data`, `ouptut`, `plots`, and `lists`.  The directories will not be overwritten if they already exist.  If you would like to choose another output location, generate a config file (`fetchsep.cfg`) with `opsep --ExportConfig` and edit the paths described there.

A configuration file may also be placed in your home directory with the name `.fetchsep`.  The configuration file does not need to be complete; you may specify only the values for which you wish to override the default.  The order of config value precidence is 1. current working directory `fetchsep.cfg`, 2. home directory `.fetchsep`, 3. fetchsep defaults.

Download Data Only
==================

FetchSEP includes many "native" experiment datasets and will continue to expand its library. For data handling, it had been built to:

* Pull data from its original source online and save to your computer in the data/ folder
* Understand file formats to extract particle fluxes
* Manage occasional problems encountered in data files, including incomplete files or missing data
* Manage different data versions
* For GOES-13, -14, and -15, choose the West-facing detector by referring to the orientation files provided by NOAA
* Output particle time series into a single, easy-to-read csv file
* Once a file has been successfully downloaded and read in, a file manager keeps track of files that are complete and do not have to be downloaded again

Users might like to download original data and convert to csv files without any of the additional processing performed in FetchSEP. This may be done with the download function:

    | python bin/download --StartDate 2025-01-01 --EndDate 2025-02-01 --Experiment GOES-18 --FluxType differential --showplot

The example downloads the daily netcdf files for GOES-18 differential flux, reads them in, plots them, and exports a single csv file of the full time series. A folder

    | output/idsep/GOES_differential

is created containing the time series file

    | fluxes_GOES-18_differential_20250101_20250131.csv



IDSEP
=====

The `idsep` code will read in a long time series and automatically identify increases above background. This is done by estimating a mean background level plus an expected level of variation (sigma). All flux less than mean + 3sigma are considered background while all points above mean + 3sigma are identified as increases.

The mean background solution and sigma are output for every timestamp. These values may be used in `opsep` to identify enhancements above background or perform background subtraction.

By assigning a set of criteria, enhancements that are most likely due to SEP events are identified and an SEP event list is output for each energy channel.

The code also outputs a file containing every single high flux point above the mean + 3sigma threshold.

Note that in fetchsep.cfg, number of sigma, the initial window used to estimate background levels, and the final sliding window used to estimate background levels can be adjusted.

    | idsep_nsigma = 3
    | init_win = 150 #days to average initial estimate of threshold
    | sliding_win = 5 #days in sliding window to calculate final threshold
    | percent_points = 0.4 #Percent of points that must be in the sliding
    |                #window to calculate the background; otherwise use
    |                #previous good value


To run:

    | python bin/idsep --StartDate 2017-01-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType integral --RemoveAbove 10 --saveplot


Note that the command above will download all GOES-13 data from 2017-01-01 to 2018-01-01 to your data/ directory then performan an analysis to identify the mean background and SEP enhancements.

An analysis of a year of GOES data will take around 20 minutes. The time period specified for analysis must exceed init_win number of days. If it is shorter, `idsep` will automatically extend the analysis time window.

For more features, run:

    | python bin/idsep --help


OPSEP
=====

The `opsep` code was previously supported at https://github.com/ktindiana/operational-sep and is now transitioned to this package going forward.

`opsep` is intended to assess each individual SEP event at a time, extracting information such as start and end times, peak fluxes, and event fluence.

**Two operational SEP event definitions are applied automatically for:**

* >10 MeV exceeds 10 pfu
* >100 MeV exceeds 1 pfu

If differential fluxes were input into opsep, it will automatically estimate >10 MeV and >100 MeV fluxes from the available energy channels.

The code will output various csv files and a json file with accompany txt files. The JSON file is in the same format as required by the SEP Scoreboard to submit forecasts.

`opsep` creates files from observations that can be directly compared to SEP model forecasts sent to the SEP Scoreboard.

For time profile SEP models, `opsep` may be used to create the JSON files that can be submitted to the SEP Scoreboard.


To run OpSEP to process individual SEP events:

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --showplot


You may add additional event definitions with the --Threshold flag.

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --Threshold "30,1;50,1" --showplot

You may search for associated flare, CME, radio, etc information by including the --Associations flag. If the SEP event is in in fetchsep/reference/SRAG_SEP_List_R11_CLEARversion.csv, the flare, etc, information will be saved to the output json and csv files.

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --Threshold "30,1;50,1" --Associations --showplot

To see how to add thresholds to differential energy channels, run:

    | python bin/opsep --help


Running `opsep` (or `idsep`) for your own time series
----------------------------------------

Users may input their own time series into `opsep` by specifying some
information in the config file:

    | ##### DELIMETER between columns of file with time series
    | user_delim = " "  #any string
    | ##### COLUMNS containing the fluxes you want to analyze
    | user_col = arr.array('i',[1,2,3,4,5,6,7,8])
    | err_col = arr.array('i',[]) #set to [] if no uncertainties. err_col only used by idsep
    | ##### ENERGY BINS associated with user file and columns
    | #For differential bins, use the format:
    | user_energy_bins = [[Elow1,Ehigh1],[Elow2,Ehigh2],etc]
    | #For integral bins, use the format:
    | user_energy_bins = [[Elow1,-1],[Elow2,-1],[Elow3,-1],etc]

For a time profile produced by a model (specify differential or integral flux as appropriate), run as:

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment user --FluxType differential --ExperimentName MyModel --UserFile my/model/timeprofile.txt --JSONType model --showplot

For a time profile produced by a satellite or experiment (specify differential or integral flux as appropriate), run as:

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment user --FluxType differential --ExperimentName MyData --UserFile my/data/timeprofile.txt --JSONType observations --showplot

Perform Background Subtraction and Identify Enhancements Above Background
=========================================================================

Background-subtraction of particle fluxes may be performed in two different ways in FetchSEP.

With OpSEP-calculated Background
--------------------------------

With `opsep`: The user may specify a specific time period to use as the background. OpSEP will calculate the mean particle flux and level of variation (sigma) for that time period. Fluxes above mean + n*sigma, where n is specified in fetchsep.cfg in the opsep_nsigma variable, are considered SEP fluxes and will be subtracted by the mean. Fluxes below mean+n*sigma are consered background and are set to zero.

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType differential --Threshold "30,1;50,1" --OPSEPSubtractBG --BGStartDate 2012-05-10 --BGEndDate 2012-05-12 --showplot

If background subtraction is performed, SEP events are automatically identified above background. An event definition with a flux threshold of 1e-6 (e.g. >10 MeV exceeds 1e-6 pfu) indicates that the SEP event was identified for all fluxes above background.

The --OPSEPEnhancement flag may be used to identify and analyze SEP events above background WITHOUT background subtraction. An event definition with a flux threshold of 1e-6 (e.g. >10 MeV exceeds 1e-6 pfu) indicates that the SEP event was identified for all fluxes above background.

GOES integral fluxes provided by NOAA already have some amount of background-subtraction. SEP enhancements above background may be evaluated without performing a background subtraction by:

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --Threshold "30,1;50,1" --OPSEPEnhancement --BGStartDate 2012-05-10 --BGEndDate 2012-05-12 --showplot

Although the background has not been subtracted, calling --IDSEPEnhancement will set values below mean + n*sigma to zero.

With IDSEP-calculated Background
--------------------------------

With `idsep`: Run `idsep` for a long timeframe (e.g. months, years) to calculate the mean background and sigma with time. Run `opsep` and use the background solution created by idsep to identify SEP enhancements above mean + n*sigma, where n is specified in fetchsep.cfg in the opsep_nsigma variable, subtract the mean background, and set background fluxes to zero. You may calculate the mean background with idsep for, e.g.,  the entire history of an experiment and keep that file around for use in background subtraction with opsep.

    | python bin/idsep --StartDate 2011-05-01 --EndDate 2013-01-01 --Experiment GOES-13 --FluxType differential --RemoveAbove 10 --showplot

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType differential --Threshold "30,1;50,1" --IDSEPSubtractBG --IDSEPPath output/idsep/GOES-13_differential/csv --showplot

If background subtraction is performed, SEP events are automatically identified above background. An event definition with a flux threshold of 1e-6 (e.g. >10 MeV exceeds 1e-6 pfu) indicates that the SEP event was identified for all fluxes above background.


The --IDSEPEnhancement flag may be used to identify and analyze SEP events above background WITHOUT background subtraction. An event definition with a flux threshold of 1e-6 (e.g. >10 MeV exceeds 1e-6 pfu) indicates that the SEP event was identified for all fluxes above background.

GOES integral fluxes provided by NOAA already have some amount of background-subtraction. SEP enhancements above background may be evaluated without performing a background subtraction by:

    | python bin/idsep --StartDate 2011-05-01 --EndDate 2013-01-01 --Experiment GOES-13 --FluxType integral --RemoveAbove 10 --showplot

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --Threshold "30,1;50,1" --IDSEPEnhancement --IDSEPPath output/idsep/GOES-13_integral/csv --showplot

Although the background has not been subtracted, calling --IDSEPEnhancement will set values below mean + n*sigma to zero.


Apply Calibrated Energy Bin Corrections to GOES Data
====================================================

Sandberg et al. (2014) and Bruno (2017) published energy bin calibrations for selected GOES satellites by comparing with IMP-8 and PAMELA, respectively. These calibrated energy bins may be applied to GOES data within FetchSEP by using the --options flag.

Bruno (2017) corrections depend on the choice of spacecraft, the A or B GOES HEPAD detectors, uncorrected versus corrected GOES fluxes, and background subtraction.

The following example applies Sandberg et al. (2014) energy channels at the lower energies (using the ones calculated for GOES-11) and Bruno (2017) at the higher energies to GOES-13 uncorrected differential fluxes with background subtracted.

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType differential --options "S14;Bruno2017;uncorrected" --Threshold "30,1;50,1" --OPSEPSubtractBG --BGStartDate 2012-05-10 --BGEndDate 2012-05-12 --showplot


Automatically generate a Processed SEP Event list
=================================================

It is possible to run both codes with a single button push to create a preliminary SEP event list.

The process described here will take a long time series, generate rough SEP event lists for each energy channel, create a batch file that is split up into quiet periods and enhanced periods, then batch run these time periods through opsep for a more careful analysis and to produce a final SEP event list and non-event list. This automated process will result in enhanced time periods contain multiple SEP events as automatic SEP event discrimination has not yet been perfected. However, it will generate lists that can be used to remove SEP event periods from data sets for GCR studies or to get a good start on creating your own curated SEP event list.

The batch file is built from a list from one of the energy channels, from the files named like SEPTimes_GOES-13_integral_10.0_to_-1.txt. The user may choose the energy channel by setting a variable in fetchsep.cfg. The energy bin must exactly match one of the energy bins in the input data. For example:

    | ref_energy_bin = [10.0,-1] #integral
    | ref_energy_bin = [11.64,23.27] #differential GOES-18
    

The example below will create a rough SEP event list for a year of GOES-13 with `idsep` then process each enhanced period and each quiet period individually with `opsep` to extract characteristics:

    | python bin/fetchsep_prepare_obs --StartDate 2017-01-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType integral --RemoveAbove 10 --IDSEPEnhancement --Threshold "30,1;50,1"

will first run `idsep` on a specified data set and identify all increases above background. Note that the first days of the dataset may not have a good background solution. Output files are created that are then used to automatically run `opsep` in batch mode to analyze each quiet and elevated period. This creates a set of json another other supporting files for each SEP event and quiet time period in the time series.

Note that manual intervention is required to get a truly good event list. The automated method is not perfect at identifying individual SEP events, but it will get you 80% of the way there. 

There will be time periods that contain multiple events. The user may edit the batch_event_list_* file which contains each individual time period and rerun the SEP analysis:
 
     | python bin/fetchsep_prepare_obs --StartDate 2017-01-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType integral --RemoveAbove 10 --IDSEPEnhancement --Threshold "30,1;50,1" --StartPoint BATCH


The CLEAR Space Weather Center of Excellence Benchmark Dataset
==============================================================
A set of curated files and scripts have been created to allow any user to run FetchSEP and generate the CLEAR Benchmark SEP dataset. The dataset spans 1986-01-01 to 2025-09-10. When run, the scripts will download all GOES data between those time frames and process each spacecraft to calculate the mean background and SEP enhancements. The `idsep` output and plots folders will contain the mean background and sigma and the `opsep` output and plots folders will contain analysis of each individual SEP event and quite time periods.

The Benchmark dataset can be produced in its entirety using the deploy_* scripts in fetchsep/references/CLEAR.

** ---> Generating the Benchmark dataset from scratch requires about 15 - 20 hours of processing time and a reliable internet connection. <--- **

The deploy scripts will perform the following steps:

* Configure your environment to use the default settings in fetchsep/reference/fetchsep_CLEAR.cfg
* Set up the CLEAR/ folder to hold the Benchmark dataset results
* Download all GOES-06 to present data from NOAA and CCMC iSWA to your data/ folder
* Calculate the mean background and level of variation (sigma) for the full time series of each GOES spacecraft (most time consuming step of 15+ hours and produces about 20 GB of data)
* Generates "batch" lists used to analyze quiet time periods and individual SEP events
* Processes each time period and calculates values for each quiet period or SEP event, generating an SEP event list for each GOES spacecraft (3 - 4 hours 3 GB of data)
* Combines the SEP parameters from only the primary GOES spacecraft at the time to create a single SEP event list from 1986-01-01 to 2025-09-10 (The CLEAR Benchmark List)

**Two Benchmark lists and supporting data are produced at the end:**

* **Operational list** - GOES integral fluxes as-is (1.7 MB)
* **Energy Bin Calibrated List** - GOES uncorrected differential fluxes with background subtracted and Sandberg et al. (2014) and Bruno 2017 effective energy bins applied (423 KB)

Create the full Benchmark dataset from scratch by running the scripts from the base fetchsep directory:

Mac:

    | ./fetchsep/reference/CLEAR/deploy_CLEAR_Mac.sh

Linux:

    | ./fetchsep/reference/CLEAR/deploy_CLEAR_Linux.sh


Windows:

    | .\\fetchsep\\reference\\CLEAR\\deploy_CLEAR_Windows.bat

If you have already created the dataset once and you have the mean background solutions (everything in the output/idsep and plots/idsep folders), you can reprocess the SEP event analysis without having to rerun everything. Say you want to add a new event definition to your analysis (e.g. >10 MeV exceeds 100 pfu), you can modify the appropriate deploy script and rerun only the SEP analysis using the LISTS flag to skip the background calculation.

Mac:

    | ./fetchsep/reference/CLEAR/deploy_CLEAR_Mac.sh LISTS

Linux:

    | ./fetchsep/reference/CLEAR/deploy_CLEAR_Linux.sh LISTS


Windows:

    | .\\fetchsep\\reference\\CLEAR\\deploy_CLEAR_Windows.bat LISTS

**Note:** The deploy_CLEAR_Mac.sh script shows each command step-by-step. This is a good primer on how to use FetchSEP. You can use these commands as an example and modify them to run other experiments, e.g. to produce a similar dataset for STEREO or SOHO, including applying thresholds to differential energy channels.

**Note:** The CLEAR batch_event_list_* files were curated by hand to ensure that each time period to be analyzed contained on a quiet period or a single SEP event. If a user would like to follow this procedue to create their own list, the bin/fetchsep_prepare_obs command described above will perform the full workflow for an experiment and make a rough draft of the batch_event_list file in the process. That batch file can be curated by the user and bin/fetchsep_prepare_obs can be run with the --StartPoint BATCH flag to reanalyze the individual quiet and SEP time periods without recalculating the mean backgrounds.


Support
=======

Do not hesitate to contact Katie Whitman at kathryn.whitman@nasa.gov for support with this code.

Credits
=======

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
