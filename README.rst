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
* Documentation: https://fetchsep.readthedocs.io.

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

The code will output various csv files and a json file with accompany txt files. The JSON file is in the same format as required by the SEP Scoreboard to submit forecasts.

`opsep` creates files from observations that can be directly compared to SEP model forecasts sent to the SEP Scoreboard.

For time profile SEP models, `opsep` may be used to create the JSON files that can be submitted to the SEP Scoreboard.


To run OpSEP to process individual SEP events:

    | opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --showplot

OR

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --showplot

Two operational SEP event definitions are applied automatically for:
>10 MeV exceeds 10 pfu
>100 MeV exceeds 1 pfu

You may add addition event definitions with the --Threshold flag.

    | python bin/opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --Threshold "30,1;50,1" --showplot

To see how to add thresholds to differential energy channels, run:

    | python bin/opsep --help


Running `opsep` for your own time series
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


Automatically generate a Processed SEP Event list
-------------------------------------------------

It is possible to run both codes with a single button push to create a preliminary SEP event list.

The example below will create a rough SEP event list for a year of GOES-13 with `idsep` then process each enhanced period and each quiet period individually with `opsep` to extract characteristics:

    | python bin/fetchsep_prepare_obs --StartDate 2017-01-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType integral --RemoveAbove 10 --Threshold "30,1;50,1"

will first run `idsep` on a specified data set and identify all increases above background. Note that the first days of the dataset may not have a good background solution. Output files are created that are then used to automatically run `opsep` in batch mode to analyze each quiet and elevated period. This creates a set of json another other supporting files for each SEP event and quiet time period in the time series.

Note that manual intervention is required to get a truly good event list. The automated method is not perfect at identifying individual SEP events, but it will get you 80% of the way there. 

 There will be time periods that contain multiple events. The user may edit the batch_event_list_* file which contains each individual time period and rerun the SEP analysis:
 
     | python bin/fetchsep_prepare_obs --StartDate 2017-01-01 --EndDate 2018-01-01 --Experiment GOES-13 --FluxType integral --RemoveAbove 10 --Threshold "30,1;50,1" --StartPoint BATCH


Support
-------

Do not hesitate to contact Katie Whitman at kathryn.whitman@nasa.gov for support with this code.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
