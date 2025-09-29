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

FetchSEP and requirements may have the best results using
python version 3.10.

By default, `idsep` and `opsep` will create the necessary output
directories in the current working directory where the command is
executed.  The directories that will be created are `data`, `ouptut`,
`plots`, and `lists`.  The directories will not be overwritten if they
already exist.  If you would like to choose another output location,
generate a config file (`fetchsep.cfg`) with `opsep --ExportConfig`
and edit the paths described there.

A configuration file may also be placed in your home directory with
the name `.fetchsep`.  The configuration file does not need to be
complete; you may specify only the values for which you wish to
override the default.  The order of config value precidence
is 1. current working directory `fetchsep.cfg`, 2. home directory
`.fetchsep`, 3. fetchsep defaults.


Run
===

To run OpSEP to process individual SEP events:

    | opsep --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --showplot

IDSEP
=====

The `idsep` code will read in a long time series and automatically identify increases above background. This is done by estimating a mean background level plus an expected level of variation (sigma). All flux less than mean + 3sigma are considered background while all points above mean + 3sigma are identified as increases.

By assigning a set of criteria, increases that are most likely due to SEP events are identified and an SEP event list is output for each energy channel. 

The code also outputs a file containing every single high flux point above the mean + 3sigma threshold.

Note that in fetchsep/utils/config.py, number of sigma, the initial window used to estimate background levels, and the final sliding window used to estimate background levels can be adjusted. 

    | idsep_nsigma = 3
    | init_win = 150 #days to average initial estimate of threshold
    | sliding_win = 27 #days in sliding window to calculate final threshold
    | percent_points = 0.9 #Percent of points that must be in the sliding
    |                #window to calculate the background; otherwise use
    |                #previous good value



OPSEP
=====

The `opsep` code was previously supported at https://github.com/ktindiana/operational-sep and is now transitioned to this package going forward. Please see the operational-sep repository for extensive documentation until the documentation in this repository can be updated.

`opsep` is intended to assess each individual SEP event at a time, extracting information such as start and end times, peak fluxes, and event fluence.

The code will output various csv files and a json file with accompany txt files. The JSON file is in the same format as required by the SEP Scoreboard to submit forecasts.

`opsep` creates files from observations that can be directly compared to SEP model forecasts sent to the SEP Scoreboard.

For time profile SEP models, `opsep` may be used to create the JSON files that can be submitted to the SEP Scoreboard.



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


    
Automatically generate a Processed SEP Event list
-------------------------------------------------

It is possible to run both codes with a single button push to create a preliminary SEP event list. 
The code:

    | fetchsep_prepare_obs

will first run `idsep` on a specified data set and identify all increases above background. Output files are created that are then used to automatically run `opsep` in batch mode to analyze each quiet and elevated period. This creates a set of json another other supporting files for each SEP event and quiet time period in the time series.

Note that manual intervention is required to get a truly good event list. The automated method is not perfect at identifying individual SEP events, but it will get you 80% of the way there. 

Support
-------

Do not hesitate to contact Katie Whitman at kathryn.whitman@nasa.gov for support with this code.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
