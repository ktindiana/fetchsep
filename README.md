# fetch-sep
Identify SEP elevations above background in a time series (idsep) and analyze events individually (opsep)

#Set Up
In fetch_sep/utils/config.py, set the location of your datapath, outpath, plotpath, listpath.
datapath will be the location that the code downloads satellite data from the internet and stores it on your computer.
Note that outpath, plotpath, and listpath each need to have subdirectories idsep and opsep. For now, you need to create them yourself.

#Run
>source env.sh to add the directory to your python path

To run OpSEP to process individual SEP events:

>python3 bin/opsep.py --StartDate 2012-05-16 --EndDate 2012-05-22 --Experiment GOES-13 --FluxType integral --showplot

## IDSEP
The idsep code will read in a long time series and automatically identify increases above background. This is done by estimating a mean background level plus an expected level of variation (sigma). All flux less than mean + 3sigma are considered background while all points above mean + 3sigma are identified as increases.

By assigning a set of criteria, increases that are most likely due to SEP events are identified and an SEP event list is output for each energy channel. 

The code also outputs a file containing every single high flux point above the mean + 3sigma threshold.

Note that in fetch-sep/utils/config.py, number of sigma, the initial window used to estimate background levels, and the final sliding window used to estimate background levels can be adjusted. 

>idsep_nsigma = 3
>init_win = 150 #days to average initial estimate of threshold
>sliding_win = 27 #days in sliding window to calculate final threshold
>percent_points = 0.9 #Percent of points that must be in the sliding
>                    #window to calculate the background; otherwise use
>                    #previous good value


##OPSEP
The OpSEP code was previously supported at https://github.com/ktindiana/operational-sep and is now transitioned to this package going forward. Please see the operational-sep repository for extensive documentation until the documentation in this repository can be updated.

OpSEP is intended to assess each individual SEP event at a time, extracting information such as start and end times, peak fluxes, and event fluence.

The code will output various csv files and a json file with accompany txt files. The JSON file is in the same format as required by the SEP Scoreboard to submit forecasts.

OpSEP creates files from observations that can be directly compared to SEP model forecasts sent to the SEP Scoreboard.

For time profile SEP models, OpSEP may be used to create the JSON files that can be submitted to the SEP Scoreboard.

##Automatically generating an SEP Event list
It is possible to run both codes with a single button push to create a preliminary SEP event list. 
The code:

>bin/prep_obs.py

will first run idsep on a specified data set and identify all increases above background. Output files are created that are then used to automatically run opsep in batch mode to analyze each quiet and elevated period. This creates a set of json another other supporting files for each SEP event and quiet time period in the time series.

Note that manual intervention is required to get a truly good event list. The automated method is not perfect at identifying individual SEP events, but it will get you 80% of the way there. 

##Support
Do not hesitate to contact Katie Whitman at kathryn.whitman@nasa.gov for support with this code.
