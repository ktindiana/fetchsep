import array as arr

#Configuration values used by idsep, opsep, read_datasets, and
#derive_background.

datapath = '/Users/kwhitman/Documents/Programs/data'
outpath = '/Users/kwhitman/Documents/Programs/FetchSEP/output'
plotpath = '/Users/kwhitman/Documents/Programs/FetchSEP/plots'
listpath = '/Users/kwhitman/Documents/Programs/FetchSEP/lists'
badval = -1 #bad data points will be set to this value; must be negative
endfac = 1. #factor multiplied by flux threshold to determine end of
            #event; SET AS 1 to get same event end definition as SWPC
            #DO NOT CHANGE FROM 1.0 UNLESS YOU ARE SURE YOU WANT TO USE
            #DIFFERENT THRESHOLDS TO DETERMINE THE START AND END
            #OF AN EVENT.
dwell_time = 3.*60.*60. #Time in seconds after which an event is determined
                    #to end. If event remains below threshold for
                    #>dwell_time, then last point above threshold
                    #is event end.
errval = "Value Not Found"  #alternative value to None to indicate value not present


#description of all possible bad values that may be found in
#the data sets that could be read into this code
#For example, so data sets may use -999 to fill in a missing
#or bad value. Some may use 1.000000e+31
#All bad values will be set to None or -1 in the code
all_badval = [-999, -1, 1.000000e+31,"NaN", "nan", "n/a", 9.999999e+06,-707, 9.999E+003, 9.999000e+05, -1.0e+31, -9.9900e+02, -100000.0, -1.00e5, -5.55555547e+28, None]

#TIME SHIFT APPLIED TO ENTIRE TIME SERIES
#float in hours !!!WILL SHIFT TIMES IN USER-INPUT FILES IF NON-ZERO!!!
time_shift = 0. # positive shifts times later, negative shifts earlier

#####UNITS
#Set the units associated with your data
#Make sure that fluence units = flux units*s*sr
#If your data set is in units other than pfu or
#MeV^-1*cm^-2*s^-1*sr^-1, the two operational thresholds
#will not be equivalent to >10 MeV, 10 pfu and >100 MeV, 1 pfu
#--energy
energy_units = "MeV"
#--integral flux
flux_units_integral = "pfu" #pfu = cm^-2*s^-1*sr^-1
fluence_units_integral = "cm^-2"
#--differential flux
flux_units_differential = "MeV^-1*cm^-2*s^-1*sr^-1"
fluence_units_differential = "MeV^-1*cm^-2"
################################



###############################################################
###### CONFIGURATION VARIABLES FOR USER TIME SERIES FILES #####
###############################################################

#Set these variables when running fetchsep with "user" as the
#experiment.

##### DELIMETER between columns (any string)
##### For whitespace separating columns, use "" or " "
user_delim = " "

##### COLUMNS containing the fluxes you want to analyze
#####(require the first (0th) column contains date in YYYY-MM-DD
##### HH:MM:SS format OR in unixtime) so flux columns start at 1 or more
user_col = arr.array('i',[1,2,3,4,5,6,7,8])
err_col = arr.array('i',[]) #set to [] if no uncertainties
                            #err_col only used by idsep

##### ENERGY BINS associated with user file and columns specified
##### above as:
#####   [[Elow1,Ehigh1],[Elow2,Ehigh2],[Elow3,Ehigh3],etc]
##### Use -1 in the second edge of the bin to specify integral channel:
#####   [[Elow1,-1],[Elow2,-1],[Elow3,-1],etc]
user_energy_bins = [[750,-1],[500,-1],[300,-1],[100,-1],\
                    [60,-1],[50,-1],[30,-1],[10,-1]]



#######################################################
############## IDSEP CONFIGURATION VARIABLES ##########
#######################################################

#Set these variables to automatically identify SEP events
#in long time series using the idsep code.

#BACKGROUND SUBTRACTION
#derive_background.py calculates the mean background plus an
#expected level of variation (sigma).
#SEP flux is flux > mean + nsigma*sigma
#Background flux is flux <= mean + nsigma*sigma
#The value of nsigma affects how the event onset is captured because it
#determines the flux level that is considered background at the beginning of the
#event
idsep_nsigma = 3
init_win = 150 #days to average initial estimate of threshold
sliding_win = 27 #days in sliding window to calculate final threshold
percent_points = 0.9 #Percent of points that must be in the sliding
                    #window to calculate the background; otherwise use
                    #previous good value

###FOR SEP IDENTIFICATION#######
##THESE VALUES ARE CURRENTLY BEING SET IN AN AUTOMATED
##WAY IN THE CODE USING THE TIME RESOLUTION OF THE DATA.
##WILL LEAVE HERE FOR NOW TO CONSIDER WHETHER TO CONTROL
##THEM HERE INSTEAD.
#nconsec = 6 #number of consecutive points that must be nonzero
#allow_miss = 1 #number of points that can be zero when checking
#               #for SEP start
#dwell_pts = 4 #number of points that can be zero after SEP starts
#              #Dwell time
#################################


#######################################################
############## OPSEP CONFIGURATION VARIABLES ##########
#######################################################

#Set these variables for the analysis of a single event using
#the opsep code.

###FOR BACKGROUND SUBTRACTION - IF USED###
#derive_background.py calculates the mean background plus an
#expected level of variation (sigma).
#SEP flux is flux > mean + nsigma*sigma
#Background flux is flux <= mean + nsigma*sigma
#The value of nsigma affects how the event onset is captured because it
#determines the flux level that is considered background at the beginning of the
#event
opsep_nsigma = 2.0
#################################


#######################################################
######### COMMONLY USED VALUES FOR REFERENCE ##########
#######################################################
##### USER COLS
#   REleASE-30 [2,3,4,5]
#   REleASE-60 [6,7,8,9]
#   REleASE-90 [10,11,12,13]
#   SPARX [1,2]
#   STAT [1,2,3]
#   STAT SOHO [1,2,3,4]
#   SEPMOD [1,2,3,4,5,6,7,8]
#   SEPCaster [1,2,3,4,5]
#   GOES-16 [1,2,3,4,5,6,7,8]
#   M-FLAMPA [1,3,5]
#   iPATH differential [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
#   Pioneer10 & 11 min [2,3]
#   Rosetta/SREM L2 [2,3,4,5,6,7,8,9,10,11,12,13,14]
#   Rosetta/SREM V0 [1,2,3,4]
#   Ulysses/COSPINKET [2,3,4,5]
#   Ulysses/COSPINLETHET [6,7]
#   Voyager1/CRS [10,11,12,13,14,15,16,17,18]
#   Voyager2/CRS [10,11,12,13,14,15,16,17,18,19]
#   STEREOA&B [8,9,10,11]
#   IMP-8/CRNC H [5,6,7,8,9]
#   IMP-8/GME [19,20,21,22,23,24,25,26,27,28,29,30]
#   IMP-8/GME for comparison with SEPEM [9,10,12,14,17,19,20,22,24,27,28]
#   IMP-8/CPME [5,6,7,8,9,10]
#   SOHO/EPHIN L2 [3,4]
#   BeppiColombo/MPO/BERM/L0 electron [[0.3,0.3]]
#   BeppiColombo/MPO/BERM/L0 electron [[9.10,13.00]]
#   SOHO/EPHIN/ele [[0.25,0.70]]

##### ENERGY BINS
#   SEPEM_H_GOES13 P3 - P7 [[4,9],[12,23],[26,38],[40,73],[100,142],[160,242]]
#   EPREM [[10,-1],[30,-1],[40,-1],[50,-1],[100,-1]]
#   SEPMOD integral [[750,-1],[500,-1],[300,-1],[100,-1],\
#                    [60,-1],[50,-1],[30,-1],[10,-1]]
#   SEPMOD differential [[1000,1000],[750,750],[500,500],[300,300],[100,100],\
#                    [60,60],[50,50],[30,30],[10,10]]
#   SEPMOD IMP8 differential [[72,72],[53,53],[36,36],[26,26],[17,17],\
#                    [8.6,8.6],[5.1,5.1],[2.6,2.6],[1.2,1.2]]
#   SPARX [[10,-1],[60,-1]]
#   STAT [[10,-1],[50,-1],[100,-1]]
#   STAT SOHO [[4.3,7.8],[7.8,25.0],[25.0,40.9],[40.9,53]]
#   ASPECS [[10,-1],[30,-1]]
#   RELeASE [[4,9],[9,15.8],[15.8,39.6],[28.2,50.1]]
#   SEPCaster [[5,-1],[15,-1],[30,-1],[60,-1],[100,-1]]
#   GOES-16 [[1,-1],[5,-1],[10,-1],[30,-1],[50,-1],[60,-1],[100,-1],[500,-1]]
#   M-FLAMPA [[10,-1],[50,-1],[100,-1]]
#   iPATH differential [[0.1,0.1],[0.14677993,0.14677993],[0.21544347,0.21544347],\
#   [0.31622777,0.31622777],[0.46415888,0.46415888],\
#   [0.68129207,0.68129207],[1,1],[1.4677993,1.4677993],\
#   [2.1544347,2.1544347],[3.1622777,3.1622777],[4.6415888,4.6415888],\
#   [6.8129207,6.8129207],[10,10],[14.677993,14.677993],\
#   [21.544347,21.544347],[31.622777,31.622777],[46.415888,46.415888],\
#   [68.129207,68.129207],[100,100],[146.77993,146.77993],\
#   [215.44347,215.44347],[316.22777,316.22777],[464.15888,464.15888],\
#   [681.29207,681.29207],[1000,1000]]
#   Pioneer10&11/CRT H [[30.55,56.47],[120.7,227.3]]
#   Rosetta/SREM L2 [[11,14],[14,17.8],[17.8,22.6],[22.6,28.7],[28.7,36.4],[36.4,46.3],[46.3,58.8],[58.8,74.8],[74.8,95],[95,120.7],[120.7,153.4],[153.4,195],[195,247.7]]
#   Rosetta/SREM V0 [[48,270],[43,86],[52,278],[76,450]]
#   Ulysses/COSPINKET [[34.1,125],[125,250],[250,2200],[2200,-1]]
#   Ulysses/COSPINLETHET [[39,70],[71,94]]
#   Voyager1/CRS [[30.000,48.000],[48.000,56.000],[74.471,83.661],[132.834,154.911],[154.911,174.866],[174.866,187.713],[187.713,220.475],[220.475,270.050],[270.050,346.034]]
#   Voyager2/CRS [[30.000,48.000],[48.000,56.000],[75.861,82.562],[130.339,154.217],[154.217,171.338],[171.338,193.643],[193.643,208.152],[208.152,245.690],[245.690,272.300],[272.300,344.010]]
#   STEREOA&B [[33.4,35.8],[35.5,40.5],[40.0,60.0],[60.0,100.0]]
#   IMP-8/CRNC H [[29.75,40.10],[40.10,50.90],[50.90,62.55],[62.55,74.50],[74.50,94.78]]
#   IMP-8/GME [[35.2,42.9],[42.9,51],[51,63.2],[63.2,81],[87,92.5],[92.5,107],[107,121],[121,154],[154,178],[178,230],[230,327],[327,485]]
#   IMP-8/GME for comparison with SEPEM [[6,7.3],[7.2,8.6],[11.1,13.6],[16.2,18.8],[24.2,28.7],[35.2,42.9],[42.9,51],[63.2,81],[92.5,107],[154,178],[178,230]]
#   IMP-8/CPME H [[4.6,15.0],[15,25],[25,48],[48,96],[96,145],[190,440]]
#   SOHO/EPHIN L3 [[25,40.9],[40.9,53]]




def about_config():
    """ About config.py
    
        VALUES SPECIFIED IN utils/config.py:
        
            :datapath: directory containing data, 'data'
            :outpath: directory for program output, 'output'
            :plotpath: directory for saving plots, 'plots'
            :listpath: directory for lists (for run_multi_sep.py)
            :badval: will set any bad data points to this value
            :endfac: multiplicative factor to define threshold for
                    end of event; threshold*endfac (default 0.85)
            :nsigma: number of sigma to define SEP versus background
                    flux in background subtraction routine
            :version: if you are running a model or data set, allows you
                    to enter a version number
            :user_col: array defining flux columns (0 is always datetime)
            :user_delim: delimeter used to separate the columns in the time
                    profile file that you will read in
            :user_energy_bins: energy bins associated with the columns
                    specified in user_col
            :energy_units: e.g. "MeV'
            :flux_units_integral: e.g. "pfu"
            :fluence_units_integral: e.g. "cm^-2"
            :flux_units_differential: e.g. "MeV^-1*cm^-2*s^-1*sr^-1" (CCMC format)
            :fluence_units_differential: e.g. "MeV^-1*cm^-2" (CCMC format)
            
            (setting the units here will make correct units on plots and in json
            file, but doesn't change operational threshold values; must be done
            accordingly by hand)
    """
