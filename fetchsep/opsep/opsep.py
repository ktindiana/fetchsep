from ..utils import classes as cl
from ..utils import read_datasets as datasets
from ..utils import config as cfg
from ..json import ccmc_json_handler as ccmc_json
from ..utils import derive_background_opsep as bgsub
from ..utils import error_check
from ..utils import tools
from ..utils import plotting_tools as plt_tools
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import math
import numpy as np
import sys
#import urllib2
import re
import calendar
import datetime
import argparse
from datetime import timedelta
import os
import wget
from calendar import monthrange
import urllib.request
import csv
from dateutil.parser import parse
import scipy.integrate
from numpy import exp
import array as arr
import pandas as pd
import scipy
from scipy import signal
from statistics import mode
from lmfit import minimize, Parameters
import array
from pprint import pprint

__version__ = "3.7"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

################ CHANGE LOG ####################
#CHANGES in 0.4: allows users to specify user model or experiment name for
#output files
#CHANGES in 0.4: user passes filename through an argument if selects "user" for
#experiment rather than setting filename at the top of the code
#Changes in 0.5: Added ability to download and plot SOHO/ERNE data for the
#>10 MeV threshold calculations. -- INCOMPLETE
#Changes in 0.5: Indicate peak flux selection on plots.
#Changes in 0.6: SEP event considered to start only after 3 consecutive points.
#   Start time set to the first of the three points.
#Changes in 0.6: Added the --TwoPeaks flag to allow users to indicate that
#   the event has an initial threshold crossing for just a few points,
#   drops below threshold, then continues the increase to the main part of
#   the event. This is phenomenalogical addition after finding a few events
#   with this problem.
#Changes in 0.7: The main function now returns flags to indicate various
#   things that may help a user identify if an event has been poorly timed.
#   These return flags are generated with the idea that this could would
#   be run in batch mode for multiple events without looking at the output.
#   The flags would help a user decide after the fact if event timing should be
#   checked more closely. The flags include:
#       - Event starts at first time point on start date (i.e. previous
#           event already ongoing)
#       - Event is less than 12 hours long (i.e. start time might indicate an
#           initial rise above threshold, but miss main event. TwoPeaks flag
#           was created to handle events like this.)
#       - >100 MeV, 1 pfu threshold start is more than 12 hours after >10 MeV,
#           10 pfu threshold start (e.g. there might be multiple events in a
#           row and the first event is only >10 MeV and second events has >100)
#       - Event ends on very last time point indicating that it ended with the
#           time range rather than dropping below the threshold
#   A flag was added with the "UMASEP" option. When this flag is selected, the
#   code finds all information for four energy channels and thresholds used by
#   UMASEP: >10 MeV, 10 pfu; >100 MeV, 1 pfu; >30 MeV, 1 pfu; >50 MeV, 1 pfu
#   The proton flux in each of these channels is reported for specific times
#   after threshold crossing (Ts) and added to the sep_values_ output file for
#   times Ts + 3, 4, 5, 6, and 7 hours
#   Returns start date of SEP event found so that know the names of the output
#   files produced by this code.
#   Added a saveplot option to automatically write plots to file with a
#   unique filename.
#Changes in 0.7: Adding SEP event onset peak, i.e. the peak associated with the
#   initial acceleration at the sun. Versions 0.6 and previous find the absolute
#   peak in the entire time period between onset and end. For lower energy
#   channels, that peak often corresponds to the ESP portion of the event.
#   Now including a calculation of onset peak.
#Changes in 1.0: Adding support for SOHO/COSTEP/EPHIN differential proton
#   flux data. The 4 energy channels span from 4.3 - 53 MeV.
#   The user must input a threshold in differential flux to calculate all
#   SEP quantities. The operational thresholds cannot be applied.
#   Add ability for user to specify a differential flux threshold to define SEP
#   event.
#Changes in 1.1: Added ability for users to input multiple thresholds. If
#   user uses differential fluxes, then may request a combination of thresholds
#   based off of differential flux bins and estimated integral fluxes.
#Changes in 1.2: Adding options to:
#       Choose corrected or uncorrected GOES fluxes.
#       Choose to apply Bruno (2017) or Sandberg et al. (2014) effective
#       energies to GOES uncorrected data.
#2020-11-13, Changes in 2.0: Complete restructuring of code. A library directory has
#   been created and some of the subroutines originally in
#   operational_sep_quantities is now in read_datasets.py. Global variables,
#   such as directory names and the information users must input to run their
#   data sets, is not in config_opsep.py.
#   A code has been written to perform a background subtraction of the SEP
#   flux using a time period specified by the user - presumably spanning a day
#   or more prior to the SEP event. The separation of SEP flux and background,
#   followed by a background subtraction of the SEP flux is in
#   derive_background.py.
#2020-11-24, Changes in 2.1: Added ability to write out json file in CCMC SEP Scoreboard
#   format. More information about those files can be found at
#   https://ccmc.gsfc.nasa.gov/challenges/sep.php#format
#   If a threshold isn't crossed, the code now calculates the maximum flux
#   during the full time period (start date to end date) and the time and saves
#   it in the onset peak and onset date. This info won't be reported in the
#   csv files, but will be reported in the json.
#2021-01-06, Changes in 2.2: DHConsultancy (Daniel Heynderickx) sent out a beta
#   version of the SEPEM RSDv3 data sets. Added here as a native data set.
#2021-01-12, Changes in 2.3: Added functionality to read in REleASE data.
#2021-01-19, Changes in 2.4: If time resolution of data set is >15 minutes,
#   relax three point requirement to exceed or fall below a threshold.
#   For finding onset time, restricted onset to fall between event start and
#   start + 18 hours or start and end time, which ever is shorter.
# 2021-01-25, Changes in 2.4.1: Changed the logic when converting differential
#   flux to integral flux (from_differential_to_integral_flux). Previously,
#   if one bin (bin[i]) had non-zero flux and the next bin (bin[i+1]) had zero
#   flux, the bin[1+1] was set to a value of 1e-15 and then an integral was
#   calculate from bin[i] to 1e-15 to get the flux contribution. This gave
#   strange results. Now, if bin[i] or bin[i+1] has a flux of zero or None, no
#   flux is added to the integral flux.
#2021-02-25, Changes in 2.5: Changing the all_integral_fluences array to
#   all_threshold_fluences. In previous versions, all_integral_fluences
#   contained fluences for integral energy channels only, i.e. >10, >100 MeV.
#   If a user input a differential channel thresold, e.g. for an energy bin
#   of 5 - 7 MeV, then this all_integral_fluxes would hold a NaN value for
#   the fluence corresponding to the threshold array element in
#   all_integral_fluences.
#   Instead, we will now have all_threshold_fluences. The fluence will
#   will correspond to the fluence in each energy channel, whether it
#   is integral or differential. If the user specifies a differential energy
#   channel (e.g. 5 - 7 MeV) with a threshold, then the corresponding element in
#   all_threshold_fluences will have the fluence in the 5 - 7 MeV bin.
#2021-02-28, Changes in 2.5.1: Adjusting the onset peak estimation to fix bugs.
#   Reorganized some of the error checking into subroutines to clean up run_all.
#2021-03-02, Changes in 2.5.2: Added resorting of bins if they are in "reverse"
#   order. For example, SEPMOD produces flux files for energies in order
#   of 1000, 750, 500, ... MeV bins. This doesn't matter for integral fluxes,
#   but from_differential_to_integral flux expects differential bins
#   in increasing energy order. Added sort_bin_order to ensure that
#   energy bins are always in increasing order of effective energies.
#!***!2021-03-03, Changes in 2.5.3: Changed threshold crossing logic.
#   Previously,required flux to exceed threshold (>) for 3 consecutive points
#   to get event start. Changed to flux >= threshold for 3 consecutive points.
#2021-05-17, changes in 2.5.4: Nothing has changed in the
#   operational_sep_quantities.py code, but the ccmc_json_handler.py
#   code was modified to v0.4 and changed the format of the outpuer
#   json file a little bit. I want to mark this change within this code
#   as well and explicitly state the date when it happened.
#2021-05-27, changes in 2.5.5: Added the NoInterp flag to allow users
#   to specify that negative flux or None values should be set to None.
#   If NoInterp is not set, the default behavior is to fill in bad data
#   points using linear interpolation in time. This may not be desired
#   for model output as it inherently changes the model predictions.
#   Zeroes are always treated as valid values and are not replaced.
#   If there are gaps in the time steps, the code does NOT try to
#   fill in the gaps. It will only perform interpolation for time
#   steps present in the input data set.
#   Additionally, changed the way in which the time resolution was
#   calculated in calculate_fluence. Previously, the time resolution
#   was determined by the time difference between the first and
#   second time points. Now the difference in time is calculated for
#   every consecutive set of time points and the most common value
#   for the difference is used as the time resolution. This is done
#   in case there are time gaps in the data set.
#2021-05-28, changes in 2.6: Discovered an unintended bug in the
#   calculate_onset_peak code. I had intended to normalize the
#   deriv variable by the changing flux value to minimize
#   the derivative for jumps and dips when the flux is already
#   elevated. Instead, I had simply normalized by the very first
#   flux value in the array.
#   deriv = smooth_flux[i][j] - smooth_flux[i][j-nwin]
#   Previously: run_deriv[i].append(deriv/smooth_flux[i][0])
#   Changed to: run_deriv[i].append(deriv/smooth_flux[i][j-nwin])
#   Now using two versions of the derivative to get the onset peak.
#   One is the derivative divided by smooth_flux[i][j-nwin], which
#   really highlights the intial rise. The second is deriv divided
#   by a constant normalization factor (the first non-zero flux in
#   the array). This is used to escape local minima.
#   I intended to always use deriv/smooth_flux[i][j-nwin], but
#   I made a typo in v0.8 that ended up dividing by a constant
#   factor. I feel that combining the two concepts works best.
#!!!!2021-07-19, changes in 3.0: Went up to the next integer version!!!
#   number because this version has reconciled all differences
#   with the CCMC SEP Scoreboard JSON format. i.e. this version
#   produces JSON files and supporting output files in exactly
#   the format required by the SEP Scoreboard. The format
#   is specified at https://ccmc.gsfc.nasa.gov/challenges/sep.php.
#       - JSON in CCMC format
#       - cleaning up code and writing a lot more small subroutines.
#       - changed calculate_fluence to multiply by 4pi to remove
#            units of sr^-1
#       - added global units variables at top of code
#       - added ability for user to set units for a user-read file
#           in config_opsep.py
#       - when threshold is not crossed, maximum flux in time period
#           is saved in peak_flux rather than onset_peak, as before
#   2021-08-06: Changing onset peak flux definition to be more
#   independent of applied threshold. For events with maxmimum flux
#   values less than 500/energy_channel, search for the onset peak
#   starting 12 hours before the flux threshold is crossed.
#   Onset peak will still only be derived if thresholds
#   are crossed.
#2021-08-30, changes in 3.1: if the user tries to apply a threshold
#   to an energy channel that isn't in the data, give an error
#   and exit. Added ability for
#   check_bin_exists to check if integral bins within energy_bins.
#   In extract_integral_fluxes, change code to save flux as -999
#   rather than 0 if a >10 or >100 MeV energy channel doesn't
#   exist in the data set. Previously stored integral fluxes
#   as zero values if couldn't find the energy bin, but this is
#   no good since 0 is a valid model prediction and a valid
#   flux value in some data sets.
#2021-9-16, changes in 3.2: add feature to specify which type of json
#   file to write if the user inputs their own experiment. JSONType flag
#   added to the inputs.
#   Fixed bug in extract_integral_fluxes that only happened when an
#   input data set did not have the >10 MeV channel.
#2021-09-25, changes in 3.3: When convert GOES corrected differential
#   flux to integral flux, perform a background subtraction on the HEPAD
#   channels, since those are not treated with the Zwickl correction.
#   Currently using very rough background estimates:
#   330 - 420 MeV: 0.001804912
#   420 - 510 MeV: 0.001014797
#   510 - 700 MeV: 0.000431988
#   >700 MeV:      0.00013735
#2021-11-16, changes in 3.4: Added support in read_datasets and
#   operational_sep_quantities to read in GOES-R differential flux
#   data (5 min averaged). Integral fluxes are not readily available
#   from the NOAA website.
#2022-02-11, change in 3.5: Added support for real time GOES-R
#   integral fluxes served by CCMC. These are the daily integral
#   fluxes provided by NOAA SWPC and archived at CCMC. They are not
#   the official science-grade integral product, which is not
#   yet available. If the user specifies GOES-16 or GOES-17,
#   the integral fluxes will come from the primary instrument.
#   It's not possible to select between the spacecraft with
#   this particular product.
#2022-03-22, 2022-05-24 changes in 3.6: Added more colors and markers
#   in the plots in run_all. (in May 2022) Commented the ad hoc background
#   subtraction added in v3.3 in "from_differential_to_integral" that I had
#   applied to GOES corrected differential HEPAD channels.
#2022-06-16, changes in 3.7: Fixed how the peak fluxes are calculated if
#   no threshold is crossed. Removed the subroutine to swap peak fluxes
#   around and made sure they were calculated correctly in the subroutines
#   where they are first calculated.
#   Added caveat for code not to calculate onset peak if duration of
#   data is less than required to get a derivative value.
#   Added capability to read in Shaowen Hu's recalibrated GOES
#   data set called SRAG1.2.
#2022-08-04, changes in 3.8: Explicity added variables,
#   doBGSub, model_name, showplot, saveplot, in call to
#   append_differential_threshold subroutine. Worked on my mac
#   without them but caused crashes on other Windows machines.
#   Added check for determine_time_resolution() to exit gracefully
#   if passed a date array of only 1 data point.
#   library/read_datasets.py extract_date_range() was updated to
#   ensure that starting point was after the specified start time.
#2022-08-22, changes in 3.9: fixed a small bug in 3.8 (two colons).
#   Modified run_all so that plots would always be generated, even
#   if thresholds were not crossed.
#2022-09-08, changes in v3.10: Changed so that OpSEP will always produce
#   plots at the end, even if no thresholds were crossed.
#   Added a new subroutine/algorithm to estimate the location of the
#   onset peak: calculate_onset_peak_from_fit(). The previous algorithm,
#   calculate_onset_peak() is still in the code, however the default
#   choice is the new one that fits a Weibull to the first hours of an
#   SEP event and uses the second derivative of the Weibull to find
#   the location of the onset peak.
#2022-09-19, changes in v3.11: GOES-16 real time integral fluxes are
#   only available starting on 2020-03-08. Added check to ensure that
#   user won't pull archived real time fluxes from earlier spacecraft,
#   which are also stored at CCMC in the same location with the same
#   filenames. ONE EXCEPTION: The HEPAD files for GOES-14 and GOES-15
#   are missing a column for the month of 2019-09, so specify GOES-16
#   to pull the real time GOES fluxes from the CCMC archive instead.
#2022-11-14, still v3.11: Updated all_program_info() to reflect new All Clear
#   logic that was coded into library/ccmc_json_handler.py > fill_json().
#   This fixes a previous bug that incorrectly resulted in 
#   all_clear_boolean = False if only a single point was above 
#   threshold (for 5 min cadence GOES data, 3 consecutive points above
#   threshold are required to define an SEP event).
#2022-12-06, changes in v3.12: Fixed some error checking in calculate_fluence
#   and from_differential_to_integral flux to account for the
#   possibility of NaN values in the flux arrays when interpolation
#   is not used. These changes were needed when testing out the
#   additional of STEREO-A and -B data as native data sets in
#   read_datasets.py v1.3. Added STEREO-A and -B to list of allowed instruments.
#2023-03-03, changes in v3.2: Restructuring to be part of fetchsep package.
#   Combing with SEPAutoID (now idsep), which identifies multiple SEP events in
#   long time series, in same package. Merging overlapping code and supporting
#   files.
#   Changed filename from operational_sep_quantities.py to opsep.py.
#2023-04-10, changes in v3.3: Changed the logic to end an event. Previously
#   used 3 points below cfg.endfac*threshold (typically endfac was set at 0.85).
#   Now the code will apply a dwell time of 3 hours to find when the flux
#   has dropped below threshold for longer than the dwell time, then choose the
#   end time as the first below below threshold. endfac should be set to 1.0 in
#   config.py to ensure the end of the event applies the same threshold as the
#   start of the event.
#2023-09-30, changes in v3.4: Changed filenames of plots to contain zulu
#   time for starting window or threshold crossing time so that they
#   are unique.
#2024-02-09, changes in v3.5: Changed onset peak fitting algorithm. Added
#   logic that identifies large and small events and handles them a little
#   differently as far as the time range for applying the fit. For lower
#   intensity events, added logic that attempts to identify when the first
#   increase above background occurred (defined as a set flux threshold,
#   currently will only be useful for integral fluxes) and extends the fit
#   time range closer to that initial increase to capture more of the onset.
#2024-03-11, changes in v3.6: If onset peak time is calculated to be after the
#   maximum flux time, then choose to set the onset peak to the same point
#   as the max flux.
#2024-04-24, changes in v3.7: The json templates and ccmc_json_handler
#   were changed to reflect CCMC's schema.
#   "model":{"flux_type":} --> "source_info":{"native_flux_type":}
#   "event_lengths": [{"threshold":}] --> "event_lengths": [{"threshold_start":}]
########################################################################

#See full program description in all_program_info() below
datapath = cfg.datapath
outpath = cfg.outpath + "/opsep"
plotpath = cfg.plotpath + "/opsep"
badval = cfg.badval #bad data points will be set to this value; must be negative

# Prepare directories
cfg.prepare_dirs()
for path in (outpath, plotpath):
    if not os.path.isdir(path):
        print("Making directory:", path)
        os.mkdir(path)

#####UNITS#####
energy_units = cfg.energy_units
flux_units_integral = cfg.flux_units_integral
fluence_units_integral = cfg.fluence_units_integral
flux_units_differential = cfg.flux_units_differential
fluence_units_differential = cfg.fluence_units_differential

######FOR USER DATA SETS######
#(expect the first (0th) column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = array.array('i', cfg.user_col)
#DELIMETER between columns; for whitespace separating columns, use " " or ""
user_delim = cfg.user_delim
#DEFINE ENERGY BINS associated with user file and columns specified above as:
user_energy_bins = cfg.user_energy_bins
############################

#FILENAME(s) containing user input fluxes - WILL BE SET THROUGH ARGUMENT
user_fname = ['tmp.txt']


def about_opsep(): #only for documentation purposes
    """ Program description for opsep.py v3.2.
    
    Formatting for Sphinx web-based documentation, which can be viewed
    within the docs/index.html directory or online at:
    https://ktindiana.github.io/operational-sep/index.html
    
    This program will calculate various useful pieces of operational
    information about SEP events from GOES-08, -09, -10, -11, -12, -13, -14, -15
    data, GOES-R (-16, -17, real time integral), SOHO/EPHIN Level 3, SOHO/EPHIN
    real time data from the REleASE website, and the SEPEM (RSDv2 and RSDv3)
    dataset.

    SEP event values are always calculated for threshold definitions:
        
        * >10 MeV exceeds 10 pfu
        * >100 MeV exceed 1 pfu

    The user may add multiple additional thresholds through the command line.
    This program will check if data is already present in a 'data' directory. If
    not, GOES or EPHIN data will be automatically downloaded from the web. SEPEM
    (RSDv2 and RSDv3) data must be downloaded by the user and unzipped inside
    the 'data' directory. Because the SEPEM data set is so large (every 5
    minutes from 1974 to 2015 for RSDv2 and to 2017 for RSDv3), the program will
    break up the data into yearly files for faster reading.
    
    --NoInterp: Data sets are checked for bad data point (negative or None value
    fluxes) and the default behavior is to fill in those bad data points by
    performing a linear interpolation with time. This choice was made to
    calculate more accurate event-intergrated fluence values from data.
    Interpolation with time is not appropriate for model predictions, as
    it will inherently change the prediction or may not be desired by
    the user for the data set. Turn off interpolation with time by
    setting the --NoInterp flag (or nointerp=True). If the interpolation
    is turned off, negative flux values will be set to None.
    Zeroes are always treated as valid values and are not replaced.
    If there are gaps in the time steps, the code does NOT try to
    fill in the gaps. It will ONLY perform interpolation for time
    steps present in the input data set. i.e. gaps in time are not
    interpolated, only time steps with negative or None flux values.

    The values calculated here are important for space radiation operations:
       
       * Onset time, i.e. time to cross thresholds
       * Onset peak intensity
       * Onset peak time
       * Maximum intensity
       * Time of maximum intensity
       * Rise time (onset to peak)
       * End time, i.e. fall below 0.85*threshold for 3 points (15 mins for GOES)
       * Duration
       * Event-integrated fluences
       * Proton fluxes at various times after threshold crossing (UMASEP option)

    UNITS: User may choose differential proton fluxes (e.g. [MeV s sr cm^2]^-1)
    or integral fluxes (e.g. [s sr cm^2]^-1 or pfu). Default units are:
    MeV, cm, s, sr
    Units are specified in config/config_opsep.py.
    The program has no internal checks or requirements on units - EXCEPT FOR
    THE THRESHOLD DEFINITIONS OF >10, 10 and >100, 1.
    If you convert those thresholds in the main program to your units,
    you should be able to generate consistent results.
    Currently there are no features to change units automatically. SEE MORE
    about units in the USER INPUT DATA section below.

    OPTIONS: User may specify various options, that currently only apply to
    GOES data:
        
        * Choose corrected or uncorrected GOES fluxes.
        * Choose to apply Bruno (2017) or Sandberg et al. (2014) effective
          energies to GOES uncorrected data.
    
    .. code-block::
    
        --options uncorrected
        --options uncorrected,S14,Bruno2017 (recommend using background subtraction)

    BACKGROUND SUBTRACTION: Users may choose to perform a background
    subtraction by specifying:
    
    .. code-block::
    
        --SubtractBG --BGStartDate YYYY-MM-DD --BGEndDate YYYY-MM-DD
    
    The user should look at the data and select an appropriate time frame
    prior to the event when the background is calm and well-defined. If
    performing background subtraction, the mean background will be subtracted
    from the fluxes in the SEP time frame (StartDate to EndDate). Plots
    showing the mean background level, the background flux only, and the
    background-subtracted SEP fluxes will be created by derive_background to
    verify the quality of the background estimation and subtraction.

    If a previous event is ongoing and the specified time period starts with a
    threshold already crossed, you may try to set the --DetectPreviousEvent
    flag. If the flux drops below threshold before the next event starts, the
    program will identify the second event. This will only work if the
    threshold is already crossed for the very first time in your specified
    time period, and if the flux drops below threshold before the next event
    starts.

    If the event has an initial increase above threshold for a few points, falls
    below threshold, then continues to increase above threshold again, you
    may try to use the --TwoPeak feature to capture the full duration of the
    event. The initial increase above threshold must be less than a day. An
    example of this scenario can be seen in >100 MeV for 2011-08-04.

    A flag was added with the "UMASEP" option. When this flag is used
    (--UMASEP), the code finds all information for four energy channels and
    thresholds used by UMASEP: >10 MeV, 10 pfu; >100 MeV, 1 pfu; >30 MeV, 1 pfu;
    >50 MeV, 1 pfu. The proton flux in each of these channels is reported for
    multiple times after threshold crossing (Ts). The applied time delays are
    as follows:
        
        * >10 MeV - Ts + 3, 4, 5, 6, 7 hours
        * >30 MeV - Ts + 3, 4, 5, 6, 7 hours
        * >50 MeV - Ts + 3, 4, 5, 6, 7 hours
        * >100 MeV - Ts + 3, 4, 5, 6, 7 hours
    
    --spase_id: If you know the appropriate spase_id for the your model or
    experiment, you may specify it here to be filled in to the json file
    for CCMC's SEP Scoreboard.
    
    Note about EVENT END TIME: Currently, the code calculates the end of the
    event as the first time that the flux drops below threshold*endfac for
    three consecutive points. The endfac variable is specified in
    config/config_opsep.py.
    Mimicking a SRAG code, endfac is set to 0.85 (85% of threshold). This
    may be changed in the future to update this to more NOAA SWPC-like logic
    to define the end of an SEP event.
    
    Note about ONSET PEAK: The algorithm to estimate the location of the
    onset peak was changed in v3.10. The previous algorithm is still in
    the code, but the new one is implemented in the overall workflow.
    The new algorithm, called calculate_onset_peak_from_fit(), was implemented
    in an attempt to make the identification of the onset peak more robust.
    Following the approach of Kahler and Ling (2017), a modified Weibull is
    fit to the time profile between the time points 6 hours prior to a threshold
    crossing out to 24 hours after the threshold crossing. The second derivative
    is taken of the fitted Weibull to find the estimated onset location. The
    final reported value is the measured maximum flux within 1 hour of the
    estimated onset peak time derived from the Weibull fit.
    

    ALL CLEAR        
    library/ccmc_json_handler.py > fill_json() contains logic to determine the 
    All Clear status (all_clear_boolean) for each energy block. For 
    >10 MeV and >100 MeV, only specific thresholds are allowed to determine 
    the All Clear status: >10 MeV, 10 pfu and >100 MeV, 1 pfu.

    For all other energy channels, the All Clear status will be filled by the first
    threshold encountered by the code. e.g. if the user runs the code for 
    >30 MeV, 1 pfu and >30 MeV, 5 pfu (--Threshold "30,1;30,5"), the all_clear_boolean
    will reflect the >30 MeV, 1 pfu status. If the user flips the call when running OpSEP 
    (e.g. --Threshold "30,5;30,1"), then all_clear_boolean will reflect >30 MeV, 5 pfu.

    
    RUN CODE FROM COMMAND LINE (put on one line), e.g.:
    
    .. code-block::
    
        python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate 2012-05-19 --Experiment GOES-13
        --FluxType integral --showplot --saveplot
        
    .. code-block::
    
        python3 operational_sep_quantities.py --StartDate "2012-05-17 01:00:00"
        --EndDate "2012-05-19 12:00:00" --Experiment GOES-13
        --FluxType integral --showplot --saveplot


    RUN CODE FROM COMMAND FOR USER DATA SET (put on one line), e.g.:
    
    .. code-block::
    
        python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate '2012-05-19 12:00:00' --Experiment user --ModelName MyModel
        --UserFile MyFluxes.txt --FluxType integral --showplot

    RUN CODE FROM COMMAND LINE AND PERFORM BACKGROUND SUBTRACTION AND APPLY
    Sandberg et al. (2014) and Bruno (2017) effective energies to the GOES bins.
    (note: cannot bg-subtract GOES integral fluxes), e.g.:
    
    .. code-block::
        
        python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate '2012-05-19 12:00:00' --Experiment GOES-13
        --FluxType differential  --showplot --options uncorrected,S14,Bruno2017
        --SubtractBG --BGStartDate 2012-05-10 --BGEndDate --2012-05-17

    RUN CODE IMPORTED INTO ANOTHER PYTHON PROGRAM, e.g.:
    
    .. code-block::
    
        import operational_sep_quantities as sep
        start_date = '2012-05-17'
        end_date = '2012-05-19 12:00:00'
        experiment = 'GOES-13'
        flux_type = 'integral'
        spase_id = ''
        model_name = '' #if experiment is user, set model_name to describe data set
        user_file = '' #if experiment is user, specify filename containing fluxes
        json_type = '' #if experiment is user, specify which type of json file
                       # should be created
        showplot = True  #Turn to False if don't want to see plots
        saveplot = False #turn to true if you want to save plots to file
        options = '' #various options: S14, Bruno2017, uncorrected
        doBGSub = False #Set true if want to perform background subtraction
        bgstart_date = "2012-05-10" #Dates used to estimate mean background if
        bgend_date = "2012-05-17"   #doBGSub is set to True
        detect_prev_event = True  #Helps if previous event causes high intensities
        two_peaks = False  #Helps if two increases above threshold in one event
        umasep = False #Set to true if want UMASEP values (see explanation above)
        threshold = '' #Add a threshold to 10,10 and 100,1: '30,1' or '4.9-7.3,0.01'
        nointerp = False #Default False; set to True to stop linear interpolatin in time

        sep_year, sep_month,sep_day, jsonfname = sep.run_all(start_date, \
            end_date, experiment, flux_type, model_name, user_file, json_type,\
            spase_id, showplot, saveplot, detect_prev_event,  \
            two_peaks, umasep, threshold, options, doBGSub, bgstart_date, \
            bgend_date,nointerp)

    Set the desired directory locations for the data and output at the beginning
    of the program in datapath and outpath. Defaults are 'data' and 'output'.

    In order to calculate the fluence, the program determines time_resolution
    (seconds) by finding the difference between every consecutive set of
    time points in the data set. The most common difference is identified as
    the time resolution. This method should find an accurate time resolution
    even if there are gaps in the time steps.
    If the time steps in the data set are truly irregular, the user will
    have to manually set the time resolution inside the subroutine
    calculate_fluence.

    OUTPUT: This program outputs 3 to 4 files, 1 per defined threshold plus
    a summary file containing all of the values calculated for each threshold.
    A file named as e.g. fluence_GOES-13_differential_gt10_2012_3_7.csv contains
    the event-integrated fluence for each energy channel using the specified
    threshold (gt10) to determine start and stop times.
    A file named as e.g. sep_values_GOES-13_differential_2012_3_7.csv contains
    start time, peak flux, etc, for each of the defined thresholds.

    The program writes to file the >10 MeV and >100 MeV time series for the
    date range input by the user. If the original data were integral fluxes,
    then the output files simply contain the >10 and >100 MeV time series from
    the input files. If the original data were differential fluxes, then the
    estimated >10 and >100 MeV fluxes are output as time series.
    
    The program also writes any time profile to file if a threshold was applied
    to it. Each energy channel is written to an independent file, named according
    to CCMC's SEP Scoreboard naming conventions. These accompany a json file
    that contains all quantities calculated for the SEP event. There is also a csv
    file that contains most of the same values recorded in the json file.
    
    Example output for values derived from SEPMOD forecasts for the 2021-05-29
    SEP event with additional threshold >10, 0.001 pfu, >100, 0.0001 pfu,
    >30 MeV, 1 pfu, >50 MeV, 1 pfu, and >60 MeV, 0.079 pfu. Some of the files
    below are only created if a threshold was crossed. A default run would produce
    only 10 and 100 MeV files, sep_values_*.csv and json file:
    
    * fluence_SEPMOD_RT_60min_integral_gt10_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt10.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt30.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt50.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt60.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt100.0_2021_5_29.csv
    * integral_fluxes_SEPMOD_RT_60min_integral_2021_5_29.csv
    * sep_values_SEPMOD_RT_60min_integral_2021_5_29.csv
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.10.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.10MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.30.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.50.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.60.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.100.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.100MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.json
    
    The json and txt files listed above would be the ones that would be appropriate
    to read into the CCMC SEP Scoreboard or to pass to the SEP validation code
    being developed in conjunction with the SEP Scoreboard. The csv files are legacy
    files, but may also be easier for some users to read.
    
    PLOTS: Prior to v3.10, plots were only generated if a threshold was crossed.
    If only a subset of the specified thresholds were crossed, then only the cases
    where a threshold was crossed would show up in the plots and the others would
    be blank spaces.
    Starting in v3.10, plots of the flux time series are ALWAYS created,
    regardless of whether any thresholds are crossed.
    

    USER INPUT DATA SETS: Users may input their own data set. For example, if an
    SEP modeler would like to feed their own intensity time series into this
    code and calculate all values in exactly the same way they were calculated
    for data, it is possible to do that. Default flux units are
    1/[MeV cm^2 s sr] or 1/[cm^2 s sr] and energy channels in MeV for the default
    thresholds to be correct. You can specify your units in the
    config/config_opsep.py file.
    
    You can use any units, as long as you are consistent with energy units in
    energy channel/bin definition and in fluxes and you MODIFY THE THRESHOLD
    VALUES TO REFLECT YOUR UNITS. If you want to use different units, but
    still have the correct operational definitions, you need to modify these
    lines in define_thresholds() below:
    
        * energy_thresholds = [10,100] #MeV; flux for particles of > this MeV
        * flux_thresholds = [10,1] #pfu; exceed this level of intensity
    
    NOTE: The first column in your flux file is assumed to be time in format
    YYYY-MM-DD HH:MM:SS. IMPORTANT FORMATTING!!
    NOTE: The flux file may contain header lines that start with a hash #,
    including blank lines.
    NOTE: Any bad or missing fluxes must be indicated by a negative value.
    NOTE: Put your flux file into the "datapath" directory. Filenames will be
    relative to this path.
    NOTE: Please use only differential or integral channels. Please do not mix
    them. You may have one integral channel in the last bin, as this is the way
    HEPAD works and the code has been written to include that HEPAD >700 MeV
    bin along with lower differential channels.
    NOTE: You must specify whether your input file represents model output or
    observations so that the correct type of JSON file may be written. Use the
    JSONType flag or json_type variable and specify "model" or "observations"
    to indicate which is the correct format. The default is set to "model"
    when run from main.

    USER VARIABLES: The user must modify the following variables in
    config/config_opsep.py:
    
        :user_col: identify columns in your file containing fluxes to analyze;
                even if your delimeter is white space, consider the date-time
                column as one single column. SET IN config/config_opsep.py.
        :user_delim: delimeter between columns, e.g. " " or ","   Use " " for
                any amount of whitespace. SET IN config/config_opsep.py.
        :user_energy_bins: define your energy bins at the top of the code in the
                variable user_energy_bins. Follow the format in the subroutine
                define_energy_bins. SET IN config/config_opsep.py.
        :user_fname: specify the name of the file containing the fluxes
                through an argument in the command line. --UserFile  The
                user_fname variable will be updated with that filename. ARGUMENT
        :time_resolution: the program determines time_resolution
                (seconds) by finding the difference between every consecutive
                set of time points in the data set. The most common difference
                is identified as the time resolution. This method should find
                an accurate time resolution even if there are gaps in the
                time steps.. AUTOMATICALLY DETERMINED.
                
    Running the code for a user-input file may look like the example below.
    Note that the --UserFile location is with respect to the "data" directory
    inside the operational-sep directory:
    
    .. code-block::
    
        python3 operational_sep_quantities.py --StartDate 2021-05-29 --EndDate 2021-06-05 --Experiment user --UserFile SEPMOD/Scoreboard/SEPMOD.20210529_000000.20210529_165133.20210529_133005_geo_integral_tseries_timestamped_60min.txt --ModelName SEPMOD_RT_60min --JSONType model --showplot --Threshold "10,0.001;100,0.0001;30,1;50,1;60,0.079" --spase_id "spase://CCMC/SimulationModel/SEPMOD" --FluxType integral
    
    VALUES SPECIFIED IN config/config_opsep.py:
    
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
    print("The docstring in about_opsep describes the full code in detail.")


def make_dirs():
    """ Make subdirectories for files written out by idsep."""

    paths = ['csv','pkl','json']
    
    for path in paths:
        check_path = os.path.join(cfg.outpath,'opsep',path)
        if not os.path.isdir(check_path):
            print('Making directory: ', check_path)
            os.makedirs(check_path)

#### TO BE INCORPORATED INTO ANALYZE CLASS IF DEEMED NECESSARY ####
#def calculate_umasep_info(energy_thresholds,flux_thresholds,dates,
#                integral_fluxes, crossing_time):
#    """ Uses the integral fluxes (either input or estimated from differential
#        channels) and all the energy and flux thresholds set in the main program
#        to calculate SEP event quantities specific to the UMASEP model.
#            Flux at threshold crossing time + 3, 4, 5, 6, 7 hours
#            
#        INPUTS:
#        
#        :energy_thresholds: (float 1xn array) - energy channels for which thresholds
#            are applied
#        :flux_thresholds: (float 1xn array) - flux thresholds that are applied
#        :dates: (datetime 1xm array) - dates associated with flux time profile
#        :integral_fluxes: (float nxm array) - fluxes for each energy channel for
#            which a threshold is applied; each is the same length as dates
#        :crossing_time: (datetime 1xn array) - threshold crossing times for each energy
#            channel for which a threshold is applied
#            
#        OUTPUTS:
#        
#        :proton_delay_times: (datetime nx5 array) - times 3, 4, 5, 6, 7 hours
#            after crossing time for n thresholds
#        :proton_flux: (float nx5 array) - value of flux at each delay time and for
#            each threshold
#        
#    """
#    nthresh = len(flux_thresholds)
#    proton_flux = []
#    delays = [datetime.timedelta(hours=3), datetime.timedelta(hours=4),
#                    datetime.timedelta(hours=5), datetime.timedelta(hours=6),
#                    datetime.timedelta(hours=7)]
#    ndelay = len(delays)
#    delay_times = []
#    proton_delay_times = []  #actual time point corresponding to flux
#
#    #Match the correct time delay with the correct threshold
#    for i in range(nthresh):
#        if crossing_time[i] == 0:
#            delay_times.append(0)
#            continue
#        all_delays = []
#        for delay in delays:
#            all_delays.append(crossing_time[i] + delay)
#            #Make sure that delayed time doesn't exceed input time range
#            if crossing_time[i] + delay > dates[len(dates)-1]:
#                sys.exit("An UMASEP delayed time (Ts+3, 4, 5, 6, 7 hrs) "
#                        "exceeded the user's input time range. Please rerun "
#                        "and extend end time.")
#
#        delay_times.append(all_delays) #all delays for a given threshold
#
#    for i in range(nthresh):
#        save_flux = [0]*ndelay
#        save_dates = [0]*ndelay
#        if delay_times[i] == 0:
#            proton_delay_times.append(0)
#            proton_flux.append(0)
#            continue
#        for k in range(ndelay):
#            save_index = -1
#            for j in range(len(dates)):
#                if dates[j] <= delay_times[i][k]:
#                    save_index = j
#
#            #GET FLUX AT DELAYED TIME WITH 10 MINUTE AVERAGE
#            #May choose to modify if input data set has something other than
#            #5 minute time cadence.
#            if save_index == -1: #should not happen, unless dates has no length
#                sys.exit("Did not find an appropriate UMASEP flux point. "
#                        "Exiting.")
#            if save_index == 0: #also should be no way for this to happen
#                save_flux[k] = (integral_fluxes[i][save_index] + \
#                            integral_fluxes[i][save_index +1])/2.
#            else:
#                save_flux[k] = (integral_fluxes[i][save_index] + \
#                            integral_fluxes[i][save_index - 1])/2.
#            save_dates[k] = dates[save_index]
#
#        proton_flux.append(save_flux)
#        proton_delay_times.append(save_dates)
#
#    return proton_delay_times, proton_flux



def load_input_data(str_startdate, str_enddate, experiment,
    flux_type, model_name, user_file, showplot, saveplot, two_peaks,
    str_thresh, options, doBGSub, str_bgstartdate, str_bgenddate,
    nointerp, spacecraft, use_bg_thresholds):
    """ Instantiate an InputData object. Load all data.
        If differential fluxes specified, estimate integral fluxes.

    INPUTS:
    
        :str_startdate: (string) - user input start date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :str_enddate: (string) - user input end date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :experiment: (string) - "GOES-08" up to "GOES-15", "SEPEM", "SEPEMv3",
            "EPHIN", "EPHIN_REleASE", or "user"
        :flux_type: (string) - "integral" or "differential" indicates the type
            of flux to read in
        :model_name: (string) - If model is "user", set model_name to describe
            your model or data set (e.g. MyModel), otherwise set to ''.
        :user_file: (string) - Default is ''. If "user" is selected for experiment,
            specify name of flux file.
        :showplot: (bool) - Set to True to show plots when run
        :saveplot: (bool) - Set to True to automatically save plots to the
            plots directory when run
        :two_peaks: (bool) - option for extending event length
        :str_thresh: (string) - user-input thresholds in the format "30,1"
            for >30 MeV exceeds 1 pfu, "4-7,0.01" for 4-7 MeV differential
            channel exceeds 0.01.  "30,1;4-7,0.01" multiple thresholds
            separated by semi-colon.
        :nointerp: (boolean) - set to true to fill in negative fluxes with None
            value rather than filling in via linear interpolation in time
        :spacecraft: (string) primary or secondary is experiment is GOES_RT
        :use_bg_thresholds: (bool) Set to true to use the threshold calculated
            by idsep

    OUTPUT:
        
        :flux_data: (Data object)
           
    """

    flux_data = cl.Data()
    
    #Load all info, including specifying desired SEP event definitions (str_thresh)
    flux_data.load_info(str_startdate, str_enddate, experiment, flux_type,
        model_name, user_file, showplot, saveplot, two_peaks,
        str_thresh, options, doBGSub, str_bgstartdate,
        str_bgenddate, nointerp, spacecraft)

    #Read in fluxes; perform any background subtraction or interpolation
    flux_data.read_in_flux_files()

    #If want to use IDSEP background and threshold outputs to define an
    #event definition for the flux to exceed background levels, this would be
    #the place to add a new event definition.
    ##### Need to write the code to read in the correct files and extract the
    #info for the time period of interest and save as event definition.

    #Extract the fluxes associated with the desired event definitions.
    #If differential fluxes, then estimate integral fluxes
    flux_data.extract_fluxes_to_evaluate()
    
    return flux_data



def calculate_event_info(flux_data):
    """ Calculate SEP event characteristics for each of the event definitions.
        Save the data to dictionaries to turn into objects, jsons, and csv files.
        
        INPUT:
        
            :input_data: (Data object) object initialized with event definitions
                and fluxes
        
        OUTPUT:
        
            
    """

    for evdef in flux_data.event_definitions:
        analyze = cl.Analyze(flux_data, evdef)
        #Calculate all SEP characteristics
        analyze.calculate_event_info(flux_data)
        #Add to the flux_data object
        flux_data.add_results(analyze)

    flux_data.get_sep_date() #Set SEP year, month, day if a SEP occurred

    return flux_data




######## MAIN PROGRAM #########
def run_all(str_startdate, str_enddate, experiment, flux_type, model_name,
    user_file, json_type, spase_id, showplot, saveplot, detect_prev_event,
    two_peaks, umasep, str_thresh, options, doBGSub, str_bgstartdate,
    str_bgenddate, nointerp=False, spacecraft='',
    use_bg_thresholds=False, location='earth', species='proton'):
    """"Runs all subroutines and gets all needed values. Takes the command line
        arguments as input. Code may be imported into other python scripts and
        run using this routine.
        
        INPUTS:
        
        :str_startdate: (string) - user input start date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :str_enddate: (string) - user input end date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :experiment: (string) - "GOES-08" up to "GOES-15", "SEPEM", "SEPEMv3",
            "EPHIN", "EPHIN_REleASE", or "user"
        :flux_type: (string) - "integral" or "differential" indicates the type
            of flux to read in
        :model_name: (string) - If model is "user", set model_name to describe
            your model or data set (e.g. MyModel), otherwise set to ''.
        :user_file: (string) - Default is ''. If "user" is selected for experiment,
            specify name of flux file.
        :spase_id: (string) - Default is ''. If you know the spase_id of you
            model or experiment, enter it here for the json file
        :showplot: (bool) - Set to True to show plots when run
        :saveplot: (bool) - Set to True to automatically save plots to the
            plots directory when run
        :detect_prev_event: (bool) - option for finding start of event
        :two_peaks: (bool) - option for extending event length
        :umasep: (boolean) - call flag to run code for specific time related
            to the UMASEP model
        :str_thresh: (string) - user-input thresholds in the format "30,1"
            for >30 MeV exceeds 1 pfu, "4-7,0.01" for 4-7 MeV differential
            channel exceeds 0.01.  "30,1;4-7,0.01" multiple thresholds
            separated by semi-colon.
        :nointerp: (boolean) - set to true to fill in negative fluxes with None
            value rather than filling in via linear interpolation in time
        :templatename: (string) optional name of user json template located in
            cfg.templatepath directory
        :spacecraft: (string) primary or secondary is experiment is GOES_RT
        :use_bg_thresholds: (bool) Set to true to use the threshold
        
        OUTPUTS:
        
        :sep_year: (integer) - year of start of SEP event or year of beginning
            of user-input time period (if no threshold crossed)
        :sep_month: (integer) - month of start of SEP event or year of beginning
            of user-input time period (if no threshold crossed)
        :sep_day: (integer) - day of start of SEP event or year of beginning
            of user-input time period (if no threshold crossed)
        :jsonfname: (string) - name of saved json file (could include issue time,
            which is the time that the program is run, so absolutely
            need to pass file name back as argument)
        :and generates multiple plots:
        
    """
    
    datasets.check_paths()
    
    #Check for empty dates
    if (str_startdate == "" or str_enddate == ""):
        sys.exit('You must enter a valid date range. Exiting.')
    
    if experiment == "GOES_RT":
        if spacecraft != "primary" and spacecraft != "secondary":
            sys.exit(f"Spacecraft must be primary or secondary. You entered {spacecraft}. Please correct and run again.")
        
    #Instantiate an InputData object to hold all of the input data
    #information and fluxes
    flux_data = load_input_data(str_startdate, str_enddate, experiment,
        flux_type, model_name, user_file, showplot, saveplot, two_peaks,
        str_thresh, options, doBGSub, str_bgstartdate, str_bgenddate,
        nointerp, spacecraft, use_bg_thresholds)

    #Calculate SEP info for each event definition and create Analyze objects.
    flux_data = calculate_event_info(flux_data)

    #Create Output object to write out results
    output_data = cl.Output(flux_data, json_type, spase_id=spase_id,
                            location=location, species=species)
    jsonfname = output_data.write_ccmc_json()
    output_data.create_csv_dict()
    output_data.plot_event_definitions()
    output_data.plot_all_fluxes()
    output_data.plot_fluence_spectra()

    print(f"A SEP occurred on {flux_data.sep_year}-{flux_data.sep_month}-{flux_data.sep_day}")

    if showplot: plt.show()
    
    return flux_data.sep_year, flux_data.sep_month, flux_data.sep_day, jsonfname



    ####NOT YET REPRODUCED IN CLASS
    #Calculate times used in UMASEP
#    umasep_times =[]
#    umasep_fluxes=[]
#    if umasep:
#        umasep_times, umasep_fluxes = calculate_umasep_info(energy_thresholds,
#                        flux_thresholds, dates, integral_fluxes, crossing_time)
