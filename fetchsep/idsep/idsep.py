from ..utils import config as cfg
from ..utils import read_datasets as datasets
from ..utils import date_handler as dateh
from ..utils import define_background_idsep as defbg
from ..utils import plotting_tools as plt_tools
from ..utils import error_check
from ..utils import tools
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
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
from astropy.time import Time
from statistics import mode
import array
import pandas as pd

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#See full program description in all_program_info() below
datapath = cfg.datapath
outpath = os.path.join(cfg.outpath, "idsep")
plotpath = os.path.join(cfg.plotpath, "idsep")

# Prepare directories
cfg.prepare_dirs()
for path in (outpath, plotpath):
    if not os.path.isdir(path):
        print("Making directory:", path)
        os.mkdir(path)

nsigma = cfg.idsep_nsigma #threshold mean + nsigma*sigma
init_win = cfg.init_win
sliding_win = cfg.sliding_win
#Values below derived in identify_sep()
nconsec = 0 #number of consecutive points that must be nonzero
allow_miss = 0 #number of points that can be zero when checking
                    #for SEP start
dwell_pts = 0 #number of points that can be missed after SEP starts
                #Dwell time



""" About idsep.py
    
    Automatically identify SEP events by estimating background
    levels across the data set and separating backround
    and SEP fluxes.
    
    Option to perform background-subtraction of SEP fluxes.
    
    Identify start and end times of SEP events.
    
    Return arrays containing background fluxes only and
    SEP fluxes only.
    
"""



def read_in_flux_files(experiment, flux_type, user_file, startdate,
        enddate, options, dointerp, is_unixtime, write_fluxes=False,
        spacecraft=""):
    """ Read in the appropriate data or user files. Trims to dates
        between start time and end time. Interpolates bad
        points with linear interpolation in time unless nointerp True.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral, differential
        :user_file: (string) - file containing user's flux time profiles
        :model_name: (string) - model name or experiment if experiment = "user"
        :startdate: (datetime) - start date of time period entered by user
        :enddate: (datetime) - end date of time period entered by user
        :options: (string array) - options that could be applied
        :dointerp: (bool) - indicates if user DOES want to do linear
            interpolation in time for negative of bad flux values
        :is_unixtime: (bool) - flag to indicate that time in first column
            of user file is in unixtime
        :write_fluxes: (bool) True writes fluxes to standard format csv file
        
        OUTPUTS:
        
        :dates: (datetime 1xm array) - times in flux time profile trimmed
            between startdate and enddate
        :fluxes: (numpy float nxm array) - fluxes for n energy channels and m
            time steps; these are background subtracted fluxes if background
            subtraction was selected.
        :energy_bins: (array nx2 for n thresholds)
        
    """
    detector= []
    
    if experiment == "GOES":
        filenames1, filenames2, filenames_orien, detector = \
                datasets.check_data(startdate, enddate, experiment, flux_type, user_file, spacecraft=spacecraft)
    else:
        filenames1, filenames2, filenames_orien = datasets.check_data(startdate,
                enddate, experiment, flux_type, user_file, spacecraft=spacecraft)
                                    
    #read in flux files
    if experiment != "user":
        #Combine integral channels for all GOES spacecraft together into
        #a long time series. Only works for integral channels, since
        #GOES differential channels differ across experiments.
        if experiment == "GOES":
            all_dates, all_fluxes, west_detector = \
                datasets.read_in_files(experiment, flux_type,
                filenames1, filenames2, filenames_orien, options, detector)
        else:
            all_dates, all_fluxes, west_detector = \
                datasets.read_in_files(experiment, flux_type,
                filenames1, filenames2, filenames_orien, options)
    
    if experiment == "user":
        all_dates, all_fluxes = datasets.read_in_user_files(filenames1,is_unixtime)
        west_detector = []
    
    #Define energy bins
    #ERNE energy bins depend on the time period of the experiment.
    #For idsep, it is acceptable to include fluxes across slightly
    #different energy channels for the purposes of SEP identification.
    if experiment == "ERNE":
        version = datasets.which_erne(startdate, enddate)
        energy_bins = datasets.define_energy_bins(version, flux_type,
                                west_detector, options)
    else:
        energy_bins = datasets.define_energy_bins(experiment, flux_type,
                                west_detector, options,spacecraft=spacecraft)


    if energy_bins == None:
        sys.exit("Could not identify energy bins for experiment " + experiment
                + " and fluxtype " + flux_type)
    
    all_fluxes, energy_bins = tools.sort_bin_order(all_fluxes, energy_bins)


    #Extract the date range specified by the user
    dates, fluxes = datasets.extract_date_range(startdate, enddate,
                            all_dates, all_fluxes)
    
    #Interpolate bad data with linear interpolation in time or set to None
    print("read_in_flux_files: Checking for bad data and performing interpolation "
        "(if not deselected).")
    fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,dointerp)
    
        
    if len(dates) <= 1:
        sys.exit("The specified start and end dates were not present in the "
                "specified input file. Exiting.")

    if write_fluxes:
        tools.write_fluxes(experiment, flux_type, options, energy_bins, dates, fluxes, "idsep",
            spacecraft=spacecraft)

    return dates, fluxes, energy_bins



def identify_sep(dates, fluxes):
    """ Identify which increases above backgrounds
        are SEP events.
        
        INPUTS:
        
        :dates: (1xn datetime array) dates for each flux point
        :fluxes: (mxn float array) m energy channels and n time points
        
        OUTPUTS:
        
        :dates: (1xn datetime array) same as in
        :fluxes_sep: (mxn float array) all points set to zero except
            those identified as SEPs
            
    """
    time_res = tools.determine_time_resolution(dates)
    print("Time resolution of the data set is: "
            + str(time_res.total_seconds()) + " seconds.")
    time_res_sec = time_res.total_seconds()
    
    #DEPENDS ON TIME RESOLUTION
    #CAN BE DIFFICULTIES IN IDENTIFYING SEP EVENTS IN VERY
    #GAPPY DATA
    time_increase = 86400/4 #86400 #Require an increase above threshold for duration
    if time_res_sec <= 60*60:
        time_increase = 3*60*60 #6*60*60 #6 hr increase for hi-res data
    global nconsec
    nconsec = max(math.ceil(time_increase/time_res_sec) + 1,3) #num points
        #3 consecutive points for Voyager
    global allow_miss
    allow_miss = max(math.ceil(nconsec/4),1)
    if nconsec == 2 or nconsec == 3: allow_miss = 0
    global dwell_pts
    dwell_pts = max(math.ceil(nconsec/2),2)
    #Rosetta? needs --> dwell_pts = max(math.ceil(nconsec),2)

    print("Requiring " + str(nconsec) + " points (" + str(nconsec*time_res_sec/(60.*60.)) + " hours) to define an onset.")
    print("Allowing " + str(allow_miss) + " points to be missed in onset definition.")
    print("Event ends after " + str(dwell_pts) + " points are below threshold (dwell time).")
                
    npts = len(dates)

    IsSPE = False
    SPEflag = False
    
    SPEstart = [[]]*len(fluxes)
    SPEend = [[]]*len(fluxes)
    SPEfluxes = [[]]*len(fluxes)
    stidx = 0
    endidx = 0
    for j in range(len(fluxes)):
        SPEfluxes[j] = [0]*npts
        for i in range(npts-nconsec):
            if fluxes[j][i] > 0 and not IsSPE:
                IsSPE = True
                
                #Check that the increase continues
                #Flux below threshold will be set to zero
                nmiss = 0
                for k in range(i, i+nconsec):
                    chk_flux = fluxes[j][k]
                    if chk_flux <= 0:
                        nmiss = nmiss + 1
                
                #too many zero points, not an SPE
                if nmiss > allow_miss: IsSPE = False
            
            #Identify the start of an SPE
            if IsSPE and not SPEflag:
                SPEflag = True
                stidx = i
                if not SPEstart[j]:
                    SPEstart[j] = [dates[i]]
                else:
                    SPEstart[j].append(dates[i])
                i = min(i+nconsec,npts-nconsec-1) #jump to end of required consecutive points
              
            #ONGOING SPE with allowed gap
            if IsSPE and SPEflag:
                if fluxes[j][i] <= 0:
                    IsSPE = False
                    end_dwell = min(i+dwell_pts,npts-nconsec-1)
                    for ii in range(i,end_dwell):
                        chk_flux = fluxes[j][ii]
                        if chk_flux > 0: IsSPE = True
                        
                if not IsSPE or i == npts-nconsec-1:
                    SPEflag = False
                    endidx = i
                    if not SPEend[j]:
                        SPEend[j] = [dates[i]]
                    else:
                        SPEend[j].append(dates[i])
                
                    #Fill in the SEP flux array with the SEP points
                    for kk in range(stidx, endidx+1):
                        SPEfluxes[j][kk] = fluxes[j][kk]
                    
    return SPEstart, SPEend, SPEfluxes



def get_bg_high(threshold, dates, fluxes):
    fluxes_bg = []
    fluxes_high = []
    for i in range(len(fluxes)):
        flux_below, flux_above = defbg.separate_with_threshold(threshold[i],\
                dates,fluxes[i])
        fluxes_bg.append(flux_below)
        fluxes_high.append(flux_above)
        
    fluxes_bg = np.array(fluxes_bg)
    fluxes_high = np.array(fluxes_high)

    return fluxes_bg, fluxes_high





def write_sep_dates(experiment, exp_name, flux_type, energy_bins, options,
                    str_stdate, str_enddate, remove_above, SEPstart, SEPend,
                    doBGSub, for_inclusive=False, spacecraft=""):
    """ Write out SEP start and end times to file.
        
        INPUTS:
        
        :experiment: (string) name of experimen or "user"
        :exp_name: (string) if experiment is user, indicates name of experiment
        :flux_type: (string) integral or differential
        :energy_bins: (nx2 float array) energy bins for each channel
        :options: (array of strings) any modifications applied to data
            (only for GOES data currently)
        :str_stdate: (string) start of time period searched for SEP events
        :str_enddate: (string) end of time period searched for SEP events
        :remove_above: (float) intial cut input by user when program called
        :SEPend: (nxm datetime array) start times for n energy channels
            and m SEP events
        :SEPend: (nxm datetime array) end times for n energy channels
            and m SEP events
        :for_inclusive: (bool) if set to true, will output end time
            1 second before the next data point starts; if end time found
            to be point on 2011-01-01 00:00:00, will make end time
            2010-12-31 23:59:59
            
        OUTPUTS:
        
        none, but writes out n files - i.e. an SEP event list for each
            energy channel
        
    """
    
    one_sec = datetime.timedelta(seconds=1)
    
    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = plt_tools.setup_modifiers(options, doBGSub, spacecraft=spacecraft)


    prename = (f"SEPTimes_{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        prename = (f"SEPTimes_{exp_name}_{flux_type}{modifier}")
    
    #####WRITE SEP DATES OUT TO FILE INSTEAD OF PRINTING##########
    for j in range(len(SEPstart)):
        fname = (f"{prename}_{energy_bins[j][0]}_to_{energy_bins[j][1]}.txt")
        fname = os.path.join(outpath,fname)
        outfile = open(fname,"w")
        outfile.write("#SEP times calculated by idsep\n")
        if experiment == "user" and exp_name != '':
            outfile.write(f"#Experiment: {exp_name}\n")
        else:
            outfile.write(f"#Experiment: {experiment}\n")
        outfile.write(f"#Flux type: {flux_type}\n")
        outfile.write(f"#Energy channel: {energy_bins[j][0]} - {energy_bins[j][1]}\n")
        outfile.write(f"#Selected options: {options}\n")
        outfile.write(f"#Searched date range: {str_stdate} to {str_enddate}\n")
        outfile.write("#User applied an initial cut of " + str(remove_above)
                + ", initial averaging window of " + str(init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(sliding_win) + " days\n")
        outfile.write("#The percentage of points in the sliding window that must be "
                "background was " + str(cfg.percent_points) + "\n")
        outfile.write("#Threshold defined as mean background + " + str(nsigma)
                + " x sigma\n")
        outfile.write("#SEPs were identified when " + str(nconsec) + " points "
                + " exceeded threshold, allowing up to " + str(allow_miss)
                + " points to be missed. The SEP event ended after " +
                str(dwell_pts) + " were below threshold (dwell time).\n")
        outfile.write("#Start Time    End Time\n")
        
        for k in range(len(SEPstart[j])):
            SEPst = SEPstart[j][k]
            SEPed = SEPend[j][k]
            if for_inclusive:
                SEPed = SEPed - one_sec
            outfile.write(str(SEPst) + " " + str(SEPed) + "\n")
            
        outfile.close()

    
    

def write_all_high_points(experiment, exp_name, flux_type, energy_bins, options,
                    str_stdate, str_enddate, remove_above, dates, fluxes_high,
                    doBGSub, for_inclusive=False, spacecraft=""):
    """ Write out SEP start and end times to file.
        
        INPUTS:
        
        :experiment: (string) name of experimen or "user"
        :exp_name: (string) if experiment is user, indicates name of experiment
        :flux_type: (string) integral or differential
        :energy_bins: (nx2 float array) energy bins for each channel
        :options: (array of strings) any modifications applied to data
            (only for GOES data currently)
        :str_stdate: (string) start of time period searched for SEP events
        :str_enddate: (string) end of time period searched for SEP events
        :remove_above: (float) intial cut input by user when program called
        :dates: (1xm datetime array) dates associated with flux points
        :fluxes_high: (nxm float array) flux values for n energy channels and
            m dates; expect all values below threshold to be set to zero and
            all values above threshold to be non-zero
        :for_inclusive: (bool) if set to true, will output end time
            1 second before the next data point starts; if end time found
            to be point on 2011-01-01 00:00:00, will make end time
            2010-12-31 23:59:59
        
            
        OUTPUTS:
        
        none, but writes out n files - i.e. a list with all high points for each
            energy channel
        
    """
    #duration of each data point
    time_res = tools.determine_time_resolution(dates)
    if for_inclusive: time_res = time_res - datetime.timedelta(seconds=1)
    
    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = plt_tools.setup_modifiers(options, doBGSub, spacecraft=spacecraft)


    prename = (f"HighPoints_{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        prename = (f"HighPoints_{exp_name}_{flux_type}{modifier}")
    
    #####WRITE SEP DATES OUT TO FILE INSTEAD OF PRINTING##########
    for j in range(len(fluxes_high)):
        fname = (f"{prename}_{energy_bins[j][0]}_to_{energy_bins[j][1]}.txt")
        fname = os.path.join(outpath, fname)
        outfile = open(fname,"w")
        outfile.write("#All high points above mean background + 3*sigma calculated by idsep\n")
        if experiment == "user" and exp_name != '':
            outfile.write(f"#Experiment: {exp_name}\n")
        else:
            outfile.write(f"#Experiment: {experiment}\n")
        outfile.write(f"#Flux type: {flux_type}\n")
        outfile.write(f"#Energy channel: {energy_bins[j][0]} - {energy_bins[j][1]}\n")
        outfile.write(f"#Selected options: {options}\n")
        outfile.write(f"#Searched date range: {str_stdate} to {str_enddate}\n")
        outfile.write("#User applied an initial cut of " + str(remove_above)
                + ", initial averaging window of " + str(init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(sliding_win) + " days\n")
        outfile.write("#The percentage of points in the sliding window that must be "
                "background was " + str(cfg.percent_points) + "\n")
        outfile.write("#Threshold defined as mean background + " + str(nsigma)
                + " x sigma\n")
        outfile.write("#High flux points were identified when flux values "
                    "exceeded mean background flux + 3*sigma "
                    "by applying a " + str(sliding_win) + " days "
                    "backward sliding window to estimate the background levels.\n")
        outfile.write("#Start Time    End Time\n")
        
        for k in range(len(fluxes_high[j])):
            if fluxes_high[j][k] > 0:
                outfile.write(str(dates[k]) + " " + str(dates[k] + time_res) + "\n")
            
        outfile.close()


def separate_sep_with_dates(dates, fluxes, SEPstart, SEPend, padding):
    """ Using a list of dates, create flux files of only background and
        only SEP.
    
    """
    fluxes_bg = []
    fluxes_sep = []
    for i in range(len(fluxes)):
        flux_below, flux_above = defbg.separate_with_dates(dates, fluxes[i], SEPstart[i], SEPend[i], padding)
        fluxes_bg.append(flux_below)
        fluxes_sep.append(flux_above)
        
    fluxes_bg = np.array(fluxes_bg)
    fluxes_sep = np.array(fluxes_sep)

    return fluxes_bg, fluxes_sep



def rough_cut(init_win, dates, fluxes, nsigma, remove_above, energy_bins, experiment,
        flux_type, exp_name, options, doBGSub, showplot, saveplot, spacecraft=""):
    """ Remove fluxes above remove_above.
        Group fluxes into time periods of init_win.
        For each time period, calculate the mean background and sigma.
        Return the mean background, sigma, and threshold values of
        mean + nsigma*sigma.
        
        
        INPUTS:
        
            :init_win: (float) number of days in the initial window used to
                calculate the first rough cut of the background. e.g. 150 days.
            :dates: (1xn datetime array)
            :fluxes: (mxn float array) fluxes for m energy bins, n dates
            :nsigma: (float) number of sigma above the mean background to set the
                threshold that separates background from enhanced flux
            :remove_above: (float) and initial cut on flux applied to all
                energy bins. Flux > remove_above is considered an enhancement.
        
        OUTPUTS:
        
            :fluxes_bg: (mxn float array) background fluxes for m energy bins, n dates
            :fluxes_high: (mxn float array) enhanced fluxes for m energy bins, n dates
            
    """
    print("Preforming an initial rough cut between the background and enhanced fluxes.")
    print(f"Init win: {init_win}, Nsigma: {nsigma}, Remove above: {remove_above}")
    ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold =\
                defbg.ndays_average_optimized(init_win, dates, fluxes, nsigma, remove_above)
    
    
    #INITIAL SEPARATION OF BG AND HIGHER THAN BG
    fluxes_bg, fluxes_high = get_bg_high(threshold,dates,fluxes)
    
    if showplot or saveplot:
        plt_tools.idsep_make_plots(str(init_win)+"days", experiment, flux_type, exp_name, options, dates, fluxes, energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold, doBGSub, showplot, saveplot, spacecraft=spacecraft)
        
        plt_tools.idsep_make_bg_sep_plot(str(init_win)+"days", experiment, flux_type, exp_name, options, dates, fluxes_bg, fluxes_high, energy_bins, doBGSub,
            showplot, saveplot, spacecraft=spacecraft)

    
    return fluxes_bg, fluxes_high


def apply_sliding_window(sliding_win, dates, fluxes_bg_in, fluxes, nsigma,
    iteration=0):
    """ Identify the background value for every day using a sliding window.
        Use the initial estimated background from rough_cut and refine
        by applying a sliding window of sliding_win days to get a background
        value and sigma for every day of the data set. Create a daily threshold
        and apply to extract background fluxes, fluxes_bg, and enhanced fluxes,
        fluxes_high.
        
        INPUTS:
        
            :sliding_win: (float) number of days for sliding window
            :dates: (1xn datetime array)
            :fluxes_bg_in: (mxn float array) background fluxes for m energy bins,
                n dates
            :fluxes: (mxn float array) all fluxes for m energy bins, n dates
            :nsigma: (float) global number of sigma above background to define an
                enhancement
                
        OUTPUTS:
        
            :fluxes_bg: (mxn float array) background fluxes for m energy bins, n dates
            :fluxes_high: (mxn float array) enhanced fluxes for m energy bins, n dates
    
    """
    ave_background, ave_sigma, threshold =\
            defbg.backward_window_background_optimized(sliding_win, dates, fluxes_bg_in,
            nsigma, iteration)
    
    for i in range(len(fluxes_bg_in)):
        if None in fluxes_bg_in[i]:
            print("None values present in second: in " + str(i))
    
    fluxes_bg, fluxes_high = get_bg_high(threshold,dates,fluxes)

    return fluxes_bg, fluxes_high, ave_background, ave_sigma, threshold

def write_sep_fluxes(dates, fluxes, fluxes_bg):
    """ Write out final SEP fluxes and bg-subtracted fluxes.
        Subtract fluxes (e.g. SEP fluxes subtracted by the mean background).
        If fluxes already at a value of zero, no background subtraction is
        performed.
        
    """
    dict = {'dates': dates}
    cols = []
    for ii in range(len(fluxes)):
        dict.update({'fluxes'+str(ii): fluxes[ii]})
        cols.append('fluxes'+str(ii))
    df = pd.DataFrame(dict)
    df_bg = pd.DataFrame(dict)
    
    defbg.write_df(df,'SEP_fluxes_FINAL')
    
    for ii in range(len(fluxes)):
        df.update({cols[ii]:fluxes[ii]})
        df_bg.update({cols[ii]:fluxes_bg[ii]})

    df[cols] = df[cols] - df_bg[cols]
    df[df[cols] < 0][cols] = 0
    
    defbg.write_df(df,'background_subtracted_SEP_fluxes_FINAL')



def make_dirs():
    """ Make subdirectories for files written out by idsep."""

    paths = ['csv','pkl']
    
    for path in paths:
        check_path = os.path.join(cfg.outpath,'idsep',path)
        if not os.path.isdir(check_path):
            print('Making directory: ', check_path)
            os.makedirs(check_path)


def run_all(str_startdate, str_enddate, experiment,
        flux_type, exp_name, user_file, is_unixtime, options, doBGSub, dointerp,
        remove_above, for_inclusive, plot_timeseries_only, showplot, saveplot,
        write_fluxes=True, spacecraft=""):
    """ Run all the steps to do background and SEP separation.
    
        INPUTS:
        
        :write_fluxes: (bool) Write fluxes to csv file after read in and processed 
            for bad points
        :spacecraft: (string) primary or secondary if exp_name = GOES_RT

    
    """
    print("TIMESTAMP: Starting idsep " + str(datetime.datetime.now()))
    startdate = dateh.str_to_datetime(str_startdate)
    enddate = dateh.str_to_datetime(str_enddate)
    eff_startdate = startdate
    
    #If the user entered a date range shorter than required for the
    #initial window used to identify the background, extend the date
    #range
    #Note that the user should consider adding up to two months prior to the dates
    #of interest because the background solution for the first dates of the
    #timeseries are not accurate
    if not plot_timeseries_only:
        #Extend timeseries to have a buffer in beginning
        #Initially included, but hard to control starting date with this
        #added to the code
#        ndays = max(27*2,sliding_win*2)
 #       eff_startdate = startdate - datetime.timedelta(days=ndays)

        #Extend timeseries to cover init_win
        diff = (enddate - startdate).days
        if diff < init_win*2:
            eff_startdate = enddate - datetime.timedelta(days=init_win*2)
        
    
    error_check.error_check_options(experiment, flux_type, options, doBGSub)
    error_check.error_check_inputs(startdate, enddate, experiment, flux_type,
        subroutine='idsep')
    datasets.check_paths()
    make_dirs()
    
    #READ IN FLUXES
    print("TIMESTAMP: Reading in flux files " + str(datetime.datetime.now()))
    dates, fluxes, energy_bins = read_in_flux_files(experiment,
        flux_type, user_file, eff_startdate, enddate, options, dointerp, is_unixtime,
        write_fluxes=write_fluxes, spacecraft=spacecraft)
            
            
    if plot_timeseries_only:
        unique_id = "FluxTimeSeries"
        plt_tools.idsep_make_timeseries_plot(unique_id, experiment, flux_type, exp_name,
        options, dates, fluxes, energy_bins, doBGSub, showplot, saveplot)
        if showplot:
            plt.show()
        sys.exit("Time series plot completed. Exiting.")
    
    #ITERATION 1: DEFINE AN INITIAL "MOVING" THRESHOLD W/DATE
    print("TIMESTAMP: Creating rough cut first guess at threshold " + str(datetime.datetime.now()))
    fluxes_bg_init, fluxes_high_init = rough_cut(init_win, dates, fluxes, nsigma,
        remove_above, energy_bins, experiment, flux_type, exp_name, options, doBGSub,
        False, saveplot, spacecraft=spacecraft)
    
    
    #ITERATE over the identification of background and SEP periods
    #This process refines the background and the identification of SEP start and
    #end times
    niter = 3 #min = 2, max 5; tests show 3 iterations gives same result as 5
    fluxes_sep = []
    fluxes_bg = []
    fluxes_high = []
    ave_background = []
    ave_sigma = []
    threshold = []
    print(f"TIMESTAMP: Starting background and SEP event identification for {niter} iterations, {datetime.datetime.now()}.")
    for iter in range(niter):
        print(f"TIMESTAMP: Performing iteration {iter}, {datetime.datetime.now()}")
        post = "_iter" + str(iter)
        if iter == niter-1:
            post += "_FINAL"
        #Separate high and low flux by applying a sliding smoothing window to the background
        #fluxes_bg_init is used to get the mean, sigma, and threshold, then fluxes is split
        #into fluxes_bg and fluxes_high
        #The mean background is sensitive to the background selection in fluxes_bg_init
        print(f"TIMESTAMP: Starting sliding window background calculation, {datetime.datetime.now()}")
        fluxes_bg, fluxes_high, ave_background, ave_sigma, threshold =\
            apply_sliding_window(sliding_win, dates, fluxes_bg_init, fluxes, nsigma,
            iteration=iter)
        print(f"TIMESTAMP: Completed sliding window background calculation, {datetime.datetime.now()}")


        if showplot or saveplot:
            plt_tools.idsep_make_plots(str(sliding_win)+"window" + post, experiment, flux_type, exp_name, options, dates, fluxes_bg_init, energy_bins, dates, ave_background, ave_sigma, dates, threshold, doBGSub, False, saveplot,
                spacecraft=spacecraft)
            
            plt_tools.idsep_make_plots(str(sliding_win)+"window_nosigma" + post,
                experiment, flux_type,
                exp_name, options, dates, fluxes_bg_init, energy_bins, dates,
                ave_background, ave_sigma, dates, threshold, doBGSub, showplot,
                saveplot, True, spacecraft=spacecraft) #disable_sigma
            
            plt_tools.idsep_make_bg_sep_plot(str(sliding_win)+"window" + post, experiment,
                flux_type, exp_name, options, dates, fluxes_bg, fluxes_high, energy_bins,
                doBGSub, showplot, saveplot, spacecraft=spacecraft)


        #Identify SEP events in full time range
        #SEP identification proves to be very robust
        SEPstart, SEPend, fluxes_sep = identify_sep(dates, fluxes_high)

        if showplot or saveplot:
            plt_tools.idsep_make_bg_sep_plot("OnlySEP"+post, experiment, flux_type,
                exp_name, options, dates, fluxes, fluxes_sep, energy_bins, doBGSub,
                showplot, saveplot, spacecraft=spacecraft)


        #Taking the estimated background flux and remove SEP periods
        padding = 2 #niter - iter #number of days on either side of SEP start and end
        if iter == niter-2: #Final round
            fluxes_bg_init, fluxes_sep_padded = separate_sep_with_dates(dates, fluxes, SEPstart, SEPend, padding)
        if iter < niter-2:
            fluxes_bg_init, fluxes_sep_padded = separate_sep_with_dates(dates, fluxes_bg, SEPstart, SEPend, padding)

    print(f"TIMESTAMP: Completed background and SEP event separation, {datetime.datetime.now()}")
    

    #Trim fluxes to the date range specified by the user
    trim_dates, trim_fluxes_high = datasets.extract_date_range(startdate,
                enddate, dates, fluxes_high)
    trim_dates, trim_fluxes = datasets.extract_date_range(startdate,
                enddate, dates, fluxes)
    trim_bg_dates, trim_fluxes_bg = datasets.extract_date_range(startdate,
                enddate, dates, fluxes_bg)
    #Expect that this solution is exactly the same as the last one in the loop above,
    #UNLESS the dates need to be trimmed down
    SEPstart, SEPend, final_fluxes_sep = identify_sep(trim_dates, trim_fluxes_high)
    
    #Write background-subtracted SEP only fluxes to file
    write_sep_fluxes(trim_dates, final_fluxes_sep, trim_fluxes_bg)
    
    #Write start and end times to file
    write_sep_dates(experiment, exp_name, flux_type, energy_bins,
                    options, str_startdate, str_enddate, remove_above,
                    SEPstart, SEPend, doBGSub, for_inclusive, spacecraft=spacecraft)
    write_all_high_points(experiment, exp_name, flux_type, energy_bins,
                    options, str_startdate, str_enddate, remove_above,
                    dates, fluxes_high, doBGSub, for_inclusive, spacecraft=spacecraft)

 
    if showplot or saveplot:
        plt_tools.idsep_make_bg_sep_plot("FINALSEP", experiment, flux_type, exp_name, options,
            trim_dates, trim_fluxes, final_fluxes_sep, energy_bins, doBGSub, showplot, saveplot, spacecraft=spacecraft)
 
        
        if showplot: plt.show()

    print("TIMESTAMP: Completed idsep " + str(datetime.datetime.now()))
