from ..utils import config as cfg
from ..utils import directories as dirs
from ..utils import parameters as fsparam
from ..utils import read_datasets as datasets
from ..utils import download as fsdl
from ..utils import date_handler as dh
from ..utils import analysis
from ..utils import define_background_idsep as defbg
from ..utils import plotting_tools as plt_tools
from ..utils import tools
from ..utils import names
from ..utils import experiments as expts
from ..json import ccmc_json_handler as ccmc_json
import datetime
from datetime import timedelta
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import pandas as pd
import copy
from pathlib import Path


__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"




#Values below derived in identify_sep()
#Define them here as global variables so that they can be filled
#in the routines in a different module
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



def get_bg_high(threshold, dates, fluxes):
    fluxes_bg = []
    fluxes_high = []
    for i in range(len(fluxes)):
        print(f"get_bg_high {i} LENGTH OF DATES {len(dates)}. LENGTH of fluxes {len(fluxes[i])}. LENGTH of threshold {len(threshold[i])}")
        flux_below, flux_above = defbg.separate_with_threshold(threshold[i],\
                dates,fluxes[i])
        fluxes_bg.append(flux_below)
        fluxes_high.append(flux_above)
        
    fluxes_bg = np.array(fluxes_bg)
    fluxes_high = np.array(fluxes_high)

    return fluxes_bg, fluxes_high



def write_sep_dates(params, energy_bins, SEPstart, SEPend):
    """ Write out SEP start and end times to file.
        
        INPUTS:
        
        :params: (FetchSEP Parameters object) 
        :energy_bins: (nx2 float array) energy bins for each channel
        :SEPend: (nxm datetime array) start times for n energy channels
            and m SEP events
        :SEPend: (nxm datetime array) end times for n energy channels
            and m SEP events
            
        OUTPUTS:
        
        none, but writes out n files - i.e. an SEP event list for each
            energy channel
        
    """
    
    one_sec = datetime.timedelta(seconds=1)
    
    name = params.idsep_subdir
    prename = (f"SEPTimes_{name}")
    zstdate = dh.time_to_zulu(params.startdate)
    zstdate = zstdate.replace(':', '')
    zenddate = dh.time_to_zulu(params.enddate)
    zenddate = zenddate.replace(':', '')
    
    #####WRITE SEP DATES OUT TO FILE##########
    for j in range(len(SEPstart)):
        fname = (f"{prename}_{zstdate}_{zenddate}_{energy_bins[j][0]}_to_{energy_bins[j][1]}.csv")
        fname = os.path.join(params.module_outpath, fname)
        outfile = open(fname,"w")
        outfile.write("#SEP times calculated by idsep\n")
        if params.experiment == "user" and params.user_name != '' and params.user_name != None:
            outfile.write(f"#Experiment: {params.user_name}\n")
        else:
            outfile.write(f"#Experiment: {params.experiment}\n")
        outfile.write(f"#Flux type: {params.flux_type}\n")
        outfile.write(f"#Energy channel: {energy_bins[j][0]} - {energy_bins[j][1]}\n")
        outfile.write(f"#Selected options: {params.options}\n")
        outfile.write(f"#Searched date range: {params.startdate} to {params.enddate}\n")
        outfile.write("#User applied an initial cut of " + str(params.remove_above)
                + ", initial averaging window of " + str(params.init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(params.sliding_win) + " days\n")
        outfile.write("#The percentage of points in the sliding window that must be "
                "background was " + str(params.percent_points) + "\n")
        outfile.write("#Threshold defined as mean background + " + str(params.idsep_nsigma)
                + " x sigma\n")
        outfile.write("#SEPs were identified when " + str(nconsec) + " points "
                + " exceeded threshold, allowing up to " + str(allow_miss)
                + " points to be missed. The SEP event ended after " +
                str(dwell_pts) + " were below threshold (dwell time).\n")
        outfile.write("#Start Time,End Time\n")
        
        for k in range(len(SEPstart[j])):
            SEPst = SEPstart[j][k]
            SEPed = SEPend[j][k]
            if params.for_inclusive:
                SEPed = SEPed - one_sec
            outfile.write(str(SEPst) + "," + str(SEPed) + "\n")
            
        outfile.close()

    

def write_all_high_points(params, energy_bins, dates, fluxes_high):
    """ Write out SEP start and end times to file.
        
        INPUTS:
        
        :params: (FetchSEP Parameters object) 
        :energy_bins: (nx2 float array) energy bins for each channel
        :dates: (1xm datetime array) dates associated with flux points
        :fluxes_high: (nxm float array) flux values for n energy channels and
            m dates; expect all values below threshold to be set to zero and
            all values above threshold to be non-zero
        
            
        OUTPUTS:
        
        none, but writes out n files - i.e. a list with all high points for each
            energy channel
        
    """
    #duration of each data point
    time_res = analysis.determine_time_resolution(dates)
    if params.for_inclusive: time_res = time_res - datetime.timedelta(seconds=1)
    
    #Additions to titles and filenames according to user-selected options
    name = params.idsep_subdir

    prename = (f"HighPoints_{name}")
    zstdate = dh.time_to_zulu(params.startdate)
    zstdate = zstdate.replace(':', '')
    zenddate = dh.time_to_zulu(params.enddate)
    zenddate = zenddate.replace(':', '')
    
    #####WRITE SEP DATES OUT TO FILE INSTEAD OF PRINTING##########
    for j in range(len(fluxes_high)):
        fname = (f"{prename}_{zstdate}_{zenddate}_{energy_bins[j][0]}_to_{energy_bins[j][1]}.csv")
        fname = os.path.join(params.module_outpath, fname)
        outfile = open(fname,"w")
        outfile.write("#All high points above mean background + 3*sigma calculated by idsep\n")
        if params.experiment == "user" and params.user_name != '' and params.user_name != None:
            outfile.write(f"#Experiment: {params.user_name}\n")
        else:
            outfile.write(f"#Experiment: {params.experiment}\n")
        outfile.write(f"#Flux type: {params.flux_type}\n")
        outfile.write(f"#Energy channel: {energy_bins[j][0]} - {energy_bins[j][1]}\n")
        outfile.write(f"#Selected options: {params.options}\n")
        outfile.write(f"#Searched date range: {params.startdate} to {params.enddate}\n")
        outfile.write("#User applied an initial cut of " + str(params.remove_above)
                + ", initial averaging window of " + str(params.init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(params.sliding_win) + " days\n")
        outfile.write("#The percentage of points in the sliding window that must be "
                "background was " + str(params.percent_points) + "\n")
        outfile.write("#Threshold defined as mean background + " + str(params.idsep_nsigma)
                + " x sigma\n")
        outfile.write("#High flux points were identified when flux values "
                    "exceeded mean background flux + 3*sigma "
                    "by applying a " + str(params.sliding_win) + " days "
                    "backward sliding window to estimate the background levels.\n")
        outfile.write("#Start Time,End Time\n")
        
        for k in range(len(fluxes_high[j])):
            if fluxes_high[j][k] > 0:
                outfile.write(str(dates[k]) + "," + str(dates[k] + time_res) + "\n")
            
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



def rough_cut(params, dates, fluxes, energy_bins):
    """ Remove fluxes above remove_above.
        Group fluxes into time periods of init_win.
        For each time period, calculate the mean background and sigma.
        Return the mean background, sigma, and threshold values of
        mean + nsigma*sigma.
        
        
        INPUTS:
        
            :dates: (1xn datetime array)
            :fluxes: (mxn float array) fluxes for m energy bins, n dates
            :energy_bins: (array) energy bins that go with fluxes 
            :params: (FetchSEP Parameters object) 
        
        OUTPUTS:
        
            :fluxes_bg: (mxn float array) background fluxes for m energy bins, n dates
            :fluxes_high: (mxn float array) enhanced fluxes for m energy bins, n dates
            
    """
    print("Preforming an initial rough cut between the background and enhanced fluxes.")
    print(f"Init win: {params.init_win}, Nsigma: {params.idsep_nsigma}, Remove above: {params.remove_above}")
    ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold =\
                defbg.ndays_average_optimized(params.init_win, dates, fluxes, energy_bins,
                params.idsep_nsigma, params.remove_above, savepath=params.module_outpath)
    
    #INITIAL SEPARATION OF BG AND HIGHER THAN BG
    fluxes_bg, fluxes_high = get_bg_high(threshold,dates,fluxes)
    
    if params.showplot or params.saveplot:
        plt_tools.idsep_make_plots(str(params.init_win)+"days", params, dates, fluxes, energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold, showplot=False)
        
        plt_tools.idsep_make_bg_sep_plot(str(params.init_win)+"days", params, dates, fluxes_bg,
            fluxes_high, energy_bins, showplot=False)

    
    return fluxes_bg, fluxes_high


######### RESUME ##########
def df_to_arrays(df, time_col='dates'):
    """ Dataframe with dates in first column and values in the others """
    dates = df[time_col].to_list()
    cols = df.columns
    vals = np.array(df[cols[1:]])
    vals = vals.T
    return dates, vals


def resume_fluxes(params, file_str='fluxes_*.csv', remove_sep=False,
    trim_to_window=False, last_date_only=False):
    """ Read in the original fluxes from the previous run. e.g.,
        fluxes_GOES-13_integral_20110806_20120531.csv
        SEP_fluxes_background-subtracted_FINAL.csv
        SEP_fluxes_FINAL.csv
        
        If remove_sep, remove the time periods with previously identified
        SEP events to recreate the background-only fluxes used to calculate
        mean and sigma.
        
        INPUTS:
        
            :remove_sep: (bool) False will return original fluxes;
                True will return fluxes with SEP periods removed (null)
            :trim_to_window: (bool) False will return the entire time series;
                True will trim from the end point back by the sliding window
                used to calculating the mean background. e.g. only the last
                5 days of the time series will be returned
            :last_date_only: (bool) will return a single output with the
                last date of the flux files generated during the previous
                run of IDSEP
                
        OUTPUTS:
        
            :dates: (datetime arr)
            :fluxes: (arr) array of flux time series for all energy channels
        
    """
    directory = Path(params.idsep_path)
    # Find all text files in this folder only
    fname = ''
    for file in directory.glob(file_str):
        fname = file

    if fname == '':
        sys.exit(f"RESUME: resume_background_fluxes: Could not find {file_str} file in {params.idsep_path}. Exiting.")
    df = pd.read_csv(fname)
    df['dates'] = pd.to_datetime(df['dates'])
    
    #Get the last date in the files generated by the previous run of IDSEP
    lastdate = df.at[len(df['dates'])-1,'dates']
    if last_date_only:
        return lastdate

    if params.startdate is None:
        sys.exit("resume_fluxes: Startdate is None. Get the resume start date by rerunning with last_date_only = True and set params.startdate. Exiting.")

    #Remove the last point in the previous run file (i.e. the startdate
    #used by resume to add more data), because the flux read
    #in for the next flux analysis will start at the previous endpoint.
    df = df.loc[df['dates'] < params.startdate]

    #If trim to the size of the sliding window
    if trim_to_window:
        firstdate = lastdate - datetime.timedelta(days=params.sliding_win)
        firstdate = datetime.datetime(firstdate.year, firstdate.month, firstdate.day)
        df = df.loc[df['dates'] >= firstdate]

    #Do not remove SEP time periods, i.e. keep original fluxes
    if not remove_sep:
        dates, fluxes = df_to_arrays(df)
        return dates, fluxes
    
    #Read in previously identified SEP time periods for each energy channel
    #and remove SEP periods from the original fluxes
    fluxes = []
    bin_keys = df.columns.to_list()
    bin_keys = bin_keys[1:]
    for key in bin_keys:
        #Convert the energy bin names to those in the SEPTimes filenames.
        #e.g. 10.0--1 --> _10.0_to_01 or 4.2-9.0 --> 4.2_to_9.0
        if '--' in key:
            bin = key.strip().split('--')
            label = f"{bin[0]}_to_-1"
        else:
            bin = key.strip().split('-')
            label = f"{bin[0]}_to_{bin[1]}"
            
        fname = ''
        for file in directory.glob('SEPTimes_*' + label + '.csv'):
            fname = file

        if fname == '':
            sys.exit(f"RESUME: resume_fluxes: Cannot find file containing SEP times for {key}. Exiting.")

        df_sep = pd.read_csv(fname, comment='#', header=None, names=['Start Time','End Time'])
        df_sep['Start Time'] = pd.to_datetime(df_sep['Start Time'])
        df_sep['End Time'] = pd.to_datetime(df_sep['End Time'])
        SEPstart = df_sep['Start Time'].to_list()
        SEPend = df_sep['End Time'].to_list()

        dates = df['dates'].to_list()
        flux_vals = df[key].to_list()

        #Remove SEP periods
        padding = 0 #number of days on either side of SEP start and end
        flux_bg, flux_sep = defbg.separate_with_dates(dates, flux_vals, SEPstart,
                SEPend, padding)
        fluxes.append(flux_bg)

    return df, dates, fluxes


def resume_derived_values(params, filename, energy_bins, trim_to_window=False):
    """ Read in the files from a previous idsep run that will act
        as a starting point for analysis of additional data in time.
        
        The resume directory will be stored in the variable
        params.idsep_path
        
        INPUT:
        
            :params: (Parameters object)
            :energy_bins: (list) energy bins of resume data, for comparison
                with previously prepared idsep results read from the resume directory 
                as a consistency check, e.g. [[10,-1],[30,-1],..]
            :trim_to_window: (bool) False will return the entire time series;
                True will trim from the end point back by the sliding window
                used to calculating the mean background. e.g. only the last
                5 days of the time series will be returned
                
    """
    bgfilename = os.path.join(params.idsep_path, filename)

    df = pd.read_csv(bgfilename)
    df['dates'] =pd.to_datetime(df['dates'])
    cols = df.columns.to_list()

    #Bruno2017 energy bins depend on which detector is the west detector
    #and which spacecraft is used. IDSEP will use the energy bins
    #for the data at the start of the dataset, so need to relax
    alt_energy_bins = []
    if params.goes_Bruno2017:
        alt_energy_bins_A, centers_A = datasets.define_energy_bins(params, ["A"])
        alt_energy_bins_B, centers_B = datasets.define_energy_bins(params,["B"])
        alt_energy_bins = alt_energy_bins_A + alt_energy_bins_B

    previous_bin_keys = cols[1:] + alt_energy_bins
    #Check that energy bins in previous idsep data and resume data are the same
    for bin in energy_bins:
        bin_key = names.energy_bin_key(bin)
        if bin_key not in previous_bin_keys:
            sys.exit(f"RESUME: prepare_resume: Energy bins do not match between new data {bin_key} and previous idsep results {previous_bin_keys} in {params.idsep_path}. Exiting.")


    lastdate = params.startdate #df.at[len(df)-1,'dates']
    #Remove the last point in the previous run file (i.e. the startdate
    #used by resume to add more data), because the flux read
    #in for the next flux analysis will start at the previous endpoint.
    df = df.loc[df['dates'] < params.startdate]
    
    #Trim to the end of the resume arrays, keeping only the duration of one sliding window
    if trim_to_window:
        firstdate = lastdate - datetime.timedelta(days=params.sliding_win)
        firstdate = datetime.datetime(firstdate.year, firstdate.month, firstdate.day)
        df = df.loc[df['dates'] >= firstdate]
        if df.empty:
            sys.exit("RESUME: read_idsep_files: The idsep file containing the mean background "
                    f"does not cover the dates required. {firstdate} to {lastdate}")

    dates, vals = df_to_arrays(df, time_col='dates')

    return df, dates, vals


def prepare_resume(params, energy_bins):
    """ Prepare the values needed to start the IDSEP analysis from a previous run. """
    #Get the background fluxes with SEP events removed and trimmed to sliding window
    df_fluxes, flux_dates, fluxes = resume_fluxes(params, file_str='fluxes_*.csv', remove_sep=True, trim_to_window=True)

    #Previous mean background and sigma
    df_means, mdates, means = resume_derived_values(params, params.idsep_fname_background, energy_bins, trim_to_window=True)
    df_sigmas, sdates, sigmas = resume_derived_values(params, params.idsep_fname_sigma, energy_bins, trim_to_window=True)
    df_thresholds, tdates, thresholds = resume_derived_values(params, params.idsep_fname_threshold, energy_bins, trim_to_window=True)
    df_kurtosis, kdates, kurtosis = resume_derived_values(params, params.idsep_fname_kurtosis, energy_bins, trim_to_window=True)

    if flux_dates != mdates:
        print(f"flux_dates {flux_dates[0]} {flux_dates[-1]} dates {mdates[0]} {mdates[-1]}")
        sys.exit("RESUME: prepare_resume: The flux dates don't match the means, sigmas, thresholds dates. Exiting.")

    resume_arrays = {
        'dates': flux_dates,
        'fluxes': fluxes,
        'mean': means,
        'sigma': sigmas,
        'threshold': thresholds,
        'kurtosis': kurtosis,
        'df_fluxes': df_fluxes,
        'df_means': df_means,
        'df_sigmas': df_sigmas,
        'df_thresholds': df_thresholds,
        'df_kurtosis': df_kurtosis
        }

    del mdates
    del sdates
    del tdates
    del kdates

    return resume_arrays


def resume_cut(params, dates, fluxes, energy_bins):
    """ Use the last mean and sigma calculated in the previous run
        as a remove_above value and perform an initial rough separation
        of SEP enhancements and background values.
        
        Start by using mean + 4*sigma for remove_above.
        
        INPUTS:
        
            :dates: (1xn datetime array)
            :fluxes: (mxn float array) fluxes for m energy bins, n dates
            :energy_bins: (array) energy bins that go with fluxes 
            :params: (FetchSEP Parameters object) 
        
        OUTPUTS:
        
            :fluxes_bg: (mxn float array) background fluxes for m energy bins, n dates
            :fluxes_high: (mxn float array) enhanced fluxes for m energy bins, n dates
            
    """
    print("RESUME: Performing an initial rough cut between the background and enhanced fluxes using values from the previous IDSEP run.")
    
    resume_arrays = prepare_resume(params, energy_bins)
    #Prepend dates and fluxes with a sliding window's worth of fluxes from the previous run
    dates = resume_arrays['dates'] + dates
    nsigma = 4
    remove_above = []
    fluxes = np.concatenate((resume_arrays['fluxes'], fluxes), axis=1)
    for i in range(len(resume_arrays['sigma'])):
        remove_above.append(resume_arrays['mean'][i][-1] + nsigma*resume_arrays['sigma'][i][-1])
        print(f"{i} {resume_arrays['dates'][0]} {resume_arrays['fluxes'][i][0:20]}")

    #Apply remove_above = mean + 4sigma for each energy channel individually
    print(f"RESUME: resume_cut: Init win: {params.init_win}, Nsigma: {nsigma}, Remove above: {remove_above}")
    ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold =\
            defbg.ndays_average_optimized(params.init_win, dates, fluxes, energy_bins,
            nsigma, remove_above)

    for i in range(len(fluxes)):
        print(f"{i} LENGTH OF DATES {len(dates)}. LENGTH of fluxes {len(fluxes)} and {len(fluxes[i])}. LENGTH of threshold {len(threshold)} and {len(threshold[i])}")
    
    #INITIAL SEPARATION OF BG AND HIGHER THAN BG
    fluxes_bg, fluxes_high = get_bg_high(threshold,dates,fluxes)
    
    if params.showplot or params.saveplot:
        plt_tools.idsep_make_plots(str(params.init_win)+"days", params, dates, fluxes, energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold, showplot=False)
        
        plt_tools.idsep_make_bg_sep_plot(str(params.init_win)+"days", params, dates, fluxes_bg, fluxes_high, energy_bins, showplot=False)
    
    return dates, fluxes, fluxes_bg, fluxes_high, resume_arrays



def combine_resume_fluxes(params, dates, fluxes):
    """ Combine original fluxes with fluxes from previous run """
    #Read in the original fluxes from the previous run
    prev_dates, prev_fluxes = resume_fluxes(params, file_str='fluxes_*.csv', remove_sep=False)
    
    #Trim current fluxes to exclude the sliding window addition that was tacked
    #onto the beginning
    trim_dates, trim_fluxes = datasets.extract_date_range(params.startdate, params.enddate, dates, fluxes)
    
    #Concatonate
    dates = prev_dates + trim_dates
    fluxes = np.concatenate((prev_fluxes, trim_fluxes), axis=1)
    
    print(f"PREVIOUS RUN {len(prev_dates)}, CURRENT RUN {len(trim_dates)}, TOTAL {len(dates)}")
    
    return dates, fluxes
####################################



def apply_sliding_window(params, dates, fluxes_bg_in, fluxes, energy_bins, iteration=0,
    is_final=False, resume_arrays={}):
    """ Identify the background value for every day using a sliding window.
        Use the initial estimated background from rough_cut and refine
        by applying a sliding window of sliding_win days to get a background
        value and sigma for every day of the data set. Create a daily threshold
        and apply to extract background fluxes, fluxes_bg, and enhanced fluxes,
        fluxes_high.
        
        INPUTS:
        
            :params: (FetchSEP Parameters object)
            :dates: (1xn datetime array)
            :fluxes_bg_in: (mxn float array) background fluxes for m energy bins,
                n dates
            :fluxes: (mxn float array) all fluxes for m energy bins, n dates
            :prev_means: (pandas DataFrame) if resuming, one sliding window of
                mean background values from previous run
            :prev_sigmas: (pandas DataFrame) if resuming, one sliding window of
                sigma from previous run
                
        OUTPUTS:
        
            :fluxes_bg: (mxn float array) background fluxes for m energy bins, n dates
            :fluxes_high: (mxn float array) enhanced fluxes for m energy bins, n dates
    
    """
    #If resuming and on the final iteration, will return entire time series
    #including values from the previous run for background, sigma, and threshold
    ave_background, ave_sigma, threshold =\
            defbg.backward_window_background_optimized(params, dates, fluxes_bg_in,
            energy_bins, iteration, is_final=is_final, resume_arrays=resume_arrays)
    
    for i in range(len(fluxes_bg_in)):
        if None in fluxes_bg_in[i]:
            print("None values present in second: in " + str(i))
 
    ########### RESUME###############
    #If RESUMING, add in all of the results from the previous IDSEP run
    #and then use the full time period to get the final high and bg separation.
    if is_final and params.idsep_resume:
        #Add in previous original dates and fluxes
        dates, fluxes = combine_resume_fluxes(params, dates, fluxes)
    ########### RESUME###############
 
    for i in range(len(fluxes)):
        print(f"apply_sliding_window {i} LENGTH OF DATES {len(dates)}. LENGTH of fluxes {len(fluxes)} and {len(fluxes[i])}. LENGTH of threshold {len(threshold)} and {len(threshold[i])}")
 
    fluxes_bg, fluxes_high = get_bg_high(threshold,dates,fluxes)

    return fluxes_bg, fluxes_high, ave_background, ave_sigma, threshold



def write_sep_fluxes(params, dates, fluxes, fluxes_bg, energy_bins):
    """ Write out final SEP fluxes and bg-subtracted fluxes.
        Subtract fluxes (e.g. SEP fluxes subtracted by the mean background).
        If fluxes already at a value of zero, no background subtraction is
        performed.
        
    """
    savepath = params.module_outpath
    dict = {'dates': dates}
    dict_bg = {'dates': dates}
    cols = []
    for ii in range(len(fluxes)):
        bin = energy_bins[ii]
        key = names.energy_bin_key(bin)
        dict.update({key: fluxes[ii]})
        dict_bg.update({key: fluxes_bg[ii]})
        cols.append(key)
    df = pd.DataFrame(dict) #original SEP fluxes with all background set to zero
    df_bg = pd.DataFrame(dict_bg)
    
    defbg.write_df(df, params.idsep_fname_sep, savepath=savepath)
    
    df[cols] = df[cols] - df_bg[cols]
    for col in cols:
        df.loc[df[col] < 0, col] = 0
    
    defbg.write_df(df, params.idsep_fname_sep_bgsub, savepath=savepath)




def run_idsep(str_startdate, str_enddate, experiment,
    flux_type=None, spacecraft=None,
    user_name=None, user_file=None,
    idsep_path=None, idsep_resume=None,
    directory_depth=None,
    is_unixtime=None, options=None, dointerp=None,
    remove_above=None, for_inclusive=None,
    kurtosis_cut=None,
    idsep_nsigma=None,
    init_win=None,
    sliding_win=None,
    percent_points=None,
    plot_timeseries_only=None,
    showplot=None, saveplot=None,
    write_fluxes=None,
    use_absolute_datapath=None,
    path_to_data=None,
    path_to_output=None,
    path_to_plots=None,
    path_to_lists=None):
    """ Run all the steps to do background and SEP separation.
    
        INPUTS:

        :str_startdate: (string) - user input start date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :str_enddate: (string) - user input end date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :experiment: (string) - "GOES-08" up to "GOES-15", "SEPEM", "SEPEMv3",
            "EPHIN", "EPHIN_REleASE", or "user"
        :flux_type: (string) - "integral" or "differential" indicates the type
            of flux to read in
        :spacecraft: (string) primary or secondary if exp_name = GOES-RT
        :user_name: (string) - If experiment is "user", set user_name to describe
            your model or data set (e.g. MyModel), otherwise set to ''.
        :user_file: (string) - Default is ''. If "user" is selected for experiment,
            specify name of flux file.
        :is_unixtime: (bool) True indicates first column in user file is in unixtime
        :directory_depth: (int) default = 2; Subdirectories for output files may be 
                supressed by choosing the directory depth.
                0 - Files output to top directories: cfg.outpath (output/), cfg.plotpath (plots/) level; 
                1 - Files output to subdirectory at module level, cfg.outpath/module (output/opsep); 
                2 - Files output to subdirectory named according to experiment and 
                options, e.g. cfg.outpath/module/subdir (output/opsep/GOES-13_integral/
        :options: (string) may specify a series of options as a semi-colon separated list. 
            uncorrected - for GOES uncorrected differential fluxes
            S14 - apply Sandberg et al. (2014) effective energies to GOES P2-P6 
                (derived for GOES uncorrected fluxes)
            Bruno2017 - apply Bruno (2017) effective energies to GOES-13
                or GOES-15 P6-P11 channels for either corrected or uncorrected
                GOES fluxes. Bruno recommends performing background subtraction. 
            If both S14 and Bruno2017 are specified for GOES-13 or GOES-15, 
            S14 bins will be applied to P2-P5 and Bruno2017 bins will be applied 
            to P6-P11 for uncorrected fluxes. e.g. "uncorrected;S14;Bruno17"
        :dointerp: (boolean) - set to true to fill in data gaps via linear interpolation in time, otherwise fill with nan values
        :remove_above: (float) Remove all flux points above a specified value. 
            Helps to exclude high values above background during the first iteration
            to estimate the background.
        :for_inclusive: (bool) Write out SEP end times such that they end 1 second 
            before the next data point begins.
        :plot_timeseries_only: (bool) True to only download the data and plot the 
            flux timeseries without calculating background and SEP events.
        :write_fluxes: (bool) True to write fluxes to csv file after read in and processed 
            for bad points (default = True)
        :path_to_data: (string) path where satellite data should be downloaded. Will default to 
            datapath listed in fetchsep.cfg if a value is not specified.
        :path_to_output: (string) path where output files should be saved. Will default to 
            outpath listed in fetchsep.cfg if a value is not specified.
        :path_to_plots: (string) path where plots should be saved. Will default to
            plotpath listed in fetchsep.cfg if a value is not specified.
        :path_to_lists: (string) path where lists should be saved. Will default to
            listpath listed in fetchsep.cfg if a value is not specified.
    
    """
    print("TIMESTAMP: Starting idsep " + str(datetime.datetime.now()))
    expts.set_config_energy_units(experiment)
    expts.set_config_flux_units(experiment)
    cfg.set_config_paths(path_to_data=path_to_data, path_to_output=path_to_output,
        path_to_plots=path_to_plots, path_to_lists=path_to_lists)

    #### SET UP EXPERIMENT VALUES #####
    params = fsparam.Parameters('idsep', str_startdate, str_enddate, experiment, idsep_resume=idsep_resume)
    params.set_values(flux_type=flux_type, spacecraft=spacecraft,
        user_name=user_name, user_file=user_file,
        idsep_path=idsep_path, idsep_resume=idsep_resume,
        is_unixtime=is_unixtime,
        options=options, dointerp=dointerp,
        showplot=showplot, saveplot=saveplot,
        directory_depth=directory_depth,
        use_absolute_datapath=use_absolute_datapath,
        remove_above=remove_above,
        for_inclusive=for_inclusive,
        kurtosis_cut=kurtosis_cut,
        idsep_nsigma=idsep_nsigma,
        init_win=init_win,
        sliding_win=sliding_win,
        percent_points=percent_points,
        write_fluxes=write_fluxes)
    #################

    ########### RESUME###############
    #If RESUME, set end date of previous IDSEP run to startdate of this one
    if idsep_resume: #False if False or None
        lastdate = resume_fluxes(params, file_str='fluxes_*.csv', last_date_only=True)
        params.startdate = lastdate
    #################################

    print(f"IDSEP START DATE {params.startdate}")
    #Once params have been set up, need to use variables from params because
    #some processing and logic have been applied when creating the Parameters object
    eff_startdate = params.startdate
    params_cp = copy.deepcopy(params)

    #If the user entered a date range shorter than required for the
    #initial window used to identify the background, extend the date range
    #Note that the user should consider adding up to two months prior to the dates
    #of interest because the background solution for the first dates of the
    #timeseries are not accurate
    if not plot_timeseries_only and not params.idsep_resume:
        #Extend timeseries to cover init_win
        diff = (params.enddate - params.startdate).days
        if diff < params.init_win*2:
            eff_startdate = params.enddate - datetime.timedelta(days=params.init_win*2)
            
    params_cp.startdate = eff_startdate

    #READ IN FLUXES
    print("TIMESTAMP: Reading in flux files at time " + str(datetime.datetime.now()))
    #DOWNLOAD AND READ IN DATA
    dl_outpath, dl_plotpath, dates, fluxes, energy_bins, energy_bin_centers =\
        fsdl.get_data(params_cp, showplot=False)

    if plot_timeseries_only:
        sys.exit("Time series plot completed. Exiting.")
    
    #ITERATION 1: DEFINE AN INITIAL "MOVING" THRESHOLD W/DATE
    print("TIMESTAMP: Creating rough cut first guess at threshold at time " + str(datetime.datetime.now()))
 
    if params.idsep_resume:
        #Here, dates and fluxes will be extended backwards by one sliding window
        dates, fluxes, fluxes_bg_init, fluxes_high_init, resume_arrays = resume_cut(params, dates, fluxes, energy_bins)
    else:
        fluxes_bg_init, fluxes_high_init = rough_cut(params, dates, fluxes, energy_bins)
    
    
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
    is_final = False
    close_plot = True
    for iter in range(niter):
        print(f"TIMESTAMP: Performing iteration {iter}, {datetime.datetime.now()}")
        post = "_iter" + str(iter)
        if iter == niter-1:
            post += "_FINAL"
            is_final = True
            close_plot = False

        #Separate high and low flux by applying a sliding smoothing window to the background
        #fluxes_bg_init is used to get the mean, sigma, and threshold, then fluxes is split
        #into fluxes_bg and fluxes_high
        #The mean background is sensitive to the background selection in fluxes_bg_init
        #If resuming, the previously calculated mean background, sigmas, thresholds, and kurtosis
        #will be read in and prepended to the final values calculated for the current date range
        print(f"TIMESTAMP: Starting sliding window background calculation, {datetime.datetime.now()}")
        fluxes_bg, fluxes_high, ave_background, ave_sigma, threshold =\
            apply_sliding_window(params, dates, fluxes_bg_init, fluxes, energy_bins,
            iteration=iter, is_final=is_final, resume_arrays=resume_arrays)
        print(f"TIMESTAMP: Completed sliding window background calculation, {datetime.datetime.now()}")

        ########### RESUME###############
        #If RESUMING and at the last step, add in all of the results from the previous IDSEP run
        if is_final and params.idsep_resume:
            #Add in previous original dates and fluxes
            dates, fluxes = combine_resume_fluxes(params, dates, fluxes)

            #If resuming, set the startdate to the beginning of the previous time period
            #for trimming and writing out files.
            params.startdate = dates[0]
        ########### RESUME###############

        if params.showplot or params.saveplot:
            plt_tools.idsep_make_plots(str(params.sliding_win)+"window" + post, params, dates,
                fluxes_bg, energy_bins, dates, ave_background, ave_sigma, dates, threshold,
                showplot=False, close_plot=close_plot)
            
            plt_tools.idsep_make_plots(str(params.sliding_win)+"window_nosigma" + post,
                params, dates, fluxes_bg, energy_bins, dates, ave_background, ave_sigma,
                dates, threshold, close_plot=close_plot, disable_sigma=True) #disable sigma
            
            plt_tools.idsep_make_bg_sep_plot(str(params.sliding_win)+"window" + post, params,
                dates, fluxes_bg, fluxes_high, energy_bins, close_plot=close_plot)


        #Identify SEP events in full time range
        #SEP identification proves to be very robust
        #This plot for niter-1 is the same as FINALSEP. Will only be different if idsep
        #automatically extended the data set to accomodate the required
        #date lengths. Don't calculate and plot if it will be redundant.
        get_sep = True
        if is_final and dates[0] == params.startdate:
            get_sep = False
        if get_sep:
            global dwell_pts #to get value from tools and print to screen
            SEPstart, SEPend, fluxes_sep = analysis.identify_sep_above_background(dates,
                fluxes_high)

            if showplot or saveplot:
                plt_tools.idsep_make_bg_sep_plot("OnlySEP"+post, params, dates, fluxes,
                    fluxes_sep, energy_bins, close_plot=close_plot)


        #Taking the estimated background flux and remove SEP periods
        padding = 2 #number of days on either side of SEP start and end
        if iter <= niter-2:
            fluxes_bg_init, fluxes_sep_padded = separate_sep_with_dates(dates, fluxes, SEPstart, SEPend, padding)
 

    print(f"TIMESTAMP: Completed background and SEP event separation, {datetime.datetime.now()}")
    
    #Trim fluxes to the date range specified by the user
    trim_dates, trim_fluxes_high = datasets.extract_date_range(params.startdate, params.enddate, dates, fluxes_high)
    trim_dates, trim_fluxes = datasets.extract_date_range(params.startdate, params.enddate, dates, fluxes)
    trim_bg_dates, trim_ave_bg = datasets.extract_date_range(params.startdate, params.enddate, dates, ave_background)
    #Expect that this solution is exactly the same as the last one in the loop above,
    #UNLESS the dates need to be trimmed down
    SEPstart, SEPend, final_fluxes_sep = analysis.identify_sep_above_background(trim_dates, trim_fluxes_high)
    
    #Write SEP only fluxes to file
    write_sep_fluxes(params, trim_dates, final_fluxes_sep, trim_ave_bg, energy_bins)
    
    #Write start and end times to file
    write_sep_dates(params, energy_bins, SEPstart, SEPend)
    write_all_high_points(params, energy_bins, dates, fluxes_high)

    outputs = {
        "idsep_subdir": params.module_subdir,
        "idsep_outpath": params.module_outpath,
        "idsep_plotpath": params.module_plotpath,
    }

    outputs.update({"config": cfg.output_config()})
    outputs.update({"parameters": params.output_parameters()})

    stdtz = dh.time_to_zulu(str_startdate).replace(":","")
    outputs_fname = os.path.join(outputs["idsep_outpath"], f"{outputs['idsep_subdir']}.{stdtz}_idsep_outputs.json")
    ccmc_json.write_json(outputs,outputs_fname)

 
    if params.showplot or params.saveplot:
        plt_tools.idsep_make_bg_sep_plot("FINALSEP", params, trim_dates, trim_fluxes,
            final_fluxes_sep, energy_bins)
 
        if params.showplot: plt.show()

    #Clean up selected variables manually
    del dates
    del fluxes
    del trim_dates
    del trim_fluxes
    del trim_fluxes_high
    del final_fluxes_sep
    del trim_bg_dates
    del trim_ave_bg
    del fluxes_high

    print("TIMESTAMP: Completed idsep " + str(datetime.datetime.now()))
    return outputs
