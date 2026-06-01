from ..utils import config as cfg
from ..utils import directories as dirs
from ..utils import read_datasets as datasets
from ..utils import download as fsdl
from ..utils import date_handler as dh
from ..utils import analysis
from ..utils import define_background_idsep as defbg
from ..utils import plotting_tools as plt_tools
from ..utils import error_check
from ..utils import tools
from ..utils import names
from ..utils import experiments as expts
import datetime
from datetime import timedelta
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import pandas as pd

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
        flux_below, flux_above = defbg.separate_with_threshold(threshold[i],\
                dates,fluxes[i])
        fluxes_bg.append(flux_below)
        fluxes_high.append(flux_above)
        
    fluxes_bg = np.array(fluxes_bg)
    fluxes_high = np.array(fluxes_high)

    return fluxes_bg, fluxes_high





def write_sep_dates(experiment, exp_name, flux_type, energy_bins, options,
    str_stdate, str_enddate, remove_above, SEPstart, SEPend,
    for_inclusive=False, spacecraft='', modifier='', savepath=''):
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
        :modifier: (string) modifier for path and filenames
        :savepath: (string) path to save output file
            
        OUTPUTS:
        
        none, but writes out n files - i.e. an SEP event list for each
            energy channel
        
    """
    
    one_sec = datetime.timedelta(seconds=1)
    
    name = names.idsep_naming_scheme(experiment, flux_type, exp_name, modifier=modifier)
    prename = (f"SEPTimes_{name}")
    zstdate = dh.time_to_zulu(str_stdate)
    zstdate = zstdate.replace(':', '')
    zenddate = dh.time_to_zulu(str_enddate)
    zenddate = zenddate.replace(':', '')
    
    #####WRITE SEP DATES OUT TO FILE##########
    for j in range(len(SEPstart)):
        fname = (f"{prename}_{zstdate}_{zenddate}_{energy_bins[j][0]}_to_{energy_bins[j][1]}.txt")
        fname = os.path.join(savepath, fname)
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
                + ", initial averaging window of " + str(cfg.init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(cfg.sliding_win) + " days\n")
        outfile.write("#The percentage of points in the sliding window that must be "
                "background was " + str(cfg.percent_points) + "\n")
        outfile.write("#Threshold defined as mean background + " + str(cfg.idsep_nsigma)
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
    for_inclusive=False, spacecraft="", modifier='', savepath=''):
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
    time_res = analysis.determine_time_resolution(dates)
    if for_inclusive: time_res = time_res - datetime.timedelta(seconds=1)
    
    #Additions to titles and filenames according to user-selected options
    name = names.idsep_naming_scheme(experiment, flux_type, exp_name, modifier=modifier)

    prename = (f"HighPoints_{name}")
    zstdate = dh.time_to_zulu(str_stdate)
    zenddate = dh.time_to_zulu(str_enddate)
    
    #####WRITE SEP DATES OUT TO FILE INSTEAD OF PRINTING##########
    for j in range(len(fluxes_high)):
        fname = (f"{prename}_{zstdate}_{zenddate}_{energy_bins[j][0]}_to_{energy_bins[j][1]}.txt")
        fname = os.path.join(savepath, fname)
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
                + ", initial averaging window of " + str(cfg.init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(cfg.sliding_win) + " days\n")
        outfile.write("#The percentage of points in the sliding window that must be "
                "background was " + str(cfg.percent_points) + "\n")
        outfile.write("#Threshold defined as mean background + " + str(cfg.idsep_nsigma)
                + " x sigma\n")
        outfile.write("#High flux points were identified when flux values "
                    "exceeded mean background flux + 3*sigma "
                    "by applying a " + str(cfg.sliding_win) + " days "
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



def rough_cut(dates, fluxes, remove_above, energy_bins, experiment,
        flux_type, exp_name, options, showplot, saveplot, spacecraft="",
        modifier='', title_mod='', id_outpath='', id_plotpath=''):
    """ Remove fluxes above remove_above.
        Group fluxes into time periods of init_win.
        For each time period, calculate the mean background and sigma.
        Return the mean background, sigma, and threshold values of
        mean + nsigma*sigma.
        
        
        INPUTS:
        
            :dates: (1xn datetime array)
            :fluxes: (mxn float array) fluxes for m energy bins, n dates
            :remove_above: (float) and initial cut on flux applied to all
                energy bins. Flux > remove_above is considered an enhancement.
            :id_outpath: (string) path to save output files
            :id_plotpath: (string) path to save plots
        
        OUTPUTS:
        
            :fluxes_bg: (mxn float array) background fluxes for m energy bins, n dates
            :fluxes_high: (mxn float array) enhanced fluxes for m energy bins, n dates
            
    """
    print("Preforming an initial rough cut between the background and enhanced fluxes.")
    print(f"Init win: {cfg.init_win}, Nsigma: {cfg.idsep_nsigma}, Remove above: {remove_above}")
    ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold =\
                defbg.ndays_average_optimized(cfg.init_win, dates, fluxes, energy_bins,
                cfg.idsep_nsigma, remove_above, savepath=id_outpath)
    
    #INITIAL SEPARATION OF BG AND HIGHER THAN BG
    fluxes_bg, fluxes_high = get_bg_high(threshold,dates,fluxes)
    
    if showplot or saveplot:
        plt_tools.idsep_make_plots(str(cfg.init_win)+"days", experiment, flux_type, exp_name, options, dates, fluxes, energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold, showplot, saveplot, spacecraft=spacecraft, modifier=modifier, title_mod=title_mod, savepath=id_plotpath)
        
        plt_tools.idsep_make_bg_sep_plot(str(cfg.init_win)+"days", experiment, flux_type, exp_name, options, dates, fluxes_bg, fluxes_high, energy_bins, showplot, saveplot, spacecraft=spacecraft,
            modifier=modifier, title_mod=title_mod, savepath=id_plotpath)

    
    return fluxes_bg, fluxes_high


def apply_sliding_window(dates, fluxes_bg_in, fluxes, energy_bins, iteration=0,
    is_final=False, modifier='', title_mod='', savepath=''):
    """ Identify the background value for every day using a sliding window.
        Use the initial estimated background from rough_cut and refine
        by applying a sliding window of sliding_win days to get a background
        value and sigma for every day of the data set. Create a daily threshold
        and apply to extract background fluxes, fluxes_bg, and enhanced fluxes,
        fluxes_high.
        
        INPUTS:
        
            :dates: (1xn datetime array)
            :fluxes_bg_in: (mxn float array) background fluxes for m energy bins,
                n dates
            :fluxes: (mxn float array) all fluxes for m energy bins, n dates
                
        OUTPUTS:
        
            :fluxes_bg: (mxn float array) background fluxes for m energy bins, n dates
            :fluxes_high: (mxn float array) enhanced fluxes for m energy bins, n dates
    
    """
    ave_background, ave_sigma, threshold =\
            defbg.backward_window_background_optimized(cfg.sliding_win, dates, fluxes_bg_in,
            energy_bins, cfg.idsep_nsigma, iteration, is_final=is_final, savepath=savepath)
    
    for i in range(len(fluxes_bg_in)):
        if None in fluxes_bg_in[i]:
            print("None values present in second: in " + str(i))
    
    fluxes_bg, fluxes_high = get_bg_high(threshold,dates,fluxes)

    return fluxes_bg, fluxes_high, ave_background, ave_sigma, threshold


def write_sep_fluxes(dates, fluxes, fluxes_bg, energy_bins, savepath=''):
    """ Write out final SEP fluxes and bg-subtracted fluxes.
        Subtract fluxes (e.g. SEP fluxes subtracted by the mean background).
        If fluxes already at a value of zero, no background subtraction is
        performed.
        
    """
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
    
    defbg.write_df(df,'SEP_fluxes_FINAL', savepath=savepath)
    
    df[cols] = df[cols] - df_bg[cols]
    for col in cols:
        df.loc[df[col] < 0, col] = 0
    
    defbg.write_df(df,'SEP_fluxes_background-subtracted_FINAL',savepath=savepath)



def run_idsep(str_startdate, str_enddate, experiment,
    flux_type=None, spacecraft="",
    exp_name=None, user_file=None,
    directory_depth=1,
    is_unixtime=False, options=None, dointerp=False,
    remove_above=999999, for_inclusive=False,
    plot_timeseries_only=False,
    showplot=False, saveplot=False,
    write_fluxes=True,
    path_to_data=None,
    path_to_output=None,
    path_to_plots=None,
    path_to_lists=None,
    kurtosis_cut=None,
    idsep_nsigma=None,
    init_win=None,
    sliding_win=None,
    percent_points=None):
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

    if flux_type == '' or flux_type == None:
        flux_type = expts.get_flux_type(experiment)
    expts.set_config_kurtosis_cut(experiment, flux_type)
    expts.set_config_energy_units(experiment)
    expts.set_config_flux_units(experiment)
    cfg.set_config_paths(path_to_data=path_to_data, path_to_output=path_to_output,
        path_to_plots=path_to_plots, path_to_lists=path_to_lists)
    cfg.configure_idsep(kurtosis_cut=kurtosis_cut, idsep_nsigma=idsep_nsigma,
        init_win=init_win, sliding_win=sliding_win, percent_points=percent_points)
    cfg.print_configured_values()

    dirs.check_paths(experiment)


    #### SET UP EXPERIMENT VALUES #####
    #If user specifies a spacecraft but isn't relevant to experiment,
    #overrides and sets spacecraft to ''
    spacecraft = expts.get_spacecraft(experiment, spacecraft)

    #Check for empty dates
    if (str_startdate == "" or str_enddate == ""):
        sys.exit('You must enter start and end dates. Exiting.')

    if experiment != 'user':
        exp_info = expts.experiment_info(experiment)

    #################

    startdate = dh.str_to_datetime(str_startdate)
    enddate = dh.str_to_datetime(str_enddate)
    eff_startdate = startdate
    
    options = options.split(";")
    #If the user entered a date range shorter than required for the
    #initial window used to identify the background, extend the date range
    #Note that the user should consider adding up to two months prior to the dates
    #of interest because the background solution for the first dates of the
    #timeseries are not accurate
    if not plot_timeseries_only:
        #Extend timeseries to cover init_win
        diff = (enddate - startdate).days
        if diff < cfg.init_win*2:
            eff_startdate = enddate - datetime.timedelta(days=cfg.init_win*2)
        
    
    error_check.error_check_options(experiment, flux_type, options, False, spacecraft=spacecraft)
    error_check.error_check_inputs(startdate, enddate, experiment, flux_type, subroutine='idsep')

    # Create directory names and prepare directories
    modifier, title_mod = names.setup_modifiers(options, spacecraft=spacecraft)
    subdir = names.idsep_naming_scheme(experiment, flux_type, exp_name, modifier=modifier)
    id_outpath = dirs.create_subdirectories(cfg.outpath, module='idsep',
        subdir=subdir, directory_depth=directory_depth)
    id_plotpath = dirs.create_subdirectories(cfg.plotpath, module='idsep',
        subdir=subdir, directory_depth=directory_depth)


    #READ IN FLUXES
    print("TIMESTAMP: Reading in flux files at time " + str(datetime.datetime.now()))
    #DOWNLOAD AND READ IN DATA
    dl_outpath, dl_plotpath, dates, fluxes, energy_bins, energy_bin_centers =\
        fsdl.get_data(eff_startdate, enddate, experiment, flux_type=flux_type,
            exp_name=exp_name, user_file=user_file, is_unixtime=is_unixtime,
            spacecraft=spacecraft, directory_depth=directory_depth, module='idsep',
            options=options, dointerp=dointerp, showplot=plot_timeseries_only, saveplot=saveplot,
            write_fluxes=write_fluxes, dl_outpath=id_outpath, dl_plotpath=id_plotpath,
            path_to_data=cfg.datapath, path_to_output=cfg.outpath,
            path_to_plots=cfg.plotpath, path_to_lists=cfg.listpath)

    if plot_timeseries_only:
        sys.exit("Time series plot completed. Exiting.")
    
    #ITERATION 1: DEFINE AN INITIAL "MOVING" THRESHOLD W/DATE
    print("TIMESTAMP: Creating rough cut first guess at threshold at time " + str(datetime.datetime.now()))
    fluxes_bg_init, fluxes_high_init = rough_cut(dates, fluxes,
        remove_above, energy_bins, experiment, flux_type, exp_name, options,
        False, saveplot, spacecraft=spacecraft, modifier=modifier, title_mod=title_mod,
        id_outpath=id_outpath, id_plotpath=id_plotpath)
    
    
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
        print(f"TIMESTAMP: Starting sliding window background calculation, {datetime.datetime.now()}")
        fluxes_bg, fluxes_high, ave_background, ave_sigma, threshold =\
            apply_sliding_window(dates, fluxes_bg_init, fluxes, energy_bins,
            iteration=iter, is_final=is_final, savepath=id_outpath)
        print(f"TIMESTAMP: Completed sliding window background calculation, {datetime.datetime.now()}")


        if showplot or saveplot:
            plt_tools.idsep_make_plots(str(cfg.sliding_win)+"window" + post, experiment, flux_type, exp_name, options, dates, fluxes_bg_init, energy_bins, dates, ave_background, ave_sigma, dates, threshold, False, saveplot, spacecraft=spacecraft, close_plot=close_plot,
                modifier=modifier, title_mod=title_mod, savepath=id_plotpath)
            
            plt_tools.idsep_make_plots(str(cfg.sliding_win)+"window_nosigma" + post,
                experiment, flux_type,
                exp_name, options, dates, fluxes_bg_init, energy_bins, dates,
                ave_background, ave_sigma, dates, threshold, showplot,
                saveplot, True, spacecraft=spacecraft, close_plot=close_plot,
                modifier=modifier, title_mod=title_mod, savepath=id_plotpath) #disable sigma
            
            plt_tools.idsep_make_bg_sep_plot(str(cfg.sliding_win)+"window" + post, experiment,
                flux_type, exp_name, options, dates, fluxes_bg, fluxes_high, energy_bins,
                showplot, saveplot, spacecraft=spacecraft, close_plot=close_plot,
                modifier=modifier, title_mod=title_mod, savepath=id_plotpath)


        #Identify SEP events in full time range
        #SEP identification proves to be very robust
        #This plot for niter-1 is the same as FINALSEP. Will only be different if idsep
        #automatically extended the data set to accomodate the required
        #date lengths. Don't calculate and plot if it will be redundant.
        get_sep = True
        if is_final and dates[0] == startdate:
            get_sep = False
        if get_sep:
            global dwell_pts #to get value from tools and print to screen
            SEPstart, SEPend, fluxes_sep = analysis.identify_sep_above_background(dates, fluxes_high)

            if showplot or saveplot:
                plt_tools.idsep_make_bg_sep_plot("OnlySEP"+post, experiment, flux_type,
                    exp_name, options, dates, fluxes, fluxes_sep, energy_bins,
                    showplot, saveplot, spacecraft=spacecraft, close_plot=close_plot,
                    modifier=modifier, title_mod=title_mod, savepath=id_plotpath)


        #Taking the estimated background flux and remove SEP periods
        padding = 2 #number of days on either side of SEP start and end
        if iter <= niter-2: #Final round
            fluxes_bg_init, fluxes_sep_padded = separate_sep_with_dates(dates, fluxes, SEPstart, SEPend, padding)
 

    print(f"TIMESTAMP: Completed background and SEP event separation, {datetime.datetime.now()}")
    

    #Trim fluxes to the date range specified by the user
    trim_dates, trim_fluxes_high = datasets.extract_date_range(startdate, enddate, dates, fluxes_high)
    trim_dates, trim_fluxes = datasets.extract_date_range(startdate, enddate, dates, fluxes)
    trim_bg_dates, trim_ave_bg = datasets.extract_date_range(startdate, enddate, dates, ave_background)
    #Expect that this solution is exactly the same as the last one in the loop above,
    #UNLESS the dates need to be trimmed down
    SEPstart, SEPend, final_fluxes_sep = analysis.identify_sep_above_background(trim_dates, trim_fluxes_high)
    
    #Write background-subtracted SEP only fluxes to file
    write_sep_fluxes(trim_dates, final_fluxes_sep, trim_ave_bg, energy_bins, savepath=id_outpath)
    
    #Write start and end times to file
    write_sep_dates(experiment, exp_name, flux_type, energy_bins,
                    options, str_startdate, str_enddate, remove_above,
                    SEPstart, SEPend, for_inclusive, spacecraft=spacecraft,
                    modifier=modifier, savepath=id_outpath)
    write_all_high_points(experiment, exp_name, flux_type, energy_bins,
                    options, str_startdate, str_enddate, remove_above,
                    dates, fluxes_high, for_inclusive, spacecraft=spacecraft,
                    modifier=modifier, savepath=id_outpath)

 
    if showplot or saveplot:
        plt_tools.idsep_make_bg_sep_plot("FINALSEP", experiment, flux_type, exp_name, options,
            trim_dates, trim_fluxes, final_fluxes_sep, energy_bins, showplot, saveplot,
            spacecraft=spacecraft, modifier=modifier, title_mod=title_mod, savepath=id_plotpath)
 
        
        if showplot: plt.show()

    print("TIMESTAMP: Completed idsep " + str(datetime.datetime.now()))
    return id_outpath, id_plotpath
