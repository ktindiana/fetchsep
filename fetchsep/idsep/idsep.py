from ..utils import config as cfg
from ..utils import read_datasets as datasets
from ..utils import date_handler as dateh
from ..utils import define_background_idsep as defbg
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

__version__ = "0.1"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#See full program description in all_program_info() below
datapath = cfg.datapath
outpath = cfg.outpath + "/idsep"
plotpath = cfg.plotpath + "/idsep"

# Prepare directories
cfg.prepare_dirs()
for path in (outpath, plotpath):
    if not os.path.isdir(path):
        print("Making directory:", path)
        os.mkdir(path)

badval = cfg.badval #bad data points will be set to this value; must be negative
nsigma = cfg.idsep_nsigma #threshold mean + nsigma*sigma
init_win = cfg.init_win
sliding_win = cfg.sliding_win
nconsec = 0 #number of consecutive points that must be nonzero
allow_miss = 0 #number of points that can be zero when checking
                    #for SEP start
dwell_pts = 0 #number of points that can be missed after SEP starts
                #Dwell time


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
#Can be list of files that are continuous in time
 #      e.g. user_fname = ["file1.txt","file2.txt"]
user_fname = ['tmp.txt']



def about_idsep():
    """ About idsep.py
        
        Automatically identify SEP events by estimating background
        levels across the data set and separating backround
        and SEP fluxes.
        
        Option to perform background-subtraction of SEP fluxes.
        
        Identify start and end times of SEP events.
        
        Return arrays containing background fluxes only and
        SEP fluxes only.
        
    """

def error_check_inputs(startdate, enddate, experiment, flux_type):
    """ Check that all of the user inputs make sense and fall within bounds.
        
        INPUTS:
        
        :startdate: (datetime) - start of time period entered by user
        :enddate: (datetime) - end of time period entered by user
        :experiment: (string) - name of experiment specifed by user
        :flux_type: (string) - integral or differential
            
        OUTPUTS:
        
        None, but system exit if error found
    """
    #CHECKS ON INPUTS
    if (enddate < startdate):
        sys.exit('End time before start time! Enter a valid date range. '
                'Exiting.')
                
    if flux_type == "":
        sys.exit('User must indicate whether input flux is integral or '
                'differential. Exiting.')

    if (experiment == "SEPEM" and flux_type == "integral"):
        sys.exit('The SEPEM (RSDv2) data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    if (experiment == "SEPEMv3" and flux_type == "integral"):
        sys.exit('The SEPEM (RSDv3) data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    if ((experiment == "EPHIN" or experiment == "EPHIN_REleASE") \
        and flux_type == "integral"):
        sys.exit('The SOHO/EPHIN data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    sepem_end_date = datetime.datetime(2015,12,31,23,55,00)
    if(experiment == "SEPEM" and (startdate > sepem_end_date or
                   enddate > sepem_end_date)):
        sys.exit('The SEPEM (RSDv2) data set only extends to '
                  + str(sepem_end_date) +
            '. Please change your requested dates. Exiting.')

    sepemv3_end_date = datetime.datetime(2017,12,31,23,55,00)
    if(experiment == "SEPEMv3" and (startdate > sepemv3_end_date or
                   enddate > sepemv3_end_date)):
        sys.exit('The SEPEM (RSDv3) data set only extends to '
                  + str(sepemv3_end_date) +
            '. Please change your requested dates. Exiting.')



def sort_bin_order(all_fluxes, energy_bins):
    """Check the order of the energy bins. Usually, bins go from
        low to high energies, but some modelers or users may
        go in reverse order. Usually expect:
        [[10,20],[20,30],[30,40]]
        But user may instead input files with fluxes in order of:
        [[30,40],[20,30],[10,20]]
        
        This subroutine will reorder the fluxes and energy_bins
        to go in increasing order. If differential fluxes were input,
        this reordering will ensure that integral fluxes are
        estimated properly.
        
        INPUTS:
        
        :all_fluxes: (float nxm array) - fluxes for n energy channels
            and m time points
        :energy_bins: (float 2xn array) - energy bins for each of the
            energy channels
            
        OUTPUTS:
        
        :sort_fluxes: (float nxm array) - same as above, but sorted so
            that the lowest energy channel is first and highest is last
        :sort_bins: (float 2xn array) - same as above, but sorted
        
    """

    nbins = len(energy_bins)
    #Rank energy bins in order of lowest to highest effective
    #energies
    eff_en = []
    for i in range(nbins):
        if energy_bins[i][1] == -1:
            eff_en.append(energy_bins[i][0])
        else:
            midpt = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
            eff_en.append(midpt)
            
    eff_en_np = np.array(eff_en)
    sort_index = np.argsort(eff_en_np) #indices in sorted order
    
    sort_fluxes = np.array(all_fluxes)
    sort_bins = []
    for i in range(nbins):
        sort_fluxes[i] = all_fluxes[sort_index[i]]
        sort_bins.append(energy_bins[sort_index[i]])
    
    return sort_fluxes, sort_bins



def read_in_flux_files(experiment, flux_type, user_file, startdate,
        enddate, options, dointerp, is_unixtime):
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
                datasets.check_data(startdate, enddate, experiment, flux_type, user_file)
    else:
        filenames1, filenames2, filenames_orien = datasets.check_data(startdate,
                                enddate, experiment, flux_type, user_file)
                                    
    #read in flux files
    if experiment != "user":
        if experiment == "GOES":
            all_dates, all_fluxes, west_detector = \
                datasets.read_in_files(experiment, flux_type, \
                filenames1, filenames2, filenames_orien, options, detector)
        else:
            all_dates, all_fluxes, west_detector = \
                datasets.read_in_files(experiment, flux_type, \
                filenames1, filenames2, filenames_orien, options)
    
    if experiment == "user":
        all_dates, all_fluxes = datasets.read_in_user_files(filenames1,is_unixtime)
        west_detector = []
    
    print("+===ALL_FLUXES=====")
    print(len(all_fluxes))
    for ii in range(len(all_fluxes)):
        print(len(all_fluxes[ii]))
    
    #Define energy bins
    energy_bins = datasets.define_energy_bins(experiment, flux_type, \
                                west_detector, options)
    
    all_fluxes, energy_bins = sort_bin_order(all_fluxes, energy_bins)


    #Extract the date range specified by the user
    dates, fluxes = datasets.extract_date_range(startdate, enddate,
                            all_dates, all_fluxes)
    
    #Interpolate bad data with linear interpolation in time or set to None
    fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,dointerp)
    
        
    if len(dates) <= 1:
        sys.exit("The specified start and end dates were not present in the "
                "specified input file. Exiting.")
    
    return dates, fluxes, energy_bins


def determine_time_resolution(dates):
    """ The time resolution is found by taking the difference between
        every consecutive data point. The most common difference is
        taken as the time resolution. Even if the data set has gaps,
        if there are enough consecutive time points in the observational
        or model output, the correct time resolution should be identified.
        
        INPUTS:
        
        :dates: (datetime 1xm array) - dates associated with flux time profile
        
        OUTPUTS:
        
        :time_resolution: (time delta object)
        
    """
    ndates = len(dates)
    time_diff = [a - b for a,b in zip(dates[1:ndates],dates[0:ndates-1])]
    time_resolution = mode(time_diff)
    return time_resolution


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
    time_res = determine_time_resolution(dates)
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
    
    print("Requiring " + str(nconsec) + " points (" + str(time_increase/(60.*60.)) + " hours) to define an onset.")
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



def make_plots(unique_id, experiment, flux_type, exp_name, options, dates, fluxes,\
         energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates,\
        threshold, doBGSub, showplot, saveplot, disable_sigma=False):
    
    #Additions to titles and filenames according to user-selected options
    modifier = ''
    title_mod = ''
    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '


    figname = experiment + '_' + flux_type + modifier \
            + '_' + 'All_Bins_' + unique_id
    if experiment == 'user' and exp_name != '':
        figname = exp_name + '_' + flux_type + modifier \
                + '_' + 'All_Bins_' + unique_id
    
    fig = plt.figure(figname,figsize=(12,8))
    plt.rcParams.update({'font.size': 16})
    ax = plt.subplot(111)
    nbins = len(energy_bins)
    ifig = 0
    for i in range(nbins):
        if (i != 0 and not i%3) or nbins < 3:
            if saveplot:
                fig.savefig(plotpath + '/' +figname + '.png')
            figname_plt = figname + str(i)
            fig = plt.figure(figname_plt,figsize=(12,8))
            ax = plt.subplot(111)
            ifig = 0
    
        ax = plt.subplot(min(nbins,3), 1, ifig+1)
        ifig = ifig + 1
        legend_label = ""
        if energy_bins[i][1] != -1:
            legend_label = str(energy_bins[i][0]) + '-' \
                           + str(energy_bins[i][1]) + ' ' + energy_units
        else:
            legend_label = '>'+ str(energy_bins[i][0]) + ' ' + energy_units

        maskfluxes = np.ma.masked_less_equal(fluxes[i], 0)
        ax.plot_date(dates,maskfluxes,'.-',label=legend_label)
    
        if not disable_sigma:
            ax.errorbar(ave_dates, ave_fluxes[i],fmt='.', yerr=[ave_sigma[i][0], ave_sigma[i][1]], label="ave " + legend_label,zorder=100)
        if disable_sigma:
            ax.errorbar(ave_dates, ave_fluxes[i],fmt='-',
                label="ave " + legend_label,zorder=100)
        
        ax.plot_date(threshold_dates,threshold[i],'-',label="threshold\n" + legend_label, zorder=200)
        
        flux_units = ''
        if flux_type == "integral": flux_units = flux_units_integral
        if flux_type == "differential": flux_units = flux_units_differential
        
        if i==0:
            plt.title(experiment + ' '+ title_mod + ' ' + unique_id)
            if experiment == 'user' and exp_name != '':
                plt.title(exp_name + ' '+ title_mod + ' ' + unique_id)
        plt.xlabel('Date')
      #  plt.ylabel('Flux [' + flux_units + ']')
        plt.ylabel(r'Flux (MeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
        plt.yscale("log")
        fig.autofmt_xdate(rotation=45)
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05), fontsize=11)
    
        if saveplot and i == nbins-1:
            fig.savefig(plotpath + '/' +figname_plt + '.png')
    
    if not showplot:
        plt.close(fig)



def make_timeseries_plot(unique_id, experiment, flux_type, exp_name,\
        options, dates, fluxes, energy_bins, doBGSub, showplot, saveplot):

    #Additions to titles and filenames according to user-selected options
    modifier = ''
    title_mod = ''
    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '


    figname = experiment + '_' + flux_type + modifier \
            + '_' + 'All_Bins_' + unique_id
    if experiment == 'user' and exp_name != '':
        figname = exp_name + '_' + flux_type + modifier \
                + '_' + 'All_Bins_' + unique_id
    
    fig = plt.figure(figname,figsize=(12,8))
    plt.rcParams.update({'font.size': 16})
    nbins = len(energy_bins)
    ifig = 0
    for i in range(nbins):
        if i!= 0 and not i%3:
            if saveplot:
                fig.savefig(plotpath + '/' +figname + '.png')
            figname_plt = figname + str(i)
            fig = plt.figure(figname_plt,figsize=(12,8))
            ax = plt.subplot(111)
            ifig = 0
    
        ax = plt.subplot(min(nbins,3), 1, ifig+1)
        ifig = ifig + 1
        legend_label = ""
        if energy_bins[i][1] != -1:
            legend_label = str(energy_bins[i][0]) + '-' \
                           + str(energy_bins[i][1]) + ' ' + energy_units
        else:
            legend_label = '>'+ str(energy_bins[i][0]) + ' ' + energy_units

        maskfluxes = np.ma.masked_less_equal(fluxes[i], 0)
        ax.plot_date(dates,maskfluxes,'.-',label=legend_label)
    
        
        flux_units = ''
        if flux_type == "integral": flux_units = flux_units_integral
        if flux_type == "differential": flux_units = flux_units_differential
        
        if i==0:
            plt.title(experiment + ' '+ title_mod + ' ' + unique_id)
            if experiment == 'user' and exp_name != '':
                plt.title(exp_name + ' '+ title_mod + ' ' + unique_id)
            #plt.ylabel('Flux [' + flux_units + ']')
            plt.ylabel(r'Flux (MeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
        
        plt.xlabel('Date')
        
        plt.yscale("log")
        fig.autofmt_xdate(rotation=45)
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05),fontsize=11)
    
        if saveplot and i == nbins-1:
            figname_plt = figname + str(i)
            fig.savefig(plotpath + '/' +figname_plt + '.png')
    
    if not showplot:
        plt.close(fig)



def make_bg_sep_plot(unique_id, experiment, flux_type, exp_name, options,\
            dates, fluxes_bg, fluxes_sep, energy_bins, doBGSub,
            showplot, saveplot):
    
    #Additions to titles and filenames according to user-selected options
    modifier = ''
    title_mod = ''
    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '


    figname = experiment + '_' + flux_type + modifier \
            + '_' + 'SEP_BG_' + unique_id
    if experiment == 'user' and exp_name != '':
        figname = exp_name + '_' + flux_type + modifier \
                + '_' + 'SEP_BG_' + unique_id
    
    fig = plt.figure(figname,figsize=(12,8))
    plt.rcParams.update({'font.size': 16})
    ax = plt.subplot(111)
    nbins = len(energy_bins)
    ifig = 0
    for i in range(nbins):
        if i != 0 and not i%3:
            if saveplot:
                fig.savefig(plotpath + '/' +figname + '.png')
            figname = figname + str(i)
            fig = plt.figure(figname,figsize=(12,8))
            ax = plt.subplot(111)
            ifig = 0
    
        ax = plt.subplot(min(nbins,3), 1, ifig+1)
        ifig = ifig + 1
        legend_label = ""
        if energy_bins[i][1] != -1:
            legend_label = str(energy_bins[i][0]) + '-' \
                           + str(energy_bins[i][1]) + ' ' + energy_units
        else:
            legend_label = '>'+ str(energy_bins[i][0]) + ' ' + energy_units

    
        maskfluxes_bg = np.ma.masked_less_equal(fluxes_bg[i], 0)
        ax.plot_date(dates,maskfluxes_bg,'.-',label="bg " + legend_label)
        maskfluxes_sep = np.ma.masked_less_equal(fluxes_sep[i], 0)
        ax.plot_date(dates,maskfluxes_sep,'.-',label="sep " + legend_label, zorder=100)
        
        
        flux_units = ''
        if flux_type == "integral": flux_units = flux_units_integral
        if flux_type == "differential": flux_units = flux_units_differential
        
        if i==0:
            plt.title(experiment + ' '+ title_mod + ' ' + unique_id)
            if experiment == 'user' and exp_name != '':
                plt.title(exp_name + ' '+ title_mod + ' ' + unique_id)
        plt.xlabel('Date')
        #plt.ylabel('Flux [' + flux_units + ']')
        plt.ylabel(r'Flux (MeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
        plt.yscale("log")
        fig.autofmt_xdate(rotation=45)
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05),fontsize=11)
    
        if saveplot and i == nbins-1:
            fig.savefig(plotpath + '/' +figname + '.png')
    
    if not showplot:
        plt.close(fig)




def make_diff_plot(unique_id, experiment, flux_type, exp_name, options, dates,\
            diff_fluxes, ave_sigma, energy_bins, doBGSub, showplot, saveplot):
    
    #Additions to titles and filenames according to user-selected options
    modifier = ''
    title_mod = ''
    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '


    figname = experiment + '_' + flux_type + modifier \
            + '_' + 'Diff_' + unique_id
    if experiment == 'user' and exp_name != '':
        figname = exp_name + '_' + flux_type + modifier \
                + '_' + 'Diff_' + unique_id
    
    fig = plt.figure(figname,figsize=(12,8))
    plt.rcParams.update({'font.size': 16})
    ax = plt.subplot(111)
    nbins = len(energy_bins)
    ifig = 0
    for i in range(nbins):
        thresh = np.multiply(ave_sigma[i][1],nsigma)
        if i != 0 and not i%3:
            if saveplot:
                fig.savefig(plotpath + '/' +figname + '.png')
            figname = figname + str(i)
            fig = plt.figure(figname,figsize=(12,8))
            ax = plt.subplot(111)
            ifig = 0
    
        ax = plt.subplot(min(3,nbins), 1, ifig+1)
        ifig = ifig + 1
        legend_label = ""
        if energy_bins[i][1] != -1:
            legend_label = str(energy_bins[i][0]) + '-' \
                           + str(energy_bins[i][1]) + ' ' + energy_units
        else:
            legend_label = '>'+ str(energy_bins[i][0]) + ' ' + energy_units

        ax.plot_date(dates,diff_fluxes[i],'.',label="diff " + legend_label)
        ax.plot_date(dates,thresh,'-',label="threshold " + legend_label, zorder=100)
        
        
        flux_units = ''
        if flux_type == "integral": flux_units = flux_units_integral
        if flux_type == "differential": flux_units = flux_units_differential
        
        if i==0:
            plt.title(experiment + ' '+ title_mod + ' ' + unique_id\
                    + "\nDiff = Flux - Mean BG")
            if experiment == 'user' and exp_name != '':
                plt.title(exp_name + ' '+ title_mod + ' ' + unique_id\
                    + "\nDiff = Flux - Mean BG")
        plt.xlabel('Date')
        #plt.ylabel('Flux [' + flux_units + ']')
        plt.ylabel(r'Flux (MeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
        fig.autofmt_xdate(rotation=45)
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05),fontsize=11)
 
 
        if saveplot and i == nbins-1:
            fig.savefig(plotpath + '/' +figname + '.png')
    
    if not showplot:
        plt.close(fig)




def get_bg_sep(threshold, dates, fluxes):
    fluxes_bg = []
    fluxes_sep = []
    for i in range(len(fluxes)):
        flux_below, flux_above = defbg.separate_with_threshold(threshold[i],\
                dates,fluxes[i])
        fluxes_bg.append(flux_below)
        fluxes_sep.append(flux_above)

    return fluxes_bg, fluxes_sep





def write_sep_dates(experiment, exp_name, flux_type, energy_bins, options,
                    str_stdate, str_enddate, remove_above, SEPstart, SEPend,
                    doBGSub, for_inclusive=False):
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
    modifier = ''
    title_mod = ''
    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '


    prename = 'SEPTimes_' + experiment + '_' + flux_type + modifier
    if experiment == 'user' and exp_name != '':
        prename = 'SEPTimes_' + exp_name + '_' + flux_type + modifier
    
    #####WRITE SEP DATES OUT TO FILE INSTEAD OF PRINTING##########
    for j in range(len(SEPstart)):
        fname = outpath + "/" + prename + '_' + str(energy_bins[j][0]) + '_to_'\
                + str(energy_bins[j][1]) + '.txt'
        outfile = open(fname,"w")
        outfile.write("#SEP times calculated by SEPAutoID\n")
        if experiment == "user" and exp_name != '':
            outfile.write("#Experiment: " + exp_name + "\n")
        else:
            outfile.write("#Experiment: " + experiment + "\n")
        outfile.write("#Flux type: " + flux_type + "\n")
        outfile.write("#Energy channel: " + str(energy_bins[j][0]) \
                    + " - " + str(energy_bins[j][1]) + "\n")
        outfile.write("#Selected options: " + str(options) + "\n")
        outfile.write("#Searched date range: " + str_stdate + " to "\
                    + str_enddate + "\n")
        outfile.write("#User applied an initial cut of " + str(remove_above)
                + ", initial averaging window of " + str(init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(sliding_win) + " days\n")
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
                SEPst = SEPst - one_sec
                SEPed = SEPed - one_sec
            outfile.write(str(SEPst) + " " + str(SEPed) + "\n")
            
        outfile.close()

    
    

def write_all_high_points(experiment, exp_name, flux_type, energy_bins, options,
                    str_stdate, str_enddate, remove_above, dates, fluxes_high,
                    doBGSub, for_inclusive=False):
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
    time_res = determine_time_resolution(dates)
    if for_inclusive: time_res = time_res - datetime.timedelta(seconds=1)
    
    #Additions to titles and filenames according to user-selected options
    modifier = ''
    title_mod = ''
    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '


    prename = 'HighPoints_' + experiment + '_' + flux_type + modifier
    if experiment == 'user' and exp_name != '':
        prename = 'HighPoints_' + exp_name + '_' + flux_type + modifier
    
    #####WRITE SEP DATES OUT TO FILE INSTEAD OF PRINTING##########
    for j in range(len(fluxes_high)):
        fname = outpath + "/" + prename + '_' + str(energy_bins[j][0]) + '_to_'\
                + str(energy_bins[j][1]) + '.txt'
        outfile = open(fname,"w")
        outfile.write("#All high points above mean background + 3*sigma calculated by SEPAutoID\n")
        if experiment == "user" and exp_name != '':
            outfile.write("#Experiment: " + exp_name + "\n")
        else:
            outfile.write("#Experiment: " + experiment + "\n")
        outfile.write("#Flux type: " + flux_type + "\n")
        outfile.write("#Energy channel: " + str(energy_bins[j][0]) \
                    + " - " + str(energy_bins[j][1]) + "\n")
        outfile.write("#Selected options: " + str(options) + "\n")
        outfile.write("#Searched date range: " + str_stdate + " to "\
                    + str_enddate + "\n")
        outfile.write("#User applied an initial cut of " + str(remove_above)
                + ", initial averaging window of " + str(init_win) + " days"
                + ", final threshold with a sliding window of "
                + str(sliding_win) + " days\n")
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




def run_all(str_startdate, str_enddate, experiment,
        flux_type, exp_name, user_file, is_unixtime, options, doBGSub, dointerp,
        remove_above, for_inclusive, plot_timeseries_only, showplot, saveplot):
    """ Run all the steps to do background and SEP separation.
    
    """
    
    startdate = dateh.str_to_datetime(str_startdate)
    enddate = dateh.str_to_datetime(str_enddate)
    
    error_check_inputs(startdate, enddate, experiment, flux_type)
        
    #READ IN FLUXES
    dates, fluxes, energy_bins = read_in_flux_files(experiment,
        flux_type, user_file, startdate, enddate, options, dointerp,is_unixtime)
            
            
    if plot_timeseries_only:
        unique_id = "FluxTimeSeries"
        make_timeseries_plot(unique_id, experiment, flux_type, exp_name,\
        options, dates, fluxes, energy_bins, doBGSub, showplot, saveplot)
        if showplot:
            plt.show()
        sys.exit("Time series plot completed. Exiting.")
    
    #DEFINE AN INITIAL "MOVING" THRESHOLD W/DATE IN FIRST ITERATION
    ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold =\
                defbg.ndays_average(init_win, dates, fluxes, nsigma, remove_above)
    
    #INITIAL SEPARATION OF BG AND HIGHER THAN BG
    fluxes_bg, fluxes_high = get_bg_sep(threshold,dates,fluxes)
    
    #ITERATE AGAIN AND GET BACKGROUND VALUE FOR EVERY DAY USING A
    #BACKWORD SMOOTHING WINDOW
    ave_background2, ave_sigma2, threshold2, diff_fluxes =\
                defbg.backward_window_background(sliding_win, dates, fluxes_bg, nsigma)
    
    diff_fluxes2 = [[]]*len(ave_background2)
    for i in range(len(ave_background2)):
        diff_fluxes2[i] = fluxes[i] - ave_background2[i]
        fluxes_nozero, bad_index = defbg.remove_zero_one(fluxes[i])
        for bad in bad_index:
            diff_fluxes2[i][bad] = None #zero fluxes not included in plot
    
    for i in range(len(fluxes_bg)):
        if None in fluxes_bg[i]:
            print("None values present in second: in " + str(i))
    fluxes_bg2, fluxes_high2 = get_bg_sep(threshold2,dates,fluxes)
    
    
    #Identify SEP events
    SEPstart, SEPend, fluxes_sep = identify_sep(dates, fluxes_high2)
    
    #Write start and end times to file
    write_sep_dates(experiment, exp_name, flux_type, energy_bins, \
                    options, str_startdate, str_enddate, remove_above, \
                    SEPstart, SEPend, doBGSub, for_inclusive)
    write_all_high_points(experiment, exp_name, flux_type, energy_bins, \
                    options, str_startdate, str_enddate, remove_above,\
                    dates, fluxes_high2, doBGSub, for_inclusive)


    if showplot or saveplot:
        make_plots(str(init_win)+"days", experiment, flux_type, exp_name, options, dates, fluxes, energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold, doBGSub, showplot, saveplot)
        
        make_bg_sep_plot(str(init_win)+"days", experiment, flux_type, exp_name, options, dates, fluxes_bg, fluxes_high, energy_bins, doBGSub,
            showplot, saveplot)
        
        make_plots(str(sliding_win)+"window", experiment, flux_type, exp_name, options, dates, fluxes_bg, energy_bins, dates, ave_background2, ave_sigma2, dates, threshold2, doBGSub, showplot, saveplot)
        
        make_plots(str(sliding_win)+"window_nosigma", experiment, flux_type, exp_name, options, dates, fluxes_bg, energy_bins, dates, ave_background2, ave_sigma2, dates, threshold2, doBGSub, showplot, saveplot,
            True) #disable_sigma

        
        make_bg_sep_plot(str(sliding_win)+"window", experiment, flux_type, exp_name, options, dates, fluxes_bg2, fluxes_high2, energy_bins, doBGSub,
            showplot, saveplot)
        
        make_bg_sep_plot("OnlySEP", experiment, flux_type, exp_name, options, dates, fluxes, fluxes_sep, energy_bins, doBGSub, showplot, saveplot)
        
        make_diff_plot("Background-Subtracted", experiment, flux_type, exp_name, options, dates,diff_fluxes2, ave_sigma2, energy_bins, doBGSub,
            showplot, saveplot)
 
        
        if showplot: plt.show()

