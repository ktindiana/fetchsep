from ..utils import config as cfg
from ..utils import read_datasets as datasets
from ..utils import date_handler as dateh
from ..utils import plotting_tools as plt_tools
from ..utils import error_check
from ..utils import tools
import datetime
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#See full program description in all_program_info() below
#datapath = cfg.datapath
outpath = cfg.outpath
plotpath = cfg.plotpath



""" About download.py
    
    Download data from internet sources into FetchSEP's
    specified data directory.
    
    File presence and completeness is managed by FetchSEP's
    data manager.
    
"""



def read_in_flux_files(experiment, flux_type, startdate,
        enddate, options, dointerp, write_fluxes=False,
        spacecraft="", user_file="", exp_name=""):
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
        filenames1, filenames2, filenames_orien, detector = datasets.check_data(startdate,
                enddate, experiment, flux_type, user_file, spacecraft=spacecraft)
    else:
        filenames1, filenames2, filenames_orien = datasets.check_data(startdate,
                enddate, experiment, flux_type, user_file, spacecraft=spacecraft)
                            
    #read in flux files
    if experiment != "user":
        #Combine integral channels for all GOES spacecraft together into
        #a long time series. Only works for integral channels, since
        #GOES differential channels differ across experiments.
        if experiment == "GOES":
            all_dates, all_fluxes, west_detector, energy_bins = \
                datasets.read_in_files(experiment, flux_type,
                filenames1, filenames2, filenames_orien, options, detector,
                spacecraft=spacecraft)
        else:
            all_dates, all_fluxes, west_detector = \
                datasets.read_in_files(experiment, flux_type,
                filenames1, filenames2, filenames_orien, options)
    
    if experiment == "user":
        all_dates, all_fluxes = datasets.read_in_user_files(filenames1,is_unixtime)
        west_detector = []
    
    if len(all_fluxes) == 0:
        sys.exit("Could not read in flux files. Check for bad date ranges or missing files.")
    
    #Define energy bins
    #ERNE energy bins depend on the time period of the experiment.
    #For idsep, it is acceptable to include fluxes across slightly
    #different energy channels for the purposes of SEP identification.
    if experiment == "ERNE":
        version = datasets.which_erne(startdate, enddate)
        energy_bins, energy_bin_centers = datasets.define_energy_bins(version, flux_type,
                                west_detector, options)
    elif experiment != "GOES":
        energy_bins, energy_bin_centers = datasets.define_energy_bins(experiment, flux_type,
                                west_detector, options, spacecraft=spacecraft)


    if energy_bins == None:
        sys.exit("Could not identify energy bins for experiment " + experiment
                + " and fluxtype " + flux_type)

    all_fluxes, energy_bins, energy_bin_centers = tools.sort_bin_order(all_fluxes, energy_bins)


    #Extract the date range specified by the user
    dates, fluxes = datasets.extract_date_range(startdate, enddate,
                            all_dates, all_fluxes)
    
    #Interpolate bad data with linear interpolation in time or set to None
    print("read_in_flux_files: Checking for bad data.")
    if dointerp:
        print("Performing interpolation with time.")
    else:
        print("Setting bad values to NaN.")
    fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,dointerp)
    
        
    if len(dates) <= 1:
        sys.exit("The specified start and end dates were not present in the "
                "specified input file. Exiting.")

    if write_fluxes:
        dir = tools.idsep_naming_scheme(experiment, flux_type, exp_name, options, spacecraft=spacecraft)
        path = os.path.join(outpath, 'idsep', dir)
        if not os.path.isdir(path):
            print("Making directory:", path)
            os.mkdir(path)
            
        tools.write_fluxes(experiment, flux_type, exp_name, options, energy_bins, dates, fluxes,
            "idsep", spacecraft=spacecraft)

    return dates, fluxes, energy_bins




def get_data(str_startdate, str_enddate, experiment,
        flux_type, options, dointerp, showplot, saveplot,
        write_fluxes=True, spacecraft="", path_to_data=''):
    """ Download data. Create an output file of all fluxes in the
        specified date range.
    
        INPUTS:
        
        :str_startdate: (str) start date in string format
        :str_enddate: (str) end date in string format
        :experiment: (str) name of experiment native to FetchSEP
        :flux_type: (str) integral or differential, mainly matters for GOES
        :options: (str) options that may be used with GOES data
        :dointerp: (bool) will interpolate bad data linearly with time
        :showplot: (bool) set True to show plots on screen
        :saveplot: (bool) set True to save plots in plotpath
        :write_fluxes: (bool) Write fluxes to csv file after read in and processed 
            for bad points
        :spacecraft: (string) primary or secondary if exp_name = GOES_RT
        :path_to_data: (string) if set to a value, will be used as the location to
            store downloaded data

    
    """
    cfg.configure_for(experiment)

    if path_to_data != '':
        cfg.set_datapath(path_to_data)
#        global datapath
#        datapath = path_to_data
#        print(f"get_data: Setting datapath to {datapath}")

    cfg.print_configured_values()

    # Prepare directories
    cfg.prepare_dirs()
    for path in (outpath, plotpath):
        if not os.path.isdir(path):
            print("Making directory:", path)
            os.mkdir(path)

    startdate = dateh.str_to_datetime(str_startdate)
    enddate = dateh.str_to_datetime(str_enddate)
    
    options = options.split(";")

    error_check.error_check_options(experiment, flux_type, options, False, spacecraft=spacecraft)
    error_check.error_check_inputs(startdate, enddate, experiment, flux_type)
    datasets.check_paths()


    #READ IN FLUXES
    dates, fluxes, energy_bins = read_in_flux_files(experiment,
        flux_type, startdate, enddate, options, dointerp,
        write_fluxes=write_fluxes, spacecraft=spacecraft)
            
            
    if showplot or saveplot:
        unique_id = "FluxTimeSeries"
        exp_name = ''
        plt_tools.idsep_make_timeseries_plot(unique_id, experiment, flux_type, exp_name,
        options, dates, fluxes, energy_bins, False, showplot, saveplot, spacecraft=spacecraft)
        if showplot:
            plt.show()
    
