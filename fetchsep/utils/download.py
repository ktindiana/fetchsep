from ..utils import config as cfg
from ..utils import directories as dirs
from ..utils import read_datasets as datasets
from ..utils import date_handler as dh
from ..utils import plotting_tools as plt_tools
from ..utils import error_check
from ..utils import tools
from ..utils import names
from ..utils import experiments as expts
from ..utils import parameters as fsparam
import datetime
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


""" About download.py
    
    Download data from internet sources into FetchSEP's
    specified data directory.
    
    File presence and completeness is managed by FetchSEP's
    data manager.
    
"""


def read_in_flux_files(params):
    """ Read in the appropriate data or user files. Trims to dates
        between start time and end time. Interpolates bad
        points with linear interpolation if dointerp True.
        
        INPUTS:
        
        :params: (FetchSEP Parameters object)
        
        OUTPUTS:
        
        :dates: (datetime 1xm array) - times in flux time profile trimmed
            between startdate and enddate
        :fluxes: (numpy float nxm array) - fluxes for n energy channels and m
            time steps; these are background subtracted fluxes if background
            subtraction was selected.
        :energy_bins: (array nx2 for n thresholds)
        
    """
    detector= []
    west_detector = []

    ##### Check if data is available and download
    if params.experiment == "GOES":
        filenames1, filenames2, filenames_orien, detector = datasets.check_data(params)
    else:
        filenames1, filenames2, filenames_orien = datasets.check_data(params)
                            
    ###### Read in flux files
    if params.experiment != "user":
        #Combine integral channels for all GOES spacecraft together into
        #a long time series. Only works for integral channels, since
        #GOES differential channels differ across experiments.
        if params.experiment == "GOES":
            all_dates, all_fluxes, west_detector, energy_bins, energy_bin_centers = \
                datasets.read_in_files(params, filenames1, filenames2, filenames_orien, detector)
        else:
            all_dates, all_fluxes, west_detector = \
                datasets.read_in_files(params, filenames1, filenames2, filenames_orien, detector)
    
    if params.experiment == "user":
        all_dates, all_fluxes = datasets.read_in_user_files(filenames1, params.is_unixtime)
 
 
    if len(all_fluxes) == 0:
        sys.exit("download: Could not read in flux files. Check for bad date ranges or missing files. Exiting.")
    
    #Define energy bins
    #ERNE energy bins depend on the time period of the experiment.
    if params.experiment == "SOHO_ERNE":
        version = datasets.which_erne(params.startdate, params.enddate)
        params_cp = copy.deepcopy(params)
        params_cp.experiment = version
        energy_bins, energy_bin_centers = datasets.define_energy_bins(params_cp, west_detector)
    elif params.experiment != "GOES":
        energy_bins, energy_bin_centers = datasets.define_energy_bins(params, west_detector=west_detector)

    if energy_bins == None:
        sys.exit("Could not identify energy bins for experiment " + params.experiment
                + " and fluxtype " + params.flux_type)

    all_fluxes, energy_bins, energy_bin_centers = tools.sort_bin_order(all_fluxes, energy_bins, energy_bin_centers=energy_bin_centers)

    #Extract the date range specified by the user
    dates, fluxes = datasets.extract_date_range(params.startdate, params.enddate, all_dates, all_fluxes)

    #Interpolate bad data with linear interpolation in time or set to None
    print("read_in_flux_files: Checking for bad data.")
    if params.do_interpolation:
        print("Performing interpolation with time.")
    else:
        print("Setting bad values to NaN.")
    fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,params.do_interpolation)
    
        
    if len(dates) <= 1:
        sys.exit("The specified start and end dates were not present in the "
                "specified input file. Exiting.")

    return dates, fluxes, energy_bins, energy_bin_centers


def load_parameters(str_startdate, str_enddate, experiment,
    flux_type=None, spacecraft=None, user_name=None, user_file=None, is_unixtime=None,
    directory_depth=None, options=None, dointerp=None,
    showplot=None, saveplot=None, write_fluxes=None,
    path_to_data=None, path_to_output=None, path_to_plots=None,
    path_to_lists=None, use_absolute_datapath=None):
    """ Create FetchSEP Parameters object """

    cfg.set_config_paths(path_to_data=path_to_data, path_to_output=path_to_output,
        path_to_plots=path_to_plots, path_to_lists=path_to_lists)

    #### SET UP EXPERIMENT VALUES #####
    params = fsparam.Parameters('download', str_startdate, str_enddate, experiment)
    params.set_values(flux_type=flux_type, spacecraft=spacecraft, user_name=user_name,
        user_file=user_file, is_unixtime=is_unixtime, options=options, dointerp=dointerp,
        showplot=showplot, saveplot=saveplot, directory_depth=directory_depth,
        write_fluxes=write_fluxes, use_absolute_datapath=use_absolute_datapath)

    return params


def get_data(params,
    showplot=None,
    saveplot=None,
    path_to_data=None,
    path_to_output=None,
    path_to_plots=None,
    path_to_lists=None,
    format='dict'):
    """ Download data. Create an output file of all fluxes in the
        specified date range.
    
        INPUTS:
        
        :str_startdate: (str) start date in string format
        :str_enddate: (str) end date in string format
        :experiment: (str) name of experiment native to FetchSEP
        :flux_type: (str) integral or differential, mainly matters for GOES
        :directory_depth: (int) default = 2; Subdirectories for output files may be 
                supressed by choosing the directory depth.
                0 - Files output to top directories: cfg.outpath (output/), cfg.plotpath (plots/) level; 
                1 - Files output to subdirectory at module level, cfg.outpath/module (output/opsep); 
                2 - Files output to subdirectory named according to experiment and 
                options, e.g. cfg.outpath/module/subdir (output/opsep/GOES-13_integral/
        :options: (str) options that may be used with GOES data
        :dointerp: (bool) will interpolate bad data linearly with time
        :module: (string) fetchsep module - idsep, opsep, download (default=download)
        :showplot: (bool) set True to show plots on screen
        :saveplot: (bool) set True to save plots in plotpath
        :write_fluxes: (bool) Write fluxes to csv file after read in and processed 
            for bad points
        :spacecraft: (string) primary or secondary if exp_name = GOES-RT
        :user_energy_bins: (2d array) [[10,-1], [30,-1],...] in case the user wants to override
            default or config energy bins
        :dl_outpath: (string) directly specify the path for output files
        :dl_plotpath: (string) directly specify the path for plots
        :path_to_data: (string) if set to a value, will be used as the location to
            store downloaded data
        :format: (string) 'dict' to get dictionary with list of dates and energy_bins 
            with fluxes in numpy array; 'dataframe' or 'df' to get timeseries in 
            pandas DataFrame.

    
    """

    cfg.set_config_paths(path_to_data=path_to_data, path_to_output=path_to_output,
        path_to_plots=path_to_plots, path_to_lists=path_to_lists)
    cfg.print_configured_values()

    if showplot == None: showplot = params.showplot
    if saveplot == None: saveplot = params.saveplot

    #READ IN FLUXES
    dates, fluxes, energy_bins, energy_bin_centers = read_in_flux_files(params)
  
    if params.write_fluxes:
        fluxes_filename = tools.write_fluxes(params, energy_bins, dates, fluxes, suffix="original_fluxes")
 
    if showplot or saveplot:
        unique_id = "FluxTimeSeries"
        plt_tools.idsep_make_timeseries_plot(unique_id, params, dates, fluxes, energy_bins)
        if showplot:
            plt.show()

    return params.module_outpath, params.module_plotpath, dates, fluxes, energy_bins, energy_bin_centers
