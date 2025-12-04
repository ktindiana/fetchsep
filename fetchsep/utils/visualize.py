import os
import sys
import matplotlib.pyplot as plt
from . import config as cfg
from . import tools
from . import plotting_tools as plt_tools
from ..json import ccmc_json_handler as ccmc_json
import pandas as pd
import numpy as np


""" Visualize outputs produced by FetchSEP. """

def col_to_bin(col):
    """ Convert string in column headers to energy bins 

        The columns are labeled with the energy bins, like:
            10--1 ([10,-1] for >10 MeV integral channel)
            9.4-15.9 ([9.4,15.9] differential channel)    
            
        INPUTS:
        
            :col: (str) like 10--1 or 9.4-15.9
            
        OUTPUTS:
        
            :energy_bin: (list) like [10,-1] or [9.4,15.9]
            
    """
    energy_bin = []
    
    #No or empty bin provided
    if col == '' or col == None:
        energy_bin = []

    #Integral channel
    elif '--1' in col:
        col = col.strip().split('--1')
        energy_bin = [float(col[0]),-1]
    
    #Differential channel
    else:
        col = col.strip().split('-')
        energy_bin = [float(col[0]), float(col[1])]
        
    return energy_bin



def read_fluxes_all_bins(filename):
    """ Read files produced by FetchSEP that have
        datetime (YYYY-MM-DD HH:MM:SS) in the first column and 
        fluxes for multiple energy bins in the next columns.
        
        The columns are labeled with the energy bins, like:
            10--1 ([10,-1] for >10 MeV integral channel)
            9.4-15.9 ([9.4,15.9] differential channel)

        INPUTS:
        
            :filename: (string) Full path to the file
            
        OUTPUTS:
        
            :df: (pandas DataFrame) time series
            :energy_bins: (list of energy bins for each flux column)
        
    """
    if not os.path.isfile(filename):
        print(f"read_fluxes_all_bins: Cannot find file {filename}. Exiting")
        sys.exit()
        
    df = pd.read_csv(filename)
    
    #The 0th column should always be datetime
    cols = df.columns.tolist()
    df[cols[0]] = pd.to_datetime(df[cols[0]])
    
    energy_bins = []
    for col in cols[1:]:
        energy_bin = col_to_bin(col)
        energy_bins.append(energy_bin)
    
    return df, energy_bins


def read_ccmc_time_series(filename):
    """ Format specified in the CCMC SEP Scoreboard json schema
        as a time series file that accompanies the json. 
        zulu time in first column and flux in second column.
        
    """
    
    if not os.path.isfile(filename):
        print(f"read_ccmc_time_series: Cannot find file {filename}. Exiting")
        sys.exit()

    dates = []
    fluxes = []
    with open(filename,'r') as file:
        for line in file:
            line = line.strip().split()
            if line == '': continue
            date = ccmc_json.zulu_to_time(line[0])
            flux = float(line[1])
            dates.append(date)
            fluxes.append(flux)
            
    return dates, fluxes


def choose_flux_units(flux_type):
    flux_units = ''
    if flux_type == 'integral':
        flux_units = cfg.flux_units_integral
    if flux_type == 'differential':
        flux_units = cfg.flux_units_differential
        
    return flux_units



def setup_labels(experiment='', flux_type='', options='',
    bgsub=None, spacecraft=''):
    
    doBGSubOPSEP = False
    doBGSubIDSEP = False
    if bgsub == 'OPSEP':
        doBGSubOPSEP = True
    if bgsub == 'IDSEP':
        doBGSubIDSEP = True

    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP)

    flux_units = choose_flux_units(flux_type)

    ylabel = f"Flux [${flux_units}$]"
    ylabel = plt_tools.make_math_label(ylabel)

    return title_mod, ylabel



def plot_fluxes_all_bins(filename, experiment='', flux_type='',
    options='', bgsub=None, spacecraft='', ylog = True):
    """ Plot files produced by FetchSEP that have
        datetime (YYYY-MM-DD HH:MM:SS) in the first column and 
        fluxes for multiple energy bins in the next columns.
        
        The columns are labeled with the energy bins, like:
            10--1 ([10,-1] for >10 MeV integral channel)
            9.4-15.9 ([9.4,15.9] differential channel)
    
        INPUTS:
        
            :filename: (string) Full path to the file
            :experiment: (string) used for plot labeling
            :flux_type: (string) integral or differential, used 
                for plot labeling
            :options: (string) options applied to data separated by 
                a semicolon, e.g. 'S14;Bruno2017;uncorrected', used
                for plot labeling
            bgsub: (string or None) 
                None if no background subtraction performed
                IDSEP if used IDSEP background
                OPSEP if used background subtraction feature in OPSEP
                Used for plot labeling
            :spacecraft: (string) primary or secondary, if relevant
            :ylog: (bool) make the y-axis logarithmic
    
        OUTPUTS:
        
            Plot of time series
    
    """
    
    df, energy_bins = read_fluxes_all_bins(filename)
    cols = df.columns.tolist()
    dates = df[cols[0]]

    title_mod, ylabel = setup_labels(experiment=experiment,
        flux_type=flux_type, options=options, bgsub=bgsub,
        spacecraft=spacecraft)

    flux_units = choose_flux_units(flux_type)

    plot_title = f"All Energy Bins\n {experiment} {title_mod} {flux_type}"

    fig = plt.figure(figsize=(13.5,8))
    ax = plt.subplot(111)
    
    for i, col in enumerate(cols[1:]):
        bin_label = tools.setup_energy_bin_label(energy_bins[i])
        fluxes = np.array(df[col].to_list())
        maskfluxes = np.ma.masked_where(fluxes <= 0, fluxes)
        plt.plot(dates, maskfluxes, label=bin_label, marker='.')
    
    plt.xlabel("Date")
    plt.ylabel(ylabel)
    plt.title(plot_title)
    if ylog: plt.yscale("log")
    plt.legend(loc="upper right")
    plt.grid(True, which="both")
    


def plot_event_definition(filename, experiment='', flux_type='', options='',
    bgsub=None, spacecraft='', energy_bin=[], threshold=np.nan,
    sep_start_time=pd.NaT, sep_end_time=pd.NaT,
    max_flux=np.nan, max_time=pd.NaT,
    onset_peak=np.nan, onset_time=pd.NaT, ylog = True):
    """ Plot any threshold crossings, start and end times, max flux, 
        and onset peak, if available.
    
        Plot is for a single energy channel in files like 
        GOES-06_integral_enhance_idsep.1986-01-01T000000Z.10.0.MeV.txt
        
        with zulu time in first column and flux in second column.
    
    """

    dates, fluxes = read_ccmc_time_series(filename)

    title_mod, ylabel = setup_labels(experiment=experiment,
        flux_type=flux_type, options=options, bgsub=bgsub,
        spacecraft=spacecraft)

    flux_units = choose_flux_units(flux_type)

    bin_label = tools.setup_energy_bin_label(energy_bin)
    plot_title = f"SEP Event Definition {bin_label}, {threshold} {flux_units}\n {experiment} {title_mod} {flux_type}"

    fig = plt.figure(figsize=(13.5,8))
    ax = plt.subplot(111)
    
    fluxes = np.array(fluxes)
    maskfluxes = np.ma.masked_where(fluxes <= 0, fluxes)
    plt.plot(dates, maskfluxes, label=bin_label, marker='.')
    
    if not pd.isnull(threshold):
        ax.axhline(threshold,color='red',linestyle=':', label="Threshold")
    if not pd.isnull(sep_start_time):
        ax.axvline(sep_start_time,color='black',linestyle=':')
    if not pd.isnull(sep_end_time):
        ax.axvline(sep_end_time,color='black',linestyle=':',label="Start, End")

    if not pd.isnull(onset_peak) and not pd.isnull(onset_time):
        ax.plot_date(onset_time,onset_peak,'o',color="black",
                label="Onset Peak")
    if not pd.isnull(max_flux) and not pd.isnull(max_time):
        ax.plot_date(max_time,max_flux,'ro',mfc='none', label="Max Flux")

    plt.xlabel("Date")
    plt.ylabel(ylabel)
    plt.title(plot_title)
    if ylog: plt.yscale("log")
    plt.legend(loc="upper right")
    plt.grid(True, which="both")


def cast_time_columns(df):
    
    cols = df.columns.tolist()
    for col in cols:
        if 'Time' in col and 'Series' not in col:
            df[col] = pd.to_datetime(df[col])
            
    return df


def sep_column_label(flux_type, energy_bin, threshold):
    """ Create the column names used in the SEP lists """

    threshold_units = choose_flux_units(flux_type)
    energy_units = cfg.energy_units
    
    if energy_bin[1] == -1:
        channel_label = f">{float(energy_bin[0])} {energy_units}"
    else:
        channel_label = f"{float(energy_bin[0])}-{float(energy_bin[1])} {energy_units}"

    threshold_label = f"{float(threshold)} {threshold_units}"

    label = f"{channel_label} {threshold_label}"
    
    return label


def find_column(cols, id):
    """ Provided a unique id containing part of the column name,
        return the desired column name.
        
        Do any filtering of column labels ahead of time so
        that the id will give a unique result.
        
    """

    sel_col = ''
    for col in cols:
        if id in col:
            sel_col = col

    return sel_col


def extract_sep_info(df, index, sep_cols):
    """ Extract SEP start and end time, onset peak, max flux. 
        sep_cols include the column labels for only one event
        definition. e.g. all cols with >10.0 MeV 10 pfu in label
    
    """
    flux_type = df['Flux Type'].values[index]
    flux_units = choose_flux_units(flux_type)

    col = find_column(sep_cols, 'SEP Start Time')
    sep_start_time = df[col].values[index]

    col = find_column(sep_cols, 'SEP End Time')
    sep_end_time = df[col].values[index]

    col = find_column(sep_cols, f"Onset Peak ({flux_units})")
    onset_peak = df[col].values[index]

    col = find_column(sep_cols, f"Onset Peak Time")
    onset_time = df[col].values[index]

    col = find_column(sep_cols, f"Max Flux ({flux_units})")
    max_flux = df[col].values[index]

    col = find_column(sep_cols, f"Max Flux Time")
    max_time = df[col].values[index]

    return sep_start_time, sep_end_time, onset_peak, onset_time, max_flux, max_time


def derive_path(sep_list, time_series_filename):
    """ OpSEP outputs filenames with same name as directory in
        the filename:
        GOES-06_integral_enhance_idsep.1986-02-03T070500Z.10.0.MeV.txt
        
        will be in the directory:
        output/opsep/GOES-06_integral_enhance_idsep/
        
    """
    path = os.path.dirname(sep_list)
    time_series_filename = time_series_filename.strip().split('.')
    subdir = time_series_filename[0]
    
    if subdir in path:
        return path
    
    if 'output' not in path:
        path = os.path.join(path, 'output')
        
    if 'opsep' not in path:
        path = os.path.join(path,'opsep')
        
    path = os.path.join(path, subdir)

    return path



def plot_list(sep_list, energy_bin, threshold):
    """ SEP event list or non-event list produced after batch running
        opsep. e.g. CLEAR Benchmark dataset or single experiment lists:
        GOES_integral_PRIMARY.1986-02-03.2025-09-10_sep_events.csv
        GOES-06_integral_enhance_idsep.1986-01-01.1994-11-30_non_events.csv
        GOES-06_integral_enhance_idsep.1986-01-01.1994-11-30_sep_events.csv
    
        sep_list is the full path and filename
        
        INPUTS:
        
            :sep_list: full path name to SEP or non-event list
            :energy_bin: (list) [10,-1]
            :threshold: (float) 10
            
        OUTPUTS:
        
            Make plots
    
    """

    if not os.path.isfile(sep_list):
        print(f"plot_list: Cannot find file {sep_list}. Exiting")
        sys.exit()

    df = pd.read_csv(sep_list)
    df = cast_time_columns(df)
    
    cols = df.columns.tolist()
    
    #step through list and make plots
    for index, row in df.iterrows():
        experiment = row['Experiment']
        flux_type = row['Flux Type']
        options = row['Options']
        if pd.isnull(options): options = ''
        options = options.strip().split(';')
        bgsub = row['Background Subtraction']

        sep_label = sep_column_label(flux_type, energy_bin, threshold)
        print(f"Plotting SEP information for {sep_label}")
        #Pull all relevant SEP columns
        sep_cols = [col for col in cols if sep_label in col]
        if len(sep_cols) == 0: continue


        #plot the event definition and any derived values if there
        #was an SEP event
        col = find_column(sep_cols, 'Flux Time Series')
        event_def_filename = row[col]
        path = derive_path(sep_list, event_def_filename)
        event_def_filename = os.path.join(path, event_def_filename)

        sep_start_time, sep_end_time, onset_peak, onset_time, max_flux, max_time = extract_sep_info(df, index, sep_cols)

        plot_event_definition(event_def_filename, experiment=experiment, flux_type=flux_type, options=options,
            bgsub=bgsub, energy_bin=energy_bin, threshold=threshold,
            sep_start_time=sep_start_time, sep_end_time=sep_end_time,
            max_flux=max_flux, max_time=max_time,
            onset_peak=onset_peak, onset_time=onset_time)


        #Plot all fluxes in the timeframe
        all_bins_filename = row['All Fluxes Time Series']
        path = derive_path(sep_list, all_bins_filename)
        all_bins_filename = os.path.join(path, all_bins_filename)
        plot_fluxes_all_bins(all_bins_filename, experiment=experiment,
            flux_type=flux_type, options=options, bgsub=bgsub)
        
        

        plt.show()
