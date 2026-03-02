from . import config as cfg
from . import goes
from . import experiments as expts
from . import date_handler as dh
import pandas as pd
import datetime
import os
import sys
import math
import numpy as np

def sort_bin_order(all_fluxes, energy_bins, energy_bin_centers=[]):
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
        :energy_bin_centers: (float 1xn array) - corresponding energy bin centers
            
        OUTPUTS:
        
        :sort_fluxes: (float nxm array) - same as above, but sorted so
            that the lowest energy channel is first and highest is last
        :sort_bins: (float 2xn array) - same as above, but sorted
        
    """

    nbins = len(energy_bins)
    #Rank energy bins in order of lowest to highest effective
    #energies
    if len(energy_bin_centers) == 0:
        for i in range(nbins):
            if energy_bins[i][1] == -1:
                energy_bin_centers.append(energy_bins[i][0])
            else:
                midpt = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
                energy_bin_centers.append(midpt)
            
    eff_energies = np.array(energy_bin_centers)
    sort_index = np.argsort(eff_energies) #indices in sorted order

    sort_fluxes = np.array(all_fluxes)
    sort_bins = []
    sort_centers = []
    for i in range(nbins):
        sort_fluxes[i] = all_fluxes[sort_index[i]]
        sort_bins.append(energy_bins[sort_index[i]])
        sort_centers.append(energy_bin_centers[sort_index[i]])
    
    return sort_fluxes, sort_bins, sort_centers


##### NAMING SCHEMA #####
def setup_modifiers(options, spacecraft="", doBGSubOPSEP=False, doBGSubIDSEP=False,
        OPSEPEnhancement=False, IDSEPEnhancement=False):
    """ Add modifier strings according to options.
    
    """
    modifier = '' #for appending to filenames
    title_mod = '' #for appending to plot titles

    doBGSub = False
    if doBGSubOPSEP == True or doBGSubIDSEP == True:
        doBGSub = True

    module = ''
    if doBGSubOPSEP or OPSEPEnhancement:
        module = 'opsep'
    if doBGSubIDSEP or IDSEPEnhancement:
        module = 'idsep'

    if "uncorrected" in options:
        modifier = modifier + '_uncor'
        title_mod = title_mod + 'uncorrected '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_B17'
        title_mod = title_mod + 'Bruno2017 '
    if spacecraft:
        modifier = modifier + '_' + spacecraft
        title_mod = title_mod + spacecraft + ' '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if IDSEPEnhancement or OPSEPEnhancement:
        modifier = modifier + '_enhance'

    if module != '':
        modifier = modifier + '_' + module
        title_mod = title_mod + ' (' + module + ')'

    return modifier, title_mod


#####Refer to these subroutines to set units everywhere
def get_energy_units():
    energy_units = cfg.energy_units
    return energy_units

def get_flux_units(flux_type):
    if flux_type == "integral": flux_units = cfg.flux_units_integral
    if flux_type == "differential": flux_units = cfg.flux_units_differential
    return flux_units

def get_flux_units_bin(energy_bin):
    if energy_bin[1] == -1:
        flux_units = cfg.flux_units_integral
    else:
        flux_units = cfg.flux_units_differential
        
    return flux_units

def get_fluence_units(flux_type):
    if flux_type == "integral": fluence_units = cfg.fluence_units_integral
    if flux_type == "differential": fluence_units = cfg.fluence_units_differential
    return fluence_units
    
def get_fluence_units_bin(energy_bin):
    if energy_bin[1] == -1:
        fluence_units = cfg.fluence_units_integral
    else:
        fluence_units = cfg.fluence_units_differential
        
    return fluence_units
#########################


def setup_energy_bin_label(energy_bin):
    """ Label for a single energy bin.
    
    """
    label = ""
    energy_units = get_energy_units()
    if energy_bin[1] != -1:
        label = (f"{energy_bin[0]}-{energy_bin[1]} {energy_units}")
    else:
        label = (f">{energy_bin[0]} {energy_units}")

    return label


def energy_bin_key(bin):
    """ Create key for dataframe or columns header in dataframe and
        csv files.
        
    """
    return f"{bin[0]}-{bin[1]}"


def idsep_naming_scheme(experiment, flux_type, exp_name, options, spacecraft=''):
    """ Create naming scheme for subfolders in output/idsep 
    
        Used in:
        SEPTimes files
        HighPoints files
        subdirectories in output
        
    """

    modifier, title_mod = setup_modifiers(options, spacecraft=spacecraft)
    name = (f"{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        name = (f"{exp_name}_{flux_type}{modifier}")
    
    return name
    

def opsep_subdir(experiment, flux_type, exp_name, options,
    spacecraft='', doBGSubOPSEP=False, doBGSubIDSEP=False, OPSEPEnhancement=False,
    IDSEPEnhancement=False):

    modifier, title_mod = setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    #Return a corresponding directory name
    dir = (f"{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        dir = (f"{exp_name}_{flux_type}{modifier}")
 
    return dir


def opsep_naming_scheme(date, suffix, experiment, flux_type, exp_name, options,
    spacecraft='', doBGSubOPSEP=False, doBGSubIDSEP=False, OPSEPEnhancement=False,
    IDSEPEnhancement=False):
    """ Create naming scheme for subfolders in output/idsep """

    tzulu = dh.time_to_zulu(date)
    tzulu = tzulu.replace(":","")

    modifier, title_mod = setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    name = (f"{tzulu}_{experiment}_{flux_type}{modifier}_{suffix}")
    if experiment == 'user' and exp_name != '':
        name = (f"{tzulu}_{exp_name}_{flux_type}{modifier}_{suffix}")

    #Return a corresponding directory name
    dir = opsep_subdir(experiment, flux_type, exp_name, options, spacecraft=spacecraft,
        doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    return name, dir


def write_fluxes(experiment, flux_type, exp_name, options, energy_bins, dates, fluxes,
    module, spacecraft="", doBGSubOPSEP=False, doBGSubIDSEP=False, OPSEPEnhancement=False,
    IDSEPEnhancement=False, suffix=''):
    """ Write dates, fluxes to a standard format csv file with datetime in the 
        first column and fluxes in the remaining column. The energy bins will 
        be indicated in the file comments.
        
        INPUTS:

        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        :options: (string array) options that may be applied to GOES data
        :dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :fluxes: (float nxm array) fluxes for n energy channels and m
            time points
        :module: (string) "idsep" or "opsep"
        
        OUTPUTS:
            
            csv file written to outpath
                 
    """
    name = ''
    stdate = dates[0].strftime("%Y%m%d")
    enddate = dates[-1].strftime("%Y%m%d")
    if module == 'idsep' or module == 'download':
        name = idsep_naming_scheme(experiment, flux_type, exp_name, options,
                    spacecraft=spacecraft)
        fname = (f"fluxes_{name}_{stdate}_{enddate}.csv")
        fname = os.path.join(cfg.outpath, module, name, fname)
    elif module == 'opsep':
        #name like the json schema
        name, dir = opsep_naming_scheme(dates[0], suffix, experiment, flux_type, exp_name,
            options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP,
            doBGSubIDSEP=doBGSubIDSEP, OPSEPEnhancement=OPSEPEnhancement,
            IDSEPEnhancement=IDSEPEnhancement)
        tzulu = dh.time_to_zulu(dates[0])
        tzulu = tzulu.replace(":","")
        fname = f"{dir}.{tzulu}.csv"
        if suffix != '':
            fname = f"{dir}.{tzulu}_{suffix}.csv"
        fname = os.path.join(cfg.outpath, module, dir, fname)
        
    keys = []
    for bin in energy_bins:
        keys.append(energy_bin_key(bin))
    
    dict = {"dates":dates}
    for i in range(len(fluxes)):
        dict.update({keys[i]:fluxes[i]})
        
    df = pd.DataFrame(dict)
    df.to_csv(fname, index=False)
    print("Wrote " + fname + " to file.")

    return fname


def create_primary_goes_sep_list(lists, prefix='GOES', path_to_data=None,
    path_to_output=None, path_to_plots=None, path_to_lists=None,):
    """ Provided a set of SEP events lists created by batch running
        opsep for multiple GOES spacecraft, create a single list by extracting 
        the SEP events from the primary GOES spacecraft, when available. 
        
        The lists are named like: 
        cfg.outpath/opsep/*/GOES-06_integral_enhance_idsep.1986-01-01.1994-11-30_sep_events.csv
        
        The current version of this subroutine is good for short periods
        of time, like the duration of SEP events, for which the primary
        satellite is unlikely to change. The primary satellite on the
        date of the Time Period Start will be identified as the GOES primary.
    
        INPUT:
        
            :lists: (list of strings) List of the full path to each sep_events
                list that will be read; or filename of a list containing lists
            :path_to_data: (string) path where satellite data should be downloaded. Will default to 
                datapath listed in fetchsep.cfg if a value is not specified.
            :path_to_output: (string) path where output files should be saved. Will default to 
                outpath listed in fetchsep.cfg if a value is not specified.
            :path_to_plots: (string) path where plots should be saved. Will default to
                plotpath listed in fetchsep.cfg if a value is not specified.
            :path_to_lists: (string) path where lists should be saved. Will default to
                listpath listed in fetchsep.cfg if a value is not specified.
    
        OUTPUT:
        
            :df_primary: (pandas dataframe) containing a compiled list of one entry
                per SEP event from only the GOES primary spacecraft
            
            outputs a file to cfg.outpath/opsep/GOES_PRIMARY.YYYY-MM-DD.YYYY-MM-DD_sep_events.csv
    
    """
    cfg.set_config_paths(path_to_data=path_to_data, path_to_output=path_to_output,
        path_to_plots=path_to_plots, path_to_lists=path_to_lists)
        
    goes_R = expts.goes_R()

    df = pd.DataFrame()

    if isinstance(lists,str):
        arr = make_lists_array(lists)
        lists = arr

    for list in lists:
        list = list.strip()
        if not os.path.isfile(list):
            print(f"create_primary_goes_list: File does not exist. Skipping. {list}")
            continue
        
        if 'GOES-RT' in list and 'primary' not in list:
            #Take GOES-RT from the primary spacecraft only, which will be
            #indicated in the filename
            print(f"create_primary_goes_list: Real-time GOES list is not for the primary spacecraft. Skipping. {list}")
            continue
        print(f"processing {list}")
        fpath = os.path.dirname(list)
        df_in = pd.read_csv(list)
        df_in['Analyzed Period Start'] = pd.to_datetime(df_in['Analyzed Period Start'])
        df_in['Analyzed Period End'] = pd.to_datetime(df_in['Analyzed Period End'])
        
        #Add relative paths
        path_cols = ["All Fluxes Time Series", "JSON", "Flux Time Series"]
        cols = df_in.columns.to_list()
        for col1 in cols:
            for col2 in path_cols:
                if col2 in col1:
                    for index, row in df_in.iterrows():
                        df_in.loc[index,col1] = os.path.join(fpath, row[col1])
        
        df = pd.concat([df, df_in], ignore_index=True)
        
    #Sort by time
    df = df.sort_values(by='Analyzed Period Start')
    df.reset_index(drop=True, inplace=True)
    
    df_primary = pd.DataFrame()
    
    for index, row in df.iterrows():
        date = row['Analyzed Period Start']
        date_end = row['Analyzed Period End']
        sc = row['Experiment']
        
        goes_primary = goes.goes_primary_lookup(date)
        goes_primary_end = goes.goes_primary_lookup(date_end)
        

        #GOES integral fluxes read by FetchSEP are labeled GOES-RT because
        #NOAA does not yet provide and archive of those files. So any
        #GOES-R+ integral fluxes will come from CCMC iSWA and are labeled GOES-RT.
        
        if sc != goes_primary and sc != "GOES-RT":
            continue
        elif sc == "GOES-RT" and goes_primary in goes_R:
            df_primary = pd.concat([df_primary, df.iloc[[index]]], ignore_index=True)
        elif sc == goes_primary and goes_primary == "GOES-RT":
            df_primary = pd.concat([df_primary, df.iloc[[index]]], ignore_index=True)
        elif sc == goes_primary:
            df_primary = pd.concat([df_primary, df.iloc[[index]]], ignore_index=True)
        elif sc == goes_primary_end:
            df_primary = pd.concat([df_primary, df.iloc[[index]]], ignore_index=True)
            
    start_date = df_primary['Analyzed Period Start'].iloc[0]
    end_date = df_primary['Analyzed Period End'].iloc[-1]
    stdate = start_date.strftime("%Y-%m-%d")
    enddate = end_date.strftime("%Y-%m-%d")
    
    fname = f"{prefix}_PRIMARY.{stdate}.{enddate}_sep_events.csv"
    fname = os.path.join(cfg.outpath,fname)
    df_primary.to_csv(fname, index=False)
    print(f"create_primary_goes_sep_list: Wrote file {fname}.")
    
    return df_primary
    
