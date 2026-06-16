from . import config as cfg
from . import experiments as expts
import pandas as pd
import datetime
import os
import sys

#Information about GOES-Primary Spacecraft
def import_goes_status():
    """ Read in GOES_Primary_Secondary_Status.csv provided by
        Kim Moreland (NOAA SWPC) on September 5, 2024.
        The dates in the file go up to September 6, 2024.
        
        Provides a date ranges and specifies the primary and secondary
        GOES satellites for protons and X-rays during these times.

        INPUTS:
        
            None
            
        OUTPUTS:
        
            df: (pandas DataFrame) dataframe containing dates and goes status

    """

    filename = 'GOES_Primary_Secondary_Status.csv'
    filename = os.path.join('fetchsep', 'reference', filename)
    if not os.path.isfile(filename):
        sys.exit(f"The GOES status file, {filename}, is not present in the "
            "fetchsep/utils/ directory. Please download from the fetchsep repository "
            "or contact kathryn.whitman@nasa.gov for assistance. "
            "https://github.com/ktindiana/fetchsep")
    
    df = pd.read_csv(filename, parse_dates=True)
    df['Start Date'] = pd.to_datetime(df['Start Date'])
    df['End Date'] = pd.to_datetime(df['End Date'])

    return df
    
    
def id_goes_spacecraft(df, status, instrument, date):
    """ Identify the primary or secondary GOES satellite for protons
        or X-rays in a certain date range.
        
        INPUTS:
        
            :status: (string) "Primary" or "Secondary"
            :instrument: (string) "Protons" or "X-rays"
            :date: (datetime) date of interest
            
        OUTPUTS:
        
            :experiment: (string) name of GOES spacecraft that satisfies
                the query
        
    """
    
    sub = df.loc[(df['Status']==status) & (df['Instrument']==instrument) & (df['Start Date'] <= date) & (df['End Date'] > date)]

    if sub.empty:
        print(f"Cannot identify the {status} {instrument} GOES satellite for {date}.")
        return None

    goes = sub['Satellite'].iloc[0]
    
    return goes
    
    
    
def goes_primary_lookup(date):
    """ Determine which was the primary GOES experiment for a given timeframe.
     
        Lookup table provided by Kim Moreland (NOAA SWPC) on Sep 6, 2024.
        
        Note that all GOES real time integral fluxes retrieved from the CCMC hapi server
        are, by definition, primary GOES fluxes extracted in real time from
        SWPC's 3-day jsons. There is no need to check for primary spacecraft if
        using that data source.
        
    """

    #Dates for which historical GOES primary/secondary are known
    table_date = datetime.datetime(2024,9,6,23,23,59)

    #Date beyond table
    if date > table_date:
        print("Do not have information about primary/secondary GOES for this date. "
                "Returning GOES-RT, but this is only valid for integral fluxes. "
                "You may check "
                "https://services.swpc.noaa.gov/json/goes/instrument-sources.json "
                "if you are looking for current GOES data. If GOES integral fluxes are "
                "desired, specify GOES-RT and primary satellite fluxes will be used "
                "by default, as CCMC downloads the SWPC 3-day jsons, which are by "
                "definition from the primary satellite. This function will be updated "
                "when a reliable, living archive of GOES primary/secondary satellites "
                "becomes available.")
        
        goes = "GOES-RT"


    if date <= datetime.datetime(2024,9,6):
        df = import_goes_status()
        
        goes = id_goes_spacecraft(df, "Primary","Protons",date)
        
    return goes


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
        
        df = pd.concat([df, df_in], ignore_index=True)
        
    #Sort by time
    df = df.sort_values(by='Analyzed Period Start')
    df.reset_index(drop=True, inplace=True)
    
    df_primary = pd.DataFrame()
    
    for index, row in df.iterrows():
        date = row['Analyzed Period Start']
        date_end = row['Analyzed Period End']
        sc = row['Experiment']
        
        goes_primary = goes_primary_lookup(date)
        goes_primary_end = goes_primary_lookup(date_end)
        

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
    
