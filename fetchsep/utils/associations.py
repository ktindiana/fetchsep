from . import config as cfg
from ..json import ccmc_json_handler as ccmc_json
import pandas as pd
import os
import sys
import string
import datetime
import wget
import urllib.request
import requests
import numpy as np
import netCDF4 as nc
import cftime
import inspect

#GOES spacecraft categories needed for GOES X-ray data
#Spacecraft in the GOES-R+ series
goes_R = ["GOES-16", "GOES-17", "GOES-18"]
#Spacecraft prior to GOES-R
goes_sc = ["GOES-08", "GOES-09","GOES-10","GOES-11",
            "GOES-12","GOES-13","GOES-14","GOES-15"]
#Spacecraft prior to GOES-08
old_goes_sc = ["GOES-05", "GOES-06", "GOES-07"]


def clean_file(input_path, output_path):
    """ Clean hidden ascii characters from files.
    
    """
    # Define allowed characters (printable + common whitespace)
    allowed = set(string.printable)

    with open(input_path, "r", errors="ignore") as infile, open(output_path, "w") as outfile:
        for line in infile:
            # Keep only allowed characters
            cleaned = "".join(ch for ch in line if ch in allowed)
            outfile.write(cleaned)


def ccmc_cme_block():
    cme = { 'start_time': pd.NaT,
            'liftoff_time': pd.NaT,
            'lat': np.nan,
            'lon': np.nan,
            'pa': np.nan,
            'half_width': np.nan,
            'speed': np.nan,
            'acceleration': np.nan,
            'height': np.nan,
            'time_at_height': {'time': pd.NaT, 'height': np.nan},
            'coordinates': '',
            'catalog': '',
            'catalog_id': '',
            'urls': [],
            'derivation_technique': {'process': '','method': '', 'measurement_type': ''}
        }
    return cme

def ccmc_flare_block():
    flare = {'last_data_time': pd.NaT,
             'start_time': pd.NaT,
             'peak_time': pd.NaT,
             'end_time': pd.NaT,
             'location': '',
             'intensity': np.nan,
             'integrated_intensity': np.nan,
             'noaa_region': np.nan,
             'peak_ratio': np.nan,
             'catalog': '',
             'catalog_id': '',
             'urls': []
        }
    return flare

def string_location(lat, lon):
    """Put lat, lon in string location format """
    
    lat = int(lat)
    lon = int(lon)
    
    ns = ''
    if lat >= 0:
        ns='N'
    else:
        lat = abs(lat)
        ns = 'S'
        
    ew = ''
    if lon >= 0:
        ew = 'W'
    else:
        lon = abs(lon)
        ew = 'E'
        
    location = f"{ns}{lat:02d}{ew}{lon:02d}"
    return location


################### STEVE JOHNSON'S SRAG SEP LIST #########################
def srag_time_columns():
    columns = ['InitiationDT', 'SEP Reference Time', 'Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio DH Start Time', 'Radio DH End Time', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'CDAW CME First Look Time', 'DONKI CME Start Time', 'DONKI CME Time at 21.5', 'P10_FluxStart', 'P10_StartDT', 'P10_OnsetMax_DT', 'P10_PeakDT', 'P10_EndDT', 'P50_FluxStart', 'P50_StartDT', 'P50_OnsetPeakDT', 'P50_PeakDT', 'P50_EndDT', 'P100_FluxStart', 'P100_StartDT', 'Onset_P100_PkDT', 'Event_P100_PkDT', 'P100_EndDT','ESP Flare Xray Start Time Deprecated', 'ESP Flare Xray Peak Time Deprecated', 'ESP Flare Xray End Time Deprecated', 'ESP Radio m_TyII Start Time', 'ESP Radio m_TyII End Time', 'ESP CDAW CME First Look Time', 'ACE CME passage', 'Sudden Imp Time']
    return columns


def srag_float_columns():
    columns = ['Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Duration Deprecated', 'Flare Xray Time To Peak Deprecated','Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Xray Time To Peak', 'AR Area', 'AR Carrington', 'Event Location From Center', 'Event Latitude', 'Event Longitude', 'Event Location from Center 2', 'Event Latitude 2', 'Event Longitude 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME Speed', 'CDAW CME Width', 'CDAW CME Mean Position Angle', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon']
    return columns


def srag_int_columns():
    columns = ['Flare Catalog ID', 'Cycle', 'Active Region', 'GOES Xray Satellite', 'GOES Protons Satellite']
    return columns


def srag_string_columns():
    columns =  ['EventType', 'Case', 'Flare Class Deprecated', 'Flare Class', 'Flare Opt', 'AR Spot Class', 'AR Mag Class', 'Event Location Source', 'Event Location Source 2', 'Radio Station', 'Radio DH Note', 'DONKI CME Feature', 'DONKI CME Catalog ID', 'ESP_CME', 'GLE Event Number', 'PRF','Comments']
    return columns


def srag_clear_columns():
    columns = ['Cycle', 'EventType', 'Case', 'Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Class Deprecated', 'Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Flare Class', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Xray Time To Peak', 'Flare Catalog ID', 'GOES Xray Satellite', 'Flare Opt', 'Active Region', 'AR Area', 'AR Spot Class', 'AR Mag Class', 'AR Carrington', 'Event Location From Center', 'Event Latitude', 'Event Longitude', 'Event Location Source', 'Event Location from Center 2', 'Event Latitude 2', 'Event Longitude 2', 'Event Location Source 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio Station', 'Radio DH Start Time', 'Radio DH End Time', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio DH Note', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME First Look Time', 'CDAW CME Speed', 'CDAW CME Width', 'CDAW CME Mean Position Angle', 'DONKI CME Start Time', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', 'DONKI CME Time at 21.5', 'DONKI CME Catalog ID', 'ESP_CME', 'GLE Event Number', 'PRF', 'Comments']
    return columns


def is_time_column(col):
    columns = srag_time_columns()

    if col in columns:
        return True
    else:
        return False


def is_float_column(col):
    columns = srag_float_columns()

    if col in columns:
        return True
    else:
        return False


def is_int_column(col):
    columns = srag_int_columns()

    if col in columns:
        return True
    else:
        return False


def is_str_column(col):
    columns =  srag_string_columns()

    if col in columns:
        return True
    else:
        return False
    

def cast_srag_list_time_columns(df):
    columns = srag_time_columns()

    for col in columns:
        df[col] = pd.to_datetime(df[col])
        
    return df
    
    

def empty_srag_associations_series():
    """ Create a pandas Series with appropriate null values 
        for the associations from the SRAG list (used in 
        the CLEAR benchmark dataset).
    
    """

    clear_col = srag_clear_columns()
    
    dict = {}
    for col in clear_col:
        val = None
        if is_int_column(col):
            val = np.nan
        elif is_float_column(col):
            val = np.nan
        elif is_str_column(col):
            val = ''
        elif is_time_column(col):
            val = pd.NaT

        dict.update({col: val})

    return dict


def empty_srag_protons_series():
    dict= {"P10_FluxStart": pd.NaT, "P10_StartDT": pd.NaT}
    return dict
    

def combine_srag_list_comments(df):
    df["Comments1"] = df["Comments1"].where(pd.notna(df["Comments1"]), '')
    df["Comments2"] = df["Comments2"].where(pd.notna(df["Comments2"]), '')
    df["Comments"] = df["Comments1"] + df["Comments2"]
    df.drop("Comments1", axis=1, inplace=True)
    df.drop("Comments2", axis=1, inplace=True)
    return df


def read_srag_list():
    """ Read in the SRAG list provided by Steve Johnson (SRAG) 
        converted from excel into text format and cleaned.
        
    """
    srag_list = 'fetchsep/reference/SRAG_SEP_List_R11_CLEARversion.csv'
    df = pd.read_csv(srag_list)
    
    #Combine two separate comments columns into one
    #df = combine_srag_list_comments(df)
    
    #Cast time columns as datetime
    df = cast_srag_list_time_columns(df)
    
    return df


def identify_associations_in_srag_list(startdate):
    """ Given a date range, identify an event in the list associated
        with that date range. 
    
        Steve Johnson's SRAG list provides the InitiationDT which is
        related to the flare, CME, or first rise of protons, as indicated
        in InitiationType.
        
        P10_FluxStart indicates when >10 Mev goes above bakground by Steve's
        estimation. 
        
        P10_StartDT indicates when >10 MeV crosses 10 pfu.
        
        INPUTS:
        
            :df: (dataframe) Steve's SRAG list
            :start_date: (datetime) start time of SEP event that want to
                find in Steve's list; try time of enhancement above background
                for >10 MeV or proton energies ~10 MeV or threshold crossing time
                
        OUTPUT:
        
            :associations: (pandas Series) Series containing the associations
                fields used for the CLEAR benchmark dataset
        
    """
    df = read_srag_list()
    
    tolerance = datetime.timedelta(hours=6)
    columns = ["SEP Reference Time"]
    
    diff = df[columns] - startdate
    
    event_index = -1 #min index for each column
    check_min_diff = np.nan
    for col in columns:
        ix = []
        diff[col] = abs(diff[col])
        min_diff = diff[col].min()
        if min_diff <= tolerance:
            ix = diff[col].index[diff[col] == min_diff]
            
        if len(ix) == 0:
            continue
        else:
            if event_index == -1:
                event_index = ix[0]
                check_min_diff = min_diff
            else:
                if event_index != ix[0]:
                    print("identify_event_in_srag_list: Found more than one possible "
                        f"match for {startdate}.")
                    if min_diff < check_min_diff:
                        event_index = ix[0]
                        check_min_diff = min_diff

    null_assoc = empty_srag_associations_series()
    null_protons = empty_srag_protons_series()
    if event_index == -1:
        print("identify_event_in_srag_list: No match was found in the SRAG list for "
            f"{startdate}.")
        #Return series with columns with appropriate null values
        return null_assoc, null_protons
    else:
        clear_col = srag_clear_columns()
        associations = df[clear_col].iloc[event_index]
        proton_info = df[columns].iloc[event_index]
        print("identify_event_in_srag_list: Found association information for input date of " f"{startdate} with an entry in the SRAG list for >10 MeV {proton_info}")
        return associations.to_dict(), proton_info.to_dict()

    return null_assoc, null_protons


def srag_to_ccmc_cme(associations):
    """ Take dictionary containing one row of the SRAG list. 
        Put the CME columns into the appropriate fields for
        the CCMC SEP Scoreboard json format.
        
        INPUT:
        
            :associations: (dict) a dictionary containing a single row of
                the SRAG list with CME information included, e.g. identified
                with identify_associations_in_srag_list
                
        OUTPUT:
            
            :cme: (dict) CME trigger block in CCMC json schema
        
    """
    
    all_cmes = []

    #DONKI CME
    if not pd.isnull(associations['DONKI CME Speed']):
        cme = ccmc_cme_block()

        cme['speed'] = associations['DONKI CME Speed']
        cme['catalog'] = 'DONKI'
        cme['derivation_technique']['process'] = 'manual'
        cme['derivation_technique']['method'] = 'SWPC_CAT'

        if not pd.isnull(associations['DONKI CME Start Time']):
            cme['start_time'] = ccmc_json.make_ccmc_zulu_time(associations['DONKI CME Start Time'])
        else:
            del cme['start_time']

        if not pd.isnull(associations['DONKI CME Catalog ID']):
            cme['catalog_id'] = associations['DONKI CME Catalog ID']
        else:
            del cme['catalog_id']
        
        if not pd.isnull(associations['DONKI CME Lat']):
            cme['lat'] = associations['DONKI CME Lat']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lat']

        if not pd.isnull(associations['DONKI CME Lon']):
            cme['lon'] = associations['DONKI CME Lon']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lon']
            del cme['coordinates']

        if not pd.isnull(associations['DONKI CME Half Width']):
            cme['half_width'] = associations['DONKI CME Half Width']
        else:
            del cme['half_width']

        if not pd.isnull(associations['DONKI CME Time at 21.5']):
            cme['time_at_height']['time'] = ccmc_json.make_ccmc_zulu_time(associations['DONKI CME Time at 21.5'])
            cme['time_at_height']['height'] = 21.5
        else:
            del cme['time_at_height']

        #Remove unused columns
        del cme['liftoff_time']
        del cme['acceleration']
        del cme['height']
        del cme['urls']
        del cme['pa']
        try:
            del cme['derivation_technique']['measurement_type']
        except:
            pass

        all_cmes.append(cme)

    
    #CDAW CME
    if not pd.isnull(associations['CDAW CME Speed']):
        cme = ccmc_cme_block()

        cme['speed'] = associations['CDAW CME Speed']
        cme['catalog'] = 'SOHO_CDAW'
        cme['derivation_technique']['process'] = 'manual'
        cme['derivation_technique']['method'] = 'Plane-of-sky'

        if not pd.isnull(associations['CDAW CME First Look Time']):
            cme['start_time'] = associations['CDAW CME First Look Time'].to_pydatetime()
        else:
            del cme['start_time']

        if not pd.isnull(associations['Event Latitude']):
            cme['lat'] = associations['Event Latitude']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lat']

        if not pd.isnull(associations['Event Longitude']):
            cme['lon'] = associations['Event Longitude']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lon']
            del cme['coordinates']

        if not pd.isnull(associations['CDAW CME Mean Position Angle']):
            cme['pa'] = associations['CDAW CME Mean Position Angle']
        else:
            del cme['pa']
    
        width = associations['CDAW CME Width']
        if pd.isnull(width):
            del cme['half_width']
        else:
            if width == '':
                width = np.nan
            elif 'Halo' in width:
                width = 360.
            elif '>' in width:
                width = width.strip().split('>')
                width = float(width[1])
            else:
                width = float(width)
            cme['half_width'] = width/2.

        #Remove unused columns
        del cme['liftoff_time']
        del cme['acceleration']
        del cme['height']
        del cme['time_at_height']
        del cme['coordinates']
        del cme['catalog_id']
        del cme['urls']
        try:
            del cme['derivation_technique']['measurement_type']
        except:
            pass

        all_cmes.append(cme)

    return all_cmes
    
    
def srag_to_ccmc_flare(associations):
    """ Take dictionary containing one row of the SRAG list. 
        Put the Flare columns into the appropriate fields for
        the CCMC SEP Scoreboard json format.
        
        INPUT:
        
            :associations: (dict) a dictionary containing a single row of
                the SRAG list with flare information included, e.g. identified
                with identify_associations_in_srag_list
                
        OUTPUT:
            
            :flare: (dict) Flare trigger block in CCMC json schema
        
    """
    all_flares = []
    
    #SWPC Xray Science Data
    if not pd.isnull(associations['Flare Xray Start Time']):
        flare = ccmc_flare_block()
        flare['start_time'] = ccmc_json.make_ccmc_zulu_time(associations['Flare Xray Start Time'])

        if not pd.isnull(associations['Flare Xray Peak Time']):
            flare['peak_time'] = ccmc_json.make_ccmc_zulu_time(associations['Flare Xray Peak Time'])
        else:
            del flare['peak_time']
            
        if not pd.isnull(associations['Flare Xray End Time']):
            flare['end_time'] = ccmc_json.make_ccmc_zulu_time(associations['Flare Xray End Time'])
        else:
            del flare['end_time']
        
        if not pd.isnull(associations['Event Latitude']) and not pd.isnull(associations['Event Longitude']):
            lat = associations['Event Latitude']
            lon = associations['Event Longitude']
            flare['location'] = string_location(lat,lon)
        else:
            del flare['location']

        if not pd.isnull(associations['Flare Magnitude']):
            flare['intensity'] = associations['Flare Magnitude']
            flare['catalog'] = 'SWPC'
        else:
            del flare['intensity']
            del flare['catalog']
        
        if not pd.isnull(associations['Flare Integrated Flux']):
            flare['integrated_intensity'] = associations['Flare Integrated Flux']
        else:
            del flare['integrated_intensity']

        if not pd.isnull(associations['Active Region']):
            flare['noaa_region'] = associations['Active Region']
        else:
            del flare['noaa_region']

        if not pd.isnull(associations['Flare Catalog ID']):
            flare['catalog_id'] = int(associations['Flare Catalog ID'])
        else:
            del flare['catalog_id']

        #Remove unused entries
        del flare['last_data_time']
        del flare['peak_ratio']
        del flare['urls']
 
        all_flares.append(flare)
 
    #SWPC X-ray operational data that is now deprecated after GOES-R launched
    if not pd.isnull(associations['Flare Xray Start Time Deprecated']):
        flare = ccmc_flare_block()
        flare['start_time'] = ccmc_json.make_ccmc_zulu_time(associations['Flare Xray Start Time Deprecated'])

        if not pd.isnull(associations['Flare Xray Peak Time Deprecated']):
            flare['peak_time'] = associations['Flare Xray Peak Time Deprecated'].to_pydatetime()
        else:
            del flare['peak_time']
            
        if not pd.isnull(associations['Flare Xray End Time Deprecated']):
            flare['end_time'] = ccmc_json.make_ccmc_zulu_time(associations['Flare Xray End Time Deprecated'])
        else:
            del flare['end_time']
        
        if not pd.isnull(associations['Event Latitude']) and not pd.isnull(associations['Event Longitude']):
            lat = associations['Event Latitude']
            lon = associations['Event Longitude']
            flare['location'] = string_location(lat,lon)
        else:
            del flare['location']

        if not pd.isnull(associations['Flare Magnitude Deprecated']):
            flare['intensity'] = associations['Flare Magnitude Deprecated']
            flare['catalog'] = 'SWPC_DEPRECATED'
        else:
            del flare['intensity']
            del flare['catalog']
        
        if not pd.isnull(associations['Flare Integrated Flux Deprecated']):
            flare['integrated_intensity'] = associations['Flare Integrated Flux Deprecated']
        else:
            del flare['integrated_intensity']

        if not pd.isnull(associations['Active Region']):
            flare['noaa_region'] = associations['Active Region']
        else:
            del flare['noaa_region']
        
        #Remove unused entries
        del flare['last_data_time']
        del flare['peak_ratio']
        del flare['catalog_id']
        del flare['urls']
 
        all_flares.append(flare)
 
    return all_flares


##########################################################
################ NOAA X-RAY SCIENCE DATA #################
##########################################################
def goes_xray_sat_info():
    """ Info to create filenames for NOAA's naming scheme for X-ray files.
    
    """

    #Selecting the science-quality X-ray fluxes updated by NOAA in 2021+
    #startdate and enddate are inclusive here
    sat_info = {"GOES-08": {"url": 'goes08',
                        "avg1m": 'sci_xrsf-l2-avg1m_g08_',
                        "flsum": 'sci_xrsf-l2-flsum_g08_',
                        "avg_version": '_v1-0-0',
                        "flsum_version": '_v1-0-0',
                        "startdate": datetime.datetime(1995,1,1),
                        "enddate": datetime.datetime(2003,6,16)},
              
              "GOES-09": {"url": 'goes09',
                        "avg1m": 'sci_xrsf-l2-avg1m_g09_',
                        "flsum": 'sci_xrsf-l2-flsum_g09_',
                        "avg_version": '_v1-0-0',
                        "flsum_version": '_v1-0-0',
                        "startdate": datetime.datetime(1996,4,1),
                        "enddate": datetime.datetime(1998,7,28)},
              
              "GOES-10":  {"url": 'goes10',
                        "avg1m": 'sci_xrsf-l2-avg1m_g10_',
                        "flsum": 'sci_xrsf-l2-flsum_g10_',
                        "avg_version": '_v1-0-0',
                        "flsum_version": '_v1-0-0',
                        "startdate": datetime.datetime(1998,7,1),
                        "enddate": datetime.datetime(2009,12,1)},
              
              "GOES-11":  {"url": 'goes11',
                        "avg1m": 'sci_xrsf-l2-avg1m_g11_',
                        "flsum": 'sci_xrsf-l2-flsum_g11_',
                        "avg_version": '_v1-0-0',
                        "flsum_version": '_v1-0-0',
                        "startdate": datetime.datetime(2006,6,1),
                        "enddate": datetime.datetime(2008,2,10)},
              
              "GOES-12":  {"url": 'goes12',
                        "avg1m": 'sci_xrsf-l2-avg1m_g12_',
                        "flsum": 'sci_xrsf-l2-flsum_g12_',
                        "avg_version": '_v1-0-0',
                        "flsum_version": '_v1-0-0',
                        "startdate": datetime.datetime(2003,1,10),
                        "enddate": datetime.datetime(2007,4,12)},

              "GOES-13":  {"url": 'goes13',
                        "avg1m": 'sci_xrsf-l2-avg1m_g13_',
                        "flsum": 'sci_xrsf-l2-flsum_g13_',
                        "avg_version": '_v2-2-1',
                        "flsum_version": '_v2-3-0',
                        "startdate": datetime.datetime(2013,6,7),
                        "enddate": datetime.datetime(2017,12,14)},
               
              "GOES-14":  {"url": 'goes14',
                        "avg1m": 'sci_xrsf-l2-avg1m_g14_',
                        "flsum": 'sci_xrsf-l2-flsum_g14_',
                        "avg_version": '_v2-2-1',
                        "flsum_version": '_v2-3-0',
                        "startdate": datetime.datetime(2009,9,19),
                        "enddate": datetime.datetime(2020,3,4)},
               
              "GOES-15":  {"url": 'goes15',
                        "avg1m": 'sci_xrsf-l2-avg1m_g15_',
                        "flsum": 'sci_xrsf-l2-flsum_g15_',
                        "avg_version": '_v2-2-1',
                        "flsum_version": '_v2-3-0',
                        "startdate": datetime.datetime(2010,4,7),
                        "enddate": datetime.datetime(2020,3,4)},

               "GOES-16":   {"url": 'goes16',
                        "avg1m": 'sci_xrsf-l2-avg1m_g16_',
                        "flsum": 'sci_xrsf-l2-flsum_g16_',
                        "avg_version": '_v2-2-0',
                        "flsum_version": '_v2-2-0',
                        "startdate": datetime.datetime(2017,2,7),
                        "enddate": datetime.datetime(2025,4,6)},

               "GOES-17":   {"url": 'goes17',
                        "avg1m": 'sci_xrsf-l2-avg1m_g17_',
                        "flsum": 'sci_xrsf-l2-flsum_g17_',
                        "avg_version": '_v2-2-0',
                        "flsum_version": '_v2-2-0',
                        "startdate": datetime.datetime(2018,6,1),
                        "enddate": datetime.datetime(2024,8,11)},

               "GOES-18":   {"url": 'goes18',
                        "avg1m": 'sci_xrsf-l2-avg1m_g18_',
                        "flsum": 'sci_xrsf-l2-flsum_g18_',
                        "avg_version": '_v2-2-0',
                        "flsum_version": '_v2-2-0',
                        "startdate": datetime.datetime(2022,9,2),
                        "enddate": datetime.datetime.now() - datetime.timedelta(hours=24)} #update when GOES-18 ends
              }

    return sat_info


def goes_xray_sat_covers_date(experiment, request_date):
    """ Check if the spacecraft provided data during the requested date. """
    sat_info = goes_xray_sat_info()

    if experiment not in sat_info.keys():
        print(f"NOAA X-ray science data not available for {experiment}.")
        return False
    
    if request_date >= sat_info[experiment]["startdate"] \
        and request_date <= sat_info[experiment]["enddate"]:
        return True
    else:
        return False


def check_goes_xray_science_data(startdate, enddate, experiment):
    """ Check that GOES X-ray data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
        Reprocessed science X-ray data is available in daily files. 
        
        The following files are downloaded here:
            - avg1m: 1-minute averages of XRS irradiances
            - flsum: flare summary, flare detection flags such as start and peak
        
        The X-ray time profile can be plotted, analyzed from avg1m.
        The flare start, peak, end times and values and integrated flux are
        save in the flsum files.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of GOES satellite
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            XRS data that span the desired time range
            (monthly files)
        
    """
    avgfiles = [] #X-ray time series
    flsumfiles = [] #flare summary with magnitude, start, peak, end, etc

    if experiment in old_goes_sc:
        print(f"check_goes_xray_science_data: {experiment} X-ray science data is not yet available from NOAA. Returning.")
        return avgfiles, flsumfiles

        
    #GOES XRS science data is stored in daily data files
    td = enddate - startdate
    NFILES = td.days #number of data files to download
    if td.seconds > 0: NFILES = NFILES + 1

    sat_info = goes_xray_sat_info()

    #GOES-08 to -15 files like:
    #https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes08/xrsf-l2-avg1m_science/1995/01/
    #https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes09/xrsf-l2-flsum_science/1996/04/

    #GOES-R like
    #https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-avg1m_science/2017/02/
    #https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-flsum_science/2017/02/


    #for every day that data is required, check if file is present or
    #needs to be downloaded.
    for i in range(NFILES):
        date = startdate + datetime.timedelta(days=i)
        year = date.year
        month = date.month
        day = date.day
        date_suffix = 'd%i%02i%02i' % (year,month,day)

        if date < sat_info[experiment]['startdate'] or date > sat_info[experiment]['enddate']:
            print(f"{experiment} X-ray science data is not available for {date}. Skipping.")
            continue

        #Filenames like:
        #sci_xrsf-l2-avg1m_g08_d19950103_v1-0-0.nc
        #sci_xrsf-l2-flsum_g08_d19950103_v1-0-0.nc
        avgfile = f"{sat_info[experiment]['avg1m']}{date_suffix}{sat_info[experiment]['avg_version']}.nc"
        flsumfile = f"{sat_info[experiment]['flsum']}{date_suffix}{sat_info[experiment]['flsum_version']}.nc"

        #Check if files already on your computer
        fullpath = os.path.join(cfg.datapath,'GOES',avgfile)
        avg_exists = os.path.isfile(fullpath)
        if avg_exists:
            avgfiles.append(os.path.join('GOES',avgfile))
        
        fullpath = os.path.join(cfg.datapath,'GOES',flsumfile)
        flsum_exists = os.path.isfile(fullpath)
        if flsum_exists:
            flsumfiles.append(os.path.join('GOES',flsumfile))

        #If already have both files
        if avg_exists and flsum_exists:
            continue
        

        #If need to download
        ###### GOES-08 to GOES-15 #####
        if experiment in goes_sc:
            avgurl = ('https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/%s/xrsf-l2-avg1m_science/%i/%02i/%s') % (sat_info[experiment]["url"],year,month,avgfile)
             
            flsumurl = ('https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/%s/xrsf-l2-flsum_science/%i/%02i/%s') % (sat_info[experiment]["url"],year,month,flsumfile)
 

        ###### GOES-R ################
        if experiment in goes_R:
            avgurl = ('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/%s/l2/data/xrsf-l2-avg1m_science/%i/%02i/%s') % (sat_info[experiment]["url"],year,month,avgfile)
             
            flsumurl = ('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/%s/l2/data/xrsf-l2-flsum_science/%i/%02i/%s') % (sat_info[experiment]["url"],year,month,flsumfile)


        #Download
        if not avg_exists:
            try:
                try:
                    urllib.request.urlopen(avgurl)
                except:
                    raise

                wget.download(avgurl, os.path.join(cfg.datapath,'GOES',avgfile))
                print(f"\ncheck_goes_science_xray_data: Downloaded {avgurl}")
                avgfiles.append(os.path.join('GOES',avgfile))
            except urllib.request.HTTPError:
                print("Cannot access GOES file at " + avgurl + ". Skipping.")



        if not flsum_exists:
            try:
                try:
                    urllib.request.urlopen(flsumurl)
                except:
                    raise

                wget.download(flsumurl, os.path.join(cfg.datapath,'GOES',flsumfile))
                print(f"\ncheck_goes_science_xray_data: Downloaded {flsumurl}")
                flsumfiles.append(os.path.join('GOES',flsumfile))
            except urllib.request.HTTPError:
                print("Cannot access GOES file at " + flsumurl + ". Skipping.")


    return avgfiles, flsumfiles


def remove_swpc_calibration(flux):
    """ For GOES-07 and previous satellites, reprocessed NOAA Xray science data has not 
        been released. NOAA provides this guidance:
            
            To get true fluxes for GOES-1 through -15 operational data, 
            users must remove the SWPC scaling factors from the data. 
            To do this, divide the short band flux by 0.85 and divide 
            the long band flux by 0.7.
            
        The long band flux is what is used to determine flare class.
        For processing here, flare fluxes already compiled in the SRAG
        list are divided by 0.7
        
        INPUTS:
        
            :flux: (float) intensity level - flux or event-integrated fluence
            
        OUTPUTS:
        
            :calibrated_flux: (float) SWPC scaling factor removed
        
    """

    return flux/0.7

def flare_class(peak_intensity):
    """ Calculate flare class from peak intensity. """

    fl_class = ''
    if pd.isnull(peak_intensity):
        return fl_class

    classes = ["A", "B", "C", "M", "X"]
    intensities = {"X": {'min': 1e-4, 'max': 1},
                   "M": {'min':1e-5, 'max': 1e-4},
                   "C": {'min':1e-6, 'max': 1e-5},
                   "B": {'min':1e-7, 'max': 1e-6},
                   "A": {'min':1e-8, 'max': 1e-7},
                   }

    for cl in classes:
        if peak_intensity >= intensities[cl]['min'] and peak_intensity < intensities[cl]['max']:
            number = peak_intensity/intensities[cl]['min']
            number = str(number)
            fl_class = f"{cl}{number[0:3]}"
            
    return fl_class


def manual_flare_calibration(df, index):
    """ Convert deprecated flare xray columns into fluxes that align better
        with science-quality Xray data by removing SWPC calibration factors.
    
    """

    flare_info = {  'catalog_id': None, #flare_id
                    'catalog': 'SWPC',
                    'start_time': pd.NaT,
                    'peak_time': pd.NaT,
                    'end_time': pd.NaT,
                    'duration': np.nan,
                    'time_to_peak': np.nan,
                    'intensity': np.nan,
                    'class': None,
                    'instrument': None,
                    'integrated_intensity': np.nan, #start time to end time
                }

    if pd.isnull(df['Flare Magnitude Deprecated'].loc[index]):
        return flare_info

    goes = df['GOES Xray Satellite'].loc[index]
    if not pd.isnull(goes):
        goes = int(goes)
        experiment = f"GOES-{goes:02d}"
    else:
        experiment = None

    flare_info['instrument'] = experiment
    flare_info['start_time'] = df['Flare Xray Start Time Deprecated'].loc[index]
    flare_info['peak_time'] = df['Flare Xray Peak Time Deprecated'].loc[index]
    flare_info['end_time'] = df['Flare Xray End Time Deprecated'].loc[index]
    flare_info['duration'] = df['Flare Duration Deprecated'].loc[index]
    flare_info['time_to_peak'] = df['Flare Xray Time To Peak Deprecated'].loc[index]
    
    fluence = df['Flare Integrated Flux Deprecated'].loc[index]
    fluence = remove_swpc_calibration(fluence)
    flare_info['integrated_intensity'] = fluence
    
    peak_intensity = df['Flare Magnitude Deprecated'].loc[index]
    peak_intensity = remove_swpc_calibration(peak_intensity)
    fl_class = flare_class(peak_intensity)
    flare_info['intensity'] = peak_intensity
    flare_info['class'] = fl_class

    return flare_info


def get_unique(arr_in):
    """ Given a list arr, return a list of unique values within arr. """
    
    try:
        arr = arr_in.tolist()
    except:
        arr = arr_in
        
    unique = []
    for val in arr:
        if val in unique:
            continue
        else:
            unique.append(val)
            
    return unique



def extract_flare_info(experiment, request_date, flare_info={}, req_flare_id=None):
    """ Extract information for a specific flare on a specific datetime.
        Enter a date within 15 minutes of flare start or end of for any time
        between the flare start and end.
        
        request_date is taken to be the flare peak time as it is the least subjective.
        
        INPUTS:
        
            :experiment: (string) GOES-08, etc
            :date: (datetime) timestamp for a date related to a specific flare,
                can be similar to flare start, peak, end. Will look compare date to
                15 minutes from the flare start and end.
                
        OUTPUTS:
        
            Information available flsum files.
            :info: (dictionary) carries flare start, peak, end, magnitude, fluence, etc.
        
    """
    sat_info = goes_xray_sat_info()

    #Define flare info
    if not flare_info:
        flare_info = {  'catalog_id': None, #flare_id
                        'catalog': 'NOAA_NCEI',
                        'start_time': pd.NaT,
                        'peak_time': pd.NaT,
                        'end_time': pd.NaT,
                        'duration': np.nan,
                        'time_to_peak': np.nan,
                        'intensity': np.nan,
                        'class': None,
                        'instrument': experiment,
                        'integrated_intensity': np.nan, #start time to end time
                    }

    if experiment in old_goes_sc:
        #NEED TO DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES UNTIL NOAA RELEASES SCIENCE DATA
        print(f"extract_flare_info: {experiment} X-ray science data is not yet available from NOAA. Returning.")
        return flare_info


    year = request_date.year
    month = request_date.month
    day = request_date.day
    date_suffix = 'd%i%02i%02i' % (year,month,day)

    #Filenames like:
    #sci_xrsf-l2-flsum_g08_d19950103_v1-0-0.nc
    flsumfile = f"{sat_info[experiment]['flsum']}{date_suffix}{sat_info[experiment]['flsum_version']}.nc"

    fullpath = os.path.join(cfg.datapath,'GOES',flsumfile)
    flsum_exists = os.path.isfile(fullpath)
    if not flsum_exists:
        print(f"extract_flare_info: Cannot find file {fullpath}. Run check_goes_xray_science_data to download.")
        return flare_info

    data = nc.Dataset(fullpath)

    #long_name: Record start time. Neglects leap seconds since 1970-01-01.
    #TIME VARIABLES STORED AS:
    #comments: Time[UTC] =1 Jan 1970 00:00:00[UTC] + time[secs] + n[secs] where
    #n = {0/number of leap secs since 1 Jan 1970} for a conversion function that {ignores/includes} leap secs.
    #comment2: Duplicate timestamps occasionally occur when there are two adjacent flares and
    #the second flare status=EVENT_START at the same time as an EVENT_END or POST_EVENT of the previous flare.
    #units: seconds since 1970-01-01 00:00:00 UTC

    flare_dates = cftime.num2pydate(data.variables["time"][:], data["time"].units)

    #Get unique flare_ids
    flare_ids = get_unique(data.variables["flare_id"][:])
    td = datetime.timedelta(minutes=15) #request_date within 15 minutes of flare
    select_idx = []
    
    if req_flare_id != None:
        print(f"extract_flare_info: User provided flare ID (catalog_id) {req_flare_id}")


    #-->If flare_id is specified
    ############################
    if not pd.isnull(req_flare_id):
        if req_flare_id in flare_ids:
            #Compile all information per flare by extracting according to req_flare_id
            select_idx = [i for i in range(len(data.variables["flare_id"][:])) if data.variables["flare_id"][i]==req_flare_id]
            if len(select_idx) > 0:
                print(f"extract_flare_info: Found flare for flare ID (catalog_id) {req_flare_id}")

    #-->else search according to date of flare peak
    ###############################################
    else:
        peak_td = datetime.timedelta(hours=24)
        for index, flid in enumerate(flare_ids):
            #Compile all information per flare by extracting according to flare_id
            idx = [i for i in range(len(data.variables["flare_id"][:])) if data.variables["flare_id"][i]==flid]
            
            #Check if request_date is associated with this flare
            flst = None
            flpk = None
            flend = None
            for ix in idx:
                if data.variables["status"][ix] == "EVENT_START":
                    flst = flare_dates[ix]
                if data.variables["status"][ix] == "EVENT_PEAK" or data.variables["status"][ix] == "EVENT_PEAK_SATURATED":
                    flpk = flare_dates[ix]
                if data.variables["status"][ix] == "EVENT_END":
                    flend = flare_dates[ix]


            #Cycle through all flares and find the one that minimizes the time
            #difference between the request_date and the found flare peak
            #-->If flare occurred in one day
            ################################
            if flst != None and flend != None:
                if (request_date >= flst-td) and (request_date <= flend+td):
                    ptd = abs(flpk - request_date)
                    if ptd < peak_td:
                        peak_td = ptd
                        select_idx = idx
                        print(f"extract_flare_info: Found {experiment} flare {flst} to {flend} for requested date {request_date}")

            #-->PRECEDENCE. Just in case there are two flares close together, give precedence
            #to the one containing requested date
            #####################################################################
                if (request_date >= flst) and (request_date <= flend):
                    ptd = abs(flpk - request_date)
                    if ptd < peak_td:
                        peak_td = ptd
                        select_idx = idx
                        print(f"extract_flare_info: Found {experiment} flare {flst} to {flend} for requested date {request_date}")


            #-->Check if the last flare in the file is the desired one and goes on to the next day
            ######################################################################################
            if flst != None and flend == None:
                #PRECEDENCE
                if index == len(flare_ids)-1:
                    flend = datetime.datetime(year,month,day,23,59,59)
                    if (request_date >= flst-td) and (request_date <= flend+td):
                        select_idx = idx
                        print(f"extract_flare_info: Found {experiment} flare {flst} to {flend} for requested date {request_date}")
                        break
                elif flpk != None:
                    #look for the requested time within flare start and peak
                    if (request_date >= flst-td) and (request_date <= flpk+td):
                        ptd = abs(flpk - request_date)
                        if ptd < peak_td:
                            peak_td = ptd
                            select_idx = idx
                            print(f"extract_flare_info: Found {experiment} flare {flst} to {flpk+td} for requested date {request_date}")
                    #PRECEDENCE
                    if (request_date >= flst) and (request_date <= flpk+td):
                        ptd = abs(flpk - request_date)
                        if ptd < peak_td:
                            peak_td = ptd
                            select_idx = idx
                            print(f"extract_flare_info: Found {experiment} flare {flst} to {flpk+td} for requested date {request_date}")
                
                

    #If no flare found
    if len(select_idx) == 0:
        print(f"extract_flare_info: {experiment} flare not found for requested date {request_date}. Returning.")
        return flare_info

    #flare_id, status, xrsb_flux, flare_class, time, integrated_flux, flare_class
    for ix in select_idx:
        if data.variables["status"][ix] == "EVENT_START":
            flare_info['start_time'] = flare_dates[ix]
            flare_info['catalog_id'] = int(data.variables["flare_id"][ix])
            
        if data.variables["status"][ix] == "EVENT_PEAK" or data.variables["status"][ix] == "EVENT_PEAK_SATURATED":
            flare_info['peak_time'] = flare_dates[ix]
            flare_info['intensity'] = float(data.variables["xrsb_flux"][ix])
            flare_info['class'] = data.variables["flare_class"][ix]
            
        if data.variables["status"][ix] == "EVENT_END":
            flare_info['end_time'] = flare_dates[ix]
            flare_info['integrated_intensity'] = float(data.variables["integrated_flux"][ix])

    if not pd.isnull(flare_info['start_time']) and not pd.isnull(flare_info['peak_time']):
        ttp = (flare_info['peak_time'] - flare_info['start_time']).total_seconds()/60. #minutes
        flare_info['time_to_peak'] = ttp

    if not pd.isnull(flare_info['start_time']) and not pd.isnull(flare_info['end_time']):
        dur = (flare_info['end_time'] - flare_info['start_time']).total_seconds()/60. #minutes
        flare_info['duration'] = dur

    return flare_info


def update_srag_list_flares(df):
    """ SRAG's SEP event list contains the original X-ray fluxes made available
        by SWPC in their various products and archived data. 
        
        New science-quality X-ray data was released for GOES-08 to GOES-15, which
        should be used instead for space weather understanding, model training,
        and validation.
    
        For each flare in the SRAG list, check if the science data has been added.
        If not, extract the flare values and add to the list.
        
        INPUTS:
        
            :df: (pandas dataframe) SRAG SEP list
            
        OUTPUTS:
        
            :df: (pandas dataframe) updated with new flare information
        
    """
    #In order of preference
    xray_sat = ['GOES-08', 'GOES-09', 'GOES-10', 'GOES-12', 'GOES-15', 'GOES-14', 'GOES-16', 'GOES-18', 'GOES-11', 'GOES-17']

    #Add a column for Flare ID at end of flare info
    idx = df.columns.get_loc('Flare Xray Time To Peak')
    df.insert(loc=idx + 1, column='Flare Catalog ID', value=[None]*len(df['Flare Xray Time To Peak']))


    for index, row in df.iterrows():
        
        #Check if a flare was provided in the SRAG list
        if pd.isnull(row['Flare Magnitude Deprecated']):
            continue

        flare_info = {}

        goes = df['GOES Xray Satellite'].loc[index]
        if not pd.isnull(goes):
            goes = int(goes)
            experiment = f"GOES-{goes:02d}"
        else:
            experiment = None

        goes_r_primary_date = datetime.datetime(2017,2,7) #GOES-R+ primary XRS from this date forward

        #Old GOES-07 and previous
        #########################
        if experiment in old_goes_sc:
            #NEED TO DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES UNTIL NOAA RELEASES SCIENCE DATA
            print(f"extract_flare_info: {experiment} X-ray science data is not yet available from NOAA. Manually removing SWPC calibration factor per NOAA guidance.")
            flare_info = manual_flare_calibration(df, index)
            
        #Ian's List
        ###########
        #Ian Richardson's list has less information and ambiguous times. If entry from his
        #list, check if the provided time was associated with a flare.
        elif 'IGR List' in row['Comments']:
            request_date = row['SEP Reference Time']
            print(f"Flare sourced from Ian's list {request_date}")
            stdate = request_date - datetime.timedelta(hours=24)
            enddate = request_date + datetime.timedelta(hours=48)

            #-->Don't know which satellite was used, so check until find one covered
            #the desired date
            flare_info = {}
            for experiment in xray_sat:
                if not goes_xray_sat_covers_date(experiment, request_date):
                    continue
                
                check_goes_xray_science_data(stdate, enddate, experiment)
                flare_info = extract_flare_info(experiment, request_date)
                if not pd.isnull(flare_info['catalog_id']):
                    #The time provided in the IGR list aren't always accurate enough to find the correct
                    #flare. Us the deprecated peak value provided to compare with the identified flare.
                    if not pd.isnull(flare_info['intensity']):
                        dep_peak = row['Flare Magnitude Deprecated']
                        sci_peak = flare_info['intensity']
                        ratio = (dep_peak/0.7)/sci_peak
                        #Check if found flare has value within +-10% of expected value after
                        #removing SWPC correction factor
                        if (ratio >= 0.9) and (ratio <= 1.1):
                            break
                        else:
                            flare_info = {}
                else:
                    flare_info = {} #so will manually apply calibration

            #-->If didn't find a measurement in the science data, then apply calibration manually
            if not flare_info and request_date < goes_r_primary_date:
                #DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES
                print(f"extract_flare_info: Manually removing SWPC calibration factor from value provided in Ian Richardson's list.")
                flare_info = manual_flare_calibration(df, index)

        #Steve's List
        #############
        #If row is from Steve's SRAG list, much more information
        elif experiment not in old_goes_sc:
            request_date = row['Flare Xray Peak Time Deprecated']
            print(f"Flare sourced from Steve's list {request_date}")
 
            stdate = request_date - datetime.timedelta(hours=24)
            enddate = request_date + datetime.timedelta(hours=48)
            
            #-->Find flare in experiment specified in Steve's list
            check_goes_xray_science_data(stdate, enddate, experiment)
            flare_info = extract_flare_info(experiment, request_date)
            
            #-->If that particular satellite didn't have the needed data, search for flare
            if not goes_xray_sat_covers_date(experiment, request_date) or pd.isnull(flare_info['catalog_id']):
                for experiment in xray_sat:
                    if not goes_xray_sat_covers_date(experiment, request_date):
                        continue
                    
                    check_goes_xray_science_data(stdate, enddate, experiment)
                    flare_info = extract_flare_info(experiment, request_date)
                    if not pd.isnull(flare_info['catalog_id']):
                        break

            #-->If didn't find a measurement in the science data, then apply calibration manually
            if pd.isnull(flare_info['catalog_id']) and request_date < goes_r_primary_date:
                #DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES
                print(f"extract_flare_info: Manually removing SWPC calibration factor from value provided in Steve Johnson's SRAG list.")
                flare_info = manual_flare_calibration(df, index)

        #If none of the above and no experiment
        #######################################
        elif pd.isnull(experiment) and request_date < goes_r_primary_date:
            #Manually DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES
            print(f"extract_flare_info: No experiment specified and date prior to {goes_r_primary_date}. Manually removing SWPC calibration factor per NOAA guidance.")
            flare_info = manual_flare_calibration(df, index)


        #Check for missing/no results
        #############################
        #Didn't find any flare information for any experiments
        if not flare_info: continue

        #If found the flare in the flare summary files, but it spans files on two days
        if not pd.isnull(flare_info['catalog_id']):
            if not pd.isnull(flare_info['start_time']) and pd.isnull(flare_info['end_time']):
                print(flare_info)
                #Look for the flare end in the file for the next day
                flare_info = extract_flare_info(flare_info['instrument'], request_date+datetime.timedelta(hours=24), flare_info=flare_info, req_flare_id=flare_info['catalog_id'])

        #If experiment covered date range, but didn't have requested flare info
        if pd.isnull(flare_info['intensity']): continue

        #Save results
        #############
        #If some flare information was found, enter in list
        df.loc[index,'Flare Xray Start Time'] = flare_info['start_time']
        df.loc[index,'Flare Xray Peak Time'] = flare_info['peak_time']
        df.loc[index,'Flare Xray End Time'] = flare_info['end_time']
        df.loc[index,'Flare Class'] = flare_info['class']
        df.loc[index,'Flare Magnitude'] = flare_info['intensity']
        df.loc[index,'Flare Integrated Flux'] = flare_info['integrated_intensity']
        df.loc[index,'Flare Duration'] = flare_info['duration']
        df.loc[index,'Flare Xray Time To Peak'] = flare_info['time_to_peak']
        df.loc[index,'Flare Catalog ID'] = flare_info['catalog_id']
        if not pd.isnull(flare_info['instrument']):
            exper = flare_info['instrument'].lstrip('GOES-')
            df.loc[index,'GOES Xray Satellite'] = int(exper)
        else:
            df.loc[index,'GOES Xray Satellite'] = None
        
    return df


##########################################################
####################### CMEs #############################
##########################################################
#Original DONKI scraping code provided by Luke Stegeman
def reject_entry(data_dict, filters):

    if filters == {}:
        return False

    #Reject plane-of-sky measurements
    if 'sky' in data_dict['measurementTechnique'] or 'Sky' in data_dict['measurementTechnique']:
        return True

    for parameter, minimum in filters.items():
        if minimum is not None:
            if data_dict[parameter] < minimum:
                return True
    return False


def convert_to_time(tzulu):
    """ Zulu string to datetime """
    if tzulu == None:
        return None
    time = datetime.datetime.strptime(tzulu, '%Y-%m-%dT%H:%MZ')
    return time


def convert_to_float(val):
    if val == None:
        return np.nan
    else:
        return float(val)


def get_cmes(start_date, end_date, minimum_speed=None, minimum_halfAngle=None):
    """ Extract the most relevant CME analysis for each CME """

    #Start/End date in YYYY-MM-DD format
    url = 'https://kauai.ccmc.gsfc.nasa.gov/DONKI/WS/get/CME?startDate=' + start_date + '&endDate=' + end_date
    try:
        response = requests.get(url, timeout=3)
        cmes_ = response.json()
    except:
        print(f"get_cmes: could not make connection for {url}")
        return []

    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    signature_args = inspect.signature(get_cmes).parameters.keys()
    arg_dict = {arg : values[arg] for arg in signature_args}
    
    filters = {}
    for arg, value in arg_dict.items():
        if not (arg in ['start_date', 'end_date']):
            if value is not None:
                filters[arg.lstrip('minimum_')] = value

    cmes = []
    for cme in cmes_:
        cme_analyses = cme['cmeAnalyses']

        if cme_analyses is None:
            continue

        #Pull the most accuracte CMEs
        cme_down_select = [cme_analysis for cme_analysis in cme_analyses if cme_analysis['isMostAccurate']]
        #Exclude plane-of-sky measurements
        cme_down_select = [cme_analysis for cme_analysis in cme_down_select if 'sky' not in cme_analysis['measurementTechnique']]

        if len(cme_down_select) == 0:
            continue

        best_cme = {}
        if len(cme_down_select) == 1:
            best_cme = cme_down_select[0]
        else:
            #Select the 0th entry - generally most relevant or recent measurement
            best_cme = cme_down_select[0]
            for cme_analysis in cme_down_select:
                #Give preference to SEPVAL2023 CME analyses
                if 'SEPVAL' in cme_analysis['note']:
                    #Prefer leading edge measurements over SH, if available
                    if cme_analysis['featureCode'] == 'LE' or cme_analysis['featureCode'] is None:
                            best_cme = cme_analysis
                elif best_cme['levelOfData'] < cme_analysis['levelOfData']:
                    #Prefer leading edge measurements over SH, if available
                    if cme_analysis['featureCode'] == 'LE':
                        if best_cme['featureCode'] is not None:
                            best_cme = cme_analysis
                else:
                    #Prefer leading edge measurements over SH, if available
                    if cme_analysis['featureCode'] == 'LE':
                        if best_cme['featureCode'] is not None:
                            best_cme = cme_analysis


        if not best_cme:
            continue
            
        cme['cmeAnalysisPrime'] = best_cme

        reject = reject_entry(cme['cmeAnalysisPrime'], filters)
        if not reject:
            cmes.append(cme)
    return cmes


def get_cme_start_time(cmes):
    cme_info = []
    for i, cme in enumerate(cmes):
        if cme['startTime'] is None:
            cme_info.append(None)
        else:
            time = convert_to_time(cme['startTime'])
            cme_info.append(time)
    return cme_info


def get_cme_header_info(cmes, quantity):
    """ Information related to CME included in header """

    cme_info = [None]*len(cmes)
    for i, cme in enumerate(cmes):
        info = cme[quantity]
        if info is None:
            info = None
        cme_info[i] = info
    return cme_info


def get_cme_info(cmes, quantity):
    """ Information from CME analysis. """
    
    cme_info = [None]*len(cmes)
    for i, cme in enumerate(cmes):
        info = cme['cmeAnalysisPrime'][quantity]
        if info is None:
            info = None
        cme_info[i] = info
    return cme_info



def get_cme(starttime):
    """ Get a specific DONKI CME identified by start time (first look time).
        
        INPUTS:
        
            :starttime: (datetime) timestamp close to the start time of the CME, 
                typically defined as the first time at which the CME appears on the
                coronagraph. For SOHO, will be the LASCO C2 timestamp.
                
        OUTPUTS:
        
            :cme: (dict) CME info save in a dictionary
    
    """
    startdate = starttime - datetime.timedelta(hours=6)
    enddate = starttime + datetime.timedelta(hours=6)
    
    start_date = startdate.strftime('%Y-%m-%d')
    end_date = enddate.strftime('%Y-%m-%d')
    
    cmes = get_cmes(start_date, end_date, minimum_speed=450, minimum_halfAngle=15)

    #CME ANALYSES
    cme_start_times = get_cme_start_time(cmes)
    cme_speeds = get_cme_info(cmes, 'speed')
    cme_half_angles = get_cme_info(cmes, 'halfAngle')
    cme_longitudes = get_cme_info(cmes, 'longitude')
    cme_latitudes = get_cme_info(cmes, 'latitude')
    cme_21_5 = get_cme_info(cmes, 'time21_5')
    cme_techniques = get_cme_info(cmes, 'measurementTechnique')
    cme_feature = get_cme_info(cmes, 'featureCode')

    #CME HEADER INFO
    cme_ids = get_cme_header_info(cmes,'activityID')
    cme_notes = get_cme_header_info(cmes, 'note')
    cme_submission_times = get_cme_header_info(cmes, 'submissionTime')
    cme_linked_events = get_cme_header_info(cmes, 'linkedEvents')
    cme_links = get_cme_header_info(cmes, 'link')

    #Find requested CME
    td = datetime.timedelta(hours=2, minutes=30) #At least within 2.5 hours
    ix = -1
    for i, time in enumerate(cme_start_times):
        diff = abs(time - starttime)
        if diff < td:
            td = diff
            ix = i
    
    if ix == -1:
        print(f"get_cme: No DONKI CME found for {starttime}.")
        return {}

    cme_select = {"DONKI CME Start Time": cme_start_times[ix],
                  "DONKI CME Speed": convert_to_float(cme_speeds[ix]),
                  "DONKI CME Half Width": convert_to_float(cme_half_angles[ix]),
                  "DONKI CME Lat": convert_to_float(cme_latitudes[ix]),
                  "DONKI CME Lon": convert_to_float(cme_longitudes[ix]),
                  "DONKI CME Time at 21.5": convert_to_time(cme_21_5[ix]),
                  "DONKI CME Catalog ID": cme_ids[ix],
                  "DONKI CME Submission Time": convert_to_time(cme_submission_times[ix]),
                  "DONKI CME Feature": cme_feature[ix],
                  "DONKI CME Note": cme_notes[ix],
                  "DONKI CME Linked Events": cme_linked_events[ix],
                  "DONKI CME Link": cme_links[ix]
                  }

    return cme_select


def update_srag_list_donki_cme(df):
    """ Add DONKI CME information to SRAG SEP list """
    
    #Add a column
    # Get the index of column 'A'
    idx = df.columns.get_loc('DONKI CME Time at 21.5')

    # Insert a new column named 'B' after 'A'
    df.insert(loc=idx + 1, column='DONKI CME Feature', value=[None]*len(df['DONKI CME Time at 21.5']))
    
    for index, row in df.iterrows():
        starttime = row['CDAW CME First Look Time']
        if pd.isnull(starttime):
            continue
        
        cme_info = get_cme(starttime)
        
        if not cme_info:
            continue

        df.loc[index,'DONKI CME Start Time'] = cme_info['DONKI CME Start Time']
        df.loc[index,'DONKI CME Speed'] = cme_info['DONKI CME Speed']
        df.loc[index,'DONKI CME Half Width'] = cme_info['DONKI CME Half Width']
        df.loc[index,'DONKI CME Lat'] = cme_info['DONKI CME Lat']
        df.loc[index,'DONKI CME Lon'] = cme_info['DONKI CME Lon']
        df.loc[index,'DONKI CME Time at 21.5'] = cme_info['DONKI CME Time at 21.5']
        df.loc[index,'DONKI CME Feature'] = cme_info['DONKI CME Feature']
        df.loc[index,'DONKI CME Catalog ID'] = cme_info['DONKI CME Catalog ID']

    return df


##########################################################
################## SRAG LIST #############################
##########################################################
def update_srag_list_associations():
    """ Update SRAG SEP event list with flare X-ray science data and
        more parameters for DONKI CMEs. 
    
        INPUTS:
            
            None
            
        OUTPUTS:
        
            :df: (pandas DataFrame) updated dataframe of SRAG list
            Write out SRAG list file to fetchsep/references/
    
    """
    df = read_srag_list()
    df = update_srag_list_flares(df)
    df = update_srag_list_donki_cme(df)

    outfname = os.path.join('fetchsep','reference','SRAG_SEP_List_R11_CLEARversion_UPDATED.csv')
    df.to_csv(outfname)
    print(f"Wrote updated SEP list to {outfname}")
