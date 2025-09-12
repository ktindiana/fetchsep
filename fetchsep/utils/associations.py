import pandas as pd
import os
import sys
import string
import datetime
import numpy as np


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
        ns = 'S'
        
    ew = ''
    if lon >= 0:
        ew = 'W'
    else:
        ew = 'E'
        
    location = f"{ns}{lat:02d}{ew}{lon:02d}"
    return location


################### STEVE JOHNSON'S SRAG SEP LIST #########################
def is_time_column(col):
    columns = ['InitiationDT', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare X-ray End Time', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio DH Start Time', 'Radio DH End Time', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'CME CDAW First Look Time', 'P10_FluxStart', 'P10_StartDT', 'P10_OnsetMax_DT', 'P10_PeakDT', 'P10_EndDT', 'P50_FluxStart', 'P50_StartDT', 'P50_OnsetPeakDT', 'P50_PeakDT', 'P50_EndDT', 'P100_FluxStart', 'P100_StartDT', 'Onset_P100_PkDT', 'Event_P100_PkDT', 'P100_EndDT','ESP Flare Xray Start Time', 'ESP Flare Xray Peak Time', 'ESP Flare X-ray End Time', 'ESP Radio m_TyII Start Time', 'ESP Radio m_TyII End Time', 'ESP CME First Look Time', 'ACE CME passage', 'Sudden Imp Time']

    if col in columns:
        return True
    else:
        return False


def is_float_column(col):
    columns = ['Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Xray Time To Peak', 'AR Area', 'AR Carrington', 'Event Location From Center', 'Event Latitude', 'Event Longitude', 'Event Location from Center 2', 'Event Latitude 2', 'Event Longitude 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME Speed', 'DONKI CME Speed', 'CME Width', 'CME Mean Position Angle']

    if col in columns:
        return True
    else:
        return False


def is_int_column(col):
    columns = ['Cycle', 'Active Region', 'GOES Xray Satellite', 'GOES Protons Satellite']

    if col in columns:
        return True
    else:
        return False


def is_str_column(col):
    columns =  ['EventType', 'Case', 'Flare Class', 'Flare Opt', 'AR Spot Class', 'AR Mag Class', 'Event Location Source', 'Event Location Source 2', 'Radio Station', 'Radio DH Note', 'ESP_CME', 'GLE Event Number', 'PRF','Comments']

    if col in columns:
        return True
    else:
        return False
    

def cast_srag_list_time_columns(df):
    columns = ['InitiationDT', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare X-ray End Time', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio DH Start Time', 'Radio DH End Time', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'CME CDAW First Look Time', 'P10_FluxStart', 'P10_StartDT', 'P10_OnsetMax_DT', 'P10_PeakDT', 'P10_EndDT', 'P50_FluxStart', 'P50_StartDT', 'P50_OnsetPeakDT', 'P50_PeakDT', 'P50_EndDT', 'P100_FluxStart', 'P100_StartDT', 'Onset_P100_PkDT', 'Event_P100_PkDT', 'P100_EndDT']

    for col in columns:
        df[col] = pd.to_datetime(df[col])
        
    return df
    
    

def srag_list_clear_columns():
    columns = ['Cycle', 'EventType', 'Case', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare X-ray End Time', 'Flare Class', 'Flare Opt', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Xray Time To Peak', 'Active Region', 'AR Area', 'AR Spot Class', 'AR Mag Class', 'AR Carrington', 'Event Location From Center', 'Event Latitude', 'Event Longitude', 'Event Location Source', 'Event Location from Center 2', 'Event Latitude 2', 'Event Longitude 2', 'Event Location Source 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio Station', 'Radio DH Start Time', 'Radio DH End Time', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio DH Note', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CME CDAW First Look Time', 'CDAW CME Speed', 'DONKI CME Speed', 'CME Width', 'CME Mean Position Angle', 'ESP_CME', 'GLE Event Number', 'PRF', 'Comments']

    return columns


def empty_srag_associations_series():
    """ Create a pandas Series with appropriate null values 
        for the associations from the SRAG list (used in 
        the CLEAR benchmark dataset).
    
    """

    clear_col = srag_list_clear_columns()
    
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
    df = combine_srag_list_comments(df)
    
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
    
    tolerance = datetime.timedelta(hours=4)
    columns = ["P10_FluxStart", "P10_StartDT"] #["P10_FluxStart", "InitiationDT", "P10_StartDT"]
    
    diff = df[columns] - startdate
    
    event_index = -1 #min index for each column
    check_min_diff = np.nan
    for col in columns:
        ix = []
        min_diff = min(abs(diff[col]))
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
        clear_col = srag_list_clear_columns()
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
    

    cme = ccmc_cme_block()
    
    if not pd.isnull(associations['CME CDAW First Look Time']):
        cme['start_time'] = associations['CME CDAW First Look Time']
    else:
        del cme['start_time']
        
    if not pd.isnull(associations['Event Latitude']):
        cme['lat'] = associations['Event Latitude']
    else:
        del cme['lat']
        
    if not pd.isnull(associations['Event Longitude']):
        cme['lon'] = associations['Event Longitude']
    else:
        del cme['lon']
    
    if not pd.isnull(associations['CME Mean Position Angle']):
        cme['pa'] = associations['CME Mean Position Angle']
    else:
        del cme['pa']
    
    width = associations['CME Width']
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


    if not pd.isnull(associations['CDAW CME Speed']):
        cme['speed'] = associations['CDAW CME Speed']
        cme['catalog'] = 'SOHO_CDAW'
        cme['derivation_technique']['process'] = 'manual'
        cme['derivation_technique']['method'] = 'Plane-of-sky'
    else:
        del cme['speed']
        del cme['catalog']
        del cme['derivation_technique']



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
    
    return cme
    
    
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
    flare = ccmc_flare_block()
    
    if not pd.isnull(associations['Flare Xray Start Time']):
        flare['start_time'] = associations['Flare Xray Start Time']
    else:
        del flare['start_time']

    if not pd.isnull(associations['Flare Xray Peak Time']):
        flare['peak_time'] = associations['Flare Xray Peak Time']
    else:
        del flare['peak_time']
        
    if not pd.isnull(associations['Flare X-ray End Time']):
        flare['end_time'] = associations['Flare X-ray End Time']
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
    
    #Remove unused entries
    del flare['last_data_time']
    del flare['peak_ratio']
    del flare['catalog_id']
    del flare['urls']
    
    return flare
