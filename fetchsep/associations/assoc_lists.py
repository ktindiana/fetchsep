from . import fetch_flare
from . import fetch_cme
from . import fetch_cycle
from ..utils import config as cfg
from ..json import ccmc_json_handler as ccmc_json
from ..utils import date_handler as dh
from ..utils import experiments as expts
import pandas as pd
import os
import sys
import string
import datetime
import numpy as np
import netCDF4 as nc


#GOES spacecraft categories needed for GOES X-ray data
#Spacecraft in the GOES-R+ series
goes_R = expts.goes_R()
#Spacecraft prior to GOES-R
goes_sc = expts.goes_sc()
#Spacecraft prior to GOES-08
old_goes_sc = expts.old_goes_sc()

#For use to convert Steve's list to a machine-readable list and remove
#hidden characters from excel
def clean_file(input_path, output_path):
    """ Clean hidden ascii characters from files. 
        Useful when converting a list from Excel to plain text.
    
    """
    # Define allowed characters (printable + common whitespace)
    allowed = set(string.printable)

    with open(input_path, "r", errors="ignore") as infile, open(output_path, "w") as outfile:
        for line in infile:
            # Keep only allowed characters
            cleaned = "".join(ch for ch in line if ch in allowed)
            outfile.write(cleaned)


def string_location(lat, lon):
    """ Put lat, lon in string location format 
    
        INPUTS:
        
            :lat: (float, int) latitude
            :lon: (float, int) longitude
            
        OUTPUTS:
        
            :location: (string) e.g. N20W35
            
    """
    
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


def numerical_location(location):
    """ Put string location into lat, lon 
 
        INPUTS:
                    
            :location: (string) e.g. N20W35
            
        OUTPUTS:
        
            :lat: (float) latitude
            :lon: (float) longitude
    
    """
    location = locations.strip()
    
    if pd.isnull(location) or location == '':
        return np.nan, np.nan
    
    if 'W' in location:
        location = location.split('W')
        lon = float(location[1])
    elif 'E' in location:
        location = location.split('E')
        lon = -float(location[1])

    if 'N' in location[0]:
        lat = float(location[0].split('N')[1])
    elif 'S' in location[0]:
        lat = -float(location[0].split('S')[1])

    return lat, lon



###########################################################################
################## GENERAL ASSOCIATION LIST FUNCTIONS #####################
###########################################################################
#Columns output in final SEP list, e.g. for the CLEAR SEP Benchmark Dataset
#All possible association information
list_output_columns = ['Solar Cycle', 'EventType', 'Case', 'Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Class Deprecated', 'Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Flare Class', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Time To Peak', 'Flare Catalog ID', 'GOES Xray Satellite', 'Flare Opt', 'Active Region', 'AR Area', 'AR Spot Class', 'AR Mag Class', 'AR Carrington', 'Event Source Location From Center', 'Event Source Latitude', 'Event Source Longitude', 'Event Source Reference', 'Event Source Location from Center 2', 'Event Source Latitude 2', 'Event Source Longitude 2', 'Event Source Reference 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio Type III Start Time', 'Radio Type III End Time', 'Radio TyIII_Imp', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio Station', 'Radio DH Start Time', 'Radio DH End Time', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio DH Note', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME First Look Time', 'CDAW CME Speed', 'CDAW CME Width', 'CDAW CME Mean Position Angle', 'CDAW CME Central Position Angle', 'DONKI CME Start Time', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', 'DONKI CME Time at 21.5', 'DONKI CME Feature', 'DONKI CME Catalog ID', 'DONKI CME Measurement Technique', 'ESP_CME', 'ESP CDAW CME First Look Time', 'ESP CDAW CME LASCO Speed', 'Shock Passage Speed', 'ACE CME Passage Time', 'ACE Bz', 'Sudden Impulse Time', 'Sudden Impulse Amplitude', 'GLE Event Number', 'PRF', 'Associations List', 'Comments']

list_time_columns = ['Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio DH Start Time', 'Radio DH End Time', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'CDAW CME First Look Time', 'DONKI CME Start Time', 'DONKI CME Time at 21.5', 'ESP Flare Xray Start Time Deprecated', 'ESP Flare Xray Peak Time Deprecated', 'ESP Flare Xray End Time Deprecated', 'ESP Radio m_TyII Start Time', 'ESP Radio m_TyII End Time', 'Radio Type III Start Time', 'Radio Type III End Time', 'ESP CDAW CME First Look Time', 'ACE CME Passage Time', 'Sudden Impulse Time']

list_float_columns = ['Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Duration Deprecated', 'Flare Time To Peak Deprecated','Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Time To Peak', 'AR Area', 'AR Carrington', 'Event Source Location From Center', 'Event Source Latitude', 'Event Source Longitude', 'Event Source Location from Center 2', 'Event Source Latitude 2', 'Event Source Longitude 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME Speed', 'CDAW CME Mean Position Angle', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', 'ESP CDAW CME LASCO Speed', 'Shock Passage Speed', 'ACE Bz', 'Sudden Impulse Amplitude']

list_int_columns = ['Flare Catalog ID', 'Solar Cycle', 'Active Region', 'GOES Xray Satellite']

list_string_columns =  ['EventType', 'Case', 'Flare Class Deprecated', 'Flare Class', 'Flare Opt', 'AR Spot Class', 'AR Mag Class', 'Event Source Reference', 'Event Source Reference 2', 'Radio Station', 'Radio DH Note', 'CDAW CME Width', 'CDAW CME Central Position Angle', 'DONKI CME Feature', 'DONKI CME Catalog ID', 'DONKI CME Measurement Technique', 'ESP_CME', 'GLE Event Number', 'PRF', 'List', 'Comments']


flare_columns = ['Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time',
    'Flare Class', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration',
    'Flare Time To Peak', 'Flare Catalog ID', 'GOES Xray Satellite']

donki_cme_columns = ['DONKI CME Start Time', 'DONKI CME Speed', 'DONKI CME Half Width',
    'DONKI CME Lat', 'DONKI CME Lon', 'DONKI CME Time at 21.5', 'DONKI CME Feature',
    'DONKI CME Measurement Technique', 'DONKI CME Catalog ID']

cdaw_cme_columns = ['CDAW CME First Look Time', 'CDAW CME Speed', 'CDAW CME Width',
    'CDAW CME Central Position Angle', 'CDAW CME Mean Position Angle']

source_loc_columns = ['Event Source Location From Center', 'Event Source Latitude',
    'Event Source Longitude', 'Event Source Reference']


def is_time_column(col):
    if col in list_time_columns:
        return True
    else:
        return False


def is_float_column(col):
    if col in list_float_columns:
        return True
    else:
        return False


def is_int_column(col):
    if col in list_int_columns:
        return True
    else:
        return False


def is_str_column(col):
    if col in list_string_columns:
        return True
    else:
        return False


def cast_value(col, val):
    """ Cast value according to its type """
    if pd.isnull(val):
        return val
    
    if is_float_column(col):
        return float(val)

    if is_int_column(col):
        return int(val)
        
    if is_str_column(val):
        return str(val)

    return val


#Dictionary with appropriate null values for the columns in list_output_col
#This associations dictionary matches the columns in Steve's SEP event list
#and that are exported out to the .csv file after each SEP event is analyzed
def empty_associations_dict():
    """ Create a pandas Series with appropriate null values 
        for the associations that will be output by FetchSEP, 
        e.g., used in the CLEAR benchmark dataset).
    
    """
    
    dict = {}
    for col in list_output_columns:
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


def add_solar_cycle(request_date, associations={}):
    """ Fill in the Solar Cycle column in associations """
    sc = fetch_cycle.get_solar_cycle(request_date)
    if not associations:
        associations = empty_associations_dict()
    
    associations['Solar Cycle'] = sc
    return associations
    

def add_farside(associations):
    """ Given an associations dictionary, check the Event Source Longitude 
        to determine if it is a farside event and add update Case column.
        
    """
    is_farside = False
    #Defer to specific event source location entries
    if not pd.isnull(associations['Event Source Longitude']):
        if associations['Event Source Longitude'] < -90 or associations['Event Source Longitude'] > 90:
            is_farside = True

    #Take advantage of DONKI CME measurements, if available
    elif not pd.isnull(associations['DONKI CME Lon']):
        if associations['DONKI CME Lon'] < -90 or associations['DONKI CME Lon'] > 90:
            is_farside = True

    if is_farside:
        case = associations['Case']
        if 'Farside' not in case:
            if pd.isnull(case) or case == '':
               associations['Case'] = 'Farside'
            else:
                associations['Case'] = f"{associations['Case']}/Farside"

    return associations


#Remove the SWPC calibration factor using the associations format
def manual_flare_calibration(df, index, remove_swpc=True):
    """ Convert deprecated flare xray columns into fluxes that align better
        with science-quality Xray data by removing SWPC calibration factors.
    
        INPUTS:
            
            :df: (pandas DataFrame) list of associations (SRAG, IGR, user)
            :index: (int) index of row that has old flare info
            :remove_swpc: (bool) manually remove SWPC 0.7 calibration factor
            
        OUTPUTS:
        
            :flare_info: (dict) contains flare info from DEPRECATED columns
                with SWPC factor removed (if True)
    
    """

    flare_info = fetch_flare.flare_info_dict() #null
    flare_info['catalog'] = 'SWPC_MANUAL'

    if pd.isnull(df['Flare Magnitude Deprecated'].loc[index]):
        return flare_info

    if 'GOES Xray Satellite' in df.columns:
        goes = df['GOES Xray Satellite'].loc[index]
        if not pd.isnull(goes):
            goes = int(goes)
            experiment = f"GOES-{goes:02d}"
        else:
            experiment = None

        flare_info['instrument'] = experiment
    
    if 'Flare Xray Start Time Deprecated' in df.columns:
        flare_info['start_time'] = df['Flare Xray Start Time Deprecated'].loc[index]

    if 'Flare Xray Peak Time Deprecated' in df.columns:
        flare_info['peak_time'] = df['Flare Xray Peak Time Deprecated'].loc[index]

    if 'Flare Xray End Time Deprecated' in df.columns:
        flare_info['end_time'] = df['Flare Xray End Time Deprecated'].loc[index]

    if 'Flare Duration Deprecated' in df.columns:
        flare_info['duration'] = df['Flare Duration Deprecated'].loc[index]

    if 'Flare Time To Peak Deprecated' in df.columns:
        flare_info['time_to_peak'] = df['Flare Time To Peak Deprecated'].loc[index]
    
    if 'Flare Integrated Flux Deprecated' in df.columns:
        fluence = df['Flare Integrated Flux Deprecated'].loc[index]
        if remove_swpc:
            fluence = fetch_flare.remove_swpc_calibration(fluence)
        flare_info['integrated_intensity'] = fluence
    
    peak_intensity = df['Flare Magnitude Deprecated'].loc[index]
    if remove_swpc:
        peak_intensity = fetch_flare.remove_swpc_calibration(peak_intensity)
    fl_class = fetch_flare.flare_class(peak_intensity)
    flare_info['intensity'] = peak_intensity
    flare_info['class'] = fl_class

    return flare_info


#Add Flare columns (e.g. if have a list with only DEPRECATED columns and values)
def add_flare_columns_to_list(df):
    """ Add columns for NOAA science flare values """
    for col in flare_columns:
        if col not in df.columns:
            df.insert(loc=len(df.columns)-1, column=col, value=[None]*len(df))

    return df


#CCMC SEP Scoreboard flare trigger JSON block converted to the associations dictionary format
def flare_ccmc_json_to_associations(flare_json, associations={}):
    """ Convert CCMC SEP Scoreboard json trigger block into an associations 
        dictionary with columns in list_output_col. 
        
        INPUTS:
        
            :flare_json: (dict) in the format of CCMC's SEP Scoreboard schema
            :associations: (dict) a dictionary in the format of empty_associations_dict() 
                containing all (or a subset of) the columns in the output_list_col variable
                used in the lists. 
                If not specified, a null associations dictionary is created and the flare 
                information is filled in. 
                If specified, the associations dictionary will be updated with the flare
                information.
                
        OUTPUTS:
        
            :associations: (dict) dictionary containing all columns in output_list_col 
                with updated flare columns.
        
    """
    if not associations:
        associations = empty_associations_dict()
    #Make sure that all values aren't null; don't want to overwrite
    flare_json = ccmc_json.clean_trigger_block(flare_json)
    if not flare_json:
        return associations

    if 'start_time' in flare_json.keys():
        associations['Flare Xray Start Time'] = dh.str_to_datetime(flare_json['start_time'])
    else:
        #If associations was provided as an input and already had values
        #in these columns, but sure to reset to null so that the
        #manually input flare and existing flare values aren't mixed.
        associations['Flare Xray Start Time'] = pd.NaT

    if 'peak_time' in flare_json.keys():
        associations['Flare Xray Peak Time'] = dh.str_to_datetime(flare_json['peak_time'])
    else:
        associations['Flare Xray Peak Time'] = pd.NaT

    if 'end_time' in flare_json.keys():
        associations['Flare Xray End Time'] = dh.str_to_datetime(flare_json['end_time'])
    else:
        associations['Flare Xray End Time'] = pd.NaT


    if 'start_time' in flare_json.keys():
        if 'end_time' in flare_json.keys():
            duration = (dh.str_to_datetime(flare_json['end_time']) - dh.str_to_datetime(flare_json['start_time'])).total_seconds()/60. #mins
            associations['Flare Duration'] = duration
        if 'peak_time' in flare_json.keys():
            ttp = (dh.str_to_datetime(flare_json['peak_time']) - dh.str_to_datetime(flare_json['start_time'])).total_seconds()/60. #mins
            associations['Flare Time To Peak'] = ttp

    if 'intensity' in flare_json.keys():
        associations['Flare Class'] = fetch_flare.flare_class(flare_json['intensity'])
        associations['Flare Magnitude'] = flare_json['intensity']
    else:
        associations['Flare Class'] = None
        associations['Flare Magnitude'] = np.nan

    if 'integrated_intensity' in flare_json.keys():
        associations['Flare Integrated Flux'] = flare_json['integrated_intensity']
    else:
        associations['Flare Integrated Flux'] = np.nan

    if 'catalog_id' in flare_json.keys():
        associations['Flare Catalog ID'] = flare_json['catalog_id']
    else:
        associations['Flare Catalog ID'] = None

    if 'location' in flare_json.keys():
        lat, lon = numerical_location(flare_json['location'])
        if not pd.isnull(lat):
            assocations['Event Source Latitude'] = lat
        if not pd.isnull(lon):
            associations['Event Source Longitude'] = lon

    return associations


#flare_info dictionary created in fetch_flare.py converted to associations dictionary format
def flare_info_to_associations(flare_info, associations={}):
    """ Convert flare_info dictionary provided by fetch_flare into the 
        dictionary with columns used in the associations lists.
        
        INPUTS:
        
            :flare_info: (dict) in the format output by fetch_flare.flare_info_dict()
            :associations: (dict) a dictionary in the format of empty_associations_dict() 
                containing all (or a subset of) the columns in the output_list_col variable. 
                If not specified, a null associations dictionary is created and the flare 
                information is filled in. 
                If specified, the associations dictionary will be updated with the flare
                information.
                
        OUTPUTS:
        
            :associations: (dict) dictionary containing all columns in output_list_col 
                with updated flare columns.
        
    """
    if not associations:
        associations = empty_associations_dict()

    associations['Flare Xray Start Time'] = flare_info['start_time']
    associations['Flare Xray Peak Time'] = flare_info['peak_time']
    associations['Flare Xray End Time'] = flare_info['end_time']
    associations['Flare Class'] = flare_info['class']
    associations['Flare Magnitude'] = flare_info['intensity']
    associations['Flare Integrated Flux'] = flare_info['integrated_intensity']
    associations['Flare Duration'] = flare_info['duration']
    associations['Flare Time To Peak'] = flare_info['time_to_peak']
    associations['Flare Catalog ID'] = flare_info['catalog_id']
    if not pd.isnull(flare_info['instrument']):
        exper = flare_info['instrument'].lstrip('GOES-')
        associations['GOES Xray Satellite'] = int(exper)
    else:
        associations['GOES Xray Satellite'] = None
    
    return associations


#Read times from a DEPRECATED column in a dataframe list, identify flares
#from NOAA science data, and fill flare columns in the dataframe list
def update_all_flares_in_list(df, window_minutes_minus=15, window_minutes_plus=15, remove_swpc_cal=True):
    """ If an SEP event list contains the original X-ray fluxes made available
        by SWPC in their various products and archived data, add new science-quality
        flare information in additional columns. If flare not found in NOAA archive,
        manually remove SWPC's 0.7 calibration factor.  
        
        The list must contain old flare information in the columns:
            'Flare Xray Start Time Deprecated'
            'Flare Magnitude Deprecated'
            'GOES Xray Satellite' - optional

        If these columns are present, they will be carried over if the SWPC calibration 
        factor is manually removed:
            'Flare Xray Peak Time Deprecated'
            'Flare Xray End Time Deprecated'
            'Flare Duration Deprecated'
            'Flare Time To Peak Deprecated'
            'Flare Integrated Flux Deprecated'
            
        New science-quality X-ray data was released for GOES-08 to GOES-15, which
        should be used instead for space weather understanding, model training,
        and validation.
    
        For each flare in the list, check if the science data is available.
        If not, manually convert the flare values by removing the SWPC calibration.
        
        INPUTS:
        
            :df: (pandas dataframe) SEP associations list
            :window_minutes_minus/_plus: (float) minutes to search before and after reference 
                time to identify flares
            :remove_swpc_cal: (bool) Set to False to use reference flare magnitude without
                removing the SWPC calibration (e.g. flare already science value)
                
        OUTPUTS:
        
            :df: (pandas dataframe) updated with new flare information
        
    """
    #In order of preference; Choose GOES-R first, if available
    xray_sat = ['GOES-16', 'GOES-18', 'GOES-17', 'GOES-19', 'GOES-08', 'GOES-09',
                'GOES-10', 'GOES-12', 'GOES-15', 'GOES-14', 'GOES-11']

    magnitude_col = 'Flare Magnitude Deprecated'

    req_cols = [magnitude_col]
#    if time_col != None: req_cols.append(time_col)
    for col in req_cols:
        if col not in df.columns:
            print(f"update_all_flares_in_list: {col} column must exist and be populated. Returning.")
            return df

    #Look for times in these columns to identify a flare listed in order of preference
    ref_time_columns = ['Flare Xray Peak Time Deprecated', 'Flare Xray Start Time Deprecated']
        
    #Check if the possible reference time columns are present in the list
    for i in range(len(ref_time_columns)-1,-1,-1):
        col = ref_time_columns[i]
        if col not in df.columns:
            ref_time_columns.pop(i)

    #If columns for the flare science data aren't present, add them
    df = add_flare_columns_to_list(df)
    for index, row in df.iterrows():
        #Check if a flare was provided in the list
        if pd.isnull(row[magnitude_col]):
            continue

        flare_info = fetch_flare.flare_info_dict()

        experiment = None
        if 'GOES Xray Satellite' in df.columns:
            goes = df['GOES Xray Satellite'].loc[index]
            if not pd.isnull(goes):
                goes = int(goes)
                experiment = f"GOES-{goes:02d}"

        goes_r_primary_date = datetime.datetime(2017,2,7) #GOES-R+ primary XRS from this date forward

        #Experiment not specified
        ##########################
        #Find a GOES satellite for the requested date
        if experiment == None:
            request_date = pd.NaT
            for ref_col in ref_time_columns:
                if pd.isnull(request_date):
                    request_date = row[ref_col]
            if pd.isnull(request_date):
                continue
                
            print(f"Flare without GOES satellite specified {request_date}")
            #-->Don't know which satellite was used, will default to primary satellite
            request_intensity = row[magnitude_col]
            if request_date < goes_r_primary_date:
                flare_info = fetch_flare.get_noaa_flare(request_date,
                    window_minutes_minus=window_minutes_minus, window_minutes_plus=window_minutes_plus,
                    flare_intensity = request_intensity, remove_swpc_cal=remove_swpc_cal)
            else:
                flare_info = fetch_flare.get_noaa_flare(request_date,
                    window_minutes_minus=window_minutes_minus, window_minutes_plus=window_minutes_plus,
                    flare_intensity = request_intensity, remove_swpc_cal=False)

            #-->If didn't find a measurement in the science data, then apply calibration manually
            #After goes_r_primary_date, SWPC provided already calibrated GOES-R X-ray fluxes
            if pd.isnull(flare_info['catalog_id']):
                remove_swpc = False
                if request_date < goes_r_primary_date:
                    #DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES
                    remove_swpc = True and remove_swpc_cal
                    print(f"update_all_flares_in_list: Manually removing SWPC calibration factor from value provided in list for flare {request_date}.")
                flare_info = manual_flare_calibration(df, index, remove_swpc=remove_swpc)
            


        #Old GOES-07 and previous
        #########################
        elif experiment in old_goes_sc:
            #NEED TO DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES UNTIL NOAA RELEASES SCIENCE DATA
            print(f"update_all_flares_in_list: {experiment} X-ray science data is not yet available from NOAA for {request_date}. Manually removing SWPC calibration factor per NOAA guidance.")
            flare_info = manual_flare_calibration(df, index, remove_swpc=remove_swpc_cal)
            
        #Experiment and other flare info specified
        ##########################################
        #If row is from e.g., Steve's SRAG list, much more information
        elif experiment not in old_goes_sc:
            request_date = pd.NaT
            for ref_col in ref_time_columns:
                if pd.isnull(request_date):
                    request_date = row[ref_col]
            if pd.isnull(request_date):
                continue

            print(f"update_all_flares_in_list: Flare {request_date} requested for {experiment}.")
 
            request_intensity = row[magnitude_col]
            if request_date < goes_r_primary_date:
                flare_info = fetch_flare.get_noaa_flare(request_date, experiment=experiment,
                    window_minutes_minus=window_minutes_minus, window_minutes_plus=window_minutes_plus,
                    flare_intensity = request_intensity, remove_swpc_cal=remove_swpc_cal)
            else:
                flare_info = fetch_flare.get_noaa_flare(request_date, experiment=experiment,
                    window_minutes_minus=window_minutes_minus, window_minutes_plus=window_minutes_plus,
                    flare_intensity = request_intensity, remove_swpc_cal=False)
 
 
            #-->If didn't find a measurement in the science data, then apply calibration manually
            if pd.isnull(flare_info['catalog_id']):
                remove_swpc = False
                if request_date < goes_r_primary_date:
                    #DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES
                    remove_swpc = True and remove_swpc_cal
                    print(f"update_all_flares_in_list: Could not find flare {request_date}. Manually removing SWPC calibration factor from values provided.")
                flare_info = manual_flare_calibration(df, index, remove_swpc=remove_swpc)


        #Check for missing/no results
        #############################
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
        df.loc[index,'Flare Time To Peak'] = flare_info['time_to_peak']
        df.loc[index,'Flare Catalog ID'] = flare_info['catalog_id']
        if not pd.isnull(flare_info['instrument']):
            exper = flare_info['instrument'].lstrip('GOES-')
            df.loc[index,'GOES Xray Satellite'] = int(exper)
        else:
            df.loc[index,'GOES Xray Satellite'] = None
        
    return df


#Add columns for DONKI CMEs
def add_donki_cme_columns(df):
    """ Add columns for DONKI CME values """
    for col in donki_cme_columns:
        if col not in df.columns:
            df.insert(loc=len(df.columns), column=col, value=[None]*len(df))

    return df


#DONKI CME format returned by the DONKI catalog API converted to associations dictionary
def donki_cme_to_associations(cme_info, associations={}):
    """ fetch_cme.get_donki_cme returns a dictionary in the same format
        and keys as used by the associations lists. Add CME values
        to a new or existing associations dictionary.
    
        INPUTS:
        
            :cme_info: (dict) with keys that match the columns names in
                the associations lists, e.g. output from fetch_cme.get_donki_cme
            :associations: (dict) a dictionary in the format of empty_associations_dict() 
                containing all (or a subset of) the columns in the output_list_col variable. 
                If not specified, a null associations dictionary is created and CME 
                information is filled in. 
                If specified, the associations dictionary will be updated with CME
                information. 
                
        OUTPUTS:
        
            :associations: (dict) dictionary with CME information added
        
    """
    if not associations:
        associations = empty_associations_dict()
    
    for col in donki_cme_columns:
        associations[col] = cme_info[col]

    return associations


#CCMC SEP Scoreboard CME trigger block converted to associations dictionary format
def cme_ccmc_json_to_associations(cme_json, associations={}, catalog='DONKI'):
    """ Convert CCMC JSON CME trigger block to the format used by the associations 
        lists. Add CME values to a new or existing associations dictionary.
        CURRENTLY ASSUMES THAT THE CME INFORMATION IS FROM THE DONKI OR CDAW CATALOGs.
        Defaults to DONKI catalog if catalog field not in the input cme_json.
    
        INPUTS:
        
            :cme_json: (dict) CCMC JSON CME trigger block for DONKI or CDAW CME
            :catalog: (string) if not specified in cme_json, may be DONKI or SOHO_CDAW
            :associations: (dict) a dictionary in the format of empty_associations_dict() 
                containing all (or a subset of) the columns in the output_list_col variable. 
                If not specified, a null associations dictionary is created and CME 
                information is filled in. 
                If specified, the associations dictionary will be updated with CME
                information. 
                
        OUTPUTS:
        
            :associations: (dict) dictionary with CME information added
        
    """
    if not associations:
        associations = empty_associations_dict()

    #Make sure that all values aren't null; don't want to overwrite
    cme_json = ccmc_json.clean_trigger_block(cme_json)
    if not cme_json:
        return associations

    if 'catalog' in cme_json.keys():
        #If not specified, use default
        if pd.isnull(cme_json['catalog']) or cme_json['catalog'] == '':
            pass
        #If specified and not DONKI or SOHO_CDAW
        elif cme_json['catalog'] != 'DONKI' and cme_json['catalog'] != 'donki' and cme_json['catalog'] != 'SOHO_CDAW' and cme_json['catalog'] != 'soho_cdaw':
            print("cme_ccmc_json_to_associations: This subroutine is only built for DONKI or CDAW CME measurements presently and will fill in the corresponding columns in the output lists. You input information from a different catalog, so these will not be entered into the output associations list. Future updates to fetchsep will expand these capabilities.")
            return associations
        #If specified and is DONKI or SOHO_CDAW
        else:
            catalog = cme_json['catalog']

    ##### CDAW CME
    if catalog == 'SOHO_CDAW' or catalog == 'soho_cdaw':
        if 'start_time' in cme_json.keys():
            associations['CDAW CME First Look Time'] = dh.str_to_datetime(cme_json['start_time'])
        else:
            #If associations was provided as an input and already had values
            #in these columns, but sure to reset to null so that the
            #manually input CME and existing CME info aren't mixed.
            associations['CDAW CME First Look Time'] = pd.NaT

        if 'speed' in cme_json.keys():
            associations['CDAW CME Speed'] = float(cme_json['speed'])
        else:
            associations['CDAW CME Speed'] = np.nan

        if 'half_width' in cme_json.keys():
            associations['CDAW CME Width'] = float(cme_json['half_width'])*2
        else:
            associations['CDAW CME Width'] = np.nan

        if 'pa' in cme_json.keys():
            associations['CDAW CME Mean Position Angle'] = float(cme_json['pa'])
        else:
            associations['CDAW CME Mean Position Angle'] = np.nan


    #####DONKI CME
    if catalog == 'DONKI' or catalog == 'donki':
        if 'start_time' in cme_json.keys():
            associations['DONKI CME Start Time'] = dh.str_to_datetime(cme_json['start_time'])
        else:
            associations['DONKI CME Start Time'] = pd.NaT

        if 'speed' in cme_json.keys():
            associations['DONKI CME Speed'] = float(cme_json['speed'])
        else:
            associations['DONKI CME Speed'] = np.nan

        if 'half_width' in cme_json.keys():
            associations['DONKI CME Half Width'] = float(cme_json['half_width'])
        else:
            associations['DONKI CME Half Width'] = np.nan

        if 'lat' in cme_json.keys():
            associations['DONKI CME Lat'] = float(cme_json['lat'])
        else:
            associations['DONKI CME Lat'] = np.nan

        if 'lon' in cme_json.keys():
            associations['DONKI CME Lon'] = float(cme_json['lon'])
        else:
            associations['DONKI CME Lon'] = np.nan

        if 'time_at_height' in cme_json.keys():
            if 'time' in cme_json['time_at_height'].keys() and 'height' in cme_json['time_at_height'].keys():
                if cme_json['time_at_height']['height'] == 21.5:
                    associations['DONKI CME Time at 21.5'] = dh.str_to_datetime(cme_json['time_at_height']['time'])
                else:
                    associations['DONKI CME Time at 21.5'] = pd.NaT

        if 'derivation_technique' in cme_json.keys():
            if 'measurement_type' in cme_json['derivation_technique'].keys():
                associations['DONKI CME Feature'] = cme_json['derivation_technique']['measurement_type']
            else:
                associations['DONKI CME Feature'] = None

            if 'method' in cme_json['derivation_technique'].keys():
                associations['DONKI CME Measurement Technique'] = cme_json['derivation_technique']['method']
            else:
                associations['DONKI CME Measurement Technique'] = None

        if 'catalog_id' in cme_json.keys():
            associations['DONKI CME Catalog ID'] = cme_json['catalog_id']
        else:
            associations['DONKI CME Catalog ID'] = None

    return associations


#Given a set of CME start times, identify DONKI CMEs and populated CME columns in a list
def update_all_donki_cmes_in_list(df, minimum_speed=None, minimum_halfAngle=None,
            feature='LE', window_minutes=60, ref_col='CDAW CME First Look Time'):
    """ Given a list with CDAW LASCO CME first look times stored in a column:
            'CDAW CME First Look Time'
            
        Find corresponding DONKI CMEs. Suggested preferred selections for larger SEP events:
            minimum speed: default > 450 km/s
            minimum half angle: default > 15 km/s
            feature: Default leading edge (LE) measurements preferred over shock (SH)
            
            :window_minutes: (int) find CME within window_minutes from requested time
            :ref_col: (string) name of column with reference CME start time; default is
                    'CDAW CME First Look Time'
                    
    """

    if ref_col not in df.columns:
        print(f"update_all_donki_cmes_in_list: {ref_col} column must exist and be populated. Returning.")
        return df

    #Add DONKI CME columns if not present
    df = add_donki_cme_columns(df)


    for index, row in df.iterrows():
        starttime = row[ref_col]
        if pd.isnull(starttime):
            continue
        
        cme_info = fetch_cme.get_donki_cme(starttime, minimum_speed=minimum_speed,
            minimum_halfAngle=minimum_halfAngle, feature=feature, format='dict',
            window_minutes=window_minutes)

        df.loc[index,'DONKI CME Start Time'] = cme_info['DONKI CME Start Time']
        df.loc[index,'DONKI CME Speed'] = cme_info['DONKI CME Speed']
        df.loc[index,'DONKI CME Half Width'] = cme_info['DONKI CME Half Width']
        df.loc[index,'DONKI CME Lat'] = cme_info['DONKI CME Lat']
        df.loc[index,'DONKI CME Lon'] = cme_info['DONKI CME Lon']
        df.loc[index,'DONKI CME Time at 21.5'] = cme_info['DONKI CME Time at 21.5']
        df.loc[index,'DONKI CME Feature'] = cme_info['DONKI CME Feature']
        df.loc[index,'DONKI CME Measurement Technique'] = cme_info['DONKI CME Measurement Technique']
        df.loc[index,'DONKI CME Catalog ID'] = cme_info['DONKI CME Catalog ID']

    return df



#Add columns for DONKI CMEs
def add_cdaw_cme_columns(df):
    """ Add columns for CDAW CME values """
    for col in cdaw_cme_columns:
        if col not in df.columns:
            df.insert(loc=len(df.columns), column=col, value=[None]*len(df))

    return df


#Given a set of CME start times, identify DONKI CMEs and populated CME columns in a list
def update_all_cdaw_cmes_in_list(df, time_col, window_minutes_minus=60, window_minutes_plus=60, speed_col=None):
    """ Given a list with reference times (e.g. flare peak times)
            
        Find corresponding CDAW CMEs from:
        https://cdaw.gsfc.nasa.gov/CME_list/UNIVERSAL/text_ver/univ_all.txt
        
    """

    req_cols = [time_col, speed_col]
    for col in req_cols:
        if pd.isnull(col): continue
        if col not in df.columns:
            print(f"update_all_cdaw_cmes_in_list: {col} column must exist and be populated. Returning.")
            return df

    #Add CDAW CME columns if not present
    df = add_cdaw_cme_columns(df)

    for index, row in df.iterrows():
        starttime = row[time_col]
        if pd.isnull(starttime):
            cme_info = fetch_cme.null_cdaw_cme(type='dict')
        else:
            if not pd.isnull(row[speed_col]):
                speed = row[speed_col]
            else:
                speed = np.nan
            
            cme_info = fetch_cme.get_cdaw_cme(starttime, window_minutes_minus=window_minutes_minus, window_minutes_plus=window_minutes_plus, speed_filter=speed)

        for col in cdaw_cme_columns:
            df.loc[index,col] = cme_info[col]

    return df


#Fill an associations dictionary with information extracted from a list dataframe
def fill_associations(df, index, list=None):
    """ Extract information from SRAG, IGR, user list and put into an associations 
        dictionary. The associations dictionary contains all possible output 
        columns in list_output_col.
        
        The input dataframe must contain all or a subset of columns named exactly
        the same. For non-null entries in df at index, fill in the associations.
        
        INPUTS:
        
            :df: (pandas DataFrame) of an associations list (SRAG, IGR, user)
            :index: (int) index of relevant associations
            :list: (string) list ID
            
        OUTPUT:
        
            :associations: (dict) dictionary with associations updated from the list
    
    """
    associations = empty_associations_dict()
    associations['Associations List'] = list
    for col in df.columns:
        if col in associations.keys():
            val = df[col].iloc[index]
            associations[col] = df[col].iloc[index]

    return associations


#Given an SEP time from OpSEP, find the same SEP event in the associations list
def identify_associations_in_list(sep_start, sep_end, list_name='srag'):
    """ Given a SEP date, identify an event in the associations list associated
        with that date. 
    
        Use the reference time in a list to find the correct association
        with an observed SEP event. 

        
        INPUTS:
        
            :sep_start: (datetime) start time of SEP event that want to
                find in list; try time of enhancement above background
                for >10 MeV or proton energies ~10 MeV or threshold crossing time
            :sep_end: (datetime) end time of SEP event that want to
                find in list; try time of enhancement above background
                for >10 MeV or proton energies ~10 MeV or threshold crossing time
            :list_name: (string) list to pull from for associated flares, CMEs, etc
                "srag" is SRAG_SEP_List.csv
                "igr_25" is IGR_25MeV_SEP_List.csv
                "cane" is Cane_SEP_List.csv
                "igr_soho" is IGR_SEPlist3_20260306.csv 
                "user" is the list the user maintains in fetchsep/reference/user_associations.csv
                
        OUTPUT:
        
            :associations: (pandas Series) Series containing the associations
                fields used for the CLEAR benchmark dataset
        
    """
    null_assoc = empty_associations_dict()

    if pd.isnull(sep_start) or pd.isnull(sep_end):
        return null_assoc

    if list_name == 'srag':
        assoc_list = SRAG_List()
    elif list_name == 'igr_25':
        assoc_list = IGR_25MeV_List()
    elif list_name == 'igr_soho':
        assoc_list = IGR_SOHO_List()
    elif list_name == 'cane':
        assoc_list = Cane_List()
    elif list_name == 'user':
        assoc_list = User_List()
    else:
        sys.exit(f"Specific list {list_name} does not exist. Exiting.")
        

    df = assoc_list.read_list()
    df = assoc_list.calculate_sep_reference_times(df)
    assoc_reference_columns = assoc_list.reference_columns

    if isinstance(sep_start, str):
        sep_start=dh.str_to_datetime(sep_start)
    if isinstance(sep_end, str):
        sep_end=dh.str_to_datetime(sep_end)


    #Check if SEP falls within known SEP start and end
    #It could be that already enhanced events may have start/end times chosen a little
    #differently than in the associations list, causing an offset in overlap.
    #Check both start and end time and select the SEP event in the associations list that
    #has the most overlap.
    #Select events in associations list that contain either the sep start or sep end
    index = df.loc[((df[assoc_reference_columns[0]] <= sep_start) & (sep_start < df[assoc_reference_columns[1]])) | ((df[assoc_reference_columns[0]] < sep_end) | (sep_end <= df[assoc_reference_columns[1]]))].index
    event_index = -1
    
    ####UNIQUELY FOUND SEP BETWEEN REFERENCE FIRST AND LAST TIMES
    #If found the desired association, save and return
    if len(index) == 1:
        event_index = index[0]

    ####CONTINUE TO SEARCH ALLOWING FOR SOME FLEXIBILITY IN TIMES
    #OPTIMIZE ON START TIME
    else:
        tolerance = datetime.timedelta(hours=6) #search horizon

        #Check if sep start falls within know SEP start and end within a certain
        #tolerated time difference; e.g. sep start - 6 hrs < date < sep end + 6 hours
        #For events that are close together, check that the SEP end time isn't before
        #the start time of the known SEP events in the list
        index = df.loc[((df[assoc_reference_columns[0]]-tolerance) <= sep_start) & (sep_start < (df[assoc_reference_columns[1]]+tolerance)) & (sep_end > df[assoc_reference_columns[0]])].index

        print(index)
        ###NO SEP FOUND WITHIN TOLERANCE
        if len(index) == 0:
            event_index = -1
        
        ####UNIQUELY FOUND SEP IN EXPANDED DATE RANGE
        elif len(index) == 1:
            event_index = index[0]
            
        ####ELSE MULTIPLE MATCHES
        #OPTIMIZE BY CHOOSING CLOSEST START TIME TO sep_start
        elif len(index) > 1:
            col = assoc_reference_columns[0]
            min_diff = tolerance
            for ix in index:
                diff = df[col].iloc[ix] - sep_start
                diff = abs(diff)
                if diff <= min_diff:
                    min_diff = diff
                    event_index = ix
 
    if event_index == -1:
        print(f"identify_associations_in_list: No match was found in {assoc_list.id} for {sep_start}.")
        #Return series with columns with appropriate null values
        return null_assoc
    else:
        output_col = list_output_columns
        associations = fill_associations(df, event_index, list=assoc_list.id)
        proton_info = df[assoc_reference_columns].iloc[event_index]
        print(f"identify_associations_in_list: Analyzed {assoc_reference_columns[0]} {proton_info[0]} to {assoc_reference_columns[1]} {proton_info[1]} found in {assoc_list.id} {sep_start} to {sep_end}.")
        #Cast types in associations to avoid Int64, etc
        for key in associations.keys():
            associations[key] = cast_value(key, associations[key])
 
        return associations

    return null_assoc



#Convert associations dictionary to CCMC SEP Scoreboard CME trigger JSON format
def associations_to_ccmc_cme(associations):
    """ Take dictionary containing one row of associations list with 
        CME information. Put the CME columns into the appropriate fields for
        the CCMC SEP Scoreboard json format.
        
        INPUT:
        
            :associations: (dict) a dictionary containing a single row of
                the list with CME information included, e.g. identified
                with identify_associations_in_list
                
        OUTPUT:
            
            :cme: (dict) CME trigger block in CCMC json schema
        
    """
    
    all_cmes = []

    #DONKI CME
    donki_cme = fetch_cme.donki_cme_to_ccmc_json(associations)
    if donki_cme:
        all_cmes.append(donki_cme)

    
    #CDAW CME
    cdaw_cme = fetch_cme.cdaw_cme_to_ccmc_json(associations)
    if cdaw_cme:
        all_cmes.append(cdaw_cme)

    return all_cmes
    

#Convert associations dictionary to CCMC SEP Scoreboard flare trigger JSON format
def associations_to_ccmc_flare(associations):
    """ Take dictionary containing one row of the associations lists and 
        put the information in the Flare columns into the appropriate fields for
        the CCMC SEP Scoreboard json format.
        
        INPUT:
        
            :associations: (dict) a dictionary containing a single row of
                the list with flare information included, e.g. identified
                with identify_associations_in_list
                
        OUTPUT:
            
            :flare: (dict) Flare trigger block in CCMC json schema
        
    """
    all_flares = []
    
    #SWPC Xray Science Data
    if not pd.isnull(associations['Flare Xray Start Time']):
        times = []
        flare = ccmc_json.ccmc_flare_block()
        times.append(associations['Flare Xray Start Time'])
        flare['start_time'] = dh.time_to_zulu(associations['Flare Xray Start Time'])

        if not pd.isnull(associations['Flare Xray Peak Time']):
            times.append(associations['Flare Xray Peak Time'])
            flare['peak_time'] = dh.time_to_zulu(associations['Flare Xray Peak Time'])
        else:
            del flare['peak_time']
            
        if not pd.isnull(associations['Flare Xray End Time']):
            times.append(associations['Flare Xray End Time'])
            flare['end_time'] = dh.time_to_zulu(associations['Flare Xray End Time'])
        else:
            del flare['end_time']

        #Calculate last_data_time
        s = pd.Series(times)
        last_time = s.max()
        if not pd.isnull(last_time):
            flare['last_data_time'] = dh.time_to_zulu(last_time)
        else:
            del flare['last_data_time']

        if not pd.isnull(associations['Event Source Latitude']) and not pd.isnull(associations['Event Source Longitude']):
            lat = associations['Event Source Latitude']
            lon = associations['Event Source Longitude']
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
        del flare['peak_ratio']
        del flare['urls']
 
        all_flares.append(flare)

    #If have the science version of the flare, don't need to deprecated in the json
    if len(all_flares) > 0:
        return all_flares

    #SWPC X-ray operational data that is now deprecated after GOES-R launched
    if not pd.isnull(associations['Flare Xray Start Time Deprecated']):
        times = []
        times.append(associations['Flare Xray Start Time Deprecated'])
        flare = ccmc_json.ccmc_flare_block()
        flare['start_time'] = dh.time_to_zulu(associations['Flare Xray Start Time Deprecated'])

        if not pd.isnull(associations['Flare Xray Peak Time Deprecated']):
            times.append(associations['Flare Xray Peak Time Deprecated'])
            flare['peak_time'] = associations['Flare Xray Peak Time Deprecated'].to_pydatetime()
        else:
            del flare['peak_time']
            
        if not pd.isnull(associations['Flare Xray End Time Deprecated']):
            times.append(associations['Flare Xray End Time Deprecated'])
            flare['end_time'] = dh.time_to_zulu(associations['Flare Xray End Time Deprecated'])
        else:
            del flare['end_time']

        #Calculate last_data_time
        s = pd.Series(times)
        last_time = s.max()
        if not pd.isnull(last_time):
            flare['last_data_time'] = dh.time_to_zulu(last_time)
        else:
            del flare['last_data_time']

        if not pd.isnull(associations['Event Source Latitude']) and not pd.isnull(associations['Event Source Longitude']):
            lat = associations['Event Source Latitude']
            lon = associations['Event Source Longitude']
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
        del flare['peak_ratio']
        del flare['catalog_id']
        del flare['urls']
 
        all_flares.append(flare)
 
    return all_flares



###########################################################################
################### STEVE JOHNSON'S SRAG SEP LIST #########################
###########################################################################
class SRAG_List:
    def __init__(self):
        """ Steve Johnson's SEP list curated for the Space Radiation Analysis Group (SRAG)
            at NASA Johnson Space Center.
            
            Steve Johnson's SRAG list provides the InitiationDT which is
            related to the flare, CME, or first rise of protons, as indicated
            in InitiationType.
            
            P10_FluxStart indicates when >10 Mev goes above background by Steve's
            estimation. SEP Reference Time column is usually this value.
            If the association was taken from Ian Richardson's list, then ther 
            SEP Refernce Time column might contain the flare start time or
            the time that protons started to rise.
            
            P10_StartDT indicates when >10 MeV crosses 10 pfu.
            
        """
        self.id = "SRAG SEP List"
        self.filename = 'SRAG_SEP_List.csv'
        self.list = pd.DataFrame()
        self.reference_columns = ['First SEP Start Time', 'Last SEP End Time']

        self.time_columns = ['InitiationDT', 'SEP Reference Time', 'Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio DH Start Time', 'Radio DH End Time', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'CDAW CME First Look Time', 'DONKI CME Start Time', 'DONKI CME Time at 21.5', 'P10_FluxStart', 'P10_StartDT', 'P10_OnsetMax_DT', 'P10_PeakDT', 'P10_EndDT', 'P50_FluxStart', 'P50_StartDT', 'P50_OnsetPeakDT', 'P50_PeakDT', 'P50_EndDT', 'P100_FluxStart', 'P100_StartDT', 'Onset_P100_PkDT', 'Event_P100_PkDT', 'P100_EndDT','ESP Flare Xray Start Time Deprecated', 'ESP Flare Xray Peak Time Deprecated', 'ESP Flare Xray End Time Deprecated', 'ESP Radio m_TyII Start Time', 'ESP Radio m_TyII End Time', 'ESP CDAW CME First Look Time', 'ACE CME Passage Time', 'Sudden Impulse Time']

        self.float_columns = ['Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Duration Deprecated', 'Flare Time To Peak Deprecated','Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Time To Peak', 'AR Area', 'AR Carrington', 'Event Source Location From Center', 'Event Source Latitude', 'Event Source Longitude', 'Event Source Location from Center 2', 'Event Source Latitude 2', 'Event Source Longitude 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME Speed', 'CDAW CME Mean Position Angle', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', 'ESP CDAW CME LASCO Speed', 'ACE Bz', 'Sudden Impulse Amplitude']
            #'CDAW CME Width' contains text

        self.int_columns = ['Flare Catalog ID', 'Solar Cycle', 'Active Region', 'GOES Xray Satellite', 'GOES Protons Satellite']

        self.string_columns =  ['EventType', 'Case', 'Flare Class Deprecated', 'Flare Class', 'Flare Opt', 'AR Spot Class', 'AR Mag Class', 'Event Source Reference', 'Event Source Reference 2', 'Radio Station', 'Radio DH Note', 'DONKI CME Feature', 'DONKI CME Catalog ID', 'DONKI CME Measurement Technique', 'ESP_CME', 'GLE Event Number', 'PRF','Comments']

        return

    def is_time_column(self,col):
        if col in self.time_columns:
            return True
        else:
            return False


    def is_float_column(self,col):
        if col in self.float_columns:
            return True
        else:
            return False


    def is_int_column(self,col):
        if col in self.int_columns:
            return True
        else:
            return False


    def is_str_column(self,col):
        if col in self.string_columns:
            return True
        else:
            return False
        

    def cast_time_columns(self,df):
        for col in self.time_columns:
            df[col] = pd.to_datetime(df[col])
            
        return df
        

    def cast_string_columns(self,df):
        for col in self.string_columns:
            df[col] = df[col].fillna('').astype(str)
            df[col] = df[col].astype(str)
            
        return df


    def cast_float_columns(self,df):
        for col in self.float_columns:
            df[col] = df[col].astype(float)
            
        return df
 
 

    def combine_comments(self, df):
        df["Comments1"] = df["Comments1"].where(pd.notna(df["Comments1"]), '')
        df["Comments2"] = df["Comments2"].where(pd.notna(df["Comments2"]), '')
        df["Comments"] = df["Comments1"] + df["Comments2"]
        df.drop("Comments1", axis=1, inplace=True)
        df.drop("Comments2", axis=1, inplace=True)
        return df


    def calculate_sep_reference_times(self, df):
        """ Populate SEP list dataframe with the first time and last time 
            that should be used to  find SEP events in list. The idea is to 
            give the biggest range of times for which the SEP was observed. 
            
            For the SRAG list, we use the very first SEP start time and the
            very last SEP end time across energy channels.
            
            Steve provides timing that should completely capture SEP start
            and peak. The provided end times may not extend as long as those
            calculated by OpSEP.
            
        """
        sep_time_columns = ['SEP Reference Time', 'P10_FluxStart', 'P100_FluxStart', 'P50_FluxStart', 'P10_EndDT', 'P50_EndDT', 'P100_EndDT', 'P10_PeakDT']
        
        assoc_reference_columns = self.reference_columns
        
        if assoc_reference_columns[0] not in df.columns:
            idx = df.columns.get_loc('SEP Reference Time')
            df.insert(loc=idx+1,column=assoc_reference_columns[0], value=[None]*len(df))
        if assoc_reference_columns[1] not in df.columns:
            idx = df.columns.get_loc(assoc_reference_columns[0])
            df.insert(loc=idx+1, column=assoc_reference_columns[1], value=[None]*len(df))

        for index, row in df.iterrows():
            first_time = pd.NaT
            last_time = pd.NaT
            
            for col in sep_time_columns:
                if pd.isnull(row[col]):
                    continue
                
                if row[col] < first_time or pd.isnull(first_time):
                    first_time = row[col]
                if row[col] > last_time or pd.isnull(last_time):
                    last_time = row[col]
    
            df.at[index, assoc_reference_columns[0]] = first_time
            df.at[index, assoc_reference_columns[1]] = last_time
        
        return df


    def read_list(self):
        """ Read in the SRAG list provided by Steve Johnson (SRAG) 
            converted from excel into text format and cleaned.
            
        """
        srag_list = os.path.join('fetchsep','reference',self.filename)
        df = pd.read_csv(srag_list)
        
        #Combine two separate comments columns into one
        if 'Comments1' in df.columns:
            df = self.combine_comments(df)

        #Cast string columns
        df = self.cast_string_columns(df)

        #Cast time columns
        df = self.cast_time_columns(df)
        
        self.list = df
        
        return df


    def write_list(self, df, filename=None):
        """ Write SRAG list to file """
        
        if filename == None:
            filename = self.filename
        
        outfname = os.path.join('fetchsep','reference',filename)
        df.to_csv(outfname, index=False)
        print(f"write_list: Wrote updated SEP list to {outfname}")


    def update_flares(self, df):
        """ Update SRAG SEP event list with flare X-ray science data and
            more parameters for DONKI CMEs. 
        
            INPUTS:
                
                None
                
            OUTPUTS:
            
                :df: (pandas DataFrame) updated dataframe of SRAG list
                Write out SRAG list file to fetchsep/references/
        
        """
        
        #Add columns to the SRAG list in a specific place
        #Add a column for Flare ID at end of flare info
        if 'Flare Catalog ID' not in df.columns:
            idx = df.columns.get_loc('Flare Time To Peak')
            df.insert(loc=idx + 1, column='Flare Catalog ID', value=[None]*len(df['Flare Time To Peak']))

        df = update_all_flares_in_list(df)
        self.list = df

        return df


    def update_cmes(self, df):
        """ Update SRAG SEP event list with flare X-ray science data and
            more parameters for DONKI CMEs. 
        
            INPUTS:
                
                None
                
            OUTPUTS:
            
                :df: (pandas DataFrame) updated dataframe of SRAG list
                Write out SRAG list file to fetchsep/references/
        
        """
        
        #Add columns to the SRAG list in a specific place
        if 'DONKI CME Feature' not in df.columns:
            #Add a column for DONKI CME Feature
            idx = df.columns.get_loc('DONKI CME Time at 21.5')
            df.insert(loc=idx + 1, column='DONKI CME Feature', value=[None]*len(df['DONKI CME Time at 21.5']))

        if 'DONKI CME Measurement Technique' not in df.columns:
            #Add a column for DONKI CME Measurement Technique
            idx = df.columns.get_loc('DONKI CME Feature')
            df.insert(loc=idx + 1, column='DONKI CME Measurement Technique', value=[None]*len(df['DONKI CME Feature']))

        df = update_all_donki_cmes_in_list(df, minimum_speed=450, minimum_halfAngle=15, feature='LE')
        self.list = df

        return df


    def update_all(self):
        """ Update flares and cmes """
        
        df = self.read_list()
        df = self.update_flares(df)
        df = self.update_cmes(df)
        
        self.write_list(df, filename='SRAG_SEP_List_UPDATED.csv')


###########################################################################
#################### IAN RICHARDSON'S SOHO/STEREO LIST ####################
###########################################################################
class IGR_SOHO_List():
    def __init__(self):
        """ A list of flare and CME associations compiled by Ian G. Richardson (University
            of Maryland, NASA GSFC). 
            https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/GQPCXZ
            
            Version of list in FetchSEP has been modified for formatting and clarity.
            
            Covers dates from Dec 2006 - Jan 2024
            
        """
        self.id = "IGR SOHO/STEREO SEP List"
        self.filename = 'IGR_SEPlist3_20260306.csv'
        self.list = pd.DataFrame()
        self.reference_columns = ["SEP Reference Start", "SEP Reference End"]

        self.time_columns = ['Solar Event Date (UT)', 'SEP Event Onset Date (UT)', 'Flare Xray Start Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'CDAW CME First Look Time', 'DONKI CME Start Time', 'DONKI CME Time at 21.5', 'SEP Reference Start', 'SEP Reference End']

        self.float_columns = ['Solar Event-B (deg)', 'I(B)', 'Solar Event-Earth (deg)', 'I(Earth)', 'Solar Event-A (deg)', 'I(A)', 'CME dA (deg)', 'CME V (km/s)', 'Flare Magnitude Deprecated', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Time To Peak', 'Event Source Longitude', 'CDAW CME Speed', 'CDAW CME Width', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', ]

        self.int_columns = ['Radio TyIII_Imp', 'Radio TyII Imp', 'Flare Catalog ID', 'GOES Xray Satellite']

        self.string_columns =  ['Flare Class Deprecated','Flare Class', 'CDAW CME Central Position Angle', 'DONKI CME Feature', 'DONKI CME Measurement Technique', 'Comments']

        return


    def is_time_column(self,col):
        if col in self.time_columns:
            return True
        else:
            return False


    def is_float_column(self,col):
        if col in self.float_columns:
            return True
        else:
            return False


    def is_int_column(self,col):
        if col in self.int_columns:
            return True
        else:
            return False


    def is_str_column(self,col):
        if col in self.string_columns:
            return True
        else:
            return False
        

    def cast_time_columns(self,df):
        for col in self.time_columns:
            df[col] = pd.to_datetime(df[col])
            
        return df
        

    def cast_string_columns(self,df):
        for col in self.string_columns:
            df[col] = df[col].fillna('').astype(str)
            df[col] = df[col].astype(str)
            
        return df


    def cast_float_columns(self,df):
        for col in self.float_columns:
            df[col] = df[col].astype(float)
            
        return df


    def calculate_sep_reference_times(self,df):
        """ Does nothing because already set up with the right information """
        return df


    def read_list(self):
        """ Read in the 25 MeV SOHO/STEREO list provided by Ian Richardson
            and on Harvard Dataverse converted in machine-readable format.
            
        """
        igr_list = os.path.join('fetchsep','reference',self.filename)
        df = pd.read_csv(igr_list)
        
        #Cast string columns
        df = self.cast_string_columns(df)

        #Cast time columns
        df = self.cast_time_columns(df)
        
        self.list = df
        
        return df


    def write_list(self, df, filename=None):
        """ Write IGR List to file """

        if filename == None:
            filename = self.filename

        outfname = os.path.join('fetchsep','reference',filename)
        df.to_csv(outfname, index=False)
        print(f"write_list: Wrote updated SEP list to {outfname}")


###########################################################################
############################ HILARY CANE'S  LIST ##########################
###########################################################################
class Cane_List():
    def __init__(self):
        """ A list of flare and CME associations compiled by Ian G. Richardson (University
            of Maryland, NASA GSFC). 
            https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/GQPCXZ
            
            Version of list in FetchSEP has been modified for formatting and clarity.
            
            Covers dates from Dec 2006 - Jan 2024
            
        """
        self.id = "Cane et al. (2010) SEP List"
        self.filename = 'Cane_SEP_List.csv'
        self.list = pd.DataFrame()
        self.reference_columns = ["SEP Reference Start", "SEP Reference End"]

        self.time_columns = ['Date', 'CME First Look Time', 'Flare Xray Start Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Radio Type III Start Time', 'Radio Type III End Time', 'CDAW CME First Look Time', 'SEP Reference Start', 'SEP Reference End']

        self.float_columns = ['SEP Start Hour', 'Maximum Proton Energy', 'Proton Intensity', 'CME Width', 'CME Speed', 'Flare Magnitude Deprecated', 'Flare Time to Peak Deprecated', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Time To Peak', 'Event Source Latitude', 'Event Source Longitude', 'CDAW CME Speed', 'CDAW CME Width']

        self.int_columns = ['Group', 'Radio TyII Imp', 'Flare Catalog ID', 'GOES Xray Satellite']

        self.string_columns =  ['Type III Time', 'Flare Start Hour', 'H-alpha Location', 'CME First Appearance Hour', 'Flare Class Deprecated','Flare Class', 'CDAW CME Central Position Angle', 'Comments']

        return


    def is_time_column(self,col):
        if col in self.time_columns:
            return True
        else:
            return False


    def is_float_column(self,col):
        if col in self.float_columns:
            return True
        else:
            return False


    def is_int_column(self,col):
        if col in self.int_columns:
            return True
        else:
            return False


    def is_str_column(self,col):
        if col in self.string_columns:
            return True
        else:
            return False
        

    def cast_time_columns(self,df):
        for col in self.time_columns:
            df[col] = pd.to_datetime(df[col])
            
        return df
        

    def cast_string_columns(self,df):
        for col in self.string_columns:
            df[col] = df[col].fillna('').astype(str)
            df[col] = df[col].astype(str)
            
        return df


    def cast_float_columns(self,df):
        for col in self.float_columns:
            df[col] = df[col].astype(float)
            
        return df

    def calculate_sep_reference_times(self,df):
        """ Does nothing because already set up with the right information """
        return df

    def read_list(self):
        """ Read in the Cane et al 2010 SEP list 
            converted from excel into text format and cleaned.
            
        """
        igr_list = os.path.join('fetchsep','reference',self.filename)
        df = pd.read_csv(igr_list)
        
        #Cast string columns
        df = self.cast_string_columns(df)

        #Cast time columns
        df = self.cast_time_columns(df)
        
        self.list = df
        
        return df


    def write_list(self, df, filename=None):
        """ Write IGR List to file """

        if filename == None:
            filename = self.filename

        outfname = os.path.join('fetchsep','reference',filename)
        df.to_csv(outfname, index=False)
        print(f"write_list: Wrote updated SEP list to {outfname}")


###########################################################################
#################### IAN RICHARDSON'S OLDER 25 MeV LIST ####################
###########################################################################
class IGR_25MeV_List():
    def __init__(self):
        """ A list of flare and CME associations compiled by Ian G. Richardson (University
            of Maryland, NASA GSFC). 
            https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/GQPCXZ
            
            Version of list in FetchSEP has been modified for formatting and clarity.
            
            Covers dates from Dec 2006 - Jan 2024
            
        """
        self.id = "IGR 25 MeV SEP List"
        self.filename = 'IGR_25MeV_SEP_List.csv'
        self.list = pd.DataFrame()
        self.reference_columns = ["SEP Reference Start", "SEP Reference End"]

        self.time_columns = ['SEP Reference Start', 'SEP Reference End', 'Flare Xray Start Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'CDAW CME First Look Time']

        self.float_columns = ['25 MEV Proton Intensity (MeV s cm^2 sr)^-1', 'GOES >10 MeV (pfu)', 'CME speed (km/s)', 'Flare Magnitude Deprecated', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Time To Peak', 'Event Source Latitude', 'Event Source Longitude', 'CDAW CME Speed', 'CDAW CME Width']

        self.int_columns = ['GLE Event Number', 'Flare Catalog ID', 'GOES Xray Satellite']

        self.string_columns =  ['Flare Class Deprecated', 'Flare Class', 'CDAW CME Central Position Angle']

        return


    def is_time_column(self,col):
        if col in self.time_columns:
            return True
        else:
            return False


    def is_float_column(self,col):
        if col in self.float_columns:
            return True
        else:
            return False


    def is_int_column(self,col):
        if col in self.int_columns:
            return True
        else:
            return False


    def is_str_column(self,col):
        if col in self.string_columns:
            return True
        else:
            return False
        

    def cast_time_columns(self,df):
        for col in self.time_columns:
            df[col] = pd.to_datetime(df[col])
            
        return df
        

    def cast_string_columns(self,df):
        for col in self.string_columns:
            df[col] = df[col].fillna('').astype(str)
            df[col] = df[col].astype(str)
            
        return df


    def cast_float_columns(self,df):
        for col in self.float_columns:
            df[col] = df[col].astype(float)
            
        return df


    def read_list(self):
        """ Read in the 25 MeV list provided by Ian Richardson from 1967 to 1997 
            converted machine-readable format and cleaned.
            
        """
        igr_list = os.path.join('fetchsep','reference',self.filename)
        df = pd.read_csv(igr_list)
        
        #Cast string columns
        df = self.cast_string_columns(df)

        #Cast time columns
        df = self.cast_time_columns(df)
        
        self.list = df
        
        return df

    def calculate_sep_reference_times(self,df):
        """ Does nothing because already set up with the right information """
        return df

    def write_list(self, df, filename=None):
        """ Write IGR List to file """

        if filename == None:
            filename = self.filename

        outfname = os.path.join('fetchsep','reference',filename)
        df.to_csv(outfname, index=False)
        print(f"write_list: Wrote updated SEP list to {outfname}")



###########################################################################
################### USER-MAINTAINED ASSOCIATION LIST ######################
###########################################################################

class User_List:
    def __init__(self, type="minimal"):
        """ A list of SEP events with their flare, CME, and other associations.
            The user maintains this list, which can be used in place of, 
            or in addition to, the SRAG list.
            
            This list may have all the columns up to the ones listed in 
            list_output_columns. The column names must be exactly the same.
            
            The initial list may be created with a minimal set of columns
            for only flare and CME information, or all columns.
            
            INPUTS:
            
                :type: (string) minimal for only flare and CME columns,
                    all for all possible association columns
            
        """
        self.id = "User List"
        self.reference_columns = ["First SEP Start Time", "Last SEP End Time"]
        self.filename = os.path.join('fetchsep','reference', 'user_associations.csv')
        #All possible columns in user list
        self.time_columns = self.reference_columns + list_time_columns
        self.float_columns = list_float_columns
        self.int_columns = list_int_columns
        self.string_columns =  ['SEP Location'] + list_string_columns

        self.create_list(type=type)

        return


    def is_time_column(self,col):
        if col in self.time_columns:
            return True
        else:
            return False


    def is_float_column(self,col):
        if col in self.float_columns:
            return True
        else:
            return False


    def is_int_column(self,col):
        if col in self.int_columns:
            return True
        else:
            return False


    def is_str_column(self,col):
        if col in self.string_columns:
            return True
        else:
            return False
        

    def cast_time_columns(self,df):
        for col in self.time_columns:
            if col in df.columns:
                df[col] = pd.to_datetime(df[col])
            
        return df
        

    def cast_string_columns(self,df):
        for col in self.string_columns:
            if col in df.columns:
                df[col] = df[col].fillna('').astype(str)
                df[col] = df[col].astype(str)
            
        return df


    def cast_float_columns(self,df):
        for col in self.float_columns:
            if col in df.columns:
                df[col] = df[col].astype(float)
            
        return df


    def read_list(self):
        """ Read in the SRAG list provided by Steve Johnson (SRAG) 
            converted from excel into text format and cleaned.
            
        """
        if not os.path.exists(self.filename):
            print(f"User_List.read_list: {self.filename} does not exist. Run create_list().")
            return pd.DataFrame()
            
        df = pd.read_csv(self.filename)
        
        #Cast string columns
        df = self.cast_string_columns(df)

        #Cast time columns
        df = self.cast_time_columns(df)
        
        return df


    def calculate_sep_reference_times(self,df):
        """ Does nothing because already set up with the right information """
        return df


    def create_list(self, type="minimal"):
        """ Create the user list if it doesn't already exist.
            If the user list already exists and is of type minimal,
            specifying type 'all' will add the rest of list_output_columns 
            to the existing list.
        
            INPUTS:
            
                :type: (string) minimal or all. minimal creates a list with
                    a subset of columns for only flares, CMEs, and event source location.
                    all creates a list with all association columns in list_output_columns.
        
        """
        columns = ['SEP Location', 'First SEP Start Time', 'Last SEP End Time']
        
        if type == "minimal":
            add_columns = ['Flare Xray Start Time', 'Flare Xray Peak Time',
                'Flare Xray End Time', 'Flare Class', 'Flare Magnitude',
                'Flare Integrated Flux', 'Flare Duration', 'Flare Time To Peak',
                'Flare Catalog ID', 'GOES Xray Satellite', 'Active Region',
                'Event Source Latitude', 'Event Source Longitude',
                'CDAW CME First Look Time', 'CDAW CME Speed',
                'CDAW CME Width', 'CDAW CME Mean Position Angle',
                'DONKI CME Start Time', 'DONKI CME Speed', 'DONKI CME Half Width',
                'DONKI CME Lat', 'DONKI CME Lon', 'DONKI CME Time at 21.5',
                'DONKI CME Catalog ID','DONKI CME Measurement Technique', 'Comments']
            columns = columns + add_columns

        elif type == "all":
            columns = columns + list_output_columns

        if os.path.exists(self.filename):
            df = pd.read_csv(self.filename)
            
            #If converting from minimal to all, add columns
            for col in columns:
                if col not in df.columns:
                    df.insert(column=col, value=[None]*len(df))
        else:
            #Create empty dataframe with specified columns
            df = pd.DataFrame(columns=columns)

        df.to_csv(self.filename, index=False)
        return df


    def write_list(self, df):
        df.to_csv(self.filename, index=False)
        return


    def add_flare(self, df, index, flare, format='dict'):
        """ Add flare information to user list.
        
            INPUTS:
            
                :df: (pandas DataFrame) user list
                :index: (int) row index where flare information should be added
                :flare: (dict) either flare_info dictionary or CCMC json flare
                    trigger block format
                :format: (string) 'dict' for flare_info, 'json' for CCMC json 
                    flare trigger block
                    
            OUTPUTS:
            
                :df: (pandas DataFrame) with flare information added to index row
        
        """

        #Empty dictionary
        if not flare:
            return df

        #Intensity is required; if missing, likely a dictionary of all null
        if pd.isnull(flare['intensity']):
            return df

        if format == 'dict':
            associations = flare_info_to_associations(flare)
        elif format == 'json':
            associations = flare_ccmc_json_to_associations(flare)
            
        for col in flare_columns:
            df.loc[index,col] = associations[col]

        return df


    def add_donki_cme(self, df, index, cme):
        """ Add DONKI information to user list. File output from 
            fetch_cme.get_donki_cme with format dict.
        
            INPUTS:
            
                :df: (pandas DataFrame) user list
                :index: (int) row index where flare information should be added
                :cme: (dict) File output from fetch_cme.get_donki_cme with format=dict
                    
            OUTPUTS:
            
                :df: (pandas DataFrame) with DONKI CME information added to index row
        
        """

        #Empty dictionary
        if not cme:
            return df

        #if missing, likely dictionary contains all null
        if pd.isnull(cme['DONKI CME Speed']):
            return df

        associations = donki_cme_to_associations(cme)
            
        for col in donki_cme_columns:
            df.loc[index,col] = associations[col]

        return df


    def add_cme_ccmc_json(self, df, index, cme, catalog='DONKI'):
        """ Add CME information to user list. File output from 
            fetch_cme.manual_cme_ccmc_json in CCMC JSON CME trigger block.
            Use catalog field to determine DONKI or CDAW. Default to 
            DONKI if catalog field is blank.
        
            INPUTS:
            
                :df: (pandas DataFrame) user list
                :index: (int) row index where flare information should be added
                :cme: (dict) File output from fetch_cme.manual_cme_ccmc_json 
                    or fetch_cme.get_donki_cme with format=json in 
                    CCMC JSON CME trigger block format
                    
            OUTPUTS:
            
                :df: (pandas DataFrame) with DONKI CME information added to index row
        
        """

        #Empty dictionary
        if not cme:
            return df

        #if missing, likely dictionary contains all null
        if 'speed' not in cme.keys():
            return df
        if pd.isnull(cme['speed']):
            return df

        if 'catalog' in cme.keys():
            if not pd.isnull(cme['catalog']) and cme['catalog'] != '':
                catalog = cme['catalog']
        
        associations = cme_ccmc_json_to_associations(cme, catalog=catalog)
        
        cme_col = []
        if catalog == 'DONKI' or catalog == 'donki':
            cme_col = donki_cme_columns
        elif catalog == 'SOHO_CDAW' or catalog == 'soho_cdaw':
            cme_col = cdaw_cme_columns
        
        for col in cme_col:
            df.loc[index,col] = associations[col]

        return df

        
    def add_source_latitude(self, df, index, lat):
        """ Add event source latitude to user list. """
        
        if not pd.isnull(lat) and lat != '':
            df.loc[index, 'Event Source Latitude'] = float(lat)

        return df


    def add_source_longitude(self, df, index, lon):
        """ Add event source longitude to user list. """

        if not pd.isnull(lon) and lon != '':
            df.loc[index, 'Event Source Longitude'] = float(lon)

        return df


    def add_sep_location(self, df, index, location):
        """ Add location where SEP was observed to user list, i.e. Earth, Mars, STEREO-A """
        if not pd.isnull(location) and location != '':
            df.loc[index,'SEP Location'] = location
        
        return df


    def add_noaa_region(self, df, index, ar):
        """ Add NOAA active region number to user list. """
        
        if not pd.isnull(ar) and ar != '':
            df.loc[index, 'Active Region'] = int(ar)
            
        return df


    def add_sep_first_time(self, df, index, first_time):
        """ Add first start time of SEP event to user list. """
        
        first_time = dh.str_to_datetime(first_time)
        
        if not pd.isnull(first_time) and first_time != '':
            df.loc[index, 'First SEP Start Time'] = first_time
            
        return df


    def add_sep_last_time(self, df, index, last_time):
        """ Add last end time of SEP event to user list. """
        
        last_time = dh.str_to_datetime(last_time)
        
        if not pd.isnull(last_time) and last_time != '':
            df.loc[index, 'Last SEP End Time'] = last_time
            
        return df


    def delete_row(self, df, index):
        """ Delete the row in the user list located at index """

        df = df.drop(index=index)
        return df
