from ..json import ccmc_json_handler as ccmc_json
from ..utils import date_handler as dh
import os
import sys
import datetime
import wget
import urllib.request
import requests
import pandas as pd
import inspect
import numpy as np

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



##########################################################
#################### DONKI CMEs ##########################
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


def get_donki_cmes(start_date, end_date, minimum_speed=None, minimum_halfAngle=None,
    feature='LE'):
    """ Extract the most relevant CME analysis for each CME.
        Preference for CME that is:
            Most Accurate
            Excludes plane-of-sky measurements
            Has SEPVAL in the comments
            Has the feature code specified by the user (usually LE)
        
        
        INPUTS:
        
            :feature: (string) LE or SH, will preferably select leading edge or
                shock measurement if choice is available. Default LE.
    
    """

    #Start/End date in YYYY-MM-DD format
    url = 'https://kauai.ccmc.gsfc.nasa.gov/DONKI/WS/get/CME?startDate=' + start_date + '&endDate=' + end_date
    try:
        response = requests.get(url, timeout=15)
        cmes_ = response.json()
    except:
        print(f"get_donki_cmes: could not make connection for {url}")
        return []

    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    signature_args = inspect.signature(get_donki_cmes).parameters.keys()
    arg_dict = {arg : values[arg] for arg in signature_args}
    
    filters = {}
    for arg, value in arg_dict.items():
        if arg in ['minimum_speed', 'minimum_halfAngle']:
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
                    #Prefer specified feature, if available
                    if cme_analysis['featureCode'] == feature or cme_analysis['featureCode'] is None:
                            best_cme = cme_analysis
                elif best_cme['levelOfData'] < cme_analysis['levelOfData']:
                    #Prefer specified feature, if available
                    if cme_analysis['featureCode'] == feature:
                        if best_cme['featureCode'] is not None:
                            best_cme = cme_analysis
                else:
                    #Prefer specified feature, if available
                    if cme_analysis['featureCode'] == feature:
                        if best_cme['featureCode'] is not None:
                            best_cme = cme_analysis


        if not best_cme:
            continue
            
        cme['cmeAnalysisPrime'] = best_cme

        reject = reject_entry(cme['cmeAnalysisPrime'], filters)
        if not reject:
            cmes.append(cme)
    return cmes


def get_donki_cme_start_time(cmes):
    cme_info = []
    for i, cme in enumerate(cmes):
        if cme['startTime'] is None:
            cme_info.append(None)
        else:
            time = convert_to_time(cme['startTime'])
            cme_info.append(time)
    return cme_info


def get_donki_cme_header_info(cmes, quantity):
    """ Information related to CME included in header """

    cme_info = [None]*len(cmes)
    for i, cme in enumerate(cmes):
        info = cme[quantity]
        if info is None:
            info = None
        cme_info[i] = info
    return cme_info


def get_donki_cme_info(cmes, quantity):
    """ Information from CME analysis. """
    
    cme_info = [None]*len(cmes)
    for i, cme in enumerate(cmes):
        info = cme['cmeAnalysisPrime'][quantity]
        if info is None:
            info = None
        cme_info[i] = info
    return cme_info


def donki_cme_to_ccmc_json(list_cme):
    """ Put DONKI CME information into the CCMC SEP Scoreboard json
        trigger block format.
        
        INPUTS:

            :cme: (dict) in the format of cme_info in get_donki_cme,
                i.e. uses keys that are the same as the columns in
                fetchsep association lists
                        
    """
    if not list_cme:
        return {}
        
    if not pd.isnull(list_cme['DONKI CME Speed']):
        cme = ccmc_json.ccmc_cme_block()

        cme['speed'] = list_cme['DONKI CME Speed']
        cme['catalog'] = 'DONKI'
        cme['derivation_technique']['process'] = 'manual'

        if not pd.isnull(list_cme['DONKI CME Measurement Technique']):
            cme['derivation_technique']['method'] = list_cme['DONKI CME Measurement Technique']
        else:
            del cme['derivation_technique']['method']

        if not pd.isnull(list_cme['DONKI CME Start Time']):
            cme['start_time'] = dh.time_to_zulu(list_cme['DONKI CME Start Time'])
        else:
            del cme['start_time']

        if not pd.isnull(list_cme['DONKI CME Catalog ID']):
            cme['catalog_id'] = list_cme['DONKI CME Catalog ID']
        else:
            del cme['catalog_id']
        
        if not pd.isnull(list_cme['DONKI CME Lat']):
            cme['lat'] = list_cme['DONKI CME Lat']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lat']

        if not pd.isnull(list_cme['DONKI CME Lon']):
            cme['lon'] = list_cme['DONKI CME Lon']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lon']
            del cme['coordinates']

        if not pd.isnull(list_cme['DONKI CME Half Width']):
            cme['half_width'] = list_cme['DONKI CME Half Width']
        else:
            del cme['half_width']

        if not pd.isnull(list_cme['DONKI CME Time at 21.5']):
            cme['time_at_height']['time'] = dh.time_to_zulu(list_cme['DONKI CME Time at 21.5'])
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

    return cme


def get_donki_cme(starttime, format='dict', feature='LE', minimum_speed=None,
    minimum_halfAngle=None):
    """ Get a specific DONKI CME identified by start time (first look time).
        
        INPUTS:
        
            :starttime: (datetime, string) timestamp close to the start time of the CME, 
                typically defined as the first time at which the CME appears on the
                coronagraph. For SOHO, will be the LASCO C2 timestamp.
            :format: (string) "dict" returns dictionary with column labels for the
                associations lists; "json" returns a CCMC JSON trigger block
            :feature: (string) LE or SH, will preferably select leading edge or
                shock measurement if choice is available. Default LE.
                
        OUTPUTS:
        
            :cme: (dict) CME info save in a dictionary
    
    """
    
    if isinstance(starttime, str):
        starttime=dh.str_to_datetime(starttime)
    
    startdate = starttime - datetime.timedelta(hours=6)
    enddate = starttime + datetime.timedelta(hours=6)
    
    start_date = startdate.strftime('%Y-%m-%dT%H:%M:%SZ')
    end_date = enddate.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    cmes = get_donki_cmes(start_date, end_date, minimum_speed=minimum_speed, minimum_halfAngle=minimum_halfAngle)

    #CME ANALYSES
    cme_start_times = get_donki_cme_start_time(cmes)
    cme_speeds = get_donki_cme_info(cmes, 'speed')
    cme_half_angles = get_donki_cme_info(cmes, 'halfAngle')
    cme_longitudes = get_donki_cme_info(cmes, 'longitude')
    cme_latitudes = get_donki_cme_info(cmes, 'latitude')
    cme_21_5 = get_donki_cme_info(cmes, 'time21_5')
    cme_techniques = get_donki_cme_info(cmes, 'measurementTechnique')
    cme_feature = get_donki_cme_info(cmes, 'featureCode')

    #CME HEADER INFO
    cme_ids = get_donki_cme_header_info(cmes,'activityID')
    cme_notes = get_donki_cme_header_info(cmes, 'note')
    cme_submission_times = get_donki_cme_header_info(cmes, 'submissionTime')
    cme_linked_events = get_donki_cme_header_info(cmes, 'linkedEvents')
    cme_links = get_donki_cme_header_info(cmes, 'link')

    #Find requested CME
    td = datetime.timedelta(hours=2, minutes=30) #At least within 2.5 hours
    ix = -1
    for i, time in enumerate(cme_start_times):
        diff = abs(time - starttime)
        if diff < td:
            td = diff
            ix = i
    
    if ix == -1:
        print(f"get_donki_cme: No DONKI CME found for {starttime}.")
        return {}

    cme_info = {"DONKI CME Start Time": cme_start_times[ix],
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
                  "DONKI CME Link": cme_links[ix],
                  "DONKI CME Measurement Technique": cme_techniques[ix]
                  }

    if format == 'json':
        cme_info = donki_cme_to_ccmc_json(cme_info)

    return cme_info


def get_donki_cme_range(start_date, end_date, minimum_speed=None, minimum_halfAngle=None, feature='LE', format='pd'):
    """ Get a specific DONKI CME identified by start time (first look time).
        Preference for CME that is:
            Most Accurate
            Excludes plane-of-sky measurements
            Has SEPVAL in the comments
            Has the feature code specified by the user (usually LE)
                    
        INPUTS:
        
            :start_date: (datetime or string) YYYY-MM-DD, part of DONKI API call
            :end_date: (datetime or string) YYYY-MMM-DD, part of DONKI API call
            :feature: (string) LE or SH, will preferably select leading edge or
                shock measurement if choice is available. Default LE.
            :format: (string) pd for pandas dataframe, dict for dictionary
                
        OUTPUTS:
        
            :cme_range: (pandas Dataframe or dictionary) All CMEs
    
    """
    if isinstance(start_date,datetime.date):
        start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    if isinstance(end_date, datetime.date):
        end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    cmes = get_donki_cmes(start_date, end_date, minimum_speed=minimum_speed, minimum_halfAngle=minimum_halfAngle)

    #CME ANALYSES
    cme_start_times = get_donki_cme_start_time(cmes)
    cme_speeds = get_donki_cme_info(cmes, 'speed')
    cme_half_angles = get_donki_cme_info(cmes, 'halfAngle')
    cme_longitudes = get_donki_cme_info(cmes, 'longitude')
    cme_latitudes = get_donki_cme_info(cmes, 'latitude')
    cme_21_5 = get_donki_cme_info(cmes, 'time21_5')
    cme_techniques = get_donki_cme_info(cmes, 'measurementTechnique')
    cme_feature = get_donki_cme_info(cmes, 'featureCode')

    #CME HEADER INFO
    cme_ids = get_donki_cme_header_info(cmes,'activityID')
    cme_notes = get_donki_cme_header_info(cmes, 'note')
    cme_submission_times = get_donki_cme_header_info(cmes, 'submissionTime')
    cme_linked_events = get_donki_cme_header_info(cmes, 'linkedEvents')
    cme_links = get_donki_cme_header_info(cmes, 'link')

    cme_speeds = [convert_to_float(val) for val in cme_speeds]
    cme_half_angles = [convert_to_float(val) for val in cme_half_angles]
    cme_latitudes  = [convert_to_float(val) for val in cme_latitudes]
    cme_longitudes  = [convert_to_float(val) for val in cme_longitudes]
    cme_21_5 = [convert_to_time(val) for val in cme_21_5]
    cme_submission_times = [convert_to_time(val) for val in cme_submission_times]

    cme_range = {"DONKI CME Start Time": cme_start_times,
                  "DONKI CME Speed": cme_speeds,
                  "DONKI CME Half Width": cme_half_angles,
                  "DONKI CME Lat": cme_latitudes,
                  "DONKI CME Lon": cme_longitudes,
                  "DONKI CME Time at 21.5": cme_21_5,
                  "DONKI CME Catalog ID": cme_ids,
                  "DONKI CME Submission Time": cme_submission_times,
                  "DONKI CME Feature": cme_feature,
                  "DONKI CME Note": cme_notes,
                  "DONKI CME Linked Events": cme_linked_events,
                  "DONKI CME Link": cme_links,
                  "DONKI CME Measurement Technique": cme_techniques
                  }

    if format == 'pd':
        cme_range = pd.DataFrame(cme_range)

    return cme_range


##########################################################
##################### CDAW CMEs ##########################
##########################################################
def cdaw_cme_to_ccmc_json(list_cme):
    """ Put CDAW CME information into the CCMC SEP Scoreboard json
        trigger block format.
        
        INPUTS:

            :list_cme: (dict) in the format of the columns in
                fetchsep association lists for CDAW CMEs
                        
    """
    if not list_cme:
        return {}
        
    if not pd.isnull(list_cme['CDAW CME Speed']):
        cme = ccmc_json.ccmc_cme_block()

        cme['speed'] = list_cme['CDAW CME Speed']
        cme['catalog'] = 'SOHO_CDAW'
        cme['derivation_technique']['process'] = 'manual'
        cme['derivation_technique']['method'] = 'Plane-of-sky'

        if not pd.isnull(list_cme['CDAW CME First Look Time']):
            cme['start_time'] = list_cme['CDAW CME First Look Time'].to_pydatetime()
        else:
            del cme['start_time']

        if not pd.isnull(list_cme['Event Source Latitude']):
            cme['lat'] = list_cme['Event Source Latitude']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lat']

        if not pd.isnull(list_cme['Event Source Longitude']):
            cme['lon'] = list_cme['Event Source Longitude']
            cme['coordinates'] = 'HEEQ'
        else:
            del cme['lon']
            del cme['coordinates']

        if not pd.isnull(list_cme['CDAW CME Mean Position Angle']):
            cme['pa'] = list_cme['CDAW CME Mean Position Angle']
        else:
            del cme['pa']
    
        width = list_cme['CDAW CME Width']
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
        del cme['catalog_id']
        del cme['urls']
        try:
            del cme['derivation_technique']['measurement_type']
        except:
            pass

    return cme


##########################################################
################## MANUAL CME INPUT ######################
##########################################################
def manual_cme_ccmc_json(cme_start_time = pd.NaT,
    cme_liftoff_time=pd.NaT,
    cme_half_width_deg=np.nan,
    cme_speed_kms=np.nan,
    cme_acceleration_kms2=np.nan,
    cme_height_rs=np.nan,
    cme_time_at_height_time=pd.NaT,
    cme_time_at_height_height=np.nan,
    cme_lat=np.nan,
    cme_lon=np.nan,
    cme_pa=np.nan,
    cme_coordinates=None,
    cme_catalog=None,
    cme_catalog_id=None,
    cme_urls=[],
    cme_derivation_process=None,
    cme_derivation_method=None,
    cme_measurement_type=None):
    """ Create a CME trigger block for the CCMC SEP Scoreboard json schema
        using manual inputs.
        
    """

    cme = {
        "start_time": dh.time_to_zulu(cme_start_time),
        "liftoff_time": dh.time_to_zulu(cme_liftoff_time),
        "lat": cme_lat,
        "lon": cme_lon,
        "pa": cme_pa,
        "half_width": cme_half_width_deg,
        "speed": cme_speed_kms,
        "acceleration": cme_acceleration_kms2,
        "height": cme_height_rs,
        "time_at_height": {
            "time": dh.time_to_zulu(cme_time_at_height_time),
            "height": cme_time_at_height_height
        },
        "coordinates": cme_coordinates,
        "catalog": cme_catalog,
        "catalog_id": cme_catalog_id,
        "urls": cme_urls,
        "derivation_technique": {
            "process": cme_derivation_process,
            "method": cme_derivation_method,
            "measurement_type": cme_measurement_type
        }
    }
    
    cme = ccmc_json.clean_trigger_block(cme)

    return cme
