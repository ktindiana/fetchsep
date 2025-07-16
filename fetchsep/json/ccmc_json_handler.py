from ..utils import config as cfg
from . import keys
import json
import calendar
import datetime
from datetime import timedelta
import copy
import zulu
from astropy import units as u
import os
import sys
import pandas as pd
#import git #GitPython package
#import process
#import re

__version__ = "2.0"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#2021-01-12 changes in v0.2: Added fluence to json file. Added threshold
#   crossing entry for end of event.
#2021-02-11 changes in 0.3: Changed logic so that if integral flux was
#   derived from differential flux, and the flux did not cross the
#   threshold applied to the integral channel, then all_clear_boolean = True.
#   Previously, the code checked if the energy bin existed and set
#   all_clear_boolean to None if not. This produced the wrong value in the
#   case where differential channels were converted to integral.
#2021-02-24, changes in v0.4: Added subroutines to read in json file and
#   convert zulu time to datetime.
#   Modified fill_json to account for change in threshold fluences array in
#   v2.5 of operational_sep_quantities.py.
#2021-05-18, changes in 0.5: Realized that there were some differences
#   between the json files produced here and the CCMC format. CCMC defines
#   the fluences field as an array allowing multiple fluences, so the fluences
#   field was changed to an array here.
#   CCMC also has "event_lengths" as an array that allows the specification
#   of multiple start and end times. Changed the field name from "event_length"
#   to "event_lengths" and made it an array. Only one entry to
#   "event_lengths" is created by this program.
#2021-07-20, changes in 1.0: Went up an integer in version number because
#   making major changes to exactly match CCMC JSON format for the SEP
#   Scoreboard. Version 1.0 of this file works with v3.0 of
#   operational_sep_quantities.py, which has been modified for the same
#   purpose.
#   Removing contacts field from json due to NASA privacy rules.
#   Changed zulu time to keep seconds.
#   If multiple thresholds applied to a single energy channel,
#   values will be saved in one energy channel block as array
#   entries (previously writing separate blocks for unique
#   energy-threshold combinations.
#   Removing any attempt to differentiate between values not
#   filled out because they are not supported by the model
#   or values not filled out because there was no forecast.
#   Will have to apply that kind of information downstream.
#   JSON fields that have not been filled in are removed
#   from the final output json file.
#2021-08-18, changes in 1.1: Added subroutines from
#   validation_json_handler.py in the validation project
#   to access entries in the json file in a variety of ways.
#   This code, coupled with keys.py, will allow users
#   to access any values in the json files produced by
#   operational_sep_quantities.py
#2021-08-30, changes in 1.2: fixed bug in identifying energy
#   bins in fill_json. In clean_json, will remove the block for
#   an energy channel if the max flux is stored as a negative
#   value.
#2021-09-16, changes in 1.3: support for json_type, introduced
#   in operational_sep_quantities.py v3.2. Allows user to indicate
#   if a user input file should be written to observation or model
#   json file in fill_json and clean_json.
#2021-09-20, changes in 1.4: Added threshold units to the all clear
#   fields in fill_json (I had previously forgotten them and the
#   fields were being left as empty strings).
#2022-06-07, changed in 1.5: Don't try to guess Spase ID. If the
#   user doesn't specify an ID and the field is empty, leave it as
#   an empty string.
#2022-11-10, changes in 1.6: The all_clear logic was changed in fill_json.
#   The code no longer checks the peak flux values to determine 
#   all_clear if a threshold isn't crossed. Specific energy channel
#   and threshold crossing combinations are now strictly enforced
#   for all_clear statues (>10 MeV, 10 pfu; >100 MeV, 1 pfu)
#   All other energy channels are allowed
#   to use any threshold crossing to determine all_clear status.
#2023-04-18, changes in 1.7: Fixed bug in fill_json that assigned
#   the wrong units to the fluence spectrum if the energy block
#   for an integral channel, but the input data consisted of
#   differential fluxes. Was wrongly assigning fluence units of
#   e.g. 1/cm^2 when should have been 1/(MeV*cm^2).
#   In the get_value_by_* subroutines, added code to extract some of the
#   trigger information.
#   Added set_json_value_by_index() which allows the user to change a value in
#   the json. In the context of this code, may be most useful to add
#   the trigger information and change the issue time.
#2023-7-06, changes in 1.8: cast max flux as float in fill_json
#   because found that models with 0 flux resulted in an error
#   when dumping json.
#2023-10-01, changes in 1.9: Added subroutines to transform
#   threshold and energy channels to string keys.
#2024-04-24, changes in 2.0: Made changes to accomodate new fields
#   from CCMC. "model": {"flux_type": } --> "source_info": {"native_flux_type":}
#   Modified clean_json to remove the options field if nothing added.
#   Changed event_lengths threshold field to threshold_start and added
#   threshold_end.

def about_ccmc_json_handler():
    """ ABOUT ccmc_json_handler.py
    
        This module places all derived information into json and supporting
        files that follow the format specified for CCMC's SEP Scoreboard.
        https://ccmc.gsfc.nasa.gov/challenges/sep.php#format
        
        Reads in templates and fills in SEP event quantities organized
        by one block for each energy channel. If multiple thresholds
        were applied to the same energy channel, the derived information
        will be stored in a single block.
        
        Template for observations: observations_template.json
        
        Template for model output: model_template.json
        
        ADDITIONALLY:
        
        ccmc_jason_handler.py combined with fetchsep/json/keys.py
        provides all of the information needed to read the JSON
        files in the CCMC SEP Scoreboard format or produced
        by operational_sep_quantities.py.
        
        The user provides a unique identifier as listed in keys.py
        or the json visual scheme document on the CCMC webpage:
        https://ccmc.gsfc.nasa.gov/challenges/sep.php#format
        The CCMC identifiers are followed as closely as possible, with some
        additions specific to the operational_sep_quantities.py code and
        this SEP model validation effort.
        
        In some cases, the user must also provide an index value or
        a threshold value to uniquely identify values that are stored
        in arrays in the json file, These include 'fluences', 'fluence_spectra',
        'events_lengths', 'threshold_crossings','probabilities'
        
        In this code, a key_chain will be requested from keys.py. The
        key_chain contains the unique location of the desired value, but
        in a list of strings and integers (if an index is needed).
        validation_jason_handler.py extracts the unique value by
        extracting each sub-dictionary in an iterative manner until
        the unique value is accessed.
        
        For example, if key_chain = ['event_lengths',1,'start_time'], the
        value will be extracted from the full json dictionary (dict) as a
        loop over the key_chain, pulling out a subdictionary each time until
        the unique value is found, effectively as follows:
        
        .. code-block:: python
        
            subdict1 = dict['sep_forecast_submission']['forecasts']['event_lengths']
            subdict2 = subdict1[1]
            desired_val = subdict2['start_time']


        ALL CLEAR        
        fill_json contains logic to determine the All Clear status (all_clear_boolean)
        for each energy block. For >10 MeV and >100 MeV, only specific thresholds are 
        allowed to determine the All Clear status: >10 MeV, 10 pfu and >100 MeV, 1 pfu.

        For all other energy channels, the All Clear status will be filled by the first
        threshold encountered by the code. e.g. if the user runs the code for 
        >30 MeV, 1 pfu and >30 MeV, 5 pfu (--Threshold "30,1;30,5"), the all_clear_boolean
        will reflect the >30 MeV, 1 pfu status. If the user flips the call when running OpSEP 
        (e.g. --Threshold "30,5;30,1"), then all_clear_boolean will reflect >30 MeV, 5 pfu.
             
        
    """

def make_ccmc_zulu_time(dt):
    """ Make a datetime string in the format YYYY-MM-DDTHH:MM:SSZ
        
        INPUTS:
        
        :dt: (datetime)
        
        OUTPUTS:
        
        :zuludate: (string) in the format YYYY-MM-DDTHH:MM:SSZ
    
    """
    if dt == '':
        return ''
    if dt == None:
        return None
    if dt is pd.NaT:
        return pd.NaT
    if dt == 0:
        return 0

    zdt = zulu.create(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    stzdt = str(zdt)
    stzdt = stzdt.split('+00:00')
    zuludate = stzdt[0] + "Z"
    return zuludate
 
 
def zulu_to_time(zt):
    """ Convert Zulu time to datetime
    
        INPUTS:
        
        :zt: (string) - date in the format "YYYY-MM-DDTHH:MM:SSZ"
        
        OUTPUTS:
        
        :dt: (datetime)
        
    """
    #Time e.g. 2014-01-08T05:05:00Z or 2012-07-12T22:25Z
    if zt == '':
        return ''
    if zt == None:
        return None
    if zt is pd.NaT:
        return pd.NaT
    if zt == 0:
        return 0
  
    strzt = zt.split('T')
    strzt[1] = strzt[1].strip('Z')
    n = strzt[1].split(':')
    stdt = strzt[0] + ' ' + strzt[1]

    if len(n) == 2:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M')
    if len(n) == 3:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M:%S')
    return dt


def find_energy_bin(lowedge, energy_bins):
    """ Identify the energy bin and return low and high edges.
    
        INPUTS:
        
        :lowedge: (float) - low edge of an energy bin
        :energy_bins: (float 2xn array) - energy bins
        
        OUTPUTS:
        
        :bin: (float 2x1 array) - single energy bin
    
    """
    bin = []
    for i in range(len(energy_bins)):
        if lowedge == energy_bins[i][0]:
            bin = energy_bins[i]

    if not bin:
        print("ccmc_json_handler: find_energy_bin could not identify "
                "requested bin.")

    return bin


def id_unique_energy_channels(energy_thresholds):
    """ There may be multiple energy channel - flux threshold
        combinations. Identify the number of unique energy
        channels, e.g.
        
            * >10 MeV, 10 pfu (operations) and
            * >10 MeV, 0.001 pfu (SEPMOD testing)
        
        INPUTS:
        
        :energy_thresholds: (float 1xn array)- array containing
            all the energy channels for which a threshold
            has been applied
        
        OUTPUTS:
        
        :len(unique): (integer) - number of unique energy channels
        :unique: (float 1xm array) - the m unique energy channels to
            which a threshold was applied
        
    """
    unique = []
    for energy in energy_thresholds:
        if energy not in unique:
            unique.append(energy)
    
    return len(unique), unique



###############FILL AND WRITE JSONS##############
def set_keys(json_type):
    """ Choose the appropriate keys for the top levels of the observation
        or model jsons. Refers to keys.py
        
        #KEYS FOR OBSERVATIONS
        obs_main = 'sep_observation_submission'
        obs_exp = 'observatory'
        obs_type = 'observations'
        obs_win = 'observation_window'

        #KEYS FOR MODELS
        model_main = 'sep_forecast_submission'
        model_exp = 'model'
        model_type = 'forecasts'
        model_win = 'prediction_window'
        
    """
    if json_type == "model":
        key = keys.model_main
        type_key = keys.model_type
        win_key = keys.model_win
        exp_key = keys.model_exp
    else:
        key = keys.obs_main
        type_key = keys.obs_type
        win_key = keys.obs_win
        exp_key = keys.obs_exp

    return key, type_key, win_key, exp_key


def forecast_json():
    """ header for forecast json """

    template = {"sep_forecast_submission": {"notes": [ { "note": "produced by https://github.com/ktindiana/fetchsep"} ],
       "model": { "short_name": "", "spase_id": ""},
       "source_info": {"native_flux_type": ""},
       "options": "",
       "issue_time": "",
       "mode": "historical",
       "forecasts": []
        }}
    print("ccmc_json_handler: forecast_json: Initilizing forecast json.")
    return template


def observation_json():
    """ header for observation json """
    
    template = {"sep_observation_submission": {"notes": [ { "note": "produced by https://github.com/ktindiana/fetchsep"} ],"observatory": { "short_name": "", "spase_id": ""}, "source_info": {"native_flux_type": ""},"options": "","issue_time": "","mode": "measurement","observations": [] }}
    print("ccmc_json_handler: observation_json: Initilizing observation json.")
    return template



def fill_json_header(json_type, issue_time, experiment,
    flux_type, options, spase_id, model_name=None, mode=None,
    spacecraft=None):
    """ Fill in top level header information in json """
    
    zissue = make_ccmc_zulu_time(issue_time)
    
    short_name = experiment
    if model_name:
        short_name = model_name
    
    #If mode not specified, guess
    if not mode or mode == '':
        if json_type == "observations":
            mode = 'measurement'
        else:
            mode = 'historical'

    template = {}
    if json_type == "observations":
        template = observation_json()
    elif json_type == "model":
        template = forecast_json()

    key, type_key, win_key, exp_key = set_keys(json_type)

    template[key][exp_key]['short_name'] = short_name
    if spacecraft and spacecraft != '':
        template[key][keys.obs_exp]['short_name'] = f"{short_name} {spacecraft}"
    template[key]['source_info']['native_flux_type'] = flux_type
    if spase_id and spase_id != "":
        template[key][exp_key]['spase_id'] = spase_id
    template[key]['mode'] = mode

    template[key]['options'] = options
    template[key]['issue_time'] = zissue

    return template


def new_block(json_type):
    """ dict containing observation or forecast block """
    block = {}
    if json_type == "observations":
        block = {"energy_channel": { "min": "", "max": "", "units": ""},
                "species": "",
                "location": "",
                "observation_window": { "start_time": "", "end_time": "" },
                "peak_intensity": { "intensity": "", "units": "", "time": ""},
                "peak_intensity_max": { "intensity": "", "units": "", "time": "" },
                "event_lengths":[ { "start_time": "",  "end_time": "", "threshold_start": "", "threshold_end": "", "threshold_units": ""  }],
                "fluences": [{"fluence": "", "units": ""}],
                "fluence_spectra": [{"start_time": "", "end_time": "","threshold_start":"", "threshold_end":"", "threshold_units":"",
                "fluence_units": "",
                "fluence_spectrum":[{"energy_min": "", "energy_max":"", "fluence": ""}]}],
                "threshold_crossings": [ { "crossing_time": "", "threshold": "", "threshold_units": "" } ],
                "all_clear": { "all_clear_boolean": "", "threshold": "", "threshold_units": ""},
                "sep_profile": ""}

    if json_type == "model":
        block = {"energy_channel": { "min": "", "max": "", "units": ""},
                "species": "",
                "location": "",
                "prediction_window": { "start_time": "", "end_time": "" },
                "peak_intensity": { "intensity": "", "units": "", "time": ""},
                "peak_intensity_max": { "intensity": "", "units": "", "time": "" },
                "event_lengths":[ { "start_time": "",  "end_time": "", "threshold_start": "", "threshold_end": "", "threshold_units": ""  }],
                "fluences": [{"fluence": "", "units": ""}],
                "fluence_spectra": [{"start_time": "", "end_time": "","threshold_start":"", "threshold_end":"", "threshold_units":"",
                "fluence_units": "",
                "fluence_spectrum":[{"energy_min": "", "energy_max":"", "fluence": ""}]}],
                "threshold_crossings": [ { "crossing_time": "", "threshold": "", "threshold_units": "" } ],
                "all_clear": { "all_clear_boolean": "", "threshold": "", "threshold_units": ""},
                "sep_profile": ""}


    return block


def fill_json_block(template, json_type, energy_channel, threshold_dict, startdate,
    enddate,
    sep_start_time, sep_end_time, onset_peak, onset_peak_time, max_flux, max_flux_time,
    flux_units, fluence, fluence_units, fluence_spectrum, fluence_spectrum_units,
    fluence_energy_bins, sep_profile, location="earth", species="proton"):
    """ Fill values in a single block for a single energy channel.
        A single block may contain multiple event definitions, which means 
        multiple Analyze objects.
        
        Checks if a block of the appropriate energy channel already
        exists in the json. If not, a block is added and values filled in.
        If so, values are appended.
        
        
        energy_channel = {'min': 10.0, 'max':-1, 'units': 'MeV'}
        threshold_dict = {'threshold': 10, 'units': 'pfu'}
        
    """
    key, type_key, win_key, exp_key = set_keys(json_type)

    #Process times into the correct format
    zst = make_ccmc_zulu_time(startdate)
    if pd.isnull(zst): zst = ""
    zend = make_ccmc_zulu_time(enddate)
    if pd.isnull(zend): zend = ""
    zodate = make_ccmc_zulu_time(onset_peak_time)
    if pd.isnull(zodate): zodate = ""
    zpdate = make_ccmc_zulu_time(max_flux_time)
    if pd.isnull(zpdate): zpdate = ""
    zct = make_ccmc_zulu_time(sep_start_time)
    if pd.isnull(zct): zct = ""
    zeet = make_ccmc_zulu_time(sep_end_time)
    if pd.isnull(zeet): zeet = ""

    ######## Evaluate All Clear ##########
    all_clear = True #Assume True and check if thresholds were crossed

    #Threshold WAS crossed
    if not pd.isnull(sep_start_time):
        #For >10 MeV and >100 MeV, only the 10 pfu and 1 pfu thresholds
        #are applied for all clear as they are operational definitions.
        #>10 MeV, 10 pfu - NOAA SWPC and SRAG operations
        #>100 MeV, 1 pfu - SRAG operations
        #For all other energy channels, the all clear will be True or
        #False with respect to the thresholds applied in the event definitions

        #>10 MeV, 10 pfu
        if energy_channel['min'] == 10. and energy_channel['max'] == -1:
            if threshold_dict['threshold'] == 10:
                all_clear = False #event!
        #>100 MeV, 1 pfu
        elif energy_channel['min'] == 100. and energy_channel['max'] == -1:
            if threshold_dict['threshold'] == 1:
                all_clear = False #event!
        else:
            #Allow other flux & threshold combinations for other energy channels
            all_clear = False

    ####### CREATE THE DICTS FOR THE INDIVIDUAL VALUES #######
    #All Clear
    all_clear_dict = {"all_clear_boolean": all_clear, "threshold": threshold_dict['threshold'],
        "threshold_units": flux_units}
    
    
    #Onset Peak Flux
    onset_dict = {"intensity":float(onset_peak),"time": zodate,
                    "units":flux_units}

    #Maximum Flux
    max_dict = {"intensity":float(max_flux),"time": zpdate,
                    "units":flux_units}

    #Start and End Times
    length_dict = { "start_time": zct,  "end_time": zeet,
                    "threshold_start": threshold_dict['threshold'],
                    "threshold_end": threshold_dict['threshold']*cfg.endfac,
                    "threshold_units": threshold_dict['threshold_units'] }

    #Fluence for only the specific energy channel
    #Fluence order same as event_lengths
    fluence_dict = {"fluence": fluence,
                    "units": fluence_units}

    #Fluence spectrum
    spectrum = []
    if len(fluence_spectrum) != 0:
        for kk in range(len(fluence_energy_bins)):
            flentry = {"energy_min": fluence_energy_bins[kk][0],
                        "energy_max": fluence_energy_bins[kk][1],
                        "fluence": fluence_spectrum[kk]}
            spectrum.append(flentry)
    
    fl_spec_dict = {"start_time": zct, "end_time": zeet,
           "threshold_start":threshold_dict['threshold'],
           "threshold_end":threshold_dict['threshold']*cfg.endfac,
           "threshold_units":threshold_dict['threshold_units'],
           "fluence_units": fluence_spectrum_units,
           "fluence_spectrum":spectrum }

    #Threshold crossings
    cross_dict = { "crossing_time": zct,
                    "threshold": threshold_dict['threshold'],
                    "threshold_units": threshold_dict['threshold_units']}


    #####CREATE (if necessary) AND FILL BLOCK ####
    ##### CREATE ######
    #Check if template contains the energy channel for this event definition
    blocks = template[key][type_key]

    #Creating a new block or adding to an existing one?
    need_new_block = True

    #Identify the index of the required block
    ix = -1
    for i, block in enumerate(blocks):
        if block['energy_channel']  == energy_channel:
            ix = i
            need_new_block = False

    #If no block found, create new one
    if ix == -1:
        n = len(blocks)
        ix = n
        template[key][type_key].append(new_block(json_type))
        template[key][type_key][n]['energy_channel'] = energy_channel


    #### FILL ######

    #If added a new block, initialize all the values
    if need_new_block:
        template[key][type_key][ix]['species'] = species
        template[key][type_key][ix]['location'] = location
        template[key][type_key][ix][win_key]['start_time'] = zst
        template[key][type_key][ix][win_key]['end_time'] = zend

        template[key][type_key][ix]['peak_intensity'].update(onset_dict)
        template[key][type_key][ix]['peak_intensity_max'].update(max_dict)
        template[key][type_key][ix]['sep_profile'] = sep_profile

        template[key][type_key][ix]['event_lengths'][0].update(length_dict)
        template[key][type_key][ix]['fluences'][0].update(fluence_dict)
        template[key][type_key][ix]['fluence_spectra'][0].update(fl_spec_dict)
        template[key][type_key][ix]['threshold_crossings'][0].update(cross_dict)
        template[key][type_key][ix]['all_clear'] = all_clear_dict


    else:
 
        template[key][type_key][ix]['event_lengths'].append(length_dict)
        template[key][type_key][ix]['fluences'].append(fluence_dict)
        template[key][type_key][ix]['fluence_spectra'].append(fl_spec_dict)
        template[key][type_key][ix]['threshold_crossings'].append(cross_dict)

        #All Clear
        #Fill in All Clear if not already filled in.
        #If there was already a forecast made for a different threshold,
        #will not replace it
        #If there are multiple thresholds applied to a block, only the first
        #one will determine the all_clear field.
        if pd.isnull(template[key][type_key][ix]['all_clear']['all_clear_boolean']) \
            or template[key][type_key][ix]['all_clear']['all_clear_boolean'] == "":
            template[key][type_key][ix]['all_clear']['all_clear_boolean'] \
                            = all_clear
            template[key][type_key][ix]['all_clear']['threshold'] \
                            = threshold_dict['threshold']
            template[key][type_key][ix]['all_clear']['threshold_units'] \
                            = threshold_dict['threshold_units']
    
        #Fill the values only if the ones there are null
        if pd.isnull(template[key][type_key][ix]['peak_intensity']['intensity'])\
            or template[key][type_key][ix]['peak_intensity']['intensity'] == "":
            template[key][type_key][ix]['peak_intensity'].update(onset_dict)

        if pd.isnull(template[key][type_key][ix]['peak_intensity_max']['intensity'])\
            or template[key][type_key][ix]['peak_intensity_max']['intensity'] == "":
            template[key][type_key][ix]['peak_intensity_max'].update(max_dict)
        
        if pd.isnull(template[key][type_key][ix]['sep_profile']) \
            or template[key][type_key][ix]['sep_profile'] == "":
            template[key][type_key][ix]['sep_profile'] = sep_profile


    return template



def clean_json(template, experiment, json_type):
    """ Remove any fields that didn't get filled in and
        were left as empty strings.
        
    """
    if experiment == "user" and json_type == "model":
        key = keys.model_main
        type_key = keys.model_type
        win_key = keys.model_win
        
        if template[key]['model']['spase_id'] == "":
            template[key]['model'].pop('spase_id', None)
       
                        
    else:
        key = keys.obs_main
        type_key = keys.obs_type
        win_key = keys.obs_win


        if template[key]['observatory']['spase_id'] == "":
            template[key]['observatory'].pop('spase_id', None)

    options = template[key]['options']
    if options == [""]:
        template[key].pop('options', None)
    elif len(options) == 0:
        template[key].pop('options', None)

    nent = len(template[key][type_key])
    for i in range(nent-1,-1,-1):
        #Onset Peak Flux
        if template[key][type_key][i]['peak_intensity']['intensity'] == "":
            template[key][type_key][i].pop('peak_intensity', None) #remove from dict
        elif pd.isnull(template[key][type_key][i]['peak_intensity']['intensity']):
            template[key][type_key][i].pop('peak_intensity', None) #remove from dict
        
        #Maximum Flux
        if template[key][type_key][i]['peak_intensity_max']['intensity'] == "":
            template[key][type_key][i].pop('peak_intensity_max', None)
        elif pd.isnull(template[key][type_key][i]['peak_intensity_max']['intensity']):
            template[key][type_key][i].pop('peak_intensity_max', None)
                    
        #Maximum flux is the one field that will always be filled in,
        #regardless of whether a threshold is crossed. If max peak is
        #set to a negative number, then remove entire entry for the
        #energy channel, because it means the energy bin didn't exist
        #in the data set.
        if template[key][type_key][i]['peak_intensity_max']['intensity'] < 0:
            template[key][type_key].pop(i)
            continue
            
        #Start and End Times, Fluence
        nev = len(template[key][type_key][i]['event_lengths'])
        for j in range(nev-1,-1,-1):
            if template[key][type_key][i]['event_lengths'][j]['start_time']\
                == "":
                template[key][type_key][i]['event_lengths'].pop(j)
                template[key][type_key][i]['fluences'].pop(j) #corresponding fluences

            elif pd.isnull(template[key][type_key][i]['event_lengths'][j]['start_time']):
                template[key][type_key][i]['event_lengths'].pop(j)
                template[key][type_key][i]['fluences'].pop(j) #corresponding fluences
        
        if len(template[key][type_key][i]['event_lengths']) == 0:
            template[key][type_key][i].pop('event_lengths', None)
        if len(template[key][type_key][i]['fluences']) == 0:
            template[key][type_key][i].pop('fluences', None)
        
        
        #Fluence spectrum
        nev = len(template[key][type_key][i]['fluence_spectra'])
        for j in range(nev-1,-1,-1):
            if template[key][type_key][i]['fluence_spectra'][j]['start_time']\
                == "":
                template[key][type_key][i]['fluence_spectra'].pop(j)
        if len(template[key][type_key][i]['fluence_spectra']) == 0:
            template[key][type_key][i].pop('fluence_spectra', None)
        
        
        #Threshold crossings
        nev = len(template[key][type_key][i]['threshold_crossings'])
        for j in range(nev-1,-1,-1):
            if template[key][type_key][i]['threshold_crossings'][j]['crossing_time'] == "":
                template[key][type_key][i]['threshold_crossings'].pop(j)
            elif pd.isnull(template[key][type_key][i]['threshold_crossings'][j]['crossing_time']):
                template[key][type_key][i]['threshold_crossings'].pop(j)
        if len(template[key][type_key][i]['threshold_crossings']) == 0:
            template[key][type_key][i].pop('threshold_crossings', None)
            
        #All Clear
        if template[key][type_key][i]['all_clear']['all_clear_boolean'] == "":
            template[key][type_key][i].pop('all_clear', None)
        
        #SEP Flux Time Profile
        if template[key][type_key][i]['sep_profile'] == "":
            template[key][type_key][i].pop('sep_profile', None)
            
    return template


def write_json(template, filename):
    """Write json template to json file. """

    with open(filename, "w") as outfile:
        json.dump(template, outfile)

    if not os.path.isfile(filename):
        print("WARNING: ccmc_json_handler: write_json could not write your " \
            "file "+ str(filename))
        return False

    print("Wrote SEP values to json file --> " + filename)
    return True


###########READ JSONS####################

def read_in_json(filename):
    """Read in json file """
    if not os.path.isfile(filename):
        print("ccmc_json_handler: could not read in file " \
                + filename + ". Exiting.")
        return None
        
    with open(filename) as f:
        info=json.load(f)

    return info

def return_main_key(injson):
    """ Return the highest level key in the json file, typically
        'sep_forecast_submission' or 'sep_observation_submission'.
        Possible values assigned in keys.py.
        
        INPUTS:
        
        :injson: (dictionary) - a complete json file for a model or
            observation
            
        OUTPUTS:
        
        :main_key: (string) - the highest level key in the json file
        
    """
    #check if json file is empty
    if not injson:
        print("return_type_key: JSON file is empty.")
        return cfg.errval
    
    key_list = list(injson.keys())
    main_key = key_list[0] #only one at top level
    return main_key


def return_type_key(injson):
    """ Return the key associated with type (model or observations)
        where the forecasts or observed information are stored in
        the json file, typically 'forecasts' or 'observations'.
        Possible values assigned in keys.py.
        
        INPUTS:
        
        :injson: (dictionary) - a complete json file for a model or
            observation
            
        OUTPUTS:
        
        :type_key: (string) - a second level key in the json file which
            should match either the obs_type or model_type keys in keys.py
        
    """
    main_key = return_main_key(injson)
    
    key_list = injson[main_key].keys()
    if keys.model_type in key_list:
        type_key = keys.model_type
    elif keys.obs_type in key_list:
        type_key = keys.obs_type
    else:
        print("return_type_key: Could not identify type key in "
            "set of keys. Should be " + keys.model_type + " or "
            + keys.obs_type + ". " + "keys: " + str(key_list))
        return cfg.errval
    
    return type_key
    

def return_nforecasts(injson):
    """ Return the number of forecasts or observations
        in the file. e.g.
        
        len(json['sep_forecast_submission']['forecasts'])
        
        INPUTS:
        
        :injson: (dictionary) - a complete json file for a model or
            observation
            
        OUTPUTS:
        
        :n: (integer) - the number of forecast or observation blocks in
            the json file
            
    """
    main_key = return_main_key(injson)
    type_key = return_type_key(injson)
    if main_key == cfg.errval or type_key == cfg.errval:
        return cfg.errval
    
    arr = injson[main_key][type_key]
    n = len(arr)
    
    return n


def switch_model_to_obs_keys(key_chain):
    """ keys.py provides the location of each entry in the
        json files, but uses default identifiers that correspond
        to the ones used by CCMC for model forecasts.
        
        operational_sep_quantities.py outputs similar jsons for
        measurements, but since they are observations and not
        forecasts, some of the fields have been relabeled.
        
        If looking at an observations json, fix the keys
        to reflect the correct names of the fields.
        
        INPUTS:
        
        :key_chain: (array, list) list of strings and integers
        
        OUTPUTS:
        
        :key_chain: (array, list) list of strings and integers
        
    """
    for k in range(len(key_chain)):
        if key_chain[k] == keys.model_exp:
            key_chain[k] = keys.obs_exp
        if key_chain[k] == keys.model_type:
            key_chain[k] = keys.obs_type
        if key_chain[k] == keys.model_win:
            key_chain[k] = keys.obs_win
            
    return key_chain



def return_json_value_by_energy(injson, value, energy_channel={}, index=0):
    """ Return the value of a specific entry in the json file listed in
        the fields under 'observations' or 'forecasts'.
        
        Identify the desired block under 'observations' or 'forecasts'
        by matching the 'energy_channel' entry. The array index
        for the desired block will be discovered.
        
        matches input energy_channel dictionary to either
        injson['sep_forecast_submission']['forecasts'][i]['energy_channel']
        injson['sep_observation_submission']['observations'][i]['energy_channel']
        
        INPUTS:
    
        :injson: (dictionary) is the dictionary created by reading in a CCMC
            SEP Scoreboard JSON file for a single SEP event or
            forecast period
            
        :value: (string) indicates value desired and may be any of the
            unique identifiers in keys.py
        
        :energy_channel: (dictionary, optional) defines which energy channel
            from which the value is desired. All observations and forecasts
            are organized by energy channel. e.g.
            {'min': 10, 'max': 10, 'units': 'MeV'}
            
        :index: (int, optional) - for json elements that may be arrays, index
            indicates which array element to choose; if not specified, will
            return the 0th entry in the array
          
          
        OUTPUTS:
        
        :sub: (varies) sub is the final unique value found be equating the
            unique identier "value" with a key_chain and then extracting
            that specific value from the json dictionary; returns None
            if the value cannot be found or the field is empty
        
        NOTE: values that are strings in zulu time will be converted to
            datetime prior to returning
            
    """
    
    #Check that the json file isn't empty
    if not injson:
        print("return_json_value: JSON file is empty.")
        return cfg.errval
        
    #Get the sequence of keys for the desired value
    key_chain = keys.get_key_chain(value, index)
    if key_chain == cfg.errval:
        return cfg.errval
    
    #discover main key, e.g. 'sep_forecast_submission' or
    #'sep_observation_submission'
    main_key = return_main_key(injson)
    if main_key == cfg.errval: return cfg.errval
    
    #DEFAULT KEY CHAINS HAVE MODEL IDENTIFIERS IN THEM. IF LOOKING AT
    #OBSERVATION JSON, REPLACE APPROPRIATE FIELDS IN key_chain
    if main_key == keys.obs_main:
        key_chain = switch_model_to_obs_keys(key_chain)

    
    ##############TOP LEVEL INFO#################
    #Check if values saved at the top level, i.e. contacts,
    #model short name, etc.
    if key_chain[0] in injson[main_key].keys():
        sub = injson[main_key][key_chain[0]]
        
        if key_chain[0] == "triggers": #array
            found = False
            for entry in sub:
                if key_chain[1] in entry:
                    found = True
                    entry = entry[key_chain[1]]
                    for key in key_chain[2:]:
                        if key not in entry.keys():
                            entry = vars.errval
                            break
                        if key in entry.keys():
                            entry = entry[key]
                            if isinstance(entry,list):
                                return entry
                    return entry
            
            if not found: return vars.errval
        
        for key in key_chain[1:]:
            if isinstance(key,int): #array index value
                if key >= len(sub): #check that index inside range
                    sub = cfg.errval
                    break
                else:
                    sub=sub[key]
            else:
                if key not in sub.keys():
                    sub = cfg.errval
                    break
                if key in sub.keys():
                    sub = sub[key]
        
        if sub == cfg.errval:
            print("return_json_value: Keys for requested value " \
                + value + " not in json file: " + str(key_chain))
        else:
            #Check if value is a zulu time and convert to datetime
            for key in key_chain:
                if isinstance(key,int): continue
                if 'time' in key:
                    sub = zulu_to_time(sub)
        
        return sub #extracted down to a final value
    
    ##############UNDER FORECASTS OR OBSERVATIONS#################
    #Discover type key, i.e. 'forecasts' or 'observations'
    type_key = return_type_key(injson)
    if type_key == cfg.errval:
        return cfg.errval
    
    #Get the key for energy channel
    energy_chan_key = keys.get_key_chain(keys.id_energy_channel) #['energy_channel']
    
    #json[main_key][type_key] returns an array of forecasts or
    #observations, generally identified uniquely via energy channel.
    #Identify desired array by energy channel.
    Ndict = len(injson[main_key][type_key])
    sub = cfg.errval #subdictionaries that will be narrowed down to single value
    for i in range(Ndict):
        if energy_chan_key[0] not in injson[main_key][type_key][i]:
            continue
        else:
            #select entry with desired energy channel
            if energy_channel == injson[main_key][type_key][i][energy_chan_key[0]]:
                #Pull out subdictionary for specific energy channel
                sub_dict = injson[main_key][type_key][i]
                 
                #check if library contatining desired value is in sub_dict
                if key_chain[0] not in sub_dict.keys():
                    continue
                #Extract the value specified by the key_chain
                else:
                    if len(key_chain) == 1:
                        return sub_dict[key_chain[0]] #We have our value!
                    else: #if dictionary
                        sub = sub_dict[key_chain[0]]
                        for key in key_chain[1:]:
                            if isinstance(key,int): #array index value
                                if key >= len(sub): #check that index inside range
                                    sub = cfg.errval
                                    break
                                else:
                                    sub=sub[key]
                            else:
                                if key not in sub.keys():
                                    sub = cfg.errval
                                    break
                                if key in sub.keys():
                                    sub = sub[key]
    
    if sub == cfg.errval:
        print("return_json_value: Keys for requested value " \
            + value + " not in json file: " + str(key_chain))
    else:
        #Check if value is a zulu time and convert to datetime
        for key in key_chain:
            if isinstance(key,int): continue
            if 'time' in key:
                sub = zulu_to_time(sub)
    
    return sub #extracted down to a final value
                    
    


def return_json_value_by_index(injson, value, channel_index=0, index=0):
    """ Return the value of a specific entry in the json file listed in
        the fields under 'observations' or 'forecasts'. Select which block
        under under 'observations' or 'forecasts' to access by specifying
        channel_index.
        
        injson['forecasts'][channel_index]
        
        Then, if the desired value within that block is an array, e.g.
        fluences or event_lengths, specify which element in the array
        should be returned using index.
        
        injson['forecasts'][channel_index]['fluences'][index]
        
        INPUTS:
    
        :injson: (dictionary) is the dictionary created by reading in a CCMC
            SEP Scoreboard JSON file for a single SEP event or
            forecast period
            
        :value: (string) indicates value desired and may be any of the
            unique identifiers in keys.py
        
        :channel_index: (int) defines which entry in the forecast or
            observation array from which the value is desired. i.e.
            injson['sep_model_submission']['forecasts'] is an array and
            channel_index specifies which forecast to choose
            
        :index: (int, optional) - for json elements UNDER 'forecasts' or
            'observations' that may be arrays, index indicates which
            array element to choose, e.g. fluences is an array and
            index indicates which element of the array to choose;
            if not specified, will return the 0th entry in the array
            
        
        OUTPUTS:
        
        :sub: (varies) sub is the final unique value found be equating the
            unique identier "value" with a key_chain and then extracting
            that specific value from the json dictionary; returns None
            if the value cannot be found or the field is empty
        
        NOTE: values that are strings in zulu time will be converted to
            datetime prior to returning
        
    """
    
    #Check that the json file isn't empty
    if not injson:
        print("return_json_value: JSON file is empty.")
        return cfg.errval
        
    #Get the sequence of keys for the desired value
    key_chain = keys.get_key_chain(value, index)
    if key_chain == cfg.errval:
        return cfg.errval
        
    #discover main key, e.g. 'sep_forecast_submission' or
    #'sep_observation_submission'
    main_key = return_main_key(injson)
    if main_key == cfg.errval: return cfg.errval
    
    #DEFAULT KEY CHAINS HAVE MODEL IDENTIFIERS IN THEM. IF LOOKING AT
    #OBSERVATION JSON, REPLACE APPROPRIATE FIELDS IN key_chain
    if main_key == keys.obs_main:
        key_chain = switch_model_to_obs_keys(key_chain)
    
    ##############TOP LEVEL INFO#################
    #Check if values saved at the top level, i.e.
    #model short name, etc.
    if key_chain[0] in injson[main_key].keys():
        sub = injson[main_key][key_chain[0]]
        
        if key_chain[0] == "triggers": #array
            found = False
            for entry in sub:
                if key_chain[1] in entry:
                    found = True
                    entry = entry[key_chain[1]]
                    for key in key_chain[2:]:
                        if key not in entry.keys():
                            entry = vars.errval
                            break
                        if key in entry.keys():
                            entry = entry[key]
                            if isinstance(entry,list):
                                return entry
                    return entry
            
            if not found: return vars.errval
        
        for key in key_chain[1:]:
            if isinstance(key,int): #array index value
                if key >= len(sub): #check that index inside range
                    sub = cfg.errval
                    break
                else:
                    sub=sub[key]
            else:
                if key not in sub.keys():
                    sub = cfg.errval
                    break
                if key in sub.keys():
                    sub = sub[key]
        
        if sub == cfg.errval:
            print("return_json_value: Keys for requested value " \
                + value + " not in json file: " + str(key_chain))
        else:
            #Check if value is a zulu time and convert to datetime
            for key in key_chain:
                if isinstance(key,int): continue
                if 'time' in key:
                    sub = zulu_to_time(sub)
        
        return sub #extracted down to a final value
    
    ##############UNDER FORECASTS OR OBSERVATIONS#################
    #Discover type key, i.e. 'forecasts' or 'observations'
    type_key = return_type_key(injson)
    if type_key == cfg.errval:
        return cfg.errval
    
    #json[main_key][type_key] returns an array of forecasts or
    #observations, generally identified uniquely via energy channel
    #and threshold. Identify desired array by specified channel_index.
    #Pull out subdictionary for specific energy channel
    sub_dict = injson[main_key][type_key][channel_index]
    #check if library contatining desired value is in sub_dict
    if key_chain[0] not in sub_dict.keys():
        return cfg.errval
    #Extract the value specified by the key_chain
    else:
        if len(key_chain) == 1:
            return sub_dict[key_chain[0]] #We have our value!
        else: #if dictionary
            sub = sub_dict[key_chain[0]]
            for key in key_chain[1:]:
                if isinstance(key,int): #array index value
                    if key >= len(sub): #check that index inside range
                        sub = cfg.errval
                        break
                    else:
                        sub=sub[key]
                else:
                    if key not in sub.keys():
                        sub = cfg.errval
                        break
                    if key in sub.keys():
                        sub = sub[key]
    
    if sub == cfg.errval:
        print("return_json_value_by_index: Keys for requested value " \
            + value + " not in json file: " + str(key_chain))
    else:
        if isinstance(sub,int):
            sub = float(sub)
        #Check if value is a zulu time and convert to datetime
        for key in key_chain:
            if isinstance(key,int): continue
            if 'time' in key:
                sub = zulu_to_time(sub)
    
    return sub #extracted down to a final value


def return_json_value_by_threshold(injson, value, energy_channel={}, threshold=0):
    """ Return the value of a specific entry in the json file listed in
        the fields under 'observations' or 'forecasts'.
        
        Select the desired energy block by matching to energy_channel.
        If the desired value is inside of an array (e.g. event_lengths,
        fluences, fluence_spectra, etc), find the correct element by
        matching the applied threshold value.
        
        Match energy_channel and threshold to pull out a unique value.
        For example, to get an event start time for a 10 pfu threshold
        applied to the >10 MeV channel:
        
        * energy_channel = {"min": 10.0, "max": -1, "units": "MeV"}
        * threshold = 10
        
        injson['sep_forecast_submission']['forecasts'][i]['event_lengths'][j]['start_time']
        
        where i was found by matching energy_channel and j was found by matching
        
        injson['sep_forecast_submission']['forecasts'][i]['event_lengths'][j]['threshold']
        
        INPUTS:
    
        :injson: (dictionary) is the dictionary created by reading in a CCMC
            SEP Scoreboard JSON file for a single SEP event or
            forecast period
            
        :value: (string) indicates value desired and may be any of the
            unique identifiers in keys.py
        
        :energy_channel: (dictionary, optional) defines which energy channel
            from which the value is desired. All observations and forecasts
            are organized by energy channel. e.g.
            {'min': 10, 'max': 10, 'units': 'MeV'}
            
        :threshold: (float) - for json elements that may be arrays, identify
            the desired entry to selecting the one associated with a
            threshold value
          
          
        OUTPUTS:
        
        :sub: (varies) sub is the final unique value found be equating the
            unique identier "value" with a key_chain and then extracting
            that specific value from the json dictionary; returns None
            if the value cannot be found or the field is empty
        
        NOTE: values that are strings in zulu time will be converted to
            datetime prior to returning
            
    """
    
    #Check that the json file isn't empty
    if not injson:
        print("return_json_value: JSON file is empty.")
        return cfg.errval
        
    #Get the sequence of keys for the desired value
    key_chain = keys.get_key_chain(value)
    if key_chain == cfg.errval:
        return cfg.errval
        
    #If channel fluence is desired, need to take a different
    #approach and match to thresholds in the event_lengths
    #entry. CCMC requirement that those two arrays match up
    fluences_key_chain = []
    if 'fluences' in key_chain:
        fluences_key_chain = key_chain
        key_chain = keys.get_key_chain(keys.id_event_length_threshold)
    
    #discover main key, e.g. 'sep_forecast_submission' or
    #'sep_observation_submission'
    main_key = return_main_key(injson)
    if main_key == cfg.errval: return cfg.errval
    
    #DEFAULT KEY CHAINS HAVE MODEL IDENTIFIERS IN THEM. IF LOOKING AT
    #OBSERVATION JSON, REPLACE APPROPRIATE FIELDS IN key_chain
    if main_key == keys.obs_main:
        key_chain = switch_model_to_obs_keys(key_chain)

    
    ##############TOP LEVEL INFO#################
    #Check if values saved at the top level, i.e. contacts,
    #model short name, etc.
    if key_chain[0] in injson[main_key].keys():
        sub = injson[main_key][key_chain[0]]
        
        if key_chain[0] == "triggers": #array
            found = False
            for entry in sub:
                if key_chain[1] in entry:
                    found = True
                    entry = entry[key_chain[1]]
                    for key in key_chain[2:]:
                        if key not in entry.keys():
                            entry = vars.errval
                            break
                        if key in entry.keys():
                            entry = entry[key]
                            if isinstance(entry,list):
                                return entry
                    return entry
            
            if not found: return vars.errval
        
        for key in key_chain[1:]:
            if isinstance(key,int): #array index value
                if key >= len(sub): #check that index inside range
                    sub = cfg.errval
                    break
                else:
                    sub=sub[key]
            else:
                if key not in sub.keys():
                    sub = cfg.errval
                    break
                if key in sub.keys():
                    sub = sub[key]
        
        if sub == cfg.errval:
            print("return_json_value: Keys for requested value " \
                + value + " not in json file: " + str(key_chain))
        else:
            if isinstance(sub,int):
                sub = float(sub)
            #Check if value is a zulu time and convert to datetime
            for key in key_chain:
                if isinstance(key,int): continue
                if 'time' in key:
                    sub = zulu_to_time(sub)
        
        return sub #extracted down to a final value
    
    ##############UNDER FORECASTS OR OBSERVATIONS#################
    #Discover type key, i.e. 'forecasts' or 'observations'
    type_key = return_type_key(injson)
    if type_key == cfg.errval:
        return cfg.errval
    
    #Get the key for energy channel
    energy_chan_key = keys.get_key_chain(keys.id_energy_channel) #['energy_channel']
    
    #json[main_key][type_key] returns an array of forecasts or
    #observations, generally identified uniquely via energy channel
    #and threshold. Identify desired array by energy channel.
    Ndict = len(injson[main_key][type_key])
    sub = cfg.errval #subdictionaries that will be narrowed down to single value
    thresh_index = -1
    channel_index = -1
    for i in range(Ndict):
        if energy_chan_key[0] not in injson[main_key][type_key][i]:
            continue
        else:
            #select entry with desired energy channel
            if energy_channel == injson[main_key][type_key][i][energy_chan_key[0]]:
                #Pull out subdictionary for specific energy channel
                sub_dict = injson[main_key][type_key][i]
                channel_index = i
                
                #check if library contatining desired value is in sub_dict
                if key_chain[0] not in sub_dict.keys():
                    continue
                #Extract the value specified by the key_chain
                else:
                    if len(key_chain) == 1:
                        return sub_dict[key_chain[0]] #We have our value!
                    else: #if dictionary
                        sub = sub_dict[key_chain[0]]
                        #Get correct threshold id
                        id_threshold = ""
                        if 'fluences' in key_chain \
                            or 'event_lengths' in key_chain:
                            id_threshold = 'threshold'
                        if 'fluence_spectra' in key_chain:
                            id_threshold = 'threshold_start'
                        if 'threshold_crossings' in key_chain:
                            id_threshold = 'threshold'
                        if 'probabilities' in key_chain:
                            id_threshold = 'threshold'
                            
                        for key in key_chain[1:]:
                            if isinstance(key,int): #array index value
                                #In this case, the index stored in the key
                                #is a dummy value. We want to identify
                                #the correct element of the array by matching
                                #the threshold
                                #sub is an array
                                for j in range(len(sub)):
                                    if id_threshold not in sub[j]:
                                        continue
                                    thresh = sub[j][id_threshold]
                                    if thresh == threshold:
                                        thresh_index = j
                                
                                if thresh_index == -1: #check that index inside range
                                    sub = cfg.errval
                                    break
                                else:
                                    sub=sub[key]
                            else:
                                if key not in sub.keys():
                                    sub = cfg.errval
                                    break
                                if key in sub.keys():
                                    sub = sub[key]
    
    #IF user requested channel fluence, then the sub currently
    #holds the value for event start time, since the event_lengths
    #and fluences arrays must match up and the threshold value is
    #only stored in the event_lengths array.
    #Get the value for the fluence
    if fluences_key_chain != [] and channel_index != -1 \
        and thresh_index != -1:
        sub = return_json_value_by_index(injson, value, channel_index,
                thresh_index)
    
    if sub == cfg.errval:
        print("return_json_value: Keys for requested value " \
            + value + " not in json file: " + str(key_chain))
    else:
        if isinstance(sub,int):
            sub = float(sub)
        #Check if value is a zulu time and convert to datetime
        for key in key_chain:
            if isinstance(key,int): continue
            if 'time' in key:
                sub = zulu_to_time(sub)
    
    return sub #extracted down to a final value


def add_trigger_array(trigger_array, injson):
    """ Add triggers to the json. The triggers must be in an array
        with the correct values.
        
        There are many different types of triggers and they are stored
        in an array in the json file. e.g.
        "triggers": [
            {"cme": {"start_time": "2023-01-10T23:12:00Z", "lat": 10.0, "lon": -66.0, "half_width": 12.0, "speed": 711.0, "time_at_height": {"time": "2023-01-11T04:09Z", "height": 21.5}, "coordinates": "HEEQ", "catalog": "DONKI", "catalog_id": "2023-01-10T23:12:00-CME-001"}},
            {"cme_simulation": {"model": "WSA-ENLIL+Cone", "simulation_completion_time": "2023-01-12T13:22:27Z", "urls": ["https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/WSA-ENLIL/23267/1"]}}
            ]
        
        INPUTS:
        
        :trigger_array: (array of dictionaries) each element of the array is a
            dictionary for a specific trigger. See CCMC's json schema for
            allowed triggers.
            
        :injson: (json) file intended for output containing all SEP
            information.
            
        OUTPUTS:
        
        :injson: (json) file with trigger information added
        
    """
    main_key = return_main_key(injson)
    injson[main_key].update({"triggers": trigger_array})
    return injson
    
    

def set_json_value_by_index(value, injson, key_id, channel_index=0, index=0):
    """ !!!HERE THERE BE DRAGONS!!!!
        Please be careful using this subroutine. It will replace the values
        inside of the original injson dictionary that you feed into the
        subroutine. It takes advantage of python's way of pointing copies of
        objects to the same place in memory.
    
        Set the value of a specific entry specified by key_id in the json file
        listed in the fields under 'observations' or 'forecasts'. Select which
        block under under 'observations' or 'forecasts' to access by specifying
        channel_index.
        
        injson['forecasts'][channel_index]
        
        Then, if the desired value within that block is an array, e.g.
        fluences or event_lengths, specify which element in the array
        should be returned using index.
        
        injson['forecasts'][channel_index]['fluences'][index]
        
        INPUTS:
    
        :injson: (dictionary) is the dictionary created by reading in a CCMC
            SEP Scoreboard JSON file for a single SEP event or
            forecast period
            
        :key_id: (string) indicates key desired and may be any of the
            unique identifiers in keys.py (this is called value in the
            get_json_value_by_index subroutine)
            
        :value: (varies) the value that you want to put in the field at
            key_id.
        
        :channel_index: (int) defines which entry (json block) in the forecast
            or observation array from which the value is desired. i.e.
            injson['sep_model_submission']['forecasts'] is an array and
            channel_index specifies which forecast to choose
            
        :index: (int, optional) - for json elements UNDER 'forecasts' or
            'observations' that may be arrays, index indicates which
            array element to choose, e.g. fluences is an array and
            index indicates which element of the array to choose;
            if not specified, will return the 0th entry in the array
            
        
        OUTPUTS:
        
        :outjson: (json dictionary) the json file with the specified field
            updated to the new value
        
    """
    
    #Check that the json file isn't empty
    if not injson:
        print("set_json_value_by_index: JSON file is empty.")
        return vars.errval
        
    #Get the sequence of keys for the desired value
    key_chain = keys.get_key_chain(key_id, index)
    if key_chain == vars.errval:
        return vars.errval
        
    #discover main key, e.g. 'sep_forecast_submission' or
    #'sep_observation_submission'
    main_key = return_main_key(injson)
    if main_key == vars.errval: return vars.errval
    
    #DEFAULT KEY CHAINS HAVE MODEL IDENTIFIERS IN THEM. IF LOOKING AT
    #OBSERVATION JSON, REPLACE APPROPRIATE FIELDS IN key_chain
    if main_key == keys.obs_main:
        key_chain = switch_model_to_obs_keys(key_chain)
    
    ##############TOP LEVEL INFO#################
    #Check if values saved at the top level, i.e.
    #model short name, etc.
    if key_chain[0] in injson[main_key].keys():
        sub = injson
        sub = sub[main_key][key_chain[0]]
        for key in key_chain[1:-1]:
            if isinstance(key,int): #array index value
                if key >= len(sub): #check that index inside range
                    return vars.errval
                else:
                    sub=sub[key]
            else:
                if key not in sub.keys():
                    return vars.errval
                else:
                    sub = sub[key]
        
        #Use last key in list to set the value
        if isinstance(key_chain[-1],int): #array index value
            if key_chain[-1] >= len(sub): #check that index inside range
                return vars.errval
            else:
                sub[key_chain[-1]] = value
        else:
            if key_chain[-1] not in sub.keys():
                return vars.errval
            else:
                sub[key_chain[-1]] = value
        
        return injson #replaced value at key_chain
    
    
    ##############UNDER FORECASTS OR OBSERVATIONS#################
    #Discover type key, i.e. 'forecasts' or 'observations'
    type_key = return_type_key(injson)
    if type_key == vars.errval:
        return vars.errval
    
    #json[main_key][type_key] returns an array of forecasts or
    #observations, generally identified uniquely via energy channel
    #and threshold. Identify desired array by specified channel_index.
    #Pull out subdictionary for specific energy channel
    sub = injson
    sub = sub[main_key][type_key][channel_index]
    #check if library contatining desired value is in sub_dict
    if key_chain[0] not in sub.keys():
        return vars.errval
    #Extract the value specified by the key_chain
    else:
        if len(key_chain) == 1:
            sub[key_chain[0]] = value #We have our value!
            return injson
        else: #if dictionary
            sub = sub[key_chain[0]]
            for key in key_chain[1:-1]:
                if isinstance(key,int): #array index value
                    if key >= len(sub): #check that index inside range
                        return vars.errval
                    else:
                        sub=sub[key]
                else:
                    if key not in sub.keys():
                        return vars.errval
                    else:
                        sub = sub[key]
    
        #Use last key in list to set the value
        if isinstance(key_chain[-1],int): #array index value
            if key_chain[-1] >= len(sub): #check that index inside range
                return vars.errval
            else:
                sub[key_chain[-1]] = value
        else:
            if key_chain[-1] not in sub.keys():
                return vars.errval
            else:
                sub[key_chain[-1]] = value
        
        return injson #replaced value at key_chain
    
    
    return injson #if get here, nothing replaced


######################################################
#### Converting json fields to strings and vice versa
def convert_string_to_units(str_units):
    """ Take units written as a string following the CCMC SEP
        Scoreboard format and convert to astropy units.
        
        Expect, e.g.
        "MeV^-1*s^-1*cm^-2*sr^-1"
        "MeV"
        "cm^-2*sr^-1"
        "pfu"
        
    """
    
    str_units = str_units.replace("*",".")
    str_units = str_units.replace("^","")
    if str_units == "pfu":
        units = u.cm**-2*u.sr**-1*u.s**-1 #pfu
    else:
        units = u.Unit(str_units)
    
    return units


def convert_units_to_string(units):
    """ Convert astropy units object to a string.
    """
    return str(units)


def energy_channel_to_key(energy_channel):
    """ Want to organize observations and forecasts according
        to energy channel to reduce uneccesary elements in loops.
        
        Turn the energy channel into a string key that can
        be used to organize a dictionary.
        
    Inputs:
    
        :energy_channel: (dict)
            {'min': 10, 'max': -1, 'units': Unit("MeV")}
    
    Output:
    
        :key: (string)
    
    """

    units = energy_channel['units']
    if isinstance(units,str):
        units = convert_string_to_units(units)
        
    str_units = convert_units_to_string(units)

    key = "min." +str(float(energy_channel['min'])) + ".max." \
        + str(float(energy_channel['max'])) + ".units." \
        + str_units
    
    return key


def threshold_to_key(threshold):
    """ Want to organize observations and forecasts according
        to energy channel and thresholds.
        
        Turn the threshold into a string key that can
        be used to organize a dictionary.
        
    Inputs:
    
        :threshold: (dict)
            {'threshold': 10, 'threshold_units': Unit("1 / (cm2 s sr)")}
    
    Output:
    
        :key: (string) e.g. "threshold.10.0.units.1 / (cm2 s sr)"
    
    """

    units = threshold['threshold_units']
    if isinstance(units,str):
        units = convert_string_to_units(units)
        
    str_units = convert_units_to_string(units)

    key = "threshold." +str(float(threshold['threshold'])) \
        + ".units." + str_units
    
    return key

