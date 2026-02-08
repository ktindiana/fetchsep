from . import flare
from . import cme
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



###########################################################################
################## GENERAL ASSOCIATION LIST FUNCTIONS #####################
###########################################################################


def manual_flare_calibration(df, index):
    """ Convert deprecated flare xray columns into fluxes that align better
        with science-quality Xray data by removing SWPC calibration factors.
    
    """

    flare_info = flare.flare_info_dict()
    flare_info['catalog'] = 'SWPC_MANUAL'

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
    fluence = flare.remove_swpc_calibration(fluence)
    flare_info['integrated_intensity'] = fluence
    
    peak_intensity = df['Flare Magnitude Deprecated'].loc[index]
    peak_intensity = flare.remove_swpc_calibration(peak_intensity)
    fl_class = flare.flare_class(peak_intensity)
    flare_info['intensity'] = peak_intensity
    flare_info['class'] = fl_class

    return flare_info



def add_flare_columns(df):
    """ Add columns for NOAA science flare values """

    columns = ['Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time',
                'Flare Class', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration',
                'Flare Xray Time To Peak', 'Flare Catalog ID', 'GOES Xray Satellite']
    for col in columns:
        if col not in df.columns:
            df.insert(column=col, value=[None]*len(df))

    return df



def update_all_flares(df):
    """ If SEP event list contains the original X-ray fluxes made available
        by SWPC in their various products and archived data, add new science-quality
        flare information in additional columns. If flare not found in NOAA archive,
        manually remove SWPC's 0.7 calibration factor.  
        
        The list must contain information in the columns:
            'Flare Xray Start Time Deprecated'
            'Flare Magnitude Deprecated'
            'GOES Xray Satellite' - optional, but best to specify

        If these columns are present, they will be carried over if the SWPC calibration 
        factor is manually removed:
            'Flare Xray Peak Time Deprecated'
            'Flare Xray End Time Deprecated'
            'Flare Duration Deprecated'
            'Flare Xray Time To Peak Deprecated'
            'Flare Integrated Flux Deprecated'
            
        New science-quality X-ray data was released for GOES-08 to GOES-15, which
        should be used instead for space weather understanding, model training,
        and validation.
    
        For each flare in the list, check if the science data is available.
        If not, manually convert the flare values by removing the SWPC calibration.
        
        INPUTS:
        
            :df: (pandas dataframe) SEP associations list
            :list_name: (string) srag, user
            
        OUTPUTS:
        
            :df: (pandas dataframe) updated with new flare information
        
    """
    #In order of preference; Choose GOES-R first, if available
    xray_sat = ['GOES-16', 'GOES-18', 'GOES-17', 'GOES-19', 'GOES-08', 'GOES-09',
                'GOES-10', 'GOES-12', 'GOES-15', 'GOES-14', 'GOES-11']

    req_cols = ['Flare Xray Start Time Deprecated', 'Flare Magnitude Deprecated']
    for col in req_cols:
        if col not in df.columns:
            print(f"update_all_flares: {col} column must exist and be populated. Returning.")
            return df

    #Look for times in these columns to identify a flare
    ref_time_columns = ['Flare Xray Peak Time Deprecated', 'Flare Xray Start Time Deprecated']

    #If columns for the flare science data aren't present, add them
    df = add_flare_columns(df)

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

        #Experiment not specified, e.g. Ian's list
        #Find a GOES satellite for the requested date
        if experiment == None:
            request_date = pd.NaT
            for ref_col in ref_time_columns:
                if pd.isnull(request_date):
                    request_date = row[ref_col]
            if pd.isnull(request_date):
                continue
                
            print(f"Flare without GOES satellite specified {request_date}")
            stdate = request_date - datetime.timedelta(hours=24)
            enddate = request_date + datetime.timedelta(hours=48)

            #-->Don't know which satellite was used, so check until find one covered
            #the desired date
            flare_info = {}
            for experiment in xray_sat:
                if not flare.goes_xray_satellite_covers_date(experiment, request_date):
                    continue
                
                flare.download_goes_xray_science_data(stdate, enddate, experiment)
                flare_info = flare.get_noaa_flare(request_date, experiment=experiment)
                if not pd.isnull(flare_info['catalog_id']):
                    #The times provided in the IGR list aren't always accurate enough to find the correct
                    #flare. Use the deprecated peak value provided to compare with the identified flare.
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
            #After goes_r_primary_date, SWPC provided GOES-R X-ray fluxes
            if not flare_info and request_date < goes_r_primary_date:
                #DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES
                print(f"update_all_flares: Manually removing SWPC calibration factor from value provided in list for flare {request_date}.")
                flare_info = manual_flare_calibration(df, index)

        #Old GOES-07 and previous
        #########################
        elif experiment in old_goes_sc:
            #NEED TO DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES UNTIL NOAA RELEASES SCIENCE DATA
            print(f"update_all_flares: {experiment} X-ray science data is not yet available from NOAA for {request_date}. Manually removing SWPC calibration factor per NOAA guidance.")
            flare_info = manual_flare_calibration(df, index)
            
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

            print(f"update_all_flares: Flare {request_date} requested for {experiment}.")
 
            stdate = request_date - datetime.timedelta(hours=24)
            enddate = request_date + datetime.timedelta(hours=48)
            
            #-->Find flare in experiment specified in Steve's list
            flare.download_goes_xray_science_data(stdate, enddate, experiment)
            flare_info = flare.get_noaa_flare(request_date, experiment=experiment)
            
            #-->If that particular satellite didn't have the needed data, search for flare
            if not flare.goes_xray_satellite_covers_date(experiment, request_date) or pd.isnull(flare_info['catalog_id']):
                print(f"update_all_flares: Could not find flare {request_date} for {experiment}. Searching other GOES.")
                for exper in xray_sat:
                    if not flare.goes_xray_satellite_covers_date(exper, request_date):
                        continue
                    
                    flare.download_goes_xray_science_data(stdate, enddate, experiment)
                    flare_info = flare.get_noaa_flare(request_date, experiment=exper)
                    if not pd.isnull(flare_info['catalog_id']):
                        print(f"update_all_flares: Found flare {request_date} requested for {exper}.")
                        break

            #-->If didn't find a measurement in the science data, then apply calibration manually
            if pd.isnull(flare_info['catalog_id']) and request_date < goes_r_primary_date:
                #DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES
                print(f"update_all_flares: Could not find flare {request_date}. Manually removing SWPC calibration factor from values provided.")
                flare_info = manual_flare_calibration(df, index)


        #Check for missing/no results
        #############################
        #Didn't find any flare information for any experiments
        if not flare_info: continue

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



def add_donki_cme_columns(df):
    """ Add columns for DONKI CME values """
    columns = ['DONKI CME Start Time', 'DONKI CME Speed', 'DONKI CME Half Width',
                'DONKI CME Lat', 'DONKI CME Lon', 'DONKI CME Time at 21.5', 'DONKI CME Feature',
                'DONKI CME Measurement Technique', 'DONKI CME Catalog ID']
    for col in columns:
        if col not in df.columns:
            df.insert(column=col, value=[None]*len(df))

    return df



def update_all_donki_cmes(df, minimum_speed=450, minimum_halfAngle=15,
            feature='LE'):
    """ Given a list with CDAW LASCO CME first look times stored in a column:
            'CDAW CME First Look Time'
            
        Find corresponding DONKI CMEs with preferred selections:
            minimum speed: default > 450 km/s
            minimum half angle: default > 15 km/s
            feature: Default leading edge (LE) measurements preferred over shock (SH)
        
    """

    req_cols = ['CDAW CME First Look Time']
    for col in req_cols:
        if col not in df.columns:
            print(f"update_all_donki_cmes: {col} column must exist and be populated. Returning.")
            return df

    #Add DONKI CME columns if not present
    df = add_donki_cme_columns(df)


    for index, row in df.iterrows():
        starttime = row['CDAW CME First Look Time']
        if pd.isnull(starttime):
            continue
        
        cme_info = cme.get_donki_cme(starttime, minimum_speed=minimum_speed,
            minimum_halfAngle=minimum_halfAngle, feature=feature, format='dict')
        
        if not cme_info:
            continue

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





#Columns output in final SEP list, e.g. for the CLEAR SEP Benchmark Dataset
list_output_columns = ['Cycle', 'EventType', 'Case', 'Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Class Deprecated', 'Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Flare Class', 'Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Xray Time To Peak', 'Flare Catalog ID', 'GOES Xray Satellite', 'Flare Opt', 'Active Region', 'AR Area', 'AR Spot Class', 'AR Mag Class', 'AR Carrington', 'Event Location From Center', 'Event Latitude', 'Event Longitude', 'Event Location Source', 'Event Location from Center 2', 'Event Latitude 2', 'Event Longitude 2', 'Event Location Source 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio Station', 'Radio DH Start Time', 'Radio DH End Time', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio DH Note', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME First Look Time', 'CDAW CME Speed', 'CDAW CME Width', 'CDAW CME Mean Position Angle', 'DONKI CME Start Time', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', 'DONKI CME Time at 21.5', 'DONKI CME Catalog ID', 'DONKI CME Measurement Technique', 'ESP_CME', 'ESP CDAW CME First Look Time', 'ESP CDAW CME LASCO Speed', 'ACE CME Passage Time', 'ACE Bz', 'Sudden Impulse Time', 'Sudden Impulse Amplitude', 'GLE Event Number', 'PRF', 'Comments']

list_time_columns = ['Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio DH Start Time', 'Radio DH End Time', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'CDAW CME First Look Time', 'DONKI CME Start Time', 'DONKI CME Time at 21.5', 'ESP Flare Xray Start Time Deprecated', 'ESP Flare Xray Peak Time Deprecated', 'ESP Flare Xray End Time Deprecated', 'ESP Radio m_TyII Start Time', 'ESP Radio m_TyII End Time', 'ESP CDAW CME First Look Time', 'ACE CME Passage Time', 'Sudden Impulse Time']

list_float_columns = ['Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Duration Deprecated', 'Flare Xray Time To Peak Deprecated','Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Xray Time To Peak', 'AR Area', 'AR Carrington', 'Event Location From Center', 'Event Latitude', 'Event Longitude', 'Event Location from Center 2', 'Event Latitude 2', 'Event Longitude 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME Speed', 'CDAW CME Mean Position Angle', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', 'ESP CDAW CME LASCO Speed', 'ACE Bz', 'Sudden Impulse Amplitude']
    #'CDAW CME Width' contains text

list_int_columns = ['Flare Catalog ID', 'Cycle', 'Active Region', 'GOES Xray Satellite']

list_string_columns =  ['EventType', 'Case', 'Flare Class Deprecated', 'Flare Class', 'Flare Opt', 'AR Spot Class', 'AR Mag Class', 'Event Location Source', 'Event Location Source 2', 'Radio Station', 'Radio DH Note', 'DONKI CME Feature', 'DONKI CME Catalog ID', 'DONKI CME Measurement Technique', 'ESP_CME', 'GLE Event Number', 'PRF','Comments']


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


def empty_associations_series():
    """ Create a pandas Series with appropriate null values 
        for the associations that will be output by FetchSEP, 
        e.g., used in the CLEAR benchmark dataset).
    
    """
    
    dict = {}
    for col in list_output_columns():
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

    
def identify_associations_in_list(startdate, list_name='srag'):
    """ Given a date range, identify an event in the list associated
        with that date range. 
    
        Use the reference time in a list to find the correct association
        with an observed SEP event. 

        
        INPUTS:
        
            :startdate: (datetime) start time of SEP event that want to
                find in list; try time of enhancement above background
                for >10 MeV or proton energies ~10 MeV or threshold crossing time
            :list_name: (string) list to pull from for associated flares, CMEs, etc
                "srag" is SRAG_SEP_List_R11_CLEARversion.csv
                "user" is the list the user maintains in fetchsep/reference/user_associations.csv
                
        OUTPUT:
        
            :associations: (pandas Series) Series containing the associations
                fields used for the CLEAR benchmark dataset
        
    """

    if list_name == 'srag':
        assoc_list = SRAGList()
    else:
        sys.exit("User list is not yet implemented. Exiting.")
        

    df = assoc_list.read_list()
    
    tolerance = datetime.timedelta(hours=6) #search horizon
    columns = assoc_list.reference_columns
    
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
                    print("identify_associations_in_list: Found more than one possible "
                        f"match for {startdate}.")
                    if min_diff < check_min_diff:
                        event_index = ix[0]
                        check_min_diff = min_diff

    null_assoc = assoc_list.empty_associations_series()
    null_protons = assoc_list.empty_protons_series()
    if event_index == -1:
        print("identify_associations_in_list: No match was found in {assoc_list.id} for {startdate}.")
        #Return series with columns with appropriate null values
        return null_assoc, null_protons
    else:
        output_col = list_output_columns()
        associations = df[output_col].iloc[event_index]
        proton_info = df[columns].iloc[event_index]
        print("identify_associations_in_list: Found association information for input date of "
            f"{startdate} with an entry in {assoc_list.id} for {proton_info}")
        return associations.to_dict(), proton_info.to_dict()

    return null_assoc, null_protons



def list_to_ccmc_cme(associations):
    """ Take dictionary containing one row of list with CME information. 
        Put the CME columns into the appropriate fields for
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
    donki_cme = cme.donki_cme_to_ccmc_json(associations)
    if donki_cme:
        all_cmes.append(donki_cme)

    
    #CDAW CME
    cdaw_cme = cme.cdaw_cme_to_ccmc_json(associations)
    if cdaw_cme:
        all_cmes.append(cdaw_cme)

    return all_cmes
    
    
def list_to_ccmc_flare(associations):
    """ Take dictionary containing one row of the SRAG list. 
        Put the Flare columns into the appropriate fields for
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
        del flare['peak_ratio']
        del flare['urls']
 
        all_flares.append(flare)
 
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
        del flare['peak_ratio']
        del flare['catalog_id']
        del flare['urls']
 
        all_flares.append(flare)
 
    return all_flares



###########################################################################
################### STEVE JOHNSON'S SRAG SEP LIST #########################
###########################################################################
class SRAGList:
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
        self.reference_columns = ["SEP Reference Time"]

        self.time_columns = ['InitiationDT', 'SEP Reference Time', 'Flare Xray Start Time Deprecated', 'Flare Xray Peak Time Deprecated', 'Flare Xray End Time Deprecated', 'Flare Xray Start Time', 'Flare Xray Peak Time', 'Flare Xray End Time', 'Radio m_TyII Start Time', 'Radio m_TyII End Time', 'Radio DH Start Time', 'Radio DH End Time', 'Radio TyIV Start Time', 'Radio TyIV End Time', 'CDAW CME First Look Time', 'DONKI CME Start Time', 'DONKI CME Time at 21.5', 'P10_FluxStart', 'P10_StartDT', 'P10_OnsetMax_DT', 'P10_PeakDT', 'P10_EndDT', 'P50_FluxStart', 'P50_StartDT', 'P50_OnsetPeakDT', 'P50_PeakDT', 'P50_EndDT', 'P100_FluxStart', 'P100_StartDT', 'Onset_P100_PkDT', 'Event_P100_PkDT', 'P100_EndDT','ESP Flare Xray Start Time Deprecated', 'ESP Flare Xray Peak Time Deprecated', 'ESP Flare Xray End Time Deprecated', 'ESP Radio m_TyII Start Time', 'ESP Radio m_TyII End Time', 'ESP CDAW CME First Look Time', 'ACE CME Passage Time', 'Sudden Impulse Time']

        self.float_columns = ['Flare Magnitude Deprecated', 'Flare Integrated Flux Deprecated', 'Flare Duration Deprecated', 'Flare Xray Time To Peak Deprecated','Flare Magnitude', 'Flare Integrated Flux', 'Flare Duration', 'Flare Xray Time To Peak', 'AR Area', 'AR Carrington', 'Event Location From Center', 'Event Latitude', 'Event Longitude', 'Event Location from Center 2', 'Event Latitude 2', 'Event Longitude 2', 'Radio Rbr245Max', 'Radio Rbr2695Max', 'Radio Rbr8800', 'Radio TyIII_Imp', 'Radio TyII Imp', 'Radio TyII Speed', 'Radio m_TyII Start Frequency', 'Radio m_TyII End Frequency', 'Radio DH Start Frequency', 'Radio DH End Frequency', 'Radio TyIV Imp', 'Radio TyIV Duration', 'CDAW CME Speed', 'CDAW CME Mean Position Angle', 'DONKI CME Speed', 'DONKI CME Half Width', 'DONKI CME Lat', 'DONKI CME Lon', 'ESP CDAW CME LASCO Speed', 'ACE Bz', 'Sudden Impulse Amplitude']
            #'CDAW CME Width' contains text

        self.int_columns = ['Flare Catalog ID', 'Cycle', 'Active Region', 'GOES Xray Satellite', 'GOES Protons Satellite']

        self.string_columns =  ['EventType', 'Case', 'Flare Class Deprecated', 'Flare Class', 'Flare Opt', 'AR Spot Class', 'AR Mag Class', 'Event Location Source', 'Event Location Source 2', 'Radio Station', 'Radio DH Note', 'DONKI CME Feature', 'DONKI CME Catalog ID', 'DONKI CME Measurement Technique', 'ESP_CME', 'GLE Event Number', 'PRF','Comments']

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


    def empty_protons_series(self):
        dict= {"P10_FluxStart": pd.NaT, "P10_StartDT": pd.NaT}
        return dict
        

    def combine_comments(self, df):
        df["Comments1"] = df["Comments1"].where(pd.notna(df["Comments1"]), '')
        df["Comments2"] = df["Comments2"].where(pd.notna(df["Comments2"]), '')
        df["Comments"] = df["Comments1"] + df["Comments2"]
        df.drop("Comments1", axis=1, inplace=True)
        df.drop("Comments2", axis=1, inplace=True)
        return df


    def read_list(self):
        """ Read in the SRAG list provided by Steve Johnson (SRAG) 
            converted from excel into text format and cleaned.
            
        """
        srag_list = 'fetchsep/reference/SRAG_SEP_List_R11_CLEARversion.csv'
        df = pd.read_csv(srag_list)
        
        #Combine two separate comments columns into one
        if 'Comments1' in df.columns:
            df = self.combine_comments(df)

        #Cast string columns
        df = self.cast_string_columns(df)

        #Cast time columns
        df = self.cast_time_columns(df)
        
        return df


    def update_flares_and_cmes(self):
        """ Update SRAG SEP event list with flare X-ray science data and
            more parameters for DONKI CMEs. 
        
            INPUTS:
                
                None
                
            OUTPUTS:
            
                :df: (pandas DataFrame) updated dataframe of SRAG list
                Write out SRAG list file to fetchsep/references/
        
        """
        df = self.read_list()
        
        #Add columns to the SRAG list in a specific place
        #Add a column for Flare ID at end of flare info
        if 'Flare Catalog ID' not in df.columns:
            idx = df.columns.get_loc('Flare Xray Time To Peak')
            df.insert(loc=idx + 1, column='Flare Catalog ID', value=[None]*len(df['Flare Xray Time To Peak']))

        if 'DONKI CME Feature' not in df.columns:
            #Add a column for DONKI CME Feature
            idx = df.columns.get_loc('DONKI CME Time at 21.5')
            df.insert(loc=idx + 1, column='DONKI CME Feature', value=[None]*len(df['DONKI CME Time at 21.5']))

        if 'DONKI CME Measurement Technique' not in df.columns:
            #Add a column for DONKI CME Measurement Technique
            idx = df.columns.get_loc('DONKI CME Feature')
            df.insert(loc=idx + 1, column='DONKI CME Measurement Technique', value=[None]*len(df['DONKI CME Feature']))

        df = update_all_flares(df)
        df = update_all_donki_cmes(df, minimum_speed=450, minimum_halfAngle=15, feature='LE')

        outfname = os.path.join('fetchsep','reference','SRAG_SEP_List_R11_CLEARversion_UPDATED.csv')
        df.to_csv(outfname, index=False)
        print(f"update_flares_and_cmes: Wrote updated SEP list to {outfname}")


