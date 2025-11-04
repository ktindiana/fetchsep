from . import config as cfg
import pandas as pd
import os
import sys
import string
import datetime
import wget
import urllib.request
import numpy as np
import netCDF4 as nc
import cftime

#GOES spacecraft categories needed for GOES X-ray data
#Spacecraft in the GOES-R+ series
goes_R = ["GOES-16", "GOES-17", "GOES-18", "GOES-19"]
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
    
    tolerance = datetime.timedelta(hours=6)
    columns = ["P10_FluxStart"] #["P10_FluxStart", "InitiationDT", "P10_StartDT"]
    
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


##########################################################
#################### X-RAY DATA ##########################
##########################################################
def goes_xray_sat_info():
    """ Info to create filenames for NOAA's naming scheme for X-ray files.
    
    """

    #Selecting the science-quality X-ray fluxes updated by NOAA in 2021+
    #startdate and enddate are inclusive here
    sat_info = {"GOES-08": {"url": 'goes08',
                        "avg1m": 'sci_xrsf-l2-avg1m_g08_',
                        "flsum": 'sci_xrsf-l2-flsum_g08_',
                        "version": '_v1-0-0',
                        "startdate": datetime.datetime(1995,1,1),
                        "enddate": datetime.datetime(2003,6,16)},
              
              "GOES-09": {"url": 'goes09',
                        "avg1m": 'sci_xrsf-l2-avg1m_g09_',
                        "flsum": 'sci_xrsf-l2-flsum_g09_',
                        "version": '_v1-0-0',
                        "startdate": datetime.datetime(1996,4,1),
                        "enddate": datetime.datetime(1998,7,28)},
              
              "GOES-10":  {"url": 'goes10',
                        "avg1m": 'sci_xrsf-l2-avg1m_g10_',
                        "flsum": 'sci_xrsf-l2-flsum_g10_',
                        "version": '_v1-0-0',
                        "startdate": datetime.datetime(1998,7,1),
                        "enddate": datetime.datetime(2009,12,1)},
              
              "GOES-11":  {"url": 'goes11',
                        "avg1m": 'sci_xrsf-l2-avg1m_g11_',
                        "flsum": 'sci_xrsf-l2-flsum_g11_',
                        "version": '_v1-0-0',
                        "startdate": datetime.datetime(2006,6,1),
                        "enddate": datetime.datetime(2008,2,10)},
              
              "GOES-12":  {"url": 'goes12',
                        "avg1m": 'sci_xrsf-l2-avg1m_g12_',
                        "flsum": 'sci_xrsf-l2-flsum_g12_',
                        "version": '_v1-0-0',
                        "startdate": datetime.datetime(2003,1,10),
                        "enddate": datetime.datetime(2007,4,12)},

              "GOES-13":  {"url": 'goes13',
                        "avg1m": 'sci_xrsf-l2-avg1m_g13_',
                        "flsum": 'sci_xrsf-l2-flsum_g13_',
                        "version": '_v2-2-1',
                        "startdate": datetime.datetime(2013,6,7),
                        "enddate": datetime.datetime(2017,12,14)},
               
              "GOES-14":  {"url": 'goes14',
                        "avg1m": 'sci_xrsf-l2-avg1m_g14_',
                        "flsum": 'sci_xrsf-l2-flsum_g14_',
                        "version": '_v2-2-1',
                        "startdate": datetime.datetime(2009,9,19),
                        "enddate": datetime.datetime(2020,3,4)},
               
              "GOES-15":  {"url": 'goes15',
                        "avg1m": 'sci_xrsf-l2-avg1m_g15_',
                        "flsum": 'sci_xrsf-l2-flsum_g15_',
                        "version": '_v2-2-1',
                        "startdate": datetime.datetime(2010,4,7),
                        "enddate": datetime.datetime(2020,3,4)},

               "GOES-16":   {"url": 'goes16',
                        "avg1m": 'sci_xrsf-l2-avg1m_g16_',
                        "flsum": 'sci_xrsf-l2-flsum_g16_',
                        "version": '_v2-2-0',
                        "startdate": datetime.datetime(2017,2,7),
                        "enddate": datetime.datetime(2025,4,6)},

               "GOES-17":   {"url": 'goes17',
                        "avg1m": 'sci_xrsf-l2-avg1m_g17_',
                        "flsum": 'sci_xrsf-l2-flsum_g17_',
                        "version": '_v2-2-0',
                        "startdate": datetime.datetime(2018,6,1),
                        "enddate": datetime.datetime(2024,8,11)},

               "GOES-18":   {"url": 'goes18',
                        "avg1m": 'sci_xrsf-l2-avg1m_g18_',
                        "flsum": 'sci_xrsf-l2-flsum_g18_',
                        "version": '_v2-2-0',
                        "startdate": datetime.datetime(2022,9,2),
                        "enddate": datetime.datetime.now() - datetime.timedelta(hours=24)} #update when GOES-18 ends
              }

    return sat_info


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
        print("{experiment} X-ray science data is not yet available from NOAA. Returning.")
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
        avgfile = f"{sat_info[experiment]['avg1m']}{date_suffix}{sat_info[experiment]['version']}.nc"
        flsumfile = f"{sat_info[experiment]['flsum']}{date_suffix}{sat_info[experiment]['version']}.nc"

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


def flare_info(experiment, request_date):
    """ Extract information for a specific flare on a specific datetime.
        Enter a date within 15 minutes of flare start or end of for any time
        between the flare start and end.
        
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

    flare_info = {} #Will contain info of associated flare, if found

    year = request_date.year
    month = request_date.month
    day = request_date.day
    date_suffix = 'd%i%02i%02i' % (year,month,day)

    #Filenames like:
    #sci_xrsf-l2-flsum_g08_d19950103_v1-0-0.nc
    flsumfile = f"{sat_info[experiment]['flsum']}{date_suffix}{sat_info[experiment]['version']}.nc"

    fullpath = os.path.join(cfg.datapath,'GOES',flsumfile)
    flsum_exists = os.path.isfile(fullpath)
    if not flsum_exists:
        print(f"flare_info: Cannot find file {fullpath}. Run check_goes_xray_science_data to download.")
        return {}

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
    for flid in flare_ids:
        #Compile all information per flare by extracting according to flare_id
        idx = [i for i in range(len(data.variables["flare_id"][:])) if data.variables["flare_id"][i]==flid]
        
        #Check if request_date is associated with this flare
        flst = None
        flend = None
        for ix in idx:
            if data.variables["status"][ix] == "EVENT_START":
                flst = flare_dates[ix]
            if data.variables["status"][ix] == "EVENT_END":
                flend = flare_dates[ix]

        if flst != None and flend != None:
            if (request_date >= flst-td) and (request_date <= flend+td):
                select_idx = idx
                print(f"Found flare {flst} to {flend}")

        #Just in case there are two flares close together, give precedence
        #to the one containing requested date
        if flst != None and flend != None:
            if (request_date >= flst) and (request_date <= flend):
                select_idx = idx
                print(f"Found flare {flst} to {flend}")

    #If no flare found
    if len(select_idx) == 0:
        print(f"Flare not found for requested datetime {request_date}. Returning.")
        return flare_info

    #Compile flare info
    flare_info = {  'catalog_id': np.nan, #flare_id
                    'catalog': 'NOAA_NCEI_science_data',
                    'start_time': None,
                    'peak_time': None,
                    'end_time': None,
                    'intensity': np.nan,
                    'class': None,
                    'instrument': experiment,
                    'integrated_intensity': np.nan, #start time to end time
                }

    #flare_id, status, xrsb_flux, flare_class, time, integrated_flux, flare_class
    for ix in select_idx:
        if data.variables["status"][ix] == "EVENT_START":
            flare_info['start_time'] = flare_dates[ix]
            flare_info['catalog_id'] = int(data.variables["flare_id"][ix])
            
        if data.variables["status"][ix] == "EVENT_PEAK":
            flare_info['peak_time'] = flare_dates[ix]
            flare_info['intensity'] = float(data.variables["xrsb_flux"][ix])
            flare_info['class'] = data.variables["flare_class"][ix]
            
        if data.variables["status"][ix] == "EVENT_END":
            flare_info['end_time'] = flare_dates[ix]
            flare_info['integrated_intensity'] = float(data.variables["integrated_flux"][ix])

    return flare_info
