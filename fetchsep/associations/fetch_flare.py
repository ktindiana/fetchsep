from ..utils import date_handler as dh
from ..utils import experiments as expts
from ..json import ccmc_json_handler as ccmc_json
from ..utils import config as cfg
from ..utils import experiments as expts
import os
import sys
import datetime
import wget
import urllib.request
import requests
import cftime
import numpy as np
import pandas as pd
import netCDF4 as nc

#GOES spacecraft categories needed for GOES X-ray data
#Spacecraft in the GOES-R+ series
goes_R = expts.goes_R()
#Spacecraft prior to GOES-R
goes_sc = expts.goes_sc()
#Spacecraft prior to GOES-08
old_goes_sc = expts.old_goes_sc()

goesX_dir = 'GOES_Xray'

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



##########################################################
################ NOAA X-RAY SCIENCE DATA #################
##########################################################
def goes_xray_satellite_info():
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
                        "avg_version": '_v2-2-1',
                        "flsum_version": '_v2-2-1',
                        "startdate": datetime.datetime(2017,2,7),
                        "enddate": datetime.datetime(2025,4,6)},

               "GOES-17":   {"url": 'goes17',
                        "avg1m": 'sci_xrsf-l2-avg1m_g17_',
                        "flsum": 'sci_xrsf-l2-flsum_g17_',
                        "avg_version": '_v2-2-1',
                        "flsum_version": '_v2-2-1',
                        "startdate": datetime.datetime(2018,6,1),
                        "enddate": datetime.datetime(2024,8,11)},

               "GOES-18":   {"url": 'goes18',
                        "avg1m": 'sci_xrsf-l2-avg1m_g18_',
                        "flsum": 'sci_xrsf-l2-flsum_g18_',
                        "avg_version": '_v2-2-1',
                        "flsum_version": '_v2-2-1',
                        "startdate": datetime.datetime(2022,9,2),
                        "enddate": datetime.datetime.now()},

               "GOES-19":   {"url": 'goes19',
                        "avg1m": 'sci_xrsf-l2-avg1m_g19_',
                        "flsum": 'sci_xrsf-l2-flsum_g19_',
                        "avg_version": '_v2-2-1',
                        "flsum_version": '_v2-2-1',
                        "startdate": datetime.datetime(2025,2,24),
                        "enddate": datetime.datetime.now()} #update when GOES-19 ends

              }

    return sat_info


def goes_xray_satellite_covers_date(experiment, request_date):
    """ Check if the spacecraft provided data during the requested date. """
    sat_info = goes_xray_satellite_info()

    if isinstance(request_date, str):
        request_date=dh.str_to_datetime(request_date)

    if experiment not in sat_info.keys():
        print(f"NOAA X-ray science data not available for {experiment}.")
        return False
    
    if request_date >= sat_info[experiment]["startdate"] \
        and request_date <= sat_info[experiment]["enddate"]:
        return True
    else:
        return False


def find_goes_xray_satellite_for_date(request_date):
    """ Find a GOES X-ray satellite with data for the requested date. """

    if isinstance(request_date, str):
        request_date=dh.str_to_datetime(request_date)

    sat_info = goes_xray_satellite_info()
    keys = list(sat_info.keys())
    
    #Search in reverse time order for satellite that observed during desired date
    for i in range(len(keys)-1,-1,-1):
        stdate = sat_info[keys[i]]["startdate"]
        enddate = sat_info[keys[i]]["enddate"]

        if request_date >= stdate and request_date <= enddate:
            satellite = keys[i]
            return satellite

    #if satellite not found
    return None


def find_all_goes_xray_satellites_for_date(request_date):
    """  Return a list of all GOES X-ray satellites that were active 
        at the requested time.
        
    """
    if isinstance(request_date, str):
        request_date=dh.str_to_datetime(request_date)

    sat_info = goes_xray_satellite_info()
    keys = list(sat_info.keys())
    
    satellites = []
    #Search in reverse time order for satellite that observed during desired date
    for i in range(len(keys)-1,-1,-1):
        stdate = sat_info[keys[i]]["startdate"]
        enddate = sat_info[keys[i]]["enddate"]

        if request_date >= stdate and request_date <= enddate:
            satellites.append(keys[i])

    return satellites




def download_goes_xray_science_data(request_date, experiment):
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
        
        :request_date: (datetime, string) desired day specified by the user
        :experiment: (string) name of GOES satellite
        
        OUTPUTS:
        
        :avg_filename: (string) name of the file containing the GOES XRS data 
        :flsum_filename: (string) name of the file containing the flare summary data
        :success: (bool) True file downloaded and/or successfully on computer
        
    """
    print(f"download_goes_xray_science_data: Checking for {experiment} data on {request_date}.")

    if isinstance(request_date, str):
        request_date=dh.str_to_datetime(request_date)

    avg_filename = None #X-ray time series
    flsum_filename = None #flare summary with magnitude, start, peak, end, etc
    success = False #flsum files were effectively downloaded?

    if experiment == None:
        print("download_goes_xray_science_data: No experiment was specified.")
        return avgfile, flsumfile, success

    if experiment in old_goes_sc:
        print(f"download_goes_xray_science_data: {experiment} X-ray science data is not yet available from NOAA. Returning.")
        return avgfile, flsumfile, success

    #Create path for GOES X-ray data if doesn't exist
    if not os.path.isdir(os.path.join(cfg.datapath,goesX_dir)):
        os.mkdirs(os.path.join(cfg.datapath,goesX_dir))

    sat_info = goes_xray_satellite_info()

    #GOES-08 to -15 files like:
    #https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes08/xrsf-l2-avg1m_science/1995/01/
    #https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes09/xrsf-l2-flsum_science/1996/04/

    #GOES-R like
    #https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-avg1m_science/2017/02/
    #https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-flsum_science/2017/02/


    year = request_date.year
    month = request_date.month
    day = request_date.day
    date_suffix = 'd%i%02i%02i' % (year,month,day)

    if request_date < sat_info[experiment]['startdate'] or request_date > sat_info[experiment]['enddate']:
        print(f"{experiment} X-ray science data is not available for {date}. Skipping.")
        return avg_filename, flsum_filename, success

    #Check if file is present or needs to be downloaded.
    #Filenames like:
    #sci_xrsf-l2-avg1m_g08_d19950103_v1-0-0.nc
    #sci_xrsf-l2-flsum_g08_d19950103_v2-1-1.nc
    avgfile = f"{sat_info[experiment]['avg1m']}{date_suffix}{sat_info[experiment]['avg_version']}.nc"
    flsumfile = f"{sat_info[experiment]['flsum']}{date_suffix}{sat_info[experiment]['flsum_version']}.nc"

    #Check if files already on your computer
    avgfullpath = os.path.join(cfg.datapath,goesX_dir,avgfile)
    avg_exists = os.path.isfile(avgfullpath)
    if avg_exists:
        avg_filename = avgfile
    
    flsumfullpath = os.path.join(cfg.datapath,goesX_dir,flsumfile)
    flsum_exists = os.path.isfile(flsumfullpath)
    if flsum_exists:
        flsum_filename = flsumfile

    #If already have both files
    if avg_exists and flsum_exists:
        success = True
        return avg_filename, flsum_filename, success
    

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
            urllib.request.urlopen(avgurl)

            wget.download(avgurl, avgfullpath)
            print(f"\ndownload_goes_xray_science_data: Downloaded {avgurl}")
            avg_filename = avgfile
        except urllib.request.HTTPError as e:
            print(f"Cannot access file at {avgurl} because {e}. Please check that selected spacecraft covers date range.")
        except socket.timeout as e:
            print(f"Cannot access file at {avgurl} because {e}.")
        except Exception as e:
            print(f"Cannot access file at {avgurl} because {e}.")



    if not flsum_exists:
        try:
            urllib.request.urlopen(flsumurl)

            wget.download(flsumurl, flsumfullpath)
            print(f"\ndownload_goes_xray_science_data: Downloaded {flsumurl}")
            flsum_filename = flsumfile
            success = True
        except urllib.request.HTTPError as e:
            print(f"Cannot access file at {flsumurl} because {e}. Please check that selected spacecraft covers date range.")
        except socket.timeout as e:
            print(f"Cannot access file at {flsumurl} because {e}.")
        except Exception as e:
            print(f"Cannot access file at {flsumurl} because {e}.")


    return avg_filename, flsum_filename, success




def download_goes_xray_science_data_range(startdate, enddate, experiment):
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
        
        :startdate: (datetime, string) start of time period specified by user
        :enddate: (datetime, string) end of time period entered by user
        :experiment: (string) name of GOES satellite
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            XRS data that span the desired time range
            (monthly files)
        
    """

    if isinstance(startdate, str):
        startdate=dh.str_to_datetime(startdate)

    if isinstance(enddate, str):
        enddate=dh.str_to_datetime(enddate)

    avgfiles = [] #X-ray time series
    flsumfiles = [] #flare summary with magnitude, start, peak, end, etc
    success = False #flsum files were effectively downloaded?

    if experiment == None:
        print("download_goes_xray_science_data: No experiment was specified.")
        return avgfiles, flsumfiles, success

    if experiment in old_goes_sc:
        print(f"download_goes_xray_science_data: {experiment} X-ray science data is not yet available from NOAA. Returning.")
        return avgfiles, flsumfiles, success

        
    #GOES XRS science data is stored in daily data files
    td = enddate - startdate
    NFILES = td.days + 1 #number of data files to download
    if td.seconds > 0: NFILES = NFILES + 1

    sat_info = goes_xray_satellite_info()

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

        if date < sat_info[experiment]['startdate'] or date > sat_info[experiment]['enddate']:
            print(f"{experiment} X-ray science data is not available for {date}. Skipping.")
            continue

        avg_filename, flsum_filename, file_success = download_goes_xray_science_data(date, experiment)

        if file_success:
            avgfiles.append(os.path.join(goesX_dir,avg_filename))
            flsumfiles.append(os.path.join(goesX_dir,flsum_filename))
    
    if len(flsumfiles) == NFILES:
        success = True

    return avgfiles, flsumfiles, success


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


def flare_info_dict():
    flare_info = {  'catalog_id': None, #flare_id
                    'catalog': '',
                    'start_time': pd.NaT,
                    'peak_time': pd.NaT,
                    'end_time': pd.NaT,
                    'duration': np.nan,
                    'time_to_peak': np.nan,
                    'intensity': np.nan,
                    'class': None,
                    'instrument': '',
                    'integrated_intensity': np.nan, #start time to end time
                }
    return flare_info




def extract_flare_from_noaa_flsum(experiment, request_date, flare_info={}, req_flare_id=None):
    """ Extract information for a specific flare on a specific datetime.
        Enter a date within 15 minutes of flare start or end of for any time
        between the flare start and end. The flare peak time is the best reference.
        
        request_date is taken to be the flare peak time as it is the least subjective.
        
        INPUTS:
        
            :experiment: (string) GOES-08, etc
            :request_date: (datetime) timestamp for a date related to a specific flare,
                can be similar to flare start, peak, end. Will look compare date to
                15 minutes from the flare start and end.
                
        OUTPUTS:
        
            Information available flsum files.
            :info: (dictionary) carries flare start, peak, end, magnitude, fluence, etc.
        
    """
    sat_info = goes_xray_satellite_info()

    #Define flare info
    if not flare_info:
        flare_info = flare_info_dict()
        flare_info['instrument'] = experiment
        flare_info['catalog'] = 'NOAA_NCEI'

    if experiment in old_goes_sc:
        #NEED TO DIVIDE XRS-B FLUXES BY 0.7 TO APPROXIMATE SCIENCE VALUES UNTIL NOAA RELEASES SCIENCE DATA
        print(f"extract_flare_from_noaa_flsum: {experiment} X-ray science data is not yet available from NOAA. Returning.")
        return flare_info


    year = request_date.year
    month = request_date.month
    day = request_date.day
    date_suffix = 'd%i%02i%02i' % (year,month,day)

    #Filenames like:
    #sci_xrsf-l2-flsum_g08_d19950103_v1-0-0.nc
    flsumfile = f"{sat_info[experiment]['flsum']}{date_suffix}{sat_info[experiment]['flsum_version']}.nc"

    fullpath = os.path.join(cfg.datapath,goesX_dir,flsumfile)
    flsum_exists = os.path.isfile(fullpath)
    if not flsum_exists:
        _avgfname, _flsumfname, file_success = download_goes_xray_science_data(request_date, experiment)
        if not file_success:
            print(f"extract_flare_from_noaa_flsum: Cannot find file {fullpath} and could not download.")
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
        print(f"extract_flare_from_noaa_flsum: User provided flare ID (catalog_id) {req_flare_id}")


    #-->If flare_id is specified
    ############################
    if not pd.isnull(req_flare_id):
        if req_flare_id in flare_ids:
            #Compile all information per flare by extracting according to req_flare_id
            select_idx = [i for i in range(len(data.variables["flare_id"][:])) if data.variables["flare_id"][i]==req_flare_id]
            if len(select_idx) > 0:
                print(f"extract_flare_from_noaa_flsum: Found flare for flare ID (catalog_id) {req_flare_id}")

    #-->else search according to date of flare peak
    ###############################################
    else:
        #starting point to find least time difference between request_date and flare peak
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
                #-->PRECEDENCE. Just in case there are two flares close together, give precedence
                #to the one containing requested date
                #####################################################################
                if (request_date >= flst) and (request_date <= flend):
                    ptd = abs(flpk - request_date)
                    if ptd < peak_td:
                        peak_td = ptd
                        select_idx = idx
                        print(f"extract_flare_from_noaa_flsum: Found {experiment} flare {flst} to {flend} for requested date {request_date}")
                #slightly looser requirements in case request_date was a bit off
                elif (request_date >= flst-td) and (request_date <= flend+td):
                    ptd = abs(flpk - request_date)
                    if ptd < peak_td:
                        peak_td = ptd
                        select_idx = idx
                        print(f"extract_flare_from_noaa_flsum: Found {experiment} flare {flst} to {flend} for requested date {request_date}")



            #-->For flares without end times: account for another flare starting before the first ended
            #or the last flare in the file is the desired one and goes on to the next day
            ######################################################################################
            if flst != None and flend == None:
                #Case of last flare in file and requested time is between the flare start and
                #the end of the day
                if index == len(flare_ids)-1:
                    dayend = datetime.datetime(year,month,day,23,59,59)
                    if (request_date >= flst-td) and (request_date <= dayend):
                        select_idx = idx
                        print(f"extract_flare_from_noaa_flsum: Found {experiment} flare {flst} to {dayend} for requested date {request_date}")
                        break
                #If the flare is not the last flare and there is a peak time
                elif flpk != None:
                    #look for the requested time within flare start and peak
                    if (request_date >= flst) and (request_date <= flpk+td):
                        ptd = abs(flpk - request_date)
                        if ptd < peak_td:
                            peak_td = ptd
                            select_idx = idx
                            print(f"extract_flare_from_noaa_flsum: Found {experiment} flare {flst} to {flpk+td} for requested date {request_date}")
                    #Slightly looser timing requirements in case the request date was a bit off
                    elif (request_date >= flst-td) and (request_date <= flpk+td):
                        ptd = abs(flpk - request_date)
                        if ptd < peak_td:
                            peak_td = ptd
                            select_idx = idx
                            print(f"extract_flare_from_noaa_flsum: Found {experiment} flare {flst} to {flpk+td} for requested date {request_date}")



            #-->Check if the first flare in the file is the desired one and starts on the previous day
            ######################################################################################
            if flst == None and flend != None:
                if index == 0:
                    #If requested date is between the start of the day and the flare end
                    dayst = datetime.datetime(year,month,day,0,0,0)
                    if (request_date >= dayst) and (request_date <= flend):
                        select_idx = idx
                        print(f"extract_flare_from_noaa_flsum: Found {experiment} flare {dayst} to {flend} for requested date {request_date}")
                        break
                    #Slightly looser timing requirements in case the request date was a bit off
                    elif flpk != None:
                        if (request_date >= dayst) and (request_date <= flend+td):
                            ptd = abs(flpk - request_date)
                            if ptd < peak_td:
                                peak_td = ptd
                                select_idx = idx
                                print(f"extract_flare_from_noaa_flsum: Found {experiment} flare {dayst} to {flend+td} for requested date {request_date}")
 
 

    #If no flare found
    if len(select_idx) == 0:
        print(f"extract_noaa_flare_info: {experiment} flare not found for requested date {request_date}. Returning.")
        return flare_info

    #flare_id, status, xrsb_flux, flare_class, time, integrated_flux, flare_class
    for ix in select_idx:
        if pd.isnull(flare_info['catalog_id']):
            flare_info['catalog_id'] = int(data.variables["flare_id"][ix])

        if data.variables["status"][ix] == "EVENT_START":
            flare_info['start_time'] = flare_dates[ix]
            
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



def flare_info_to_ccmc_json(flare_info):
    """ Convert the flare_info dictionary into the CCMC trigger format """

    ccmc_flare = ccmc_json.ccmc_flare_block()
    
    times = [flare_info['start_time'], flare_info['peak_time'], flare_info['end_time']]
    #Calculate last_data_time
    s = pd.Series(times)
    last_time = s.max()
    if not pd.isnull(last_time):
        ccmc_flare['last_data_time'] = dh.time_to_zulu(last_time)
    else:
        del ccmc_flare['last_data_time']

    ccmc_flare['start_time'] = dh.time_to_zulu(flare_info['start_time'])
    ccmc_flare['peak_time'] = dh.time_to_zulu(flare_info['peak_time'])
    ccmc_flare['end_time'] = dh.time_to_zulu(flare_info['end_time'])
    ccmc_flare['intensity'] = flare_info['intensity']
    ccmc_flare['integrated_intensity'] = flare_info['integrated_intensity']
    ccmc_flare['catalog'] = flare_info['catalog']
    ccmc_flare['catalog_id'] = flare_info['catalog_id']

    ccmc_flare = ccmc_json.clean_trigger_block(ccmc_flare)

    return ccmc_flare


def get_noaa_flare(request_time, experiment=None, format='dict'):
    """ Pull flare data from NOAA using the peak time to specify
        the flare. The peak is the most robust time to use, but
        can input any time between the start and end time of
        a flare.
        
        INPUTS:
        
            :request_time: (datetime, string) best if flare peak
            :experiment: (string) GOES-18, GOES-08, etc
            :format: (string) "dict" for flare_info format used by the lists (default); 
                "json" to return CCMC SEP Scoreboard trigger block
        
    """
    flare_info = {}

    if isinstance(request_time, str):
        request_time=dh.str_to_datetime(request_time)

    #identify all GOES X-ray satellites that cover time
    experiments = find_all_goes_xray_satellites_for_date(request_time)

    if experiment in old_goes_sc:
        print(f"get_noaa_flare: {experiment} X-ray science data is not yet available from NOAA.")
        return flare_info
    #try to find a satellite
    elif experiment == None:
        if len(experiments) == 0:
            print(f"get_noaa_flare: No GOES-08 to present satellites are available that cover your requested date.")
            return flare_info

    #Download flsum files from NOAA
    search_start = request_time-datetime.timedelta(hours=24)
    search_end = request_time+datetime.timedelta(hours=24)
    
    _xray, _flsum, success = download_goes_xray_science_data_range(search_start, search_end, experiment)
    if not success:
        for exper in experiments:
            _xray, _flsum, success = download_goes_xray_science_data_range(search_start, search_end, exper)
            if success:
                print(f"get_noaa_flare: NOAA files were not available for requested {experiment} at {request_time}. Using data from {exper} instead.")
                experiment = exper
                break

    if experiment == None:
        print(f"get_noaa_flare: Could not find a satellite that covers your requested date.")
        return flare_info

    #Extract flare info for specified flare time
    flare_info = extract_flare_from_noaa_flsum(experiment, request_time)

    #If found the flare in the flare summary files, but it spans files on two days
    if not pd.isnull(flare_info['catalog_id']):
        if not pd.isnull(flare_info['start_time']) and pd.isnull(flare_info['end_time']):

            #Look for the flare end in the file for the next day
            flare_info = extract_flare_from_noaa_flsum(flare_info['instrument'], request_time+datetime.timedelta(hours=24), flare_info=flare_info, req_flare_id=flare_info['catalog_id'])

        if pd.isnull(flare_info['start_time']) and not pd.isnull(flare_info['end_time']):

            #Look for the flare start in the file for the previous day
            flare_info = extract_flare_from_noaa_flsum(flare_info['instrument'], request_time-datetime.timedelta(hours=24), flare_info=flare_info, req_flare_id=flare_info['catalog_id'])


    if format == 'json':
        flare_info = flare_info_to_ccmc_json(flare_info)

    return flare_info


##########################################################
################# LMSAL FLARE CATALOG ####################
##########################################################

def get_lmsal_flares(start_date, end_date):
    """ Extract flares from the LMSAL HEK API between dates. """

    flares = []

    if isinstance(start_date, datetime.date):
        start_date = start_date.strftime('%Y-%m-%dT%H:%M:%S')
    if isinstance(end_date, datetime.date):
        end_date = end_date.strftime('%Y-%m-%dT%H:%M:%S')

    url = 'https://www.lmsal.com/hek/her?cosec=2&cmd=search&type=column&event_type=fl&event_starttime=' + start_date + '&event_endtime=' + end_date + '&event_coordsys=helioprojective&x1=-1200&x2=1200&y1=-1200&y2=1200'
    try:
        response = requests.get(url, timeout=10)
        selected = response.json()
        print(selected)
    except:
        print(f"get_lmsal_flares: could not make connection for {url}")
        return []
    
#
#    flare_info = {  'catalog_id': None, #flare_id
#                    'catalog': '',
#                    'start_time': pd.NaT,
#                    'peak_time': pd.NaT,
#                    'end_time': pd.NaT,
#                    'duration': np.nan,
#                    'time_to_peak': np.nan,
#                    'intensity': np.nan,
#                    'class': None,
#                    'instrument': '',
#                    'integrated_intensity': np.nan, #start time to end time
#                }

#"event_starttime",  "fl_fluence", "ar_noaaclass"
