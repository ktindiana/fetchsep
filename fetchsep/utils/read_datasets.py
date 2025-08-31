from . import config as cfg
from . import tools
import pandas as pd
import re
import calendar
import datetime
import argparse
from datetime import timedelta
import os
import wget
from calendar import monthrange
import urllib.request
import csv
from dateutil.parser import parse
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import netCDF4
import requests
from bs4 import BeautifulSoup
import tarfile
import ssl
import subprocess
import gzip
import shutil

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#2021-01-06, Changes in 0.2: Added SEPEMv3 data set and made changes to
#   necessary subroutines to accomodate the new data set.
#2021-02-25, Changes 0.3: Changed GOES-13, 14, 15 S14 option to include
#   S14 corrections to channels P6 and P7. Had previously only
#   applied S14 to P2 - P5 for those experiments.
#2021-09-24, Changes in 0.4: added a global_var called time_shift
#   which allows users to shift the times in user-input files by
#   time_shift number of hours. Changed in read_in_user_files and
#   added convert_decimal_hour.
#2021-11-16, Changes in 0.5: added support for GOES-16 and GOES-17
#   differential fluxes.
#2022-02-11, Changes in 0.6: Added support for GOES-R (16&17) primary
#   integral fluxes served by CCMC. These are the real time fluxes
#   from NOAA archived on the CCMC website. These are not the
#   official NOAA L2 integral fluxes. Those are not yet
#   available.
#2022-02-18, Changes in 0.7: Added checking for data/GOES-R
#   directory and will make if not present.
#2022-03-23, Changes in 0.8: Added ability to download and read GOES-16
#   SEP event file on NOAA's website in check_goesR_data. Modified
#   read_in_goesR data since the special file contains 30 days
#   of data with slightly different variable names.
#2022-05-20. Changes in 0.9: Changed SOHO/EPHIN L3 data from 30 minute
#   to 10 min data.
#2022-06-16, changes in 1.0: Added Shaowen Hu's recalibrated GOES
#   data set as a native data set in the code (SRAG1.2)
#2022-08-04, changes in 1.1: in extract_date_range, abjusted
#   the trimming so that the selected time range starts either
#   on or one point AFTER the specified start. Previously,
#   the point right before the specified start was included.
#2022-09-19, changes in 1.2: GOES-14 and GOES-15 hepad files from
#   2019-09-01 forward are missing a column. Added code to
#   read_in_goes() to change the expected columns for later dates.
#2022-11-20, changes in 1.3: Added STEREO-A and B to native data sets.
#2023-02-09, changes in 1.4: NOAA SWPC moved the location of the historical
#   GOES-15 and previous data. Updated the url in check_goes_data().
#   Updated check_goesR_data() to account for two different version
#   numbers possible in the differential files. Updated read_in_goesR()
#   to account for the different keys used to extract the flux
#   values in the different versions.
#2023-03-2, changes in 1.5: Combining with the read_datasets file
#   used by SEPAutoID (idsep).
#3023-06-05, changes in 1.6: NOAA changed the GOES-R data to version 3
#   in May of 2022. Added ability to grab v3 data in check_goesR.
#   Fixed bug that grabbed temperature uncorrected proton fluxes in
#   read_in_goesR() and added v3 file formatting.
#2023-06-19, changes in 1.7: NOAA added a v3-0-1 format for files
#   starting in April 2023. Rewrote check_goesR to be more versatile.
#   Added checking that include v3-0-1 in read_in_goesR.
#2024-09-09: Added GOES v3-0-2 format, which NOAA began using
#   exclusively on 2023-10-19.

ssl._create_default_https_context = ssl._create_unverified_context

datapath = cfg.datapath
outpath = cfg.outpath
plotpath = cfg.plotpath
badval = cfg.badval #bad data points will be set to this value; must be negative
user_col = cfg.user_col
user_delim = cfg.user_delim
user_energy_bins = cfg.user_energy_bins

#Spacecraft in the GOES-R+ series
goes_R = ["GOES-16", "GOES-17", "GOES-18", "GOES-19"]
#Spacecraft prior to GOES-R
goes_sc = ["GOES-08", "GOES-09","GOES-10","GOES-11",
            "GOES-12","GOES-13","GOES-14","GOES-15"]
#Spacecraft prior to GOES-08
old_goes_sc = ["GOES-05", "GOES-06", "GOES-07"]

def about_read_datasets():
    """ About read_datasets.py
        
        Subroutines that are required to read in the data sets
        native to this code or the user-specified data set.
        
        Reads GOES-08 to GOES-15, GOES-R+ (-16, -17, -18),
        SOHO/ERNE data, SOHO/EPHIN Level 3 data,
        SOHO/EPHIN data from the REleASE website, SRAG's CalGOES,
        STEREO-A, STEREO-B, SEPEM RDSv2 and SEPEM RDSv3
        (if the user downloads and unzips the
        files into the data directory).
        
        When possible, data is pulled from online databases and
        saved on your computer in the directory specified by
        datapath in global_vars.py. The default path is "data".
        
        Data files will be stored in subdirectories named as, e.g.:
        
        * data/GOES/
        * data/SEPEM/  (user must make this directory and put data inside)
        * data/SEPEMv3/ (user must make this directory and put data inside)
        * data/EPHIN
        * data/EPHIN_REleASE
        
        Users may make their own directories in data for their
        own files,e.g.:
        
        * data/SEPMOD/
        * data/MyDataSet/
        
    """

def check_paths():
    """Check that the paths that hold the data and output exist. If not, create.
    """
    print('Checking that paths exist: ' + datapath + ' and ' + outpath)
    if not os.path.isdir(datapath):
        print('check_paths: Directory containing fluxes, ' + datapath +
        ', does not exist. Creating.')
        os.mkdir(datapath);

    if not os.path.isdir(os.path.join(cfg.datapath, 'GOES')):
        print('check_paths: Directory containing GOES fluxes does not exist. Creating ' + datapath + '/GOES')
        os.mkdir(os.path.join(cfg.datapath, 'GOES'));

    if not os.path.isdir(os.path.join(cfg.datapath,'GOES_RT')):
        print('check_paths: Directory containing GOES_RT fluxes does not exist. Creating ' + datapath + '/GOES_RT')
        os.mkdir(os.path.join(cfg.datapath, 'GOES_RT'));

    if not os.path.isdir(os.path.join(cfg.datapath,'SEPEM')):
        print('check_paths: Directory containing SEPEM fluxes does not exist. Creating ' + datapath + '/SEPEM')
        os.mkdir(os.path.join(cfg.datapath,'SEPEM'));

    if not os.path.isdir(os.path.join(cfg.datapath,'SEPEMv3')):
        print('check_paths: Directory containing SEPEMv3 fluxes does not exist. Creating ' + datapath + '/SEPEMv3')
        os.mkdir(os.path.join(cfg.datapath,'SEPEMv3'));

    if not os.path.isdir(os.path.join(cfg.datapath, 'EPHIN')):
        print('check_paths: Directory containing EPHIN fluxes does not exist. Creating ' + datapath + '/EPHIN')
        os.mkdir(os.path.join(cfg.datapath, 'EPHIN'));

    if not os.path.isdir(os.path.join(cfg.datapath,'ERNE')):
        print('check_paths: Directory containing ERNE fluxes does not exist. Creating ' + datapath + '/ERNE')
        os.mkdir(os.path.join(cfg.datapath,'ERNE'));
    if not os.path.isdir(os.path.join(cfg.datapath,'ERNE','export.srl.utu.fi')):
        print('check_paths: Directory containing ERNE export fluxes does not exist. Creating ' + datapath + '/ERNE/export.srl.utu.fi')
        os.mkdir(os.path.join(cfg.datapath,'ERNE','export.srl.utu.fi'));

    if not os.path.isdir(os.path.join(cfg.datapath,'CalGOES')):
        print('check_paths: Directory containing CalGOES fluxes does not exist. Creating ' + cfg.datapath + '/CalGOES')
        os.mkdir(os.path.join(cfg.datapath,'CalGOES'));

    if not os.path.isdir(os.path.join(cfg.datapath,'STEREO-A')):
        print('check_paths: Directory containing STEREO-A fluxes does not exist. Creating ' + cfg.datapath +
        '/STEREO-A')
        os.mkdir(os.path.join(cfg.datapath,'STEREO-A'));
    if not os.path.isdir(os.path.join(cfg.datapath,'STEREO-A','LET')):
        os.mkdir(os.path.join(cfg.datapath,'STEREO-A','LET'));
    if not os.path.isdir(os.path.join(cfg.datapath,'STEREO-A','HET')):
        os.mkdir(os.path.join(cfg.datapath,'STEREO-A','HET'));

    if not os.path.isdir(os.path.join(cfg.datapath, 'STEREO-B')):
        print('check_paths: Directory containing STEREO-B fluxes does not exist. Creating ' + cfg.datapath +
        '/STEREO-B')
        os.mkdir(os.path.join(cfg.datapath,'STEREO-B'));
    if not os.path.isdir(os.path.join(cfg.datapath, 'STEREO-B','LET')):
        os.mkdir(os.path.join(cfg.datapath,'STEREO-B','LET'));
    if not os.path.isdir(os.path.join(cfg.datapath, 'STEREO-B','HET')):
        os.mkdir(os.path.join(cfg.datapath,'STEREO-B','HET'));
        
    if not os.path.isdir(os.path.join(cfg.datapath, 'ACE')):
        print('check_paths: Directory containing ACE fluxes does not exist. Creating ' + cfg.datapath +
        '/ACE')
        os.mkdir(os.path.join(cfg.datapath,'ACE'));
    if not os.path.isdir(os.path.join(cfg.datapath, 'ACE','SIS')):
        os.mkdir(os.path.join(cfg.datapath,'ACE','SIS'));
    if not os.path.isdir(os.path.join(cfg.datapath, 'ACE','EPAM')):
        os.mkdir(os.path.join(cfg.datapath,'ACE','EPAM'));
        
    if not os.path.isdir(cfg.outpath):
        print('check_paths: Directory to store output information does not exist. Creating ' + cfg.outpath)
        os.mkdir(cfg.outpath);
    if not os.path.isdir(cfg.plotpath):
        print('check_paths: Directory to store plots does not exist. Creating ' + cfg.plotpath)
        os.mkdir(cfg.plotpath);
    if not os.path.isdir(cfg.templatepath):
        print('check_paths: Directory to store user json templates for opsep does not exist. Creating.')
        os.mkdir(cfg.templatepath);


def read_data_manager():
    """ Create data management file that indicates whether a file is complete
        or not.
        
    """
    fname = os.path.join(cfg.datapath,"fetchsep_data_manager.csv")
    exists = os.path.isfile(fname)
    
    if not exists:
        df = pd.DataFrame(columns = ['File', 'Experiment', 'FluxType', 'Start Time', 'Cadence', 'Resolution', 'Complete'])
        df.to_csv(fname)
    else:
        df = pd.read_csv(fname)
        
    return df
    

def add_to_data_manager(df, experiment, flux_type, filename, start_time, cadence, resolution, complete):
    """ Add a new line to the data manager dataframe

    """
    vals = [filename, experiment, flux_type, start_time, cadence, resolution, complete]
    dfadd = pd.DataFrame([vals],columns = ['File', 'Experiment', 'FluxType', 'Start Time', 'Cadence', 'Resolution', 'Complete'])

    df = pd.concat([df,dfadd])
                
    return df


def check_completeness(experiment, flux_type, filename, df=pd.DataFrame):
    """ Check a data manager file to see if data has been recorded as 
        complete or incomplete. 
        
    """
    unmanaged = ['SEPEM', 'SEPEMv3', 'EPHIN_REleASE', 'CalGOES']
    if experiment in unmanaged:
        print(f"file_completeness: {experiment} file completeness is managed by the "
            "user. Please ensure you have the most up-to-date data.")
        return True

    if experiment == 'ERNE':
        print(f"file_completeness: {experiment} uses variable time periods and cannot "
            "be checked for completeness. Double check data/ERNE/*.dates to ensure you "
            "have the most up-to-date data.")
        return True
        
    if experiment == 'IMP8_CPME':
        print(f"file_completeness: {experiment} uses variable time periods and cannot "
            "be checked for completeness. Double check data files to ensure you "
            "have the most up-to-date data.")
        return True

    if df.empty:
        read_data_manager()

    #Check if the file of interest has already been determined to be complete/incomplete
    #Read in datafile containing a list of all files downloaded by fetchsep
    #with indicated completeness
    sub = df.loc[(df['File'] == filename)]
    
    #For testing and debugging. Check if files get added more than once.
    if len(sub) > 1:
        print(f"file_completeness: Multiple entries for {filename} in the fetchsep data manager file. "
            "This should not happen. Exiting.")
        return

    complete = None
    if not sub.empty:
        complete = sub['Complete'].iloc[0]
    
    return complete


def write_data_manager(df):
    """ Write data manager dataframe to file. """
    df.to_csv(os.path.join(cfg.datapath,"fetchsep_data_manager.csv"), index=False)



def file_completeness(df, experiment, flux_type, filename, dates):
    """Depending on when a file is downloaded, the file may not contain all the data.
        i.e. a yearly file downloaded prior to the end of the year and likewise
        for a monthly or daily file.
        
        FetchSEP must be able to determine whether a file on the user's computer is
        complete or needs to be downloaded again to get the final data.
        
        This subroutine is relevant to data products that are downloaded from
        the internet by FetchSEP.
        
        SEPEM, SEPEMv3, EPHIN_REleASE, CalGOES data must be downloaded by the user
        manually, so they are not managed here.  
        
        INPUT:
            
            :df: (pandas DataFrame) contains the data manager information about
                file completeness, read in with read_data_manager
            :experiment: (string) any of the experiments that can be natively run 
                by fetchsep
            :flux_type: (string) integral or differential
            :filename: (string) full path and filename
            :dates: (array, list) dates associated with the file
            
        OUTPUT:
        
            :df: (pandas DataFrame) updated data manager dataframe
    
    """
    
    complete = check_completeness(experiment, flux_type, filename, df=df)
    if complete:
        return df

    #If dates happens to be empty
    if not dates:
        print(f"file_completeness: The dates array passed for {filename} is empty. Returning with no action.")
        return df

    #If complete is False or None

    #Note that EPHIN must match the res = '5min' variable set in check_ephin_data()
    manager = { 'GOES': {'cadence': 'month',
                        'resolution': datetime.timedelta(minutes=5)},
                'GOES_R': {'cadence': 'day',
                        'resolution': datetime.timedelta(minutes=5)},
                'GOES_RT': {'cadence': 'day',
                        'resolution': datetime.timedelta(minutes=5)},
                'EPHIN': {'cadence': 'year',
                        'resolution': datetime.timedelta(minutes=5)},
                'ERNE': {'cadence': 'variable',
                        'resolution': datetime.timedelta(minutes=1)},
                'STEREO HET': {'cadence': 'month',
                        'resolution': datetime.timedelta(minutes=1)},
                'STEREO LET': {'cadence': 'day',
                        'resolution': datetime.timedelta(minutes=1)},
                'ACE_SIS': {'cadence': 'day',
                        'resolution': datetime.timedelta(minutes=5)},
                'ACE_EPAM_electrons': {'cadence': 'day',
                        'resolution': datetime.timedelta(minutes=5)},
                'IMP8_CPME': {'cadence': 'variable',
                        'resolution': datetime.timedelta(seconds=330)}
    }

    key = experiment
    if experiment in old_goes_sc: key = 'GOES'
    if experiment in goes_sc: key = 'GOES'
    if experiment in goes_R: key = 'GOES_R'
    if 'STEREO' in experiment:
        if 'HET' in filename: key = 'STEREO HET'
        if 'LET' in filename: key = 'STEREO LET'
    
    cadence = manager[key]['cadence']
    resolution = manager[key]['resolution']
 
    year = dates[0].year
    month = dates[0].month
    day = dates[0].day

    #Assume monthly files go from 1st to 30/31st
    #Assume yearly files go Jan - Dec
    #Identify expected last timestampe
    if cadence == 'year':
        last_timestamp = datetime.datetime(year+1,1,1) - resolution

    if cadence == 'month':
        month += 1
        if month == 13:
            month = 1
            year = year + 1
        
        last_timestamp = datetime.datetime(year=year,month=month,day=1) - resolution

    if cadence == 'day':
        last_timestamp = datetime.datetime(year=year, month=month, day=day) \
                        + datetime.timedelta(days=1) - resolution

        #Sometimes missing timestamps in STEREO LET historical data and will never be complete.
        #Consider complete if data goes up to minute 58
        if key == 'STEREO LET':
            last_timestamp = datetime.datetime(year=year, month=month, day=day) \
                        + datetime.timedelta(days=1) - 2*resolution

    #Check for completeness
    complete = None
    if dates[-1] >= last_timestamp:
        complete = True
    else:
        complete = False
    
    #Update data manager dataframe
    index = df[df['File'] == filename].index.values
    if index.size == 0:
        df = add_to_data_manager(df, experiment, flux_type, filename, dates[0], cadence, resolution, complete)
    else:
        df.at[index[0],'Complete'] = complete

    print(f"file_completeness: {filename} determined to be complete {complete}.")

    return df

                        
                        

def make_yearly_files(filename):
    """ Convert a large data set into yearly files.
        
        INPUTS:
        
        :filename: (string) filename of the SEPEM or
            SEPEMv3 full data file (1974 - 2015 or 2017)
            
        OUTPUTS:
        
        No output except that yearly files are created
        with _YYYY.csv appended at the end.
    """
    print('Breaking up the SEPEM data into yearly data files. (This could '
            + 'take a while, but you will not have to do it again.)')
    fnamebase = filename.replace('.csv','')  #if csv file
    fnamebase = fnamebase.replace('.txt','')  #if txt file

    with open(datapath + '/' + filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        has_header = csv.Sniffer().has_header(csvfile.readline())
        if has_header:
            next(readCSV)  # Skip single header row.
        ncol = len(next(readCSV))
        csvfile.seek(0) #back to beginning of file
        if has_header:
            header = csvfile.readline()  # save header row.

        check_year = 0
        for row in readCSV:
            date = datetime.datetime.strptime(row[0][0:19],
                                            "%Y-%m-%d %H:%M:%S")
            year = date.year
            if check_year != year:
                if check_year != 0:
                    outfile.close()
                    outfname = fnamebase + '_' + str(year) + '.csv'
                    outfile = open(datapath + '/' + outfname,'w+')
                    check_year = year

                if check_year == 0:
                    outfname = fnamebase + '_' + str(year) + '.csv'
                    outfile = open(datapath + '/' + outfname,'w+')
                    if has_header:
                        outfile.write(header)
                    check_year = year

            outfile.write(','.join(row))
            outfile.write('\n')

    outfile.close()
    csvfile.close()
    return


def check_sepem_data(startdate, enddate, experiment, flux_type):
    """Check if SEPEM data is present on the computer. Break into yearly
        files if needed. Return SEPEM filenames for analysis.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SEPEM
            data that span the desired time range (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    filenames1 = []  #SEPEM, eps, or epead

    year = styear

    if experiment == 'SEPEM':
        basenm = 'SEPEM_H_reference'
        dir = experiment

    if experiment == 'SEPEMv3':
        basenm = 'SEPEM_RDS_v3_H'
        dir = experiment

    while (year <= endyear):
        fname = basenm + '_' + str(year) + '.csv'
        exists = os.path.isfile(datapath + '/' + dir + '/' + fname)
        if exists:
            filenames1.append(dir + '/' + fname)
            year = year + 1
        if not exists:
            full_exists = os.path.isfile(datapath + '/' + dir + '/' + \
                                 '/' + basenm + '.txt')
            if not full_exists:
                if experiment == 'SEPEM':
                    sys.exit("Please download and unzip the RSDv2 data set."
                        " You may download the file at"
                        " http://sepem.eu/help/SEPEM_RDS_v2-00.zip for full "
                        "fluxes or http://sepem.eu/help/SEPEM_RDS_v2-00.zip "
                        "for ESA background-subtracted fluxes.")
                if experiment == 'SEPEMv3':
                    sys.exit('Please contact DH Consultancy for the SEPEM '
                            'RDSv3 data set. Unzip and put SEPEM_RDS_V3_H.txt '
                            'in the data/SEPEMv3 folder.')
            if full_exists:
                #Break up SEPEM data set into yearly files
                print('The SEPEM (RSDv2 and RDSv3) is more tractable when '
                        'breaking into yearly data files. '
                        'Producing yearly files.')
                make_yearly_files(dir + '/' + basenm + '.txt')
                year = styear

    return filenames1
    
    
def check_calgoes_data(startdate, enddate, experiment, flux_type):
    """Check if SRAG1.2 (CalGOES) data is present on the computer. Break into yearly
        files if needed. Return SRAG filenames for analysis.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SEPEM
            data that span the desired time range (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    filenames1 = []  #SEPEM, eps, or epead

    year = styear

    if experiment == 'CalGOES':
        basenm = 'srag12'
        dir = experiment


    while (year <= endyear):
        fname = basenm + '_' + str(year) + '.dat'
        exists = os.path.isfile(datapath + '/' + dir + '/' + fname)
        if exists:
            filenames1.append(dir + '/' + fname)
            year = year + 1
        if not exists:
            sys.exit('Please contact Shaowen Hu (shaowen.hu-1@nasa.gov) '
                    'for his recalibrated data set. Unzip and put the .dat files '
                    'in the data/CalGOES folder.')

    return filenames1



def check_old_goes_data(startdate, enddate, experiment, flux_type):
    """Check that GOES data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
            (monthly files)

        
    """
    
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #eps
    filenames2 = []  #hepad
    filenames_orien = []  #orientation flag (not needed)

    #GOES data is stored in monthly data files
    get_years = []
    get_months = []
    test_year = styear
    test_month = stmonth
    test_date = datetime.datetime(year=test_year, month=test_month, day=1)
    while (test_date < enddate):
        get_years.append(test_year)
        get_months.append(test_month)
        test_month = test_month + 1
        if (test_month > 12):
            test_month = 1
            test_year = test_year + 1
        test_date = datetime.datetime(year=test_year, month=test_month, day=1)

    NFILES = len(get_months)  #number of data files to download

    #Set correct file prefix for data files
    if experiment == "GOES-05": #No HEPAD
        prefix1 = 'g05_eps_5m_3s_'
        prefix2 = ''
        satellite = 'goes05'

    if experiment == "GOES-06":
        prefix1 = 'g06_eps_5m_'
        prefix2 = 'g06_hepad_5m_'
        satellite = 'goes06'
    
    if experiment == "GOES-07": #No HEPAD
        prefix1 = 'g07_eps_5m_'
        prefix2 = ''
        satellite = 'goes07'


    #for every month that data is required, check if file is present or
    #needs to be downloaded.
    for i in range(NFILES):
        year = get_years[i]
        month = get_months[i]
        last_day = calendar.monthrange(year,month)[1]
        date = datetime.datetime(year=year,month=month, day=1)
        date_suffix = '%i%02i01_%i%02i%02i' % (year,month,year,month,
                        last_day)
        fname1 = prefix1 + date_suffix + '.csv'
        fullpath1 = os.path.join(cfg.datapath, 'GOES', fname1)
        exists1 = os.path.isfile(fullpath1)

        complete = None
        if exists1:
            complete = check_completeness(experiment, flux_type, fullpath1, df=df)
        
        if not exists1 or not complete: #download file if not found on your computer
            url = ('https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/' +
                '%i/%02i/%s/csv/%s' % (year,month,satellite,fname1))
            print('Downloading GOES data: ' + url)
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(fullpath1):
                    os.remove(fullpath1) # if exist, remove it directly
                
                wget.download(url, fullpath1)
            except urllib.request.HTTPError:
                print("Cannot access file at " + url +
                ". Please check that selected spacecraft covers date range. Skipping.")
                fname1 = None

        fname2 = None
        if prefix2 != '':
            fname2 = prefix2 + date_suffix + '.csv'
            fullpath2 = os.path.join(cfg.datapath,'GOES',fname2)
            exists2 = os.path.isfile(fullpath2)

            complete = None
            if exists2:
                complete = check_completeness(experiment, flux_type, fullpath2, df=df)

            if not exists2 or not complete: #download file if not found on your computer
                url = ('https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/' +
                   '%i/%02i/%s/csv/%s' % (year,month,satellite,fname2))
                print('Downloading GOES data: ' + url)
                try:
                    urllib.request.urlopen(url)
                    
                    if os.path.exists(fullpath2):
                        os.remove(fullpath2) # if exist, remove it directly
                    
                    wget.download(url, fullpath2)
                except urllib.request.HTTPError:
                    print("Cannot access file at " + url +
                   ". Please check that selected spacecraft covers date range. Skipping.")
                    fname2 = None

        if fname1 == None:
            filenames1.append(None)
        else:
            filenames1.append(os.path.join('GOES', fname1))
        if fname2 == None:
            filenames2.append(None)
        else:
            filenames2.append(os.path.join('GOES', fname2))

    return filenames1, filenames2, date



def check_goes_data(startdate, enddate, experiment, flux_type):
    """Check that GOES data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
            (monthly files)
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
        
    """
    #Have encountered some bad files and better to use a different
    #instrument
    #g14 20160801 bad because one day (last day) missing from orientation file
    #Only goes up to the 30th
    #Per communication with Juan Rodriguez, GOES-13 was always in the same
    #orientation afer 2010.
    #The orientation file is missing for 20130501_20130531. Per the previous
    #month's orientation file, B was the westward-facing detector.
#    g13_bad = ['g13_epead_cpflux_5m_20130501_20130531.csv',
#                'g13_epead_cpflux_5m_20141201_20141231.csv',
#                'g13_epead_cpflux_5m_20160801_20160831.csv',
#                'g13_epead_p17ew_5m_20160801_20160831.csv',
#                'g13_epead_p17ew_5m_20141201_20141231.csv']
#    g14_bad = ['g14_epead_cpflux_5m_20160801_20160831.csv',
#                'g14_epead_p17ew_5m_20160801_20160831.csv']
#    g15_bad = ['g15_epead_cpflux_5m_20160801_20160831.csv',
#                'g15_epead_p17ew_5m_20160801_20160831.csv',
#                'g15_epead_cpflux_5m_20190901_20190930.csv',
#                'g15_epead_cpflux_5m_20191001_20191031.csv']

    g13_bad = []
    g14_bad = []
    g15_bad = []
    
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, eps, or epead
    filenames2 = []  #hepad
    filenames_orien = []  #orientation flag for GOES-13+

    #GOES data is stored in monthly data files
    get_years = []
    get_months = []
    test_year = styear
    test_month = stmonth
    test_date = datetime.datetime(year=test_year, month=test_month, day=1)
    while (test_date < enddate):
        get_years.append(test_year)
        get_months.append(test_month)
        test_month = test_month + 1
        if (test_month > 12):
            test_month = 1
            test_year = test_year + 1
        test_date = datetime.datetime(year=test_year, month=test_month, day=1)

    NFILES = len(get_months)  #number of data files to download

    #Set correct file prefix for data files
    if experiment == "GOES-08":
        prefix1 = 'g08_eps_5m_'
        prefix2 = 'g08_hepad_5m_'
        satellite = 'goes08'

    if experiment == "GOES-09":
        prefix1 = 'g09_eps_5m_'
        prefix2 = 'g09_hepad_5m_'
        satellite = 'goes09'

    if experiment == "GOES-10":
        prefix1 = 'g10_eps_5m_'
        prefix2 = 'g10_hepad_5m_'
        satellite = 'goes10'

    if experiment == "GOES-11":
        prefix1 = 'g11_eps_5m_'
        prefix2 = 'g11_hepad_5m_'
        satellite = 'goes11'

    if experiment == "GOES-12":
        prefix1 = 'g12_eps_5m_'
        prefix2 = 'g12_hepad_5m_'
        satellite = 'goes12'

    if experiment == "GOES-13":
        prefix2 = 'g13_hepad_ap_5m_'
        prefix_orien = 'g13_epead_orientation_flag_1m_'
        satellite = 'goes13'
        if flux_type == "differential":
            prefix1 = 'g13_epead_p17ew_5m_'
        if flux_type == "integral":
            prefix1 = 'g13_epead_cpflux_5m_'

    if experiment == "GOES-14":
        prefix2 = 'g14_hepad_ap_5m_'
        prefix_orien = 'g14_epead_orientation_flag_1m_'
        satellite = 'goes14'
        if flux_type == "differential":
            prefix1 = 'g14_epead_p17ew_5m_'
        if flux_type == "integral":
            prefix1 = 'g14_epead_cpflux_5m_'

    if experiment == "GOES-15":
        prefix2 = 'g15_hepad_ap_5m_'
        prefix_orien = 'g15_epead_orientation_flag_1m_'
        satellite = 'goes15'
        if flux_type == "differential":
            prefix1 = 'g15_epead_p17ew_5m_'
        if flux_type == "integral":
            prefix1 = 'g15_epead_cpflux_5m_'

    #for every month that data is required, check if file is present or
    #needs to be downloaded.
    for i in range(NFILES):
        year = get_years[i]
        month = get_months[i]
        last_day = calendar.monthrange(year,month)[1]
        date = datetime.datetime(year=year,month=month, day=1)
        date_suffix = '%i%02i01_%i%02i%02i' % (year,month,year,month,
                        last_day)
        fname1 = prefix1 + date_suffix + '.csv'
        fullpath1 = os.path.join(cfg.datapath, 'GOES', fname1)
        exists1 = os.path.isfile(fullpath1)
        fname2 = prefix2 + date_suffix + '.csv'
        fullpath2 = os.path.join(cfg.datapath,'GOES',fname2)
        exists2 = os.path.isfile(fullpath2)
        if (experiment == "GOES-13" or experiment == "GOES-14"
            or experiment == "GOES-15"):
            fname_orien = prefix_orien + date_suffix + '_v1.0.0.csv'
            fullpath_orien = os.path.join(cfg.datapath, 'GOES', fname_orien)
            exists_orien = os.path.isfile(fullpath_orien)


        if fname1 in g13_bad or fname1 in g14_bad or fname1 in g15_bad:
            print("File set is known to be bad (" \
                    + fname1 + "). Skipping. ")
            fname1 = None

        complete = None
        if exists1:
            complete = check_completeness(experiment, flux_type, fullpath1, df=df)
        
        if not exists1 or not complete: #download file if not found on your computer
            url = ('https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/' +
                '%i/%02i/%s/csv/%s' % (year,month,satellite,fname1))
            print('Downloading GOES data: ' + url)
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(fullpath1):
                    os.remove(fullpath1) # if exist, remove it directly
                
                wget.download(url, fullpath1)
            except urllib.request.HTTPError:
                print("Cannot access file at " + url +
                ". Please check that selected spacecraft covers date range. Skipping.")
                fname1 = None

        complete = None
        if exists2:
            complete = check_completeness(experiment, flux_type, fullpath2, df=df)

        if not exists2 or not complete: #download file if not found on your computer
            url = ('https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/' +
               '%i/%02i/%s/csv/%s' % (year,month,satellite,fname2))
            print('Downloading GOES data: ' + url)
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(fullpath2):
                    os.remove(fullpath2) # if exist, remove it directly
                
                wget.download(url, fullpath2)
            except urllib.request.HTTPError:
                print("Cannot access file at " + url +
               ". Please check that selected spacecraft covers date range. Skipping.")
                fname2 = None

        if (experiment == "GOES-13" or experiment == "GOES-14"
            or experiment == "GOES-15"):
            if not exists_orien or not complete: #download file if not found on your computer
                url = ('https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/' +
                   '%i/%02i/%s/csv/%s' % (year,month,satellite,fname_orien))
                print('Downloading GOES data: ' + url)
                try:
                    urllib.request.urlopen(url)
                
                    if os.path.exists(fullpath_orien):
                        os.remove(fullpath_orien) # if exist, remove it directly
                    
                    wget.download(url, fullpath_orien)
                except urllib.request.HTTPError:
                    print("Cannot access orientation file at "
                        + url + ". Please check that selected "
                        + "spacecraft covers date range. Skipping.")
                    fname_orien = None

        if fname1 == None:
            filenames1.append(None)
        else:
            filenames1.append(os.path.join('GOES', fname1))
        if fname2 == None:
            filenames2.append(None)
        else:
            filenames2.append(os.path.join('GOES', fname2))
        if (experiment == "GOES-13" or experiment == "GOES-14"
            or experiment == "GOES-15"):
            if fname_orien == None:
                filenames_orien.append(None)
            else:
                filenames_orien.append(os.path.join('GOES', fname_orien))

    return filenames1, filenames2, filenames_orien, date




def check_goesR_data(startdate, enddate, experiment, flux_type):
    """Check that GOES data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
        GOES files are saved daily in cdf format.
        
        GOES differential fluxes are the official science-grade product
        provided by NOAA SWPC.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
            (monthly files)
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
        
    """
    if flux_type == "integral":
        sys.exit("check_goesR_data: This subroutine is only valid for GOES-R+ "
                "differential fluxes. Please set the flux_type to differential "
                "and try again.")
 
    g16_last_date = datetime.datetime(2025,4,6) #there is a file on the 7th, but has a problem
 
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #GOES-R
    filenames2 = []  #place holder
    filenames_orien = []  #place holder
    
    #SPECIAL FILE FOR 2017-09-10 SEP EVENTS
    if experiment == "GOES-16" and styear == 2017:
        fname1 = 'se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc'
        fullpath1 = os.path.join(cfg.datapath,'GOES',fname1)
        exists1 = os.path.isfile(fullpath1)
        if not exists1:
            url=('https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/goesr/solar_proton_events/sgps_sep2017_event_data/%s' % (fname1))
            try:
                urllib.request.urlopen(url)
                wget.download(url, fullpath1)
            except urllib.request.HTTPError:
                sys.exit("Cannot access SEP event file at " + url +
               ". Please check that the url is still active. Skipping.")
                

        filenames1.append(os.path.join('GOES',fname1))
        return filenames1, filenames2, filenames_orien, startdate
    
    
    #GOES-R data is stored in daily data files
    td = enddate - startdate
    NFILES = td.days #number of data files to download
    if td.seconds > 0: NFILES = NFILES + 1

    if experiment == "GOES-16":
        prefix = 'sci_sgps-l2-avg5m_g16_'
        satellite = 'goes16'

    if experiment == "GOES-17":
        prefix = 'sci_sgps-l2-avg5m_g17_'
        satellite = 'goes17'
        
    if experiment == "GOES-18": #2022-09-13 forward
        prefix = 'sci_sgps-l2-avg5m_g18_'
        satellite = 'goes18'


    #for every day that data is required, check if file is present or
    #needs to be downloaded.
    for i in range(NFILES):
        date = startdate + datetime.timedelta(days=i)
        year = date.year
        month = date.month
        day = date.day
        date_suffix = 'd%i%02i%02i' % (year,month,day)
 
        if experiment == "GOES-16":
            if date > g16_last_date:
                print(f"check_goesR_data: Requested {date}. "
                    f"Last available date for GOES-16 is {g16_last_date}. Continuing.")
 
        #GOES-R differential data has three possible version numbers
        file_ext = ['_v1-0-1.nc', '_v2-0-0.nc', '_v3-0-0.nc', '_v3-0-1.nc', '_v3-0-2.nc']
        
        foundfile = None
        for ext in file_ext:
            fname_data = prefix + date_suffix + ext
            fullpath = os.path.join(datapath,'GOES',fname_data)
            exists = os.path.isfile(fullpath)
            if exists:
                #Check if the file is listed as complete
                complete = check_completeness(experiment, flux_type, fullpath, df=df)
                if complete: foundfile = fname_data
            
        #Try versions
        if foundfile == None:
            for ext in file_ext:
                fname_data = prefix + date_suffix + ext
                fullpath = os.path.join(datapath,'GOES',fname_data)
                url=('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/%s/l2/data/sgps-l2-avg5m/%i/%02i/%s' % (satellite,year,month,fname_data))
                try:
                    urllib.request.urlopen(url)
                
                    if os.path.exists(fullpath):
                        os.remove(fullpath) # if exist, remove it directly
                        print(f"check_goesR_data: Removed existing file to avoid creating duplicates during download {fullpath}")
          
                    wget.download(url, os.path.join(datapath,'GOES',fname_data))
                    foundfile = fname_data
                    print(f"\ncheck_goesR_data: Downloaded {url}")
                    break
                except urllib.request.HTTPError:
                    foundfile = None

        if foundfile == None:
            print("Cannot access GOES file at " + url +
               ". Tried file versions " + str(file_ext) + ". Please check that selected spacecraft covers date range. Skipping.")
            
        if foundfile == None:
            filenames1.append(None)
        else:
            filenames1.append(os.path.join('GOES', foundfile))
            
        
    return filenames1, filenames2, filenames_orien, date



def rerequest(url, tries=0):
    """
    Runs requests.get() until a response is received to avoid chokepoints

    Parameters
    ----------
    url : string
        URL where desired data resides
    tries [untouched by user] : int
        Number of attempts
    """

    # end the requesting process if tries > 5
    if tries > 5:
        raise Exception(url + " refuses to respond after 5 attempts. Exiting.")
    try:
        output = requests.get(url, timeout=10.0)
        return output
    except requests.exceptions.Timeout as e:
        return rerequest(url, tries + 1)


def check_goes_RTdata(startdate, enddate, experiment, flux_type,
    spacecraft="primary"):
    """Check that GOES Real Time data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
        GOES real time integral files are saved daily in txt format.
        
        The CCMC HAPI server returns real time integral fluxes for GOES
        downloaded from the SWPC realtime 3-day jsons starting from
        2010-04-14 00:00:00.
        
        These fluxes are by default from the primary GOES satellite.
        
        The GOES-16 & 17 & 18 integral fluxes are currently only available from
        the real time product plotted by SWPC on a daily basis and archived at CCMC.
        These fluxes are not the official science-grade integral flux product from
        NOAA, as these are not yet available. Also, only the PRIMARY spacecraft
        integral fluxes are available and it isn't possible to choose
        between GOES-16 or GOES-17 or GOES-18 for the integral fluxes.
        
        When the official integral fluxes are available for GOES-R, support
        to select the spacecraft will be added to this package.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        :spacecraft: 9string) 
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
            (monthly files)
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
        
    """
    if flux_type == "differential":
        sys.exit("check_goes_RTdata: This subroutine is only valid for GOES real time "
                "integral fluxes. Please set the FluxType (flux_type) to integral and try "
                "again.")

    if spacecraft != "primary" and spacecraft != "secondary":
        sys.exit(f"check_goes_RTdata: spacecraft must be primary or secondary. You specified {spacecraft}. Please correct and try again.")

    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #GOES_RT files
    filenames2 = []  #place holder
    filenames_orien = []  #place holder

    #Choose to download daily data files
    td = enddate - startdate
    NFILES = td.days #number of data files to download
    if td.seconds > 0: NFILES = NFILES + 1

    #pulls primary spacecraft fluxes
    #fname like 20230501_Gp_part_5m.txt
    prefix = '_Gp_part_5m'


    #for every day that data is required, check if file is present or
    #needs to be downloaded.
    date = startdate
    for i in range(NFILES):
        date = startdate + datetime.timedelta(days=i)
        date2 = date + datetime.timedelta(days=1)
        year = date.year
        month = date.month
        day = date.day
        date_suffix = '%i%02i%02i' % (year,month,day)
 
        #Previously pulled a txt file archived by CCMC, but no longer
        #available, do using their HAPI API to query iSWA.
        fname1 = date_suffix + prefix + '_' + spacecraft + '.csv'
        fullpath1 = os.path.join(datapath, 'GOES_RT', fname1)
        exists1 = os.path.isfile(fullpath1)
        
        complete = False
        if exists1:
            #Check if the file is listed as complete
            complete = check_completeness(experiment, flux_type, fullpath1, df=df)

        if not exists1 or not complete:
            #https://iswa.gsfc.nasa.gov/IswaSystemWebApp/hapi/data?id=goesp_part_flux_P5M&time.min=2023-05-23T00:00:00.0Z&time.max=2023-05-24T00:00:00.0Z&format=csv
            #Note that iSWA will return an empty csv file if make a request into the future or the
            #data isn't present yet. This can cause idsep or opsep to think the data is
            #already present on the computer and leave it as an empty file.

            if spacecraft == "primary":
                url=('https://iswa.gsfc.nasa.gov/IswaSystemWebApp/hapi/data?id=goesp_part_flux_P5M&time.min=%i-%02i-%02iT00:00:00.0Z&time.max=%i-%02i-%02iT00:00:00.0Z&format=csv' % (year,month,day,date2.year,date2.month,date2.day))
            elif spacecraft == "secondary":
                url=('https://iswa.gsfc.nasa.gov/IswaSystemWebApp/hapi/data?id=goess_part_flux_P5M&time.min=%i-%02i-%02iT00:00:00.0Z&time.max=%i-%02i-%02iT00:00:00.0Z&format=csv' % (year,month,day,date2.year,date2.month,date2.day))
            
            print("Trying to download url: " + url)
            response = rerequest(url)
            if response.status_code == 200:
                data = response.text
                fileout = open(os.path.join(datapath, 'GOES_RT', fname1),'w')
                fileout.write(data)
                fileout.close()
            else:
                print(f'Failed to retrieve data. HTTP Status code: {response.status_code}. Skipping.')
                fname1 = None

        if fname1 == None:
            filenames1.append(None)
        else:
            filenames1.append(os.path.join('GOES_RT', fname1))

    return filenames1, filenames2, filenames_orien, date



def check_preferential_goes_data(startdate, enddate, experiment, flux_type, spacecraft="primary"):
    """ If the user specifies only "GOES" for the experiment, then
        read in any of the GOES experiments and string them together.
        This will allow the creation of a long flux timeseries
        with multiple GOES satellites included.
        
        Data availability are searched for spacecraft in an order
        of preference, given by:
        ["GOES-13","GOES-15","GOES-11", "GOES-14",
        "GOES-09", "GOES-08","GOES-16","GOES-17","GOES-18"]
        
        The first spacecraft encountered with data for the requested
        time period is selected.
        
        ONLY ALLOW FOR INTEGRAL FLUXES, OTHERWISE RUN INTO PROBLEM
        OF MISMATCHING ENERGY BINS.
        
    """
    #Search GOES data in a search order that emphasizes GOES-15 and
    #GOES-13 for most recent data, then looks at GOES-16 and GOES-17
    #Then continues backwards in time to find a GOES that covers
    #the requested time range
    #GOES_RT goes back to 2010 or 2011 on CCMC servers, where they
    #archived the primary and secondary integral real time streams.
    #Best to use the NOAA archived integral fluxes for those
    #time periods. However, for GOES-16 forward, the CCMC GOES_RT
    #real time stream is the only accessible archive of GOES
    #integral fluxes (as of June 2025).
    
    goes = []
    if flux_type == "integral":
        if startdate >= datetime.datetime(year=2020,month=3,day=8):
            goes = ["GOES_RT"]
        else:
            goes = ["GOES-13","GOES-15","GOES-11", "GOES-14",
                    "GOES-09", "GOES-08", "GOES-07", "GOES-06",
                    "GOES-05"]

    if flux_type == "differential":
        if startdate >= datetime.datetime(year=2020,month=3,day=8):
            goes = ["GOES-16","GOES-18", "GOES-19", "GOES-17"]
                
        else:
            goes = ["GOES-13","GOES-15","GOES-11", "GOES-14",
                    "GOES-09", "GOES-08", "GOES-07", "GOES-06",
                    "GOES-05"]
            
            
    filenames1_all = []
    filenames2_all = []
    filenames_orien_all = []
    detector = []
    
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    #Write dates and experiment used for those dates
    #This is meant to allow the user to know which
    #spacecraft is available for which time period
    outfname = os.path.join(outpath, "idsep", "goes_experiments_dates.txt")
    outfile = open(outfname,'w')
    
    
    #GOES data is stored in monthly data files
    #GOES-R is stored in daily files, but use monthly cadence
    get_years = []
    get_months = []
    test_year = styear
    test_month = stmonth
    test_date = datetime.datetime(year=test_year, month=test_month, day=1)
    while (test_date < enddate):
        get_years.append(test_year)
        get_months.append(test_month)
        test_month = test_month + 1
        if (test_month > 12):
            test_month = 1
            test_year = test_year + 1
        test_date = datetime.datetime(year=test_year, month=test_month, day=1)
    
    
    #Read in data one month at a time
    
    for k in range(len(get_months)):
        year = get_years[k]
        month = get_months[k]
        last_day = calendar.monthrange(year,month)[1]
        
        if year == endyear and month == endmonth:
            last_date = endday
        
        stdate = datetime.datetime(year=year,month=month, day=1)
        enddt = datetime.datetime(year=year,month=month, day=last_day)
        print(f"======Date: {stdate} to {enddt}")
        
        for i in range(len(goes)):
            print(f"Testing {goes[i]}")
            
            filenames1 = []
            filenames2 = []
            filenames_orien = []

            if goes[i] in old_goes_sc:
                filenames1, filenames2 = \
                     check_old_goes_data(stdate, enddt, goes[i], flux_type)

            if goes[i] in goes_sc:
                filenames1, filenames2, filenames_orien, date = \
                     check_goes_data(stdate, enddt, goes[i], flux_type)
                                     
            if goes[i]=="GOES_RT" and flux_type == "integral":
                filenames1, filenames2, filenames_orien, date = \
                 check_goes_RTdata(stdate, enddt, goes[i], flux_type, spacecraft=spacecraft)
                 
            if goes[i] in goes_R and flux_type == "differential":
                filenames1, filenames2, filenames_orien, date = \
                     check_goesR_data(stdate, enddt, goes[i], flux_type, spacecraft=spacecraft)


            if filenames1 != [] and filenames1 != None:
                filenames1_all.extend(filenames1)
                detector.extend([goes[i]]*len(filenames1))
                
                outfile.write(str(year)+","+str(month)+"," +goes[i]+"\n")
                
                if filenames2 == []:
                    filenames2_all.extend([None]*len(filenames1))
                else:
                    filenames2_all.extend(filenames2)
                
                if filenames_orien == []:
                    filenames_orien_all.extend([None]*len(filenames1))
                else:
                    filenames_orien_all.extend(filenames_orien)

                print(f"Downloaded {goes[i]}")
                break #don't test other spacecraft if found one
    
    outfile.close()
    
    return filenames1_all, filenames2_all, filenames_orien_all, detector





def identify_which_goes_spacecraft(startdate, enddate, spacecraft="primary"):
    """ For the provided start and end dates, return a list of GOES spacecraft
        which were the primary or secondary and the associated dates.
        This subroutine uses the file:
        fetchsep/reference/GOES_Primary_Secondary_Status.csv provided by 
        NOAA SWPC (Kim Moreland) for GOES-06 (1986) to September 6, 2024. 
        
        INPUTS:
        
            :startdate: (datetime) first date of date range
            :enddate: (datetime last date of date range
            :spacecraft: (string) primary or secondary
            
        OUTPUTS:
        
            :starttimes: (list of datetimes) arr of start times
            :endtimes: (list of datetimes) arr of end times
            :goes: (list of strings) associated primary or secondary GOES for those
                date ranges [GOES-11, GOES-13, etc]
        
    """
    df_goes = pd.read_csv('fetchsep/reference/GOES_Primary_Secondary_Status.csv',parse_dates=['Start Date', 'End Date'])
    df_goes = df_goes.loc[df_goes['Instrument']=='Protons']
    
    if spacecraft == 'primary':
        df_goes = df_goes.loc[df_goes['Status']=='Primary']
    elif spacecraft == 'secondary':
        df_goes = df_goes.loc[df_goes['Status']=='Secondary']
    else:
        sys.exit('identify_which_goes_spacecraft: spacecraft must be primary or secondary. You entered ' + spacecraft)
        
    firstdate = min(df_goes['Start Date'])
    lastdate = max(df_goes['End Date'])
    
    #Check that primary or secondary are specified between the dates requested
    if enddate <= firstdate or startdate >= lastdate:
        print(f"identify_which_goes_spacecraft: Primary and secondary GOES dates are available between: {firstdate} and {lastdate}. Please revise your date request: {stardate} to {enddate}. Returning empty arrays.")
        return [], []
    
    #Extract the entries between the dates of interest
    df_goes = df_goes.loc[(df_goes['Start Date'] < enddate) & (df_goes['End Date'] >= startdate)]

    #Replace the start and end dates with the ones requested by the user to get
    #date ranges for each primary/secondary spacecraft
    df_goes.loc[df_goes['Start Date'] < startdate, 'Start Date'] = startdate
    df_goes.loc[df_goes['End Date'] > enddate, 'End Date'] = enddate

    print(f"{spacecraft} spacecraft between {startdate} - {enddate}")
    print(df_goes)

    goes = df_goes['Satellite'].to_list()
    starttimes = df_goes['Start Date'].to_list()
    endtimes = df_goes['End Date'].to_list()
    
    return starttimes, endtimes, goes
    



def check_all_goes_data(startdate, enddate, experiment, flux_type, spacecraft="primary"):
    """ If the user specifies only "GOES" for the experiment, then
        read in the primary/seconday GOES experiments and string them together.
        This will allow the creation of a long flux timeseries
        with multiple GOES satellites included.
        
        Data are provided for the primary or secondary spacecraft, as
        specified.
        
        The first spacecraft encountered with data for the requested
        time period is selected.
        
        ONLY ALLOWED FOR INTEGRAL FLUXES, OTHERWISE RUN INTO PROBLEM
        OF MISMATCHING ENERGY BINS.
        
    """
    #Identify the primary or secondary GOES spacecraft during the time
    #period of interest
    starttimes, endtimes, goes = identify_which_goes_spacecraft(startdate,
            enddate, spacecraft=spacecraft)

    filenames1_all = []
    filenames2_all = []
    filenames_orien_all = []
    detector = []
    
    #Write dates and experiment used for those dates
    #This is meant to allow the user to know which
    #spacecraft is available for which time period
    outfname = os.path.join(outpath, "idsep", "goes_experiments_dates.txt")
    outfile = open(outfname,'w')

    #Check for and download (if needed) the data for the primary/secondary
    #spacecraft on a monthly basis.
    for startdate, enddate, goes_i in zip(starttimes, endtimes, goes):
        filenames1 = []
        filenames2 = []
        filenames_orien = []
        
        styear = startdate.year
        stmonth = startdate.month
        stday = startdate.day
        endyear = enddate.year
        endmonth = enddate.month
        endday = enddate.day
        
        print(f"======Date: {startdate} to {enddate}")

        if goes_i in old_goes_sc:
            filenames1, filenames2, date = \
                 check_old_goes_data(startdate, enddate, goes_i, flux_type)

        if goes_i in goes_sc:
            filenames1, filenames2, filenames_orien, date = \
                 check_goes_data(startdate, enddate, goes_i, flux_type)

        if (goes_i in goes_R) and flux_type == "differential":
            filenames1, filenames2, filenames_orien, date = \
             check_goesR_data(startdate, enddate, goes_i, flux_type, spacecraft=spacecraft)

        if (goes_i=="GOES_RT" or goes_i in goes_R) and flux_type == "integral":
            filenames1, filenames2, filenames_orien, date = check_goes_RTdata(startdate,
                enddate, goes_i, flux_type, spacecraft=spacecraft)
            goes_i = "GOES_RT"

        if len(filenames1) != 0 and filenames1 != None:
            filenames1_all.extend(filenames1)
            detector.extend([goes_i]*len(filenames1))
            if len(filenames2) == 0:
                filenames2_all.extend([None]*len(filenames1))
            else:
                filenames2_all.extend(filenames2)
            
            if len(filenames_orien) == 0:
                filenames_orien_all.extend([None]*len(filenames1))
            else:
                filenames_orien_all.extend(filenames_orien)


            #GOES data is stored in monthly data files
            #Write out goes selections to output file
            get_years = []
            get_months = []
            test_year = styear
            test_month = stmonth
            test_date = datetime.datetime(year=test_year, month=test_month, day=1)
            while (test_date <= enddate):
                get_years.append(test_year)
                get_months.append(test_month)
                test_month = test_month + 1
                if (test_month > 12):
                    test_month = 1
                    test_year = test_year + 1
                test_date = datetime.datetime(year=test_year, month=test_month, day=1)
        
            #Write out the year, month, GOES spacecraft
            for k in range(len(get_months)):
                year = get_years[k]
                month = get_months[k]
                outfile.write(str(year)+","+str(month)+"," +goes_i+"\n")

    
    outfile.close()
    
    return filenames1_all, filenames2_all, filenames_orien_all, detector


    

def check_ephin_data(startdate, enddate, experiment, flux_type):
    """Check for SOHO/COSTEP/EPHIN data on your computer. If not there,
        download from http://ulysses.physik.uni-kiel.de/costep/level3/l3i/
        5 minute data will be downloaded. Intensities are in units of
        (cm^2 s sr mev/nuc)^-1
        First available date is 1995 12 8 (DOY = 342).
        The files are available in daily or yearly format.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SOHO
            EPHIN Level 3 data that span the desired time range
            (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, EPHIN, eps, or epead

    Nyr = endyear - styear + 1
    for year in range(styear, endyear+1):
        fname = str(year) + '.l3i'
        res = '5min'
        
        svfile = os.path.join(cfg.datapath,'EPHIN',fname)
        exists = os.path.isfile(svfile)
        
        complete = False
        if exists:
            #Check if the file is listed as complete
            complete = check_completeness(experiment, flux_type, svfile, df=df)
        
        if not exists or not complete: #download file if not found on your computer
            url = ('http://ulysses.physik.uni-kiel.de/costep/level3/l3i/%s/%s'
                    % (res,fname))
            print('Downloading EPHIN data: ' + url)
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(svfile):
                    os.remove(svfile) # if exist, remove it directly
  
                wget.download(url, svfile)
            except urllib.request.HTTPError:
                sys.exit("Cannot access EPHIN file at " + url +
               ". Please check that selected spacecraft covers date range.")
               
        filenames1.append(os.path.join('EPHIN', fname))
        
    return filenames1


def check_ephin_release_data(startdate, enddate, experiment, flux_type):
    """Check for SOHO/COSTEP/EPHIN data on your computer provided by the
        HESPERIA collaboration on the website
        https://www.hesperia.astro.noa.gr/index.php/results/real-time-prediction-tools/data-retrieval-tool
        (one long URL). Monthly files of this data set are provided in the
        public git containing these codes. Otherwise, users may go to the URL
        above, download the data sets, and use the file naming convention
        adopted here:
            
        data/EPHIN_REleASE/HESPERIA_SOHO_PROTON_YYYY.txt
        
        Intensities are in units of (cm^2 s sr mev/nuc)^-1
        First available date is 1995 12 8 (DOY = 342).
        The files are saved in yearly format.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SOHO
            EPHIN data from the REleASE website that span the desired
            time range (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, EPHIN, eps, or epead

    Nyr = endyear - styear + 1
    for year in range(styear, endyear+1):
        fname = 'HESPERIA_SOHO_PROTON_' + str(year) + '.txt'

        exists = os.path.isfile(datapath + '/EPHIN_REleASE/' + fname)
        if not exists: #download file if not found on your computer
                sys.exit("Cannot access EPHIN file " + fname +
               ". Please check that selected spacecraft covers date range.")
       
        filenames1.append('EPHIN_REleASE/' + fname)
        
    return filenames1


def check_erne_data(startdate, enddate, experiment, flux_type):
    """Check for SOHO/ERNE data on your computer. If not there,
        download from https://export.srl.utu.fi.
        Intensities are in units of (cm^2 s sr mev/nuc)^-1
        First available date is 1996 05 07.
        The files are available in randomly grouped dates.
        
        .dates files on the site indicate which dates are in each
        .tgz data file. Will download all lightweight .dates files
        to find the correct data files to download.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SOHO
            EPHIN Level 3 data that span the desired time range
            (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    ernepath = os.path.join(datapath,'ERNE')
    
    #Download all the .dates files each time
    #Should always work if user has an internet connection,
    #but allow for case where user doesn't have internet access.
    #Don't want to stop workflow if user has files on the computer.
    allfiles = os.listdir(ernepath)
    datefiles = [x for x in allfiles if 'dates' in x]
    try:
        r = requests.get('https://export.srl.utu.fi')
        soup = BeautifulSoup(r.text, 'html.parser')
        links = soup.find_all('a')
        checkfiles = [x.text for x in links if 'dates' in x.text]
        for fl in checkfiles:
            if fl in datefiles:
                continue #don't need to download again
            url = 'https://export.srl.utu.fi/%s' % (fl)
            localfl = os.path.join(ernepath,fl)
            try:
                urllib.request.urlopen(url)
                wget.download(url, localfl)
                datefiles.append(fl) #add to datefiles if not present
            except:
                print('Cannot download ' + url)
    except:
        pass


    #Array of filenames that contain the data requested by the User
    filenames1 = []

    for file in datefiles:
        fnm = file.replace('.dates','.tgz') #data filename only
        localfnm = os.path.join(ernepath,fnm) #data full path
        with open(os.path.join(ernepath,file),'r') as datefile:
            lines = datefile.readlines()
            firstdate = datetime.datetime.strptime(lines[0].strip(),"%Y.%m.%d")
            lastdate = datetime.datetime.strptime(lines[-1].strip(),"%Y.%m.%d")
            
            #Subdirectory that will hold all the files that come out of
            #the ERNE data gzip file
            datedir = lines[0].strip() + '_' + lines[-1].strip()
            unzipdir = os.path.join(ernepath,datedir)

            #Are the user requested dates contained within this file?
            if (startdate >= firstdate and startdate < lastdate) or \
                (enddate > firstdate and enddate <= lastdate) or \
                (firstdate >= startdate and firstdate <= enddate):

                exists = os.path.isfile(localfnm)

                if not exists: #download file if not found on your computer
                    url = ('https://export.srl.utu.fi/%s' % (fnm))
                    print('Downloading ERNE data: ' + url)
                    try:
                        urllib.request.urlopen(url)
                        wget.download(url, localfnm)
                        #Decompress gzip file
                        tar = tarfile.open(name=localfnm,mode='r:gz')
                        tar.extractall(path=unzipdir)
                        
                    except urllib.request.HTTPError:
                        sys.exit("Cannot access EPHIN file at " + url +
                       ". Please check that selected spacecraft covers date range.")

                dfiles = os.listdir(os.path.join(unzipdir,'export.src'))
                dfiles = [os.path.join(unzipdir,'export.src',x) for x in dfiles if (('HED' in x or 'LED' in x) and ('SL2' in x))]
                #HED, LET files
                filenames1.extend(dfiles) #save filenames

    return filenames1



def check_stereo_data(startdate, enddate, experiment, flux_type):
    """ Check for 1 min STEREO data on your computer provided on:
        https://izw1.caltech.edu/STEREO/Public/LET_public.html
        https://izw1.caltech.edu/STEREO/Public/HET_public.html            
        
        
        Intensities are in units of (cm^2 s sr mev/nuc)^-1
        First available date is 2006 12 1.
        The files are saved in daily files.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the data
            that span the desired time range (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = [] #LET  
    filenames2 = [] #HET

    td = enddate - startdate
    #LET STEREO data is stored in daily data files
    NFILESd = td.days + 1
    #HET STEREO data is stored in monthly data files
    styear = startdate.year
    stmonth = startdate.month
    endyear = enddate.year
    endmonth = enddate.month
    numyr = endyear - styear
    if numyr == 0:
        NFILESm = endmonth - stmonth + 1
    if numyr > 0:
        NFILESm = (12 - stmonth) + endmonth + (numyr-1)*12 + 1


    #pulls primary spacecraft fluxes
    if experiment == 'STEREO-A':
        let_url_prefix = 'https://izw1.caltech.edu/STEREO/DATA/Level1/Public/ahead/1Minute/'
        let_prefix = 'H_summed_ahead_'
        het_url_prefix = 'https://izw1.caltech.edu/STEREO/DATA/HET/Ahead/1minute/'
        het_prefix = 'AeH'

    if experiment == 'STEREO-B':
        let_url_prefix = 'https://izw1.caltech.edu/STEREO/DATA/Level1/Public/behind/1Minute/'
        let_prefix = 'H_summed_behind_'
        het_url_prefix = 'https://izw1.caltech.edu/STEREO/DATA/HET/Behind/1minute/'
        het_prefix = 'BeH'


    #for every day that data is required, check if file is present or
    #needs to be downloaded.
    #LET
    #https://izw1.caltech.edu/STEREO/DATA/Level1/Public/behind/1Minute/2006/Summed/H/H_summed_behind_2006_317_level1_11.txt
    #HET
    #https://izw1.caltech.edu/STEREO/DATA/HET/Behind/1minute/BeH06Dec.1m
 
    #LET data
    for i in range(NFILESd):
        date = startdate + datetime.timedelta(days=i)
        year = date.year
        month = date.month
        day = date.day
        doy = date.timetuple().tm_yday
        strmonth = date.strftime("%b")
        stryr = str(year)

        #LET
        fname1 = '%s%i_%03i_level1_11.txt' % (let_prefix,year,doy)
        fullpath1 = os.path.join(cfg.datapath,experiment,'LET',fname1)
        exists1 = os.path.isfile(fullpath1)

        complete = False
        if exists1:
            #Check if the file is listed as complete
            complete = check_completeness(experiment, flux_type, fullpath1, df=df)
            if complete:
                filenames1.append(os.path.join(experiment,'LET',fname1))

        if not exists1 or not complete:
            url=(let_url_prefix + '%i/Summed/H/%s' % (year,fname1))
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(fullpath1):
                    os.remove(fullpath1) # if exist, remove it directly
  
                wget.download(url, fullpath1)
                print("Downloaded file --> " + fullpath1)
                filenames1.append(os.path.join(experiment,'LET',fname1))
            except urllib.request.HTTPError:
                sys.exit("Cannot access " + experiment + " file at " + url +
                ". Please check that selected spacecraft covers date range.")


    #HET data - monthly
    year = startdate.year
    month = startdate.month
    nyr = 0
    for i in range(NFILESm):
        if month == 13:
            month = 1
            nyr = nyr + 1
            
        year = styear + nyr
        date = datetime.datetime(year=year,month=month,day=1)
        strmonth = date.strftime("%b")
        stryr = str(year)

        #HET
        fname2 = '%s%s%s.1m' % (het_prefix,stryr[2:4],strmonth)
        fullpath2 = os.path.join(cfg.datapath,experiment,'HET',fname2)
        exists2 = os.path.isfile(fullpath2)

        complete = False
        if exists2:
            #Check if the file is listed as complete
            complete = check_completeness(experiment, flux_type, fullpath2, df=df)
            if complete:
                filenames2.append(os.path.join(experiment,'HET',fname2))

        if not exists2 or not complete:
            url=het_url_prefix + fname2
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(fullpath2):
                    os.remove(fullpath2) # if exist, remove it directly
  
                wget.download(url,fullpath2)
                print("Downloaded file --> " + fullpath2)
                filenames2.append(os.path.join(experiment,'HET',fname2))
            except urllib.request.HTTPError:
                sys.exit("Cannot access " + experiment + " file at " + url +
                    ". Please check that selected spacecraft covers date range.")

        month += 1

    return filenames1, filenames2



def check_ace_sis_data(startdate, enddate, experiment, flux_type):
    """Check for ACE/SIS data on your computer. If not there,
        download from https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/
        5 minute data will be downloaded. Only >30 and >60 MeV fluxes.
        Intensities are in units of p/cs2-sec-ster
        First available date is 2001-08-07.
        The files are available in daily format.
        e.g. 20010807_ace_sis_5m.txt
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SOHO
            EPHIN Level 3 data that span the desired time range
            (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    startdt = datetime.datetime(year=styear,month=stmonth,day=stday)
    enddt = datetime.datetime(year=endyear,month=endmonth,day=endday)
    Ndays = int((enddt - startdt)/datetime.timedelta(hours=24)) + 1

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []

    for i in range(Ndays):
        getday = startdt + i*datetime.timedelta(hours=24)
        #20010807_ace_sis_5m.txt
        fname = getday.strftime("%Y%m%d") + "_ace_sis_5m.txt"
        
        svfile = os.path.join(cfg.datapath,'ACE','SIS',fname)
        exists = os.path.isfile(svfile)
        
        complete = False
        if exists:
            #Check if the file is listed as complete
            complete = check_completeness(experiment, flux_type, svfile, df=df)
        
        if not exists or not complete: #download file if not found on your computer
            url = ('https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/%s'
                    % (fname))
            print('Downloading ACE/SIS integral data: ' + url)
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(svfile):
                    os.remove(svfile) # if exist, remove it directly
  
                wget.download(url, svfile)
            except urllib.request.HTTPError:
                sys.exit("Cannot access ACE/SIS file at " + url +
               ". Please check that selected spacecraft covers date range.")
               
        filenames1.append(os.path.join('ACE', 'SIS', fname))
        
    return filenames1


def check_ace_epam_electrons_data(startdate, enddate, experiment, flux_type):
    """Check for ACE/EPAM data on your computer. If not there,
        download from https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/
        5 minute data will be downloaded. Electron and proton data. Only take
        the energetic electrons.
        # Units: Differential Flux particles/cm2-s-ster-MeV
        First available date is 2001-08-07.
        The files are available in daily format.
        e.g. 20010807_ace_epam_5m.txt
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SOHO
            EPHIN Level 3 data that span the desired time range
            (yearly files)
        
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    startdt = datetime.datetime(year=styear,month=stmonth,day=stday)
    enddt = datetime.datetime(year=endyear,month=endmonth,day=endday)
    Ndays = int((enddt - startdt)/datetime.timedelta(hours=24)) + 1

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []

    for i in range(Ndays):
        getday = startdt + i*datetime.timedelta(hours=24)
        #20010807_ace_sis_5m.txt
        fname = getday.strftime("%Y%m%d") + "_ace_epam_5m.txt"
        
        svfile = os.path.join(cfg.datapath,'ACE','EPAM',fname)
        exists = os.path.isfile(svfile)
        
        complete = False
        if exists:
            #Check if the file is listed as complete
            complete = check_completeness(experiment, flux_type, svfile, df=df)
        
        if not exists or not complete: #download file if not found on your computer
            url = ('https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/%s'
                    % (fname))
            print('Downloading ACE/EPAM data: ' + url)
            try:
                urllib.request.urlopen(url)
                
                if os.path.exists(svfile):
                    os.remove(svfile) # if exist, remove it directly
  
                wget.download(url, svfile)
            except urllib.request.HTTPError:
                sys.exit("Cannot access ACE/EPAM file at " + url +
               ". Please check that selected spacecraft covers date range.")
               
        filenames1.append(os.path.join('ACE', 'EPAM', fname))
        
    return filenames1



def check_imp8_cpme_data(startdate, enddate, experiment, flux_type):
    """Check for IMP-8/CPME data on your computer. If not there,
        download from http://sd-www.jhuapl.edu/IMP/data/imp8/cpme/cpme_330s/protons/
        330s data will be downloaded. 
        IMP-8 CPME: 330-sec. Avg. Proton Intensities & Uncertainties [no./(cm^2-sc-ster-MeV)] 
        The files are compiled in regular doy groupings.
        e.g. h_330s_1989_240_263.txt.gz
        
        First data available: h_330s_1974_048_071.txt.gz   (Feb 17, 1974) 
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the SOHO
            EPHIN Level 3 data that span the desired time range
            (yearly files)
        
    """
    doy_groupings = [[1,23],[24,47],[48,71],[72,95],[96,119],[120,143],
                    [144,167],[168,191],[192,215],[216,239],[240,263],
                    [264,287],[288,311],[312,336],[336,365]]

    leap_year_ref = 1980

    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    stdoy = startdate.timetuple().tm_yday
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day
    enddoy = enddate.timetuple().tm_yday
    
    Nyr = (endyear - styear) + 1

    df = read_data_manager() #file completeness record

    #Array of filenames that contain the data requested by the User
    filenames1 = []

    for i in range(Nyr):
        year = styear + i
        is_leap_year = ((year - leap_year_ref)%4 == 0)

        stix = 0
        endix = 0
        if year != styear:
            stix = 0
        if year != endyear:
            endix = len(doy_groupings) - 1
        if year == styear:
            for j,grp in enumerate(doy_groupings):
                if (stdoy >= grp[0]) and (stdoy <= grp[1]):
                    stix = j
        if year == endyear:
            for j,grp in enumerate(doy_groupings):
                if (enddoy >= grp[0]) and (enddoy <= grp[1]):
                    endix = j

        for k in range(stix,endix+1):
            grp = doy_groupings[k]
            dy1 = grp[0]
            dy2 = grp[1]
            if is_leap_year and dy2 == 365: dy2 = 366
            #h_330s_1974_048_071.txt.gz
            fname = f"h_330s_{year}_{dy1:03d}_{dy2:03d}.txt"
            gzfname = f"{fname}.gz"
        
            gzsvfile = os.path.join(cfg.datapath,'IMP8','CPME',gzfname)
            svfile = os.path.join(cfg.datapath,'IMP8','CPME',fname)
            exists = os.path.isfile(svfile)
        
            complete = False
            if exists:
                #Check if the file is listed as complete
                complete = check_completeness(experiment, flux_type, svfile, df=df)
            
            if not exists or not complete: #download file if not found on your computer
                url = ('http://sd-www.jhuapl.edu/IMP/data/imp8/cpme/cpme_330s/protons/%s/%s'
                        % (year, gzfname))
                print('Downloading IMP-8/CPME data: ' + url)
                try:
                    urllib.request.urlopen(url)
                    
                    if os.path.exists(gzsvfile):
                        os.remove(gzsvfile) # if exist, remove it directly
                    if os.path.exists(svfile):
                        os.remove(svfile) # if exist, remove it directly
      
                    wget.download(url, gzsvfile)
                    #Decompress gzip file
                    # Open the gzipped file in binary read mode ('rb')
                    with gzip.open(gzsvfile, 'rb') as f_in:
                        # Open the output file in binary write mode ('wb')
                        with open(svfile, 'wb') as f_out:
                            # Copy the decompressed data from the input to the output file
                            shutil.copyfileobj(f_in, f_out)
                    
                except urllib.request.HTTPError:
                    sys.exit("Cannot access IMP-8/CPME file at " + url +
                   ". Please check that selected spacecraft covers date range.")
                   
            filenames1.append(svfile)
        
    return filenames1




def check_data(startdate, enddate, experiment, flux_type, user_file,
    spacecraft="primary"):
    """Check that the files containing the data are in the data directory. If
        the files for the requested dates aren't present, they will be
        downloaded from the NOAA website. For SEPEM (RSDv2) data, if missing,
        the program prints the URL from which the data can be downloaded and
        unzipped manually.
        The RSDv2 data set is very large and takes a long time to read as a
        single file. This program will generate files containing fluxes for
        each year for faster reading.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        :user_file: (string) name of file containing user-input data
            (if applicable)
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the data that
            span the desired time range (monthly files)
        :filenames2: (string array) if GOES, files containing HEPAD data
        :filenames_orien: (string array if GOES, files containing
            satellite orientation
        
    """
    print('Checking that the requested data is present on your computer.')
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    user_fname = [user_file]

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, eps, or epead
    filenames2 = []  #hepad
    filenames_orien = []  #orientation flag for GOES-13+

    #If user wants to use own input file (filename defined as input)
    if experiment == "user":
        nuser = len(user_fname)
        for i in range(nuser):
            exists = os.path.isfile(user_fname[i])
            if exists:
                filenames1.append(user_fname[i])
            if not exists:
                sys.exit("You have selected to read a user input "
                "file with filename " + user_fname[i]
                + ". This file is not found! Exiting.")

        return filenames1, filenames2, filenames_orien

    #SEPEM data set is continuous, but too long; prefer yearly files
    #Check if user has yearly files; if not:
        #check if user has original SEPEM, then create yearly files
        #Otherwise alert user to download data set and try again
    if (experiment == "SEPEM" or experiment == "SEPEMv3"):
        filenames1 = check_sepem_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien
        
    if experiment == "CalGOES":
        filenames1 = check_calgoes_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    #Try all GOES experiments
    if experiment == "GOES":
        filenames1, filenames2, filenames_orien, detector = \
         check_all_goes_data(startdate, enddate, experiment, flux_type,
         spacecraft=spacecraft)
        return filenames1, filenames2, filenames_orien, detector

    #Specific GOES-07 and previous
    if experiment in old_goes_sc:
        filenames1, filenames2, date = check_old_goes_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    #Specific GOES prior to GOES-R, e.g. GOES-13 or GOES-08
    if experiment in goes_sc:
        filenames1, filenames2, filenames_orien, date =\
            check_goes_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    #FILE FOR 2017-09 SEP events, but have to have this file.
 #   if (experiment == "GOES-16" or experiment == "GOES-17") and flux_type == "differential"\
 #       and enddate < datetime.datetime(2020,12,1):
 #       filenames1 = ['GOES-R/' + \
 #                   "se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc"]
 #       filenames2 = []
 #       filenames_orien =[]
 #       return filenames1, filenames2, filenames_orien


    if (experiment in goes_R) and flux_type == "differential":
        filenames1, filenames2, filenames_orien, date =\
            check_goesR_data(startdate,enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien
        
    if (experiment == "GOES_RT") and flux_type == "integral":
        filenames1, filenames2, filenames_orien, date =\
            check_goes_RTdata(startdate,enddate, experiment, flux_type, spacecraft=spacecraft)
        return filenames1, filenames2, filenames_orien
        

    if experiment == "EPHIN":
        filenames1 = check_ephin_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment == "EPHIN_REleASE":
        filenames1 = check_ephin_release_data(startdate, enddate,
            experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment == "ERNE":
        filenames1 = check_erne_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if 'STEREO' in experiment:
        filenames1, filenames2 = check_stereo_data(startdate,
            enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment == "ACE_SIS":
        filenames1 = check_ace_sis_data(startdate,
            enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment == "ACE_EPAM_electrons":
        filenames1 = check_ace_epam_electrons_data(startdate,
            enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment == "IMP8_CPME":
        filenames1 = check_imp8_cpme_data(startdate,
            enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien


    return filenames1, filenames2, filenames_orien


def find_goes_data_dimensions(filename):
    """ Input open csv file of GOES data. Identifies the start of the data by
        searching for the string 'data:', then returns the number of header
        rows and data rows present in the file.
        
        INPUTS:
        
        :filename: (string) name of GOES file containing data
        
        OUTPUTS:
        
        :nhead: (int) number of header rows
        :nrow: (int) number of rows of data
        
    """
    with open(datapath + '/' + filename) as csvfile:
        #GOES data has very large headers; figure out where the data
        #starts inside the file and skip the required number of lines
        nhead = 0
        for line in csvfile:
            nhead = nhead + 1
            if 'data:' in line: #location of line before column headers
                break
        nhead = nhead + 1 #account for column header line
        csvfile.readline() #proceed to column header line
        readCSV = csv.reader(csvfile, delimiter=',')
        nrow = len(csvfile.readlines()) #count lines of data
        #print('\nThere are ' + str(nhead) + ' header lines and '
        #    + str(nrow) + ' rows of data in ' + filename)

    csvfile.close()
    return nhead, nrow


def get_west_detector(filename, dates):
    """ For GOES-13+, identify which detector is facing west from the
        orientation flag files. Get an orientation for each data point.
        EPEAD orientation flag. 0: A/W faces East and B/E faces West.
        1: A/W faces West and B/E faces East. 2: yaw-flip in progress.
        
        INPUTS:
        
        :filename: (string) file containing satellite orientation
        :dates: (datetime 1xn array) n dates during user specified
            time period
            
        OUTPUTS:
        
        :west_detector: (string 1xn array) detector (A or B) identified
            as facing westward for each time point
       
    """
    nhead, nrow = find_goes_data_dimensions(filename)
    orien_dates = []
    orientation = []

    #Per communication with Juan Rodriguez, GOES-13 was always in the same
    #orientation afer 2010.
    #The orientation file is missing for 20130501_20130531. Per the previous
    #month's orientation file, B was the westward-facing detector.
    if dates[0] == datetime.datetime(2013,5,1) and dates[-1] == datetime.datetime(2013,5,31,23,55,0):
        return ["B"]*len(dates)

    with open(datapath + '/' + filename) as orienfile:
        #GOES data has very large headers; figure out where the data
        #starts inside the file and skip the required number of lines
        readCSV = csv.reader(orienfile, delimiter=',')
        for k in range(nhead):
            next(readCSV)  #to line with column headers

        for row in readCSV:
            date = datetime.datetime.strptime(row[0][0:19],
                                            "%Y-%m-%d %H:%M:%S")
            orien_dates.append(date)
            orientation.append(float(row[1]))

    orienfile.close()

    #orientation data is in 1 minute intervals while flux data is in 5
    #minute intervals. Identify W facing detector for flux times.
    #Assume that time stamps match every 5 minutes.
    date_index = 0
    west_detector = []
    for i in range(nrow):
        if orien_dates[i] == dates[date_index]:
            if orientation[i] == 0:
                west_detector.append("B")
            if orientation[i] == 1:
                west_detector.append("A") 
            if orientation[i] == 2:
                west_detector.append("Flip")
            date_index = date_index + 1
            if date_index == len(dates):
                break

#    print('There were ' + str(len(dates)) + ' input dates and there are ' +
#            str(len(west_detector)) + ' detector orientations. '
#            f"A: {west_detector.count('A')}, B: {west_detector.count('B')}")
    return west_detector


def read_in_sepem(experiment, flux_type, filenames1):
    """ Read in SEPEM data files from the computer.
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) "differential"
        :filenames1: (string array) names of files containing
            SEPEM data in desired time range (yearly files)
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
    """
    NFILES = len(filenames1)
    for i in range(NFILES):
        print('Reading in file ' + datapath + '/' + filenames1[i])
        with open(datapath + '/' + filenames1[i]) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            has_header = csv.Sniffer().has_header(csvfile.readline())
            if has_header:
                next(readCSV)  # Skip single header row.
            ncol = len(next(readCSV))
            csvfile.seek(0) #back to beginning of file
            if has_header:
                next(readCSV)  # Skip header row.
            nrow = len(csvfile.readlines())
            #print('There are ' + str(ncol) + ' columns and ' + str(nrow) +
            #    ' rows of data in ' + filenames1[i])

            #Define arrays that hold dates and fluxes
            dates = []
            fluxes = np.zeros(shape=(ncol-1,nrow))

            csvfile.seek(0) #back to beginning of file
            if has_header:
                next(readCSV)  # Skip header row.

            count = 0
            for row in readCSV:
                date = datetime.datetime.strptime(row[0][0:19],
                                                "%Y-%m-%d %H:%M:%S")
                dates.append(date)
                for j in range(1,ncol):
                    flux = float(row[j])
                    if flux < 0:
                        flux = badval
                    fluxes[j-1][count] = flux
                count = count + 1
        #If reading in multiple files, then combine all data into one array
        #SEPEM currently only has one file, but making generalized
        if i==0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates

    return all_dates, all_fluxes


def read_in_calgoes(experiment, filenames1):
    """ Read in the data set created by Shaowen Hu (NASA JSC SRAG)
        that generates a consistent and recalibrated GOES data set.
    """
    
    fluxes = []
    all_dates = []

    for filenm in filenames1:
        print("reading filename " + filenm)
        with open(datapath + '/' + filenm) as file:
            if len(fluxes) == 0:
                #Get number of columns in file
                #Fill in first row so can have correct format array
                firstline = file.readline().strip().split(",")
                ncol = len(firstline) - 1
                for i in range(ncol):
                    fluxes.append([float(firstline[i+1])])
                date = datetime.datetime.strptime(firstline[0][0:19],"%Y-%m-%d %H:%M:%S")
                all_dates.append(date)

            for line in file:
                if line == "": continue
                
                row = line.strip().split(",")
                if row[0] == "0": continue
                for j in range(ncol):
                    fluxes[j].append(float(row[j+1]))
                
                date = datetime.datetime.strptime(row[0][0:19],"%Y-%m-%d %H:%M:%S")
                all_dates.append(date)

    all_fluxes = np.array(fluxes)
    
    return all_dates, all_fluxes



def read_in_old_goes(experiment, flux_type, filenames1, filenames2, options):
    """Read in GOES data from your computer. User may specify option to choose
        corrected or uncorrected GOES fluxes.
        For GOES-07 and previous which did not have HEPAD.
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []
    west_detector = [] #place holder, will be filled if needed

    #Read in file that identified data files as complete
    df = read_data_manager()

    if (experiment == "GOES-05"):
        if flux_type == "differential":
            #NO CORRECTED FLUXES PRODUCED
            print("Note that GOES-05 NCEI data files contain only UNCORRECTED fluxes.")
            columns = [5,7,9,11,13,15] #eps uncorrected
            hepad_columns = []
        if flux_type == "integral":
            #NO INTEGRAL CHANNELS PRODUCED
            print("GOES-05 files do not contain integral proton fluxes.")
            columns = [] #eps
            hepad_columns = []


    if (experiment == "GOES-06"):
        if flux_type == "differential":
            if "corrected" in options or "uncorrected" not in options:
                #CORRECTED CHANNELS; DEFAULT
                columns = [16,17,18,19,20,21] #eps
                hepad_columns = [1,2,3,4]
            if "uncorrected" in options:
                #UNCORRECTED channels
                columns = [3,4,5,6,7,8] #eps
                hepad_columns = [1,2,3,4]
        if flux_type == "integral":
            #ONLY CORRECTED CHANNELS AVAILABLE
            columns = [23,24,25,26,27,28] #eps
            hepad_columns = [4]

    if (experiment == "GOES-07"):
        if flux_type == "differential":
            if "corrected" in options or "uncorrected" not in options:
                #CORRECTED CHANNELS; DEFAULT
                columns = [16,17,18,19,20,21] #eps
                hepad_columns = []
            if "uncorrected" in options:
                #UNCORRECTED channels
                columns = [3,4,5,6,7,8] #eps
                hepad_columns = []
        if flux_type == "integral":
            #ONLY CORRECTED CHANNELS AVAILABLE
            columns = [23,24,25,26,27,28] #eps
            hepad_columns = []



    ncol = len(columns)
    nhcol = len(hepad_columns)
    totcol = ncol + nhcol

    if totcol == 0:
        return all_dates, all_fluxes, west_detector

    #Read in fluxes from files
    fluxes = []
    for i in range(NFILES):
        file_dates = []
        dates = []
        if filenames1[i] != None:
            #FIRST set of files for lower energy eps or epead
            nhead, nrow = find_goes_data_dimensions(filenames1[i])
            fluxes = np.zeros(shape=(totcol,nrow))
            
            
            print('Reading in file ' + datapath + '/' + filenames1[i])
            fullpath1 = os.path.join(cfg.datapath, filenames1[i])
            with open(fullpath1) as csvfile:
                #GOES data has very large headers; figure out where the data
                #starts inside the file and skip the required number of lines
                readCSV = csv.reader(csvfile, delimiter=',')
                for k in range(nhead):
                    next(readCSV)  #to start of data

                #Get dates first; need dates to identify spacecraft orientation
                #for GOES-13+
                for row in readCSV:
                    date = datetime.datetime.strptime(row[0][0:19],
                                                    "%Y-%m-%d %H:%M:%S")
                    dates.append(date)
                    file_dates.append(date) #only dates in current file


                #Go back and get fluxes
                count = 0
                csvfile.seek(0)
                for k in range(nhead):
                    next(readCSV)  #to start of data
                for row in readCSV:
                    for j in range(ncol):
                        flux = float(row[columns[j]])
                        if flux < 0:
                            flux = badval
                        fluxes[j][count] = flux
                    count = count + 1
            csvfile.close()

            #Update file completeness
            df = file_completeness(df, experiment, flux_type, fullpath1, file_dates)


        #SECOND set of files for higher energy hepad
        if filenames2:
            if filenames2[i] != None:
                nhead, nrow = find_goes_data_dimensions(filenames2[i])
                fullpath2 = os.path.join(cfg.datapath,filenames2[i])
                print('Reading in file ' + fullpath2)
                with open(fullpath2) as csvfile:
                    readCSV = csv.reader(csvfile, delimiter=',')
                    for k in range(nhead):
                        next(readCSV)  #to start of data

                    count = 0
                    for row in readCSV:
                        for j in range(nhcol):
                            #Sometimes GOES datafiles have an incomplete row,
                            #e.g. g13_hepad_ap_5m_20101201_20101231.csv
                            #If so, rather than crash, set the flux to badval
                            try:
                                flux = float(row[hepad_columns[j]])
                            except:
                                flux = badval
                            if flux < 0:
                                flux = badval
                            fluxes[ncol+j][count] = flux
                        count = count + 1
                csvfile.close()

                #Update file completeness
                df = file_completeness(df, experiment, flux_type, fullpath2, file_dates)

        if len(fluxes) == 0: continue
        #If reading in multiple files, then combine all data into one array
        if len(all_fluxes) == 0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates
    
    if all_dates == []:
        print("read_in_old_goes: Did not find the data you were looking for.")
   
    write_data_manager(df)
    
    return all_dates, all_fluxes, west_detector




def read_in_goes(experiment, flux_type, filenames1, filenames2,
                filenames_orien, options):
    """Read in GOES data from your computer. User may specify option to choose
        corrected or uncorrected GOES fluxes.
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    #Have encountered some bad files and better to use a different
    #instrument
    
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []
    west_detector = [] #place holder, will be filled if needed

    #Read in file that identified data files as complete
    df = read_data_manager()

    if (experiment == "GOES-08" or experiment == "GOES-09" or
        experiment == "GOES-10" or experiment == "GOES-11" or
        experiment == "GOES-12"):
        if flux_type == "differential":
            if "corrected" in options or "uncorrected" not in options:
                #CORRECTED CHANNELS; DEFAULT
                columns = [18,19,20,21,22,23] #eps
                hepad_columns = [1,2,3,4]
            if "uncorrected" in options:
                #UNCORRECTED channels
                columns = [5,6,7,8,9,10] #eps
                hepad_columns = [1,2,3,4]
        if flux_type == "integral":
            #ONLY CORRECTED CHANNELS AVAILABLE
            columns = [25,26,27,28,29,30] #eps
            hepad_columns = [4]
   
    if (experiment == "GOES-13" or experiment == "GOES-14" or
        experiment == "GOES-15"):
        if flux_type == "differential":
            if "corrected" in options or "uncorrected" not in options:
                #CORRECTED CHANNELS
                columns = [16,24,32,40,48,56] #epead, default A detector
                columnsB = [12,20,28,36,44,52] #epead, B detector
                hepad_columns = [9,12,15,18]
            if "uncorrected" in options:
                #UNCORRECTED CHANNELS
                columns = [15,23,31,39,47,55] #epead, default A detector
                columnsB = [11,19,27,35,43,51] #epead, B detector
                hepad_columns = [9,12,15,18]
        if flux_type == "integral":
            #ONLY CORRECTED CHANNELS AVAILABLE
            columns = [18,20,22,24,26,28] #epead, default A detector
            columnsB = [4,6,8,10,12,14] #epead, B detector
            hepad_columns = [18]


    ncol = len(columns)
    nhcol = len(hepad_columns)
    totcol = ncol + nhcol

    #Read in fluxes from files
    for i in range(NFILES):
        file_dates = []
        dates = []
        fluxes = []
        if filenames1[i] != None:
            #FIRST set of files for lower energy eps or epead
            nhead, nrow = find_goes_data_dimensions(filenames1[i])
            fluxes = np.zeros(shape=(totcol,nrow))
            
            #GOES-14 and GOES-15 files between 2019-09-01 to 2020-03-07
            #have one column missing in the hepad file. Check the date
            #in the filenames1 and apply the correct hepad columns for
            #those specific files.
            #filename1 e.g. 15_epead_cpflux_5m_20200301_20200331.csv
            if (experiment == "GOES-14" or experiment == "GOES-15"):
                fsplit = filenames1[i].split("/") #break up path
                fnm = fsplit[-1].split(".csv") #last entry is the filename
                fnm = fnm[0].split("_")
                dt = datetime.datetime.strptime(fnm[-2],"%Y%m%d")
                if dt >= datetime.datetime(year=2019,month=9,day=1):
                    if flux_type == "differential":
                        hepad_columns = [8,11,14,17]
                    if flux_type == "integral":
                        hepad_columns = [17]
                
            
            print('Reading in file ' + datapath + '/' + filenames1[i])
            fullpath1 = os.path.join(cfg.datapath, filenames1[i])
            with open(fullpath1) as csvfile:
                #GOES data has very large headers; figure out where the data
                #starts inside the file and skip the required number of lines
                readCSV = csv.reader(csvfile, delimiter=',')
                for k in range(nhead):
                    next(readCSV)  #to start of data

                #Get dates first; need dates to identify spacecraft orientation
                #for GOES-13+
                for row in readCSV:
                    date = datetime.datetime.strptime(row[0][0:19],
                                                    "%Y-%m-%d %H:%M:%S")
                    dates.append(date)
                    file_dates.append(date) #only dates in current file

                if (experiment == "GOES-13" or experiment == "GOES-14"
                    or experiment == "GOES-15"):
                    west_detector = []
                    if filenames_orien[i] != None:
                        west_detector = get_west_detector(filenames_orien[i], dates)


                #Go back and get fluxes
                ix = 0
                csvfile.seek(0)
                for k in range(nhead):
                    next(readCSV)  #to start of data
                for row in readCSV:
                    #Get west detector
                    #Some files from NOAA have an orientation file that is missing
                    #the last day, e.g. only goes to the 30th when there are
                    #31 days in the month. To use this data from this month, assume
                    #the west-facing detector remains the same on this missing day.
                    if len(west_detector) == 0:
                        west_det = None
                    elif ix > len(west_detector)-1:
                        west_det = west_detector[-1]
                    else:
                        west_det = west_detector[ix]
                
                    for j in range(ncol):
                        flux = float(row[columns[j]])
                        #Account for orientation
                        if (experiment == "GOES-13" or experiment == "GOES-14"
                            or experiment == "GOES-15"):
                            if west_det == None:
                                flux = badval
                            elif west_det == "B":
                                flux = float(row[columnsB[j]])
                            elif west_det == "Flip":
                                flux = badval
                        if flux < 0:
                            flux = badval
                        fluxes[j][ix] = flux
                    ix = ix + 1
            csvfile.close()

            #Update file completeness
            df = file_completeness(df, experiment, flux_type, fullpath1, file_dates)

        #SECOND set of files for higher energy hepad
        if filenames2[i] != None:
            nhead, nrow = find_goes_data_dimensions(filenames2[i])
            
            if len(fluxes) == 0:
                fluxes = np.zeros(shape=(totcol,nrow))
            
            fullpath2 = os.path.join(cfg.datapath,filenames2[i])
            print('Reading in file ' + fullpath2)
            with open(fullpath2) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                for k in range(nhead):
                    next(readCSV)  #to start of data

                count = 0
                for row in readCSV:
                    for j in range(nhcol):
                        #Sometimes GOES datafiles have an incomplete row,
                        #e.g. g13_hepad_ap_5m_20101201_20101231.csv
                        #If so, rather than crash, set the flux to badval
                        try:
                            flux = float(row[hepad_columns[j]])
                        except:
                            flux = badval
                        if flux < 0:
                            flux = badval
                        fluxes[ncol+j][count] = flux
                    count = count + 1
            csvfile.close()

            #Update file completeness
            df = file_completeness(df, experiment, flux_type, fullpath2, file_dates)
    
        if len(dates) == 0:
            continue

        #If reading in multiple files, then combine all data into one array
        if len(all_fluxes) == 0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates
    
    if all_dates == []:
        print(f"read_in_goes: Did not find the data you were looking for. {filenames1}")
   
    write_data_manager(df)
    
    return all_dates, all_fluxes, west_detector



def read_in_goesR(experiment, flux_type, filenames1):
    """Read in GOES-R+ data from your computer.
        Appears that only differential channels + one >500 MeV
        integral channel are available in the files.
        
        Reading in Level 2 data.
        Flux fill value = -1e+31
        Differential flux units = protons/(cm^2 sr keV s)
        
        Time stamp is seconds since 2000-01-01 12:00:00
        
        +X and -X are stored in the same array as the 0th and 1st
        array entry... figuring out which are which
        Typicaly, -X faces West and +X faces East
        
        YawFlipFlag indicates if the spacecraft is flipped.
        Assume a value of 1 indicates the spacecraft is flipped.
        
        https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/docs/GOES-16_SEISS_SGPS_L1b_Provisional_Maturity_ReadMe.pdf
        "There are two SGPS sensor units mounted on each GOES-R series spacecraft, facing in the spacecraft -X and +X directions. When the spacecraft is not in the yaw-flipped configuration SGPS-X faces west and SGPS+X faces east. Each SGPS unit has three solid-state (silicon detector) telescopes T1, T2, and T3 for measuring 1-25, 25-80, and 80-500 MeV protons, respectively. All three telescopes have the same look direction (i.e., +X or -X). T1 and T2 have 60o (full cone angle) fields of view, and T3 has a 90o field of view. Each unit measures 1-500 MeV proton fluxes in 13 logarithmically spaced differential channels (P1-P10) and >500 proton flux in a single integral channel (P11). The L1b data product is one-second cadence fluxes. The channels generally register counts above backgrounds only during solar energetic particle events, except for P11 which measures galactic cosmic rays in the absence of a solar particle event."
        
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) the files containing the GOES-R
            netcdf data that span the desired time range
        
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    ndiff_chan = 13 #5
    conversion = 1000. #keV/MeV
    
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []
    west_detector = [] #place holder, will be filled if needed

    #Read in file that identified data files as complete
    df = read_data_manager()

    #GOES-R times are wrt reference time of 2000-01-01 12:00:00
    ref_date = datetime.datetime(year=2000, month=1, day=1, hour=12)
    
    #Read in fluxes from files
    for i in range(NFILES):
        file_dates = []
        if filenames1 != None:
            print(filenames1[i])
            fullpath = os.path.join(cfg.datapath, filenames1[i])
            infile = os.path.expanduser(fullpath)
            data = netCDF4.Dataset(infile)
            
            if "v3-0" in filenames1[i]:
                ntstep = len(data.variables["time"])
            else:
                ntstep = len(data.variables["L2_SciData_TimeStamp"])
            
            #13 differential channels, one integral channel
            #5 minute time steps
            fluxes = np.zeros(shape=(ndiff_chan+1,ntstep))
            for j in range(ntstep):
                if "v3-0" in filenames1[i]:
                    time_sec = float(data.variables["time"][j].data)
                else:
                    time_sec = float(data.variables["L2_SciData_TimeStamp"][j].data)
                td = datetime.timedelta(seconds=time_sec)
                date = ref_date + td
                all_dates.append(date)
                file_dates.append(date)
                
                #Orientation flag
                if "v3-0" in filenames1[i]:
                    flip_flag = data.variables["yaw_flip_flag"][j]
                else:
                    flip_flag = data.variables["YawFlipFlag"][j]
                idx = 0
                if flip_flag == 1:
                    idx = 1
                #if flip_flag > 1:
                #    idx = None #exclude these points because in process of flip
                
                #Extract the 13 differential channels
                for k in range(ndiff_chan):
                    ###TEMP###
                    #kk = k + 8
                    #[288 time step, 2 +/-X, 13 energy chan]
                    
                    flux = data.variables["AvgDiffProtonFlux"][j][idx][k]
                    if flux < 0:
                        flux = badval
                    fluxes[k][j] = flux*conversion
                

                flux = data.variables["AvgIntProtonFlux"][j][idx]
                if flux < 0:
                    flux = badval
                fluxes[-1][j] = flux

            #Update file completeness
            df = file_completeness(df, experiment, flux_type, fullpath, file_dates)

        if len(file_dates) == 0:
            continue
        
        if len(all_fluxes) == 0:
            all_fluxes = fluxes
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)

    write_data_manager(df)

    return all_dates, all_fluxes, west_detector


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
    if zt == 0:
        return 0
        
    if 'Z' not in zt or 'T' not in zt:
        print("zulu_to_time: Time not in proper format. Returning None.")
        return None
    
    strzt = zt.split('T')
    strzt[1] = strzt[1].strip('Z')
    n = strzt[1].split(':')
    stdt = strzt[0] + ' ' + strzt[1]

    if len(n) == 2:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M')
    if len(n) == 3:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M:%S')
    return dt



def read_in_goes_RT(experiment, flux_type, filenames1):
    """ Read in GOES_RT real time data from your computer.
        Read in the NOAA SWPC real time integral flux files from
        the 1 day json files for the primary GOES spacecraft.
        These files are archived on the CCMC website.
        
        iSWA may have time gaps that exceed 1 day.
        
        The >500 MeV fluxes are in the last row of the csv file for 
        GOES-R time periods.
        
        Reading in Level 2 data.
        Flux fill value = -100000.0
        Units: Particles = Protons/cm2-s-sr
        Units: Electrons = Electrons/cm2-s-sr
        
        Time stamp is YYYYY M D HHMM: 2022  1  20  0000
        
        No selection for +/-X direction. Assume always westward facing.
        
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) the files containing the GOES-R
            netcdf data that span the desired time range
        
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    n_chan = 8 #>1, >5, >10, >30, >50, >100, >60, >500 possible
    
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []
    west_detector = [] #place holder, will be filled if needed

    #Read in file that identified data files as complete
    df = read_data_manager()
    
    df_data = pd.DataFrame()
    cols_to_drop = [7,8,9,12,13] #Not proton fluxes
    #Read in fluxes from files
    for i in range(NFILES):
        file_dates = []
        if filenames1[i] == None:
            continue

        fullpath = os.path.join(datapath, filenames1[i])
        try:
            df_in = pd.read_csv(fullpath, header=None)
        except:
            #Sometimes files may be empty if there is a data gap > 1 day
            print(f"read_in_GOES_RT: Could not open {fullpath}. Skipping.")
            continue

        df_in[0] = df_in[0].str.replace('T',' ')
        df_in[0] = df_in[0].str.replace('Z','')
        df_in[0] = pd.to_datetime(df_in[0])
 
        #Update file completeness
        file_dates = df_in[0].to_list()
        df = file_completeness(df, experiment, flux_type, fullpath, file_dates)
 
        for col in cols_to_drop:
            if col in df_in.columns:
                df_in = df_in.drop(col,axis=1)
                        

        #Replace bad data with badval
        df_in = df_in.replace(to_replace=-100000,value=badval)

        df_data = pd.concat([df_data,df_in],ignore_index=True)
        
    if not df_data.empty:
    #These final dates and fluxes may have time gaps that were not filled in
        all_dates = df_data[0].to_list()
        cols = df_data.columns.to_list()
        all_fluxes = df_data[cols[1:]].values.T

        write_data_manager(df)

    return all_dates, all_fluxes, west_detector


def read_in_all_goes(experiment, flux_type, filenames1, filenames2,
                filenames_orien, options, detector, spacecraft=""):
    """ If the user specified "GOES" as the experiment, then
        read in the files from any of the GOES spacecraft that
        are available in the date range.
        
        ONLY ALLOW FOR INTEGRAL FLUXES, OTHERWISE RUN INTO PROBLEM
        OF ENERGY BINS. (Building in capabilities to read in
        differential fluxes, but won't have the right energy bins and
        won't work with differential channels at the moment.)
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
        :detector: (string array) for each file, the list contains
            the detector associated with that file
            
        OUTPUTS:
        
        :all_dates_out: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes_out: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    for i in range(len(filenames1)):
        print(f"{detector[i]} {filenames1[i]}")
    
    #Choose to align the fluxes with the goes_energy_bins
#    if flux_type == "integral":
#        nbins = 7
#    if flux_type == "differential":
#        nbins = 10 #all GOES differential have up to 10 channels
    
    #Get the energy bins for the first detector in the list
    ref_energy_bins = []
    
    all_dates_out = []
    all_fluxes_out = []
    west_detector = [] #place holder, not needed
    all_fluxes_out = []

    sc_starttimes = []
    sc_endtimes = []
    sc_detector = []

    #Read in each file one at a time because they may be different
    #detectors for each file
    for i in range(len(filenames1)):
        #All GOES data
        print(f"read_in_all_goes: Reading in files for {detector[i]}, File: {filenames1[i]}")
        if detector[i] in old_goes_sc:
            all_dates, all_fluxes, west_detector =\
                read_in_old_goes(detector[i], flux_type, [filenames1[i]],
                    [filenames2[i]], options)

        elif detector[i] in goes_sc:
            all_dates, all_fluxes, west_detector =\
                read_in_goes(detector[i], flux_type, [filenames1[i]],
                [filenames2[i]], [filenames_orien[i]], options)
            
        elif detector[i] in goes_R and flux_type == "differential":
            all_dates, all_fluxes, west_detector =\
                read_in_goesR(detector[i],flux_type, [filenames1[i]])
            
        elif (detector[i] == "GOES_RT" or detector[i] in goes_R) and flux_type == "integral":
            all_dates, all_fluxes, west_detector =\
                read_in_goes_RT(detector[i],flux_type, [filenames1[i]])

        if len(all_fluxes) == 0: continue
        #Energy bins for this detector
        energy_bins, energy_bin_centers = define_energy_bins(detector[i], flux_type,
                        west_detector, options,spacecraft=spacecraft)
        all_fluxes, energy_bins, energy_bin_centers = tools.sort_bin_order(all_fluxes,
                        energy_bins, energy_bin_centers)

        #Use the energy bins from the first spacecraft to define the energy
        #bins for the whole time series.
        #If the fluxes are integral, then all GOES-06+ will have the same
        #integral channels.
        #If the fluxes are differential, different spacecraft don't have the
        #same energy channels, so use with care. Will use those of the first
        #spacecraft to define the energy bins. This will affect estimated
        #>10 MeV, etc, fluxes. Better to create estimated integral channels
        #from each spacecraft, then string together, if that is the desired product.
        if len(ref_energy_bins) == 0:
            ref_energy_bins = energy_bins
            nbins = len(ref_energy_bins)
            #Target energy bins will be those in goes_energy_bins
            for k in range(nbins):
                all_fluxes_out.append([])

        #Easier to do some manipulations as a dataframe
        dict= {'Dates': all_dates}
        for kk in range(len(all_fluxes)):
            dict.update({(f"Flux{kk}"): all_fluxes[kk]})
        df_flux = pd.DataFrame(dict)

        ####PRIMARY/SECONDARY
        #Detector variable gives which detector is being read in
        #with the the filename and this was chosen to be primary or
        #secondary already in the check stage. Need to figure out how
        #to get the start and end times at exactly the transition
        #times between primary spacecraft
        startdate = all_dates[0]
        enddate = all_dates[len(all_dates)-1]
        starttimes, endtimes, goes = identify_which_goes_spacecraft(startdate,
                enddate, spacecraft=spacecraft)

        #Check the primary or secondary date ranges and fill in the lists
        #accordingly
        sc_dates = []
        sc_fluxes = []
        for j in range(len(starttimes)):
            start = starttimes[j]
            end = endtimes[j]
            det = goes[j]
            if det in goes_R and flux_type == "integral":
                det = "GOES_RT"
            if start in sc_starttimes: continue
            if det != detector[i]: continue
            sub = df_flux.loc[(df_flux['Dates'] >= start) & (df_flux['Dates'] < end)]
            if sub.empty: continue
            print(f"read_in_all_goes: Selecting fluxes for {goes[j]} between {start} and {end} for the {spacecraft} spacecraft.")
            #Extract only the dates for the primary/secondary satellite in this file
            sc_dates = sub['Dates'].to_list()
            for kk in range(len(sub.columns)-1):
                sc_fluxes.append(sub[(f"Flux{kk}")].to_list())
            #Record that these were already extracted
            sc_starttimes.append(start)
            sc_endtimes.append(end)
            sc_detector.append(det)
            break
        
        #Extend the date and flux arrays to combine the different spacecraft
        all_dates_out.extend(sc_dates)
        for k in range(len(ref_energy_bins)):
            try:
                idx = energy_bins.index(ref_energy_bins[k])
                all_fluxes_out[k].extend(sc_fluxes[idx])
            except:
                all_fluxes_out[k].extend([cfg.badval]*len(sc_dates))

    print("Combined fluxes for: ")
    for k in range(len(sc_detector)):
        print(f"{sc_detector[k]} from {sc_starttimes[k]} to {sc_endtimes[k]}")
    all_fluxes_out = np.array(all_fluxes_out)
    
    return all_dates_out, all_fluxes_out, west_detector, ref_energy_bins




def read_in_ephin(experiment, flux_type, filenames1):
    """ Read in EPHIN files from your computer.
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) names of files containing
            EPHIN data in desired time range (yearly files)
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
    
    """
    NFILES = len(filenames1)

    #Read in file that identified data files as complete
    df = read_data_manager()

    datecols = [0,1,2,4,5] #yr, mth, dy, hr, min
    fluxcols = [8,9,10,11]
    ncol= len(fluxcols)

    for i in range(NFILES):
        fullpath = os.path.join(cfg.datapath, filenames1[i])
        print('Reading in file ' + fullpath)
        with open(fullpath) as csvfile:
            #Count header lines indicated by hash #
            nhead = 0
            for line in csvfile:
                line = line.lstrip()
                if line[0] == "#":
                    nhead = nhead + 1
                else:
                    break
            #number of lines containing data
            nrow = len(csvfile.readlines())+1

            #Define arrays that hold dates and fluxes
            dates = []
            fluxes = np.zeros(shape=(ncol,nrow))

            csvfile.seek(0) #back to beginning of file
            for k in range(nhead):
                csvfile.readline()  # Skip header rows.

            count = 0
            for line in csvfile:
                if line == '': continue
                if line[0] == "#": continue

                row = line.split()
                yr = int(row[datecols[0]])
                mnth = int(row[datecols[1]])
                dy = int(row[datecols[2]])
                hr = int(row[datecols[3]])
                min = int(row[datecols[4]])
                date = datetime.datetime(yr,mnth,dy,hr,min,0,0)
                dates.append(date)
                for j in range(ncol):
                    flux = float(row[fluxcols[j]])
                    if flux < 0:
                        flux = badval
                       # print(f"SETTING EPHIN FLUX TO BADVAL {badval} for {j}, {count}")
                    fluxes[j][count] = flux
                count = count + 1

        #Update file completeness
        df = file_completeness(df, experiment, flux_type, fullpath, dates)

        #If reading in multiple files, then combine all data into one array
        if i==0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates

    write_data_manager(df)

    return all_dates, all_fluxes


def read_in_ephin_release(experiment, flux_type, filenames1):
    """ Read in EPHIN files from your computer.
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) names of files containing
            EPHIN REleASE data in desired time range (yearly files)
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
    
    """
    NFILES = len(filenames1)

    datecols = [0] #yr, mth, dy, hr, min
    fluxcols = [1,2,3,4]
    ncol= len(fluxcols)

    for i in range(NFILES):
        print('Reading in file ' + datapath + '/' + filenames1[i])
        with open(datapath + '/' + filenames1[i]) as csvfile:
            #Count header lines indicated by hash #
            nhead = 0
            for line in csvfile:
                line = line.lstrip()
                if line == '':
                    nhead = nhead + 1
                elif line[0] == "#":
                    nhead = nhead + 1
                else:
                    break
            #number of lines containing data
            nrow = len(csvfile.readlines())+1

            #Define arrays that hold dates and fluxes
            dates = []
            fluxes = np.zeros(shape=(ncol,nrow))

            csvfile.seek(0) #back to beginning of file
            for k in range(nhead):
                csvfile.readline()  # Skip header rows.

            count = 0
            for line in csvfile:
                if line == '': continue
                if line[0] == "#": continue

                row = line.split(';')
                date = datetime.datetime.strptime(row[0][0:19],
                                                "%Y-%m-%d %H:%M:%S")
                dates.append(date)
                for j in range(ncol):
                    flux = float(row[fluxcols[j]])
                    if flux < 0:
                        flux = badval
                    fluxes[j][count] = flux
                count = count + 1

        #If reading in multiple files, then combine all data into one array
        if i==0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates

    return all_dates, all_fluxes



def read_in_erne(experiment, flux_type, filenames1):
    """ Read in ERNE files from your computer.
        Each ERNE gzip file (erne-1996.05.07-1996.05.31-23.tgz) has
        a subfolder export.src which contains HED files names as
        HEDYYDOY.SL2, which contains count rates and intensities.
        LEDYYDOY.SL2, contains 10 - 13 MeV channel
        
        More information about the data formats and caveats:
        https://export.srl.utu.fi/export_data_description.txt
        
        Note that timestamps in the files are of the format
        TAI time stamp (seconds since 1.1.1958) and will not be used.

        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) names of files containing
            ERNE data in desired time range (randomly dated files)
            This list contains both LED and HED files in order of
            however os.listdir lists them.
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
    
    """
    #filenames1 contains a list of both HED and LED files in somewhat
    #random order. Need to divide the low and high energy fluxes and
    #align by date.
    HEDfiles = sorted([x for x in filenames1 if 'HED' in x])
    LEDfiles = sorted([x for x in filenames1 if 'LED' in x])
    
    if len(HEDfiles) != len(LEDfiles):
        print("LED FILES")
        print(LEDfiles)
        print("HED FILES")
        print(HEDfiles)
        sys.exit("read_datasets: read_in_erne: There are not the same number of "
                "LED and HED files. Cannot align. Exiting.")



    NFILES = len(HEDfiles)

    datecols = [3] #seconds since 1958-01-01
    fluxcols = [x for x in range(4,14)] #cols 4 - 13 for HED, LED
    fluxcols.extend(fluxcols) #2x for LED and HED
    ncol= len(fluxcols) #HED and LED each has 10 channels

    for i in range(NFILES):
        print('Reading in file ' + LEDfiles[i])
        with open(LEDfiles[i], 'r') as ledfile, open(HEDfiles[i], 'r') as hedfile:
            #Count header lines indicated by hash #
            nhead_led = 0
            for line in ledfile:
                line = line.lstrip()
                if line[0] == "#":
                    nhead_led = nhead_led + 1
                else:
                    break
            #number of lines containing data
            nrow_led = len(ledfile.readlines())+1
            
            nhead_hed = 0
            for line in hedfile:
                line = line.lstrip()
                if line[0] == "#":
                    nhead_hed = nhead_hed + 1
                else:
                    break
            #number of lines containing data
            nrow_hed = len(hedfile.readlines())+1

            if nrow_led != nrow_hed:
                print("# LED data lines = " + str(nrow_led) +", " + LEDfiles[i])
                print("# HED data lines = " + str(nrow_hed) +", " + HEDfiles[i])
                sys.exit("read_datasets: read_in_erne: There are not the same number of "
                        "lines in LED and HED files. Cannot align. Exiting.")


            #Define arrays that hold dates and fluxes
            dates = []
            fluxes = np.zeros(shape=(ncol,nrow_led))

            ledfile.seek(0) #back to beginning of file
            hedfile.seek(0) #back to beginning of file
            for k in range(nhead_led):
                ledfile.readline()  # Skip header rows and get to data.

            for k in range(nhead_hed):
                hedfile.readline()  # Skip header rows and get to data.

            count = 0 #counts rows, i.e. dates
            for k in range(nrow_led):
                line_led = ledfile.readline()
                line_hed = hedfile.readline()
                
                row_led = line_led.split()
                row_hed = line_hed.split()
                
                ts_led = int(row_led[datecols[0]])
                ts_hed = int(row_hed[datecols[0]])
                
                if ts_led != ts_hed:
                    print("Timestamps in LED and HED files on the same line "
                        "do not match. Data out of order. Look for extra blank lines "
                        "or files that do not correspond to same time period.")
                    print("LED file: " + LEDfiles[i] + ", TIMESTAMP: " + str(ts_led))
                    print("HED file: " + HEDfiles[i] + ", TIMESTAMP: " + str(ts_hed))
                    sys.exit("Exiting")
                
                date = datetime.datetime(1958,1,1) + datetime.timedelta(seconds=ts_led)
 
                dates.append(date)
                for j in range(ncol):
                    if j < ncol/2:
                        flux = float(row_led[fluxcols[j]])
                    else:
                        flux = float(row_hed[fluxcols[j]])
                    if flux < 0:
                        flux = badval
                    fluxes[j][count] = flux
                count = count + 1

        #If reading in multiple files, then combine all data into one array
        if i==0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates

    return all_dates, all_fluxes




def read_in_stereo(experiment, flux_type, filenames1, filenames2):
    """ Read in STEREO-A or STEREO-B LET (summed) and HET files
        and combine together to make a time profile across
        the energy ranged covered by the STEREO spacecraft,
        1 - 100 MeV.
        
        A complication of this paticular data set is that the lower
        energy data are stored in daily files while the higher
        energy data are stored in yearly files. Need to take care
        to match up the data points in time.

    """

    NFILESL = len(filenames1) #LET, daily
    NFILESH = len(filenames2) #HET, yearly

    #Read in file that identified data files as complete
    df = read_data_manager()

    print(f"read_in_stereo: Reading in {NFILESL} LET files going from {filenames1[0]} to {filenames1[-1]}.")
    print(f"read_in_stereo: Reading in {NFILESH} HET files going from {filenames2[0]} to {filenames2[-1]}.")

    datecolsL = [0,1,2,3,4] #yr, doy (frac), hour, min, sec - not used
    fluxcolsL = [7,8,9] #1.8-3.6, 4-6, 6-10, 10-15 <-- always empty
    ncolL = len(fluxcolsL)
    colL = ['fluxL' + str(kk) for kk in range(ncolL) ]
    colL = ['dates'] + colL
    dfL = pd.DataFrame(columns=colL)

    datecolsH = [1,2,3,4] #yr, doy, dy, hr, min - not used
    fluxcolsH = [11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]
    ncolH = len(fluxcolsH)
    colH = ['fluxH' + str(kk) for kk in range(ncolH) ]
    colH = ['dates'] + colH
    dfH = pd.DataFrame(columns=colH)
    
    #READ IN LET
    for i in range(NFILESL):
        fullpathL = os.path.join(cfg.datapath,filenames1[i])
        print(f"{datetime.datetime.now()} Reading in file {fullpathL}")
        with open(fullpathL) as infile:
            #Count header lines up until "BEGIN DATA"
            #Count remaining lines of data
            nhead = 0
            nrowL = 0
            is_header = True
            for line in infile:
                line = line.lstrip()
                if is_header:
                    nhead = nhead + 1
                
                if 'BEGIN DATA' in line:
                    #nhead = nhead + 1
                    is_header = False
                
                if not is_header and line != '': #empty line at end of file
                    nrowL = nrowL + 1
            
            #Define arrays that hold dates and fluxes
            datesL = []
 
            infile.seek(0) #back to beginning of file
            for k in range(nhead):
                infile.readline()  # Skip header rows.

            count = 0
            for line in infile:
                rowL = [] #row for dataframe
                line = line.lstrip()
                if line == '': continue
                if line[0] == "#": continue

                row = line.split()
                year = int(row[0])
                doy = math.floor(float(row[1]))
                hour = int(row[2])
                minute = int(row[3])
                sec = int(row[4])
                
                str_date = str(year) + " " + str(doy) + " " + str(hour) + " " \
                    + str(minute)
                date = datetime.datetime.strptime(str_date,"%Y %j %H %M")
                
                datesL.append(date)
                rowL.append(date)
                
                for j in range(ncolL):
                    flux = float(row[fluxcolsL[j]])
                    if flux < 0:
                        flux = badval

                    rowL.append(flux)
                
                dfLrow = pd.DataFrame([rowL], columns=colL)
                dfL = pd.concat([dfL,dfLrow],ignore_index=True)
                count = count + 1

        #Update file completeness
        df = file_completeness(df, experiment, flux_type, fullpathL, datesL)

    print(f"{datetime.datetime.now()} read_in_stereo: Finished reading LET data.")
    write_data_manager(df)

    #READ IN HET
    for i in range(NFILESH):
        fullpathH = os.path.join(datapath,filenames2[i])
        print(f"{datetime.datetime.now()} Reading in file {fullpathH}")
        with open(fullpathH) as infile:
            #Count header lines up until "BEGIN DATA"
            #Count remaining lines of data
            nhead = 0
            nrowH = 0
            is_header = True
            for line in infile:
                line = line.lstrip()
                if is_header:
                    nhead = nhead + 1
                
                if '#End' in line:
                   # nhead = nhead + 1
                    is_header = False
                
                if not is_header and line != '': #empty line at end of file
                    nrowH = nrowH + 1
            
            #Define arrays that hold dates and fluxes
            datesH = []

            infile.seek(0) #back to beginning of file
            for k in range(nhead):
                infile.readline()  # Skip header rows.

            count = 0
            for line in infile:
                rowH = []
                line = line.lstrip()
                if line == '': continue
                if line[0] == "#": continue

                row = line.split()
                
                #Date 0 2021 Jan 1 0000
                year = int(row[1])
                str_month = row[2] #str
                day = int(row[3])
                str_hr_min = row[4]

                str_date = str(year) + " " + str_month + " " + str(day) + " " + str_hr_min
                date = datetime.datetime.strptime(str_date, '%Y %b %d %H%M') #2001 Jan 1 1545
 
                datesH.append(date)
                rowH.append(date)
                
                for j in range(ncolH):
                    flux = float(row[fluxcolsH[j]])
                    if flux < 0:
                        flux = badval

                    rowH.append(flux)

                dfHrow = pd.DataFrame([rowH], columns=colH)
                dfH = pd.concat([dfH,dfHrow],ignore_index=True)
                count = count + 1
        
        #Update file completeness
        df = file_completeness(df, experiment, flux_type, fullpathH, datesH)

    print(f"{datetime.datetime.now()} read_in_stereo: Finished reading HET data.")
    write_data_manager(df)
    
    #Now we have the LET and HET data for every minute in different arrays.
    #The HET rows must be appended at the end of the LET rows (to go up
    #in energy bin)
    #There is no guarantee that they are aligned in time.
    #Need to find the appropriate alignment and fill any mismatched elements
    #with None values
    #Easiest approach is to trim both to the minimum overlapping date range
    #and hopefully this will include the dates requested by the user.
    #This approach may need updating.
    first_date = max(dfL['dates'].min(),dfH['dates'].min())
    last_date = min(dfL['dates'].max(),dfH['dates'].max())

    #Trim to the overlapping date range
    dfL = dfL.loc[(dfL['dates'] >= first_date) & (dfL['dates'] <= last_date)]
    dfH = dfH.loc[(dfH['dates'] >= first_date) & (dfH['dates'] <= last_date)]
    
    print(f"read_in_stereo: There may be missing timestamps in the LET and HET data sets. "
        f"Make a time point for every minute between {first_date} and {last_date}, then fill "
        "in with bad value.")
    badL = [badval]*ncolL
    badH = [badval]*ncolH
    nmins = math.ceil((last_date - first_date).total_seconds()/60.) + 1
    print(f"There are {len(dfL)} LET time points and {len(dfH)} HET time points between "
        f"{first_date} and {last_date}. There should be {nmins} minutes.")
        
    isgoodL = (len(dfL) == nmins)
    isgoodH = (len(dfH) == nmins)
    if isgoodL and isgoodH:
        pass
    else:
        for ii in range(nmins):
            time = first_date + datetime.timedelta(minutes = ii)

            if not isgoodL:
                have_timeL = dfL['dates'].isin([time]).any()
                if not have_timeL:
                    badrowL = [time] + badL
                    dfLbad = pd.DataFrame([badrowL],columns=colL)
                    dfL = pd.concat([dfL,dfLbad],ignore_index=True)
                    print(f"{datetime.datetime.now()} read_in_stereo: Time point missing in LET at {time}. Filled with {badval}.")

            if not isgoodH:
                have_timeH = dfH['dates'].isin([time]).any()
                if not have_timeH:
                    badrowH = [time] + badH
                    dfHbad = pd.DataFrame([badrowH], columns=colH)
                    dfH = pd.concat([dfH,dfHbad],ignore_index=True)
                    print(f"{datetime.datetime.now()} read_in_stereo: Time point missing in HET at {time}. Filled with {badval}.")


    #After bad times have been filled in, sort dataframes to be in time order again
    #and remove any duplicated time points
    dfL = dfL.sort_values(by=['dates'], ascending=True)
    dfL.drop_duplicates(subset='dates', inplace=True)
    dfL = dfL.reset_index(drop = True)
 
    dfH = dfH.sort_values(by=['dates'], ascending=True)
    dfH.drop_duplicates(subset='dates', inplace=True)
    dfH = dfH.reset_index(drop = True)

    #Check if the timestamps across the two data sets are the same
    if not dfL['dates'].equals(dfH['dates']):
        dfL.to_csv("output/idsep/STEREO_LET.csv")
        dfH.to_csv("output/idsep/STEREO_HET.csv")
        sys.exit("read_in_stereo: The LET and HET data contain different timestamps even after trimming "
            "and filling in missing time points (1 minute cadence). Exiting.")

    df_all = dfL.merge(dfH, on='dates')
    allcol = df_all.columns.to_list()

    all_dates = df_all['dates'].to_list()
    all_fluxes = df_all[allcol[1:]].values.T

    print(f"{datetime.datetime.now()} read_in_stereo: Finished reading STEREO data.")

    return all_dates, all_fluxes



def read_in_ace_sis(experiment, flux_type, filenames1):
    """ Read in ACE/SIS integral >30, >60 MeV data from your computer.
        
        # Units: proton flux p/cs2-sec-ster
        # Status(S): 0 = nominal data, 1 to 8 = bad data record, 9 = no data
        # Missing data values: -1.00e+05
        # Source: ACE Satellite - Solar Isotope Spectrometer
        #
        # 5-minute averaged Real-time Integral Flux of High-energy Solar Protons
        
                          Modified Seconds
        # UT Date   Time   Julian  of the      ---- Integral Proton Flux ----
        # YR MO DA  HHMM     Day     Day       S    > 10 MeV    S    > 30 MeV
        #--------------------------------------------------------------------
        2001 08 07  0000    52128       0      0    7.97e-01    0    5.57e-01        
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral
        :filenames1: (string array) the files containing the data
        
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    n_chan = 2 #>30, >60
    
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []

    #Read in file that identified data files as complete
    df = read_data_manager()

    #Read in fluxes from files
    for i in range(NFILES):
        file_dates = []
        if filenames1[i] == None:
            continue

        fullpath = os.path.join(datapath, filenames1[i])
        if not os.path.isfile(fullpath):
            print(f"read_in_ace_sis: Cannot read {fullpath}. Skipping.")
            continue
        with open(fullpath, 'r') as file:
            for line in file:
                if ":Data" in line: continue
                if ":Created" in line: continue
                if "#" in line: continue
                line = line.strip().split()

                #Date
                year = int(line[0])
                month = int(line[1])
                day = int(line[2])
                time = line[3]
                hr = int(time[0:2])
                min = int(time[2:4])
                date = datetime.datetime(year = year, month=month, day=day, hour=hr, minute=min)
                file_dates.append(date)
                all_dates.append(date)

                flx30 = float(line[7])
                flx60 = float(line[9])
                if flx30 < 0: flx30 = badval
                if flx60 < 0: flx60 = badval
            
                all_fluxes.append([flx30, flx60])
                


        df = file_completeness(df, experiment, flux_type, fullpath, file_dates)

    print(f"{datetime.datetime.now()} read_in_ace_sis: Finished reading ACE/SIS integral data.")
    write_data_manager(df)
    
    all_fluxes = np.array(all_fluxes).T

    return all_dates, all_fluxes


def read_in_ace_epam_electrons(experiment, flux_type, filenames1):
    """ Read in ACE/EPAM energetic electrons 175-315 keV from your computer.
        
        # Units: Differential Flux particles/cm2-s-ster-MeV
        # Units: Anisotropy Index 0.0 - 2.0
        # Status(S): 0 = nominal data, 1 to 8 = bad data record, 9 = no data
        # Missing data values: -1.00e+05, index = -1.00
        # Source: ACE Satellite - Electron, Proton, and Alpha Monitor
        # Note: 7/26/00 Corrected Differential Flux ranges.
        #
        #                      5-minute averaged Real-time Differential Electron and Proton Flux 
        # 
        #                Modified Seconds ---------------------------- Differential Flux --------------------------- 
        # UT Date   Time  Julian  of the  ----- Electron -----   ------------------- Protons keV -------------------   Anis.
        # YR MO DA  HHMM    Day    Day    S    38-53   175-315   S   47-65    112-187   310-580   761-1220 060-1910   Ratio
        #-------------------------------------------------------------------------------------------------------------------
        2001 08 07  0000   52128       0  0  7.22e+02  1.60e+01  0  1.28e+03  1.29e+02  9.20e+00  9.54e-01  2.19e-01   0.47
                
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) differential
        :filenames1: (string array) the files containing data
        
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    n_chan = 1 #175-315 keV
    
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []

    #Read in file that identified data files as complete
    df = read_data_manager()
    
    #Read in fluxes from files
    for i in range(NFILES):
        file_dates = []
        if filenames1[i] == None:
            continue

        fullpath = os.path.join(datapath, filenames1[i])
        if not os.path.isfile(fullpath):
            print(f"read_in_ace_epam_electrons: Cannot read {fullpath}. Skipping.")
            continue
        with open(fullpath, 'r') as file:
            for line in file:
                if ":Data" in line: continue
                if ":Created" in line: continue
                if "#" in line: continue
                line = line.strip().split()

                #Date
                year = int(line[0])
                month = int(line[1])
                day = int(line[2])
                time = line[3]
                hr = int(time[0:2])
                min = int(time[2:4])
                date = datetime.datetime(year = year, month=month, day=day, hour=hr, minute=min)
                file_dates.append(date)
                all_dates.append(date)

                flx = float(line[8])
                if flx < 0: flx = badval
            
                all_fluxes.append(flx)
                


        df = file_completeness(df, experiment, flux_type, fullpath, file_dates)

    print(f"{datetime.datetime.now()} read_in_ace_epam_electrons: Finished reading ACE/EPAM energetic electron data.")
    write_data_manager(df)
    
    all_fluxes = np.array([all_fluxes])

    return all_dates, all_fluxes


def read_in_imp8_cpme(experiment, flux_type, filenames1):
    """ Read in IMP-8/CPME protons from your computer.
        
        IMP-8 CPME: 330-sec. Avg. Proton Intensities & Uncertainties [no./(cm^2-sc-ster-MeV)]                                   
        year doy hr mn sc     dec_year  dec_doy  x_gse_km  y_gse_km  z_gse_km  y_gsm_km  
        z_gsm_km     p1_fx    p1_fxu     p2_fx    p2_fxu     p3_fx    p3_fxu     p4_fx    
        p4_fxu     p5_fx    p5_fxu     p7_fx    p7_fxu     p8_fx    p8_fxu     p9_fx    
        p9_fxu    p10_fx   p10_fxu    p11_fx   p11_fxu
        
        NOTE NO P6

        P5 = 4.60 - 15.0 MeV
        P7 = 15.0 - 25.0 MeV
        P8 = 25.0 - 48.0 MeV
        P9 = 48.0 - 96.0 MeV
        P10 = 96.0 - 145.0 MeV
        P11 = 145.0 - 440.0 MeV
                       
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) differential
        :filenames1: (string array) the files containing data
        
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    n_chan = 6 #P5 - P11

    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []

    #Read in file that identified data files as complete
    df = read_data_manager()
    
    #Read in fluxes from files
    for i in range(NFILES):
        file_dates = []
        if filenames1[i] == None:
            continue

        if not os.path.isfile(filenames1[i]):
            print(f"read_in_imp8_cpme: Cannot read {filenames1[i]}. Skipping.")
            continue
        with open(filenames1[i], 'r') as file:
            print(f"read_in_imp8_cpme: Reading {filenames1[i]}.")
            for line in file:
                #decimal year in column 5
                #P5 - P11 (no P6), columns 20, 22, 24, 26, 28, 30
                if "IMP-8" in line: continue
                if "year" in line: continue
                if "#" in line: continue
                line = line.strip().split()

                #Date
                frac_year = float(line[5])
                date = convert_frac_year(frac_year)
                file_dates.append(date)
                all_dates.append(date)

                flx = [float(line[20]), float(line[22]), float(line[24]), float(line[26]),
                    float(line[28]), float(line[30])]
                for jj in range(len(flx)):
                    if flx[jj] < 0: flx[jj] = badval
            
                all_fluxes.append(flx)
                


        df = file_completeness(df, experiment, flux_type, filenames1[i], file_dates)

    print(f"{datetime.datetime.now()} read_in_imp8_cpme: Finished reading IMP-8/CPME data.")
    write_data_manager(df)
    
    all_fluxes = np.array(all_fluxes).T

    return all_dates, all_fluxes



def read_in_files(experiment, flux_type, filenames1, filenames2,
                filenames_orien, options, detector=[], spacecraft=""):
    """ Read in the appropriate data files with the correct format.
        Return an array with dates and fluxes. Bad flux values (any
        negative flux) are set to -1. Format is defined to work with
        the files downloaded directly from NOAA or the RSDv2 (SEPEM)
        website as is.
        
        The fluxes output for the GOES-13+ satellites are always from
        the westward-facing detector (A or B) by referring to the
        orientation flags provided in the associated orientation file.
        Data taken during a yaw flip (orientation flag = 2) are
        excluded and fluxes are set to -1.
        
        Note that the EPS detectors on GOES-08 and -12 face westward.
        The EPS detector on GOES-10 faces eastward. GOES-11 is a
        spinning satellite.
        
        INPUTS:
        
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        :user_file: (string) name of file containing user-input data
            (if applicable)
        :filenames1: (string array) the files containing the data that
            span the desired time range
        :filenames2: (string array) if GOES, files containing HEPAD data
        :filenames_orien: (string array if GOES, files containing
            satellite orientation
        :options: (string array) options that may be applied to GOES data
        :detector: (string array) if experiment = "GOES",
            will contain which specific spacecraft is needed
            (will be used instead of experiment)
        :spacecraft: (string) FOR GOES ONLY. "primary" or "secondary" 
            to specify which GOES spacecraft to read in.

        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
        :west_detector: (string 1xm array) if GOES, indicates which detector
            is westward facing for every time point
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
       
    """
    print('Reading in data files for ' + experiment + '.')
    all_dates = []
    all_fluxes = []
    west_detector = []

    if experiment == "SEPEM" or experiment == "SEPEMv3":
        all_dates, all_fluxes = read_in_sepem(experiment, flux_type, filenames1)

    elif experiment == "CalGOES":
        all_dates, all_fluxes = read_in_calgoes(experiment, filenames1)

    #All GOES data
    elif experiment == "GOES":
        all_dates, all_fluxes, west_detector, energy_bins = read_in_all_goes(experiment,
                    flux_type, filenames1, filenames2, filenames_orien,
                    options, detector, spacecraft=spacecraft)
        energy_bin_centers = calculate_geometric_means(energy_bins)
        return all_dates, all_fluxes, west_detector, energy_bins, energy_bin_centers
                    
    elif experiment in old_goes_sc:
        all_dates, all_fluxes, west_detector =\
            read_in_old_goes(experiment, flux_type, filenames1, filenames2, options)

    elif experiment in goes_sc:
        all_dates, all_fluxes, west_detector =\
            read_in_goes(experiment, flux_type, filenames1,
                filenames2, filenames_orien, options)
        
    elif experiment in goes_R and flux_type == "differential":
        all_dates, all_fluxes, west_detector =\
            read_in_goesR(experiment,flux_type, filenames1)
        
    elif (experiment == "GOES_RT") and flux_type == "integral":
        all_dates, all_fluxes, west_detector =\
            read_in_goes_RT(experiment,flux_type, filenames1)

    elif experiment == "EPHIN":
        all_dates, all_fluxes = read_in_ephin(experiment, flux_type, filenames1)

    elif experiment == "EPHIN_REleASE":
        all_dates, all_fluxes = read_in_ephin_release(experiment, flux_type,
                    filenames1)

    elif experiment == "ERNE":
        all_dates, all_fluxes = read_in_erne(experiment, flux_type, filenames1)

    elif "STEREO" in experiment:
        all_dates, all_fluxes = read_in_stereo(experiment, flux_type,
                    filenames1, filenames2)

    elif experiment == "ACE_SIS":
        all_dates, all_fluxes = read_in_ace_sis(experiment, flux_type, filenames1)
        
    elif experiment == "ACE_EPAM_electrons":
        all_dates, all_fluxes = read_in_ace_epam_electrons(experiment, flux_type, filenames1)
        
    elif experiment == "IMP8_CPME":
        all_dates, all_fluxes = read_in_imp8_cpme(experiment, flux_type, filenames1)

    return all_dates, all_fluxes, west_detector


def convert_decimal_hour(decimal_hours):
    """ Convert an hour in decimals to number of hours, minutes,
        and seconds.
    """
    hours = math.floor(decimal_hours)
    
    frac = decimal_hours - float(hours)
    decimal_minutes = frac*60.
    
    minutes = math.floor(decimal_minutes)
    frac = decimal_minutes - float(minutes)
    
    decimal_seconds = frac*60.
    seconds = round(decimal_seconds)
    
    return hours, minutes, seconds


def convert_frac_year(frac_year):
    """ Convert a fractional year date to datetime """
    
    year = int(frac_year)
    frac = frac_year - year
    thisyear = datetime.datetime(year=year, month=1, day=1)
    nextyear = datetime.datetime(year=year+1, month=1, day=1)
    tot_sec = (nextyear - thisyear).total_seconds()
    frac_sec = tot_sec * frac
    
    date = thisyear + datetime.timedelta(seconds=frac_sec)

    return date

def read_in_user_files(filenames1, delim='', flux_col=[], is_unixtime=False):
    """ Read in file containing flux time profile information that
        was specified by the user.
        The first column MUST contain the date in YYYY-MM-DD HH:MM:SS
        format. The remaining flux columns to be read in are specified
        by the user in the variable user_col at the very beginning of
        this program.
        
        The date column should always be considered column 0, even if
        you used whitespace as your delimeter. The code will consider
        the date format YYYY-MM-DD HH:MM:SS as one column even though
        it contains whitespace.
        
        Any number of header lines are allowed, but they must be
        indicated by # at the very beginning, including empty lines.
        Be sure to add the energy bins associated with your flux
        columns in the subroutine define_energy_bins under the "user"
        if statement.
        
        Allow user to add a time shift to files using the global_var
        time_shift to specify the time shift in hours. A negative value
        shifts earlier, a positive value shifts later.
        
        INPUTS:
    
        :filenames1: (string array) the user files containing the data that span the desired time range
        :is_unixtime: (bool) True indicates the dates in the first
            column are in unixtime
                
        OUTPUTS:
       
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and
            m time points
       
    """
    print('Reading in user-specified files.')
    global user_delim
    global user_col

    if delim != "":
        user_delim = delim
        
    if flux_col:
        user_col = flux_col


    if cfg.time_shift != 0:
        print("!!!!!!!Shifting times by time_shift in global_vars.py: "
            + str(cfg.time_shift) + " hours. Set to zero if do "
            "not want to shift.")
    NFILES = len(filenames1)
    ncol = len(user_col) #include column for date
    for i in range(NFILES):
        print('Reading in ' + filenames1[i])
        with open(filenames1[i]) as csvfile:
            #Count header lines indicated by hash #
            nhead = 0
            for line in csvfile:
                line = line.lstrip()
                if line == '' or line == '\n':
                    nhead = nhead + 1
                elif line[0] == "#" or line[0] == '\"':
                    nhead = nhead + 1
                else:
                    break
            #number of lines containing data
            nrow = len(csvfile.readlines())+1

            #Define arrays that hold dates and fluxes
            dates = []
            fluxes = np.zeros(shape=(ncol,nrow))

            csvfile.seek(0) #back to beginning of file
            for k in range(nhead):
                csvfile.readline()  # Skip header rows.

            user_col_mod = []
            for j in range(len(user_col)):
                if (user_delim == " " or user_delim == "") and not is_unixtime:
                    #date takes two columns if separated by whitespace
                    #adjust the user input columns to account for this
                    user_col_mod.append(user_col[j] + 1)
                else:
                    user_col_mod.append(user_col[j])

            count = 0
            for line in csvfile:
                if line == " " or line == "":
                    continue
                if not is_unixtime:
                    if user_delim == " " or user_delim == "":
                        row = line.split()
                        str_date = row[0][0:10] + ' ' + row[1][0:8]
                    if user_delim != " " and user_delim != "":
                        row = line.split(user_delim)
                        str_date = row[0][0:19]
                    date = datetime.datetime.strptime(str_date,
                                "%Y-%m-%d %H:%M:%S")
                
                if is_unixtime:
                    if user_delim == " " or user_delim == "":
                        row = line.split()
                    else: row = line.split(user_delim)
                    utime = int(row[0])
                    date = datetime.datetime.utcfromtimestamp(utime)
                    
                #apply a time shift to user data with the variable
                #set in config.py
                if cfg.time_shift != 0:
                    hours, minutes, seconds = convert_decimal_hour(cfg.time_shift)
                    date = date + datetime.timedelta(hours=hours,
                        minutes=minutes, seconds=seconds)
                    
                dates.append(date)
                
                for j in range(len(user_col_mod)):
                    #print("Read in flux for column " + str(user_col[j]) + ': '\
                    #    + str(date)) #+ ' ' + row[user_col[j]])
                    #print(row)
                    if user_col_mod[j] >= len(row):
                        sys.exit("read_datasets: read_in_user_files: "
                            "Something is wrong with reading in the "
                            "user files (mismatch in number of "
                            "columns). Did you set the correct "
                            "information in config.py, "
                            "including the delimeter?")
                    row[user_col_mod[j]] = row[user_col_mod[j]].rstrip()
                    if row[user_col_mod[j]] == 'n/a': #REleASE
                        flux = None
                    else:
                        flux = float(row[user_col_mod[j]])
                        if flux < 0:
                            flux = badval
                    fluxes[j][count] = flux
                count = count + 1

        #Remove dates that have None values
        for k in range(len(dates)-1,-1,-1):
            if None in fluxes[:,k] or np.isnan(np.sum(fluxes[:,k])):
                del dates[k]
                fluxes = np.delete(fluxes, k, 1)

        #If reading in multiple files, then combine all data into one array
        if i==0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates

    return all_dates, all_fluxes


def extract_date_range(startdate,enddate,all_dates,all_fluxes):
    """ Extract fluxes only for the dates in the range specified by
        the user.
        
        INPUTS:
        
        :startdate: (datetime) start of desired time period
        :enddate: (datetime) end of desired time period
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and
            m time points
        
        OUTPUTS:
        
        :dates: (datetime 1xp array) dates for p time points within
            the date range specified by startdate and enddate
        :fluxes: (float nxp array) flux time profiles for n energy channels
            and p time points
        
    """
    if pd.isnull(startdate) or pd.isnull(enddate):
        return [], []

    #print('Extracting fluxes for dates: ' + str(startdate) + ' to '
    #    + str(enddate))
    ndates = len(all_dates)
    nst = 0
    nend = 0
    for i in range(ndates):
        if all_dates[i] <= startdate:
            nst = i
        if all_dates[i] <= enddate:
            nend = i
    if all_dates[nst] < startdate:
        nst = nst + 1 #move one step past the start time if no
                    #time point on exactly the start time
    nst = min(nst,ndates) #make sure nst is not past the length of the array
    nend = min(nend+1, ndates)  #add 1 to include nend in date range
    
    if nst == nend and nend < ndates:
        nend = nend + 1 #grab at least one data point at index nst
    
    dates = all_dates[nst:nend]
    fluxes = all_fluxes[:,nst:nend]

    return dates, fluxes


def do_interpolation(i,dates,flux):
    """ If bad fluxes (flux < 0) are found in the data, find the first prior
        data point and the first following data point that have good flux values.
        Perform linear interpolation in time:
        
        F(t) = F1 + (t - t1)*(F2 - F1)/(t2 - t1)
        
        This subroutine does the calculation for a single instance of bad data
        that corresponds to array index i.
        
        INPUTS:
        
        :i: (integer) index of time point to interpolate
        :dates: (datetime 1xp array) dates for p time points within
            the date range specified by startdate and enddate
        :fluxes: (float 1xp array) flux time profiles for n energy channels
            and p time points
            
        OUTPUTS:
        
        :interp_flux: (float 1xp array) flux time profile with any negative
            flux values replaced with linear interpolation with time
       
    """
    ndates = len(dates)
    preflux = badval
    postflux = badval
    
#    print("ndates: " + str(ndates) + ", i: " + str(i))

    #If first point is bad point, use the next good point to fill gap
    if i == 0:
        for j in range(i,ndates-1):
            if flux[j] != badval and not pd.isnull(flux[j]):
                postflux = flux[j]
                postdate = dates[j]
#                print('First point in array is bad. The first good value after '
#                    'the gap is on '+ str(dates[j]) + ' with value '
#                    + str(flux[j]))
                break
        preflux = postflux

    #If last point is bad point, use the first prior good point to fill gap
    if i == ndates - 1:
        for j in range(i,-1,-1):
            if flux[j] != badval and not pd.isnull(flux[j]):
                preflux = flux[j]
                predate = dates[j]
 #               print('Last point in the array is bad. The first good value '
 #                   'previous to gap is on '+ str(dates[j]) + ' with value '
 #                   + str(flux[j]))
                break
        postflux = preflux

    #Within the flux array
    if i != 0 and i != ndates-1:
        #search for first previous good value prior to the gap
        for j in range(i,-1,-1):
            if flux[j] != badval and not pd.isnull(flux[j]):
                preflux = flux[j]
                predate = dates[j]
#                print('The first good value previous to gap is on '
#                    + str(dates[j]) + ' with value ' + str(flux[j]))
                break
            if j == 0:
#                print('There is a data gap at the beginning of the '
#                        'selected time period. Program cannot estimate '
#                        f'flux in data gap. Setting to {badval}.')
                preflux = badval
                predate = None

        #search for first previous good value after to the gap
        for j in range(i,ndates-1):
            if flux[j] != badval and not pd.isnull(flux[j]):
                postflux = flux[j]
                postdate = dates[j]
#                print('The first good value after to gap is on '
#                    + str(dates[j]) + ' with value ' + str(flux[j]))
                break
            
            if j == ndates-2 and (flux[j] == badval or pd.isnull(flux[j])):
                if flux[ndates-1] != badval and not pd.isnull(flux[ndates-1]):
                    postflux = flux[ndates-1]
                    postdate = dates[ndates-1]
                else:
                    postflux = preflux
                    postdate = predate
#                    print(' Bad values continue to the end of the data set. '
#                        'Using the first good value previous to gap on '
#                        + str(postdate) + ' with value ' + str(postflux))

            
    if preflux == badval or postflux == badval:
        interp_flux = badval
#        print(f'do_interpolation could not interpolate flux at {i}. Setting to {badval}.')
    
    else:
        if preflux == postflux:
            interp_flux = preflux
        if preflux != postflux:
            interp_flux = preflux + (dates[i] - predate).total_seconds()\
                 *(postflux - preflux)/(postdate - predate).total_seconds()
#        print('Filling gap at time ' + str(dates[i])
#                + ' with interpolated flux ' + str(interp_flux))
    return interp_flux


def check_for_bad_data(dates,fluxes,energy_bins,dointerp=True):
    """ Search the data for bad values (flux < 0) and fill the
        missing data with an estimate flux found by performing a
        linear interpolation with time, using the good flux values
        immediately surrounding the data gap.
        
        INPUTS:
        
        :dates: (datetime 1xp array) dates for p time points within
            the date range specified by startdate and enddate
        :fluxes: (float nxp array) flux time profiles for n energy
            channels and p time points
        :energy_bins: (float nx2 array) energy bins associated with
            fluxes
        :dointerp: (bool) Set True to perform linear interpolation in
            time, otherwise will fill bad data points with None values
            
        OUTPUT:
        
        :fluxes: (float 1xp array) flux time profile with any negative
            flux values replaced with linear interpolated of None values
       
    """
    if dointerp:
        print('Checking for bad data values and filling with linear '
              'interpolation with time.')
    else:
        print('Checking for bad data values and filling with NaN values. ')

    ndates = len(dates)
    nbins = len(energy_bins)

    for j in range(ndates):  #flux at each time
        for i in range(nbins):
            if pd.isnull(fluxes[i,j]): #bad data
                #estimate flux with interpolation in time
                if dointerp:
                    interp_flux = do_interpolation(j,dates,fluxes[i,:])
                    fluxes[i,j] = interp_flux
#                    print('There is null data for time ' + str(dates[j])
#                        + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
#                        + str(energy_bins[i][1]) + '.'
#                        + ' Filling in missing value with linear '
#                        + 'interpolation in time. ' + str(interp_flux))
                else:
#                    print('There is bad data for time ' + str(dates[j])
#                            + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
#                            + str(energy_bins[i][1]) + '.'
#                            + ' Filling in missing value with NaN ')
                    fluxes[i,j] = np.nan

            elif fluxes[i,j] < 0:
                #estimate flux with interpolation in time
                if dointerp:
                    interp_flux = do_interpolation(j,dates,fluxes[i,:])
                    fluxes[i,j] = interp_flux
#                    print('There is negative data for time ' + str(dates[j])
#                            + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
#                            + str(energy_bins[i][1]) + '.'
#                            + ' Filling in missing value with linear '
#                            + 'interpolation in time. ' + str(interp_flux))
                else:
#                    print('There is bad data for time ' + str(dates[j])
#                            + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
#                            + str(energy_bins[i][1]) + '.'
#                            + ' Filling in missing value with NaN ')
                    fluxes[i,j] = np.nan #results in NaN value in np array

    
    print('Finished checking for bad data.')
    print()
    return fluxes


def which_erne(startdate, enddate):
    """ ERNE energy bins depend on the dates of the data.
        Identify the appropriate calibration.
        
        For now, requiring the user to choose dates that do not
        cross time periods with different energy bins.
        f40, f40brk, and f50 actually have all the same energy bins,
        so will require that users select a time period contained within:
        1995-12-02 to 2000-04-18 for f10
        2000-04-19 to present for f40
        
        From export_data_description.txt:
        
        Data versions
        -------------
        Since launch the on-board system has been updated several times. Some of the
        changes have also affected the telemetry data, these are deemed 'format
        versions'. The relevant formats are:

        f10     2 Dec 1995     Original launch format
        f40    19 Apr 2000     Major update of the on-board program
        f40brk 21 Nov 2000     HED S1X H2 E-amplifier breakdown at 00:15:44.833
        f50     3 Jul 2001     On-board SW updated to (partially) fix the HED S1X
                       E-amplifier breakdown.
    
        The format f20 has never existed and f30 was a short lived test version with
        no user level data. The version designation f40brk is unofficial and not
        used outside this document.

        NOTE: all the f40brk HED data are left out from this data set. The
            corresponding files are provided but are empty.

        NOTE: The particles through the broken S1X H2 detector are left out from
            the f50 HED PHA data.
    
        INPUTS:
            
            :dates: (datetime list) Dates for all read in data
            
        OUTPUTS:
        
            :version: (string) ERNEf10, ERNEf40, ERNEf40brk, ERNEf50
                Data version indicating which bins to choose.
                f40, f40brk, and f50 are all the same bins, but dividing
                up to match the documentation
                
    """

    date_f10 = datetime.datetime(1995,12,2)
    date_f40 = datetime.datetime(2000,4,19)
    date_f40brk = datetime.datetime(2000,11,21)
    date_f50 = datetime.datetime(2001,7,3)
    
    #SIMPLE SOLUTION TO GET IT WORK FOR THE MOMENT, MUST EXPAND ON THIS
    if startdate >= date_f10 and enddate < date_f40:
        return "ERNEf10"

    if startdate >= date_f40:
        return "ERNEf40"

    if startdate < date_f40 and enddate >= date_f40:
        sys.exit("ERNE data has different two different sets of energy channels "
                "for the date ranges specified. Please choose all data to be "
                "contained within these date ranges:\n" \
                + str(date_f10) + " to " +  str(date_f40) + "\n" \
                + str(date_f40) + " to present")


def calculate_geometric_means(energy_bins):
    """ Define the bin centers as the geometric means:
        center = sqrt(Elow*Ehigh)
        
        This approach is typically used to define the bin center of 
        proton experiments if not better known.
        
        This should only be applied to differential channels. 
        If a channel is an integral channel, will return the bottom
        energy of the channel.
        
    """
    centers = []
    for bin in energy_bins:
        if bin[1] == -1:
            centers.append(bin[0])
        else:
            centers.append(math.sqrt(bin[0]*bin[1]))
        
    return centers


def define_energy_bins(experiment, flux_type, west_detector, options,
    spacecraft="primary", user_bins=[]):
    """ Define the energy bins for the selected spacecraft or data set.
        If the user inputs their own file, they must set the
        user_energy_bins variable in config/config_opsep.py.
        User may select options to apply Sandberg et al. (2014)
        effective energies for GOES EPS by specifying "S14" and/or
        apply Bruno (2017) effective energies for GOES-13 or -15 P6, P7
        and HEPAD by specifying "Bruno2017"
        
        INPUTS:
        
        :experiment: (string) name of experiment or "user"
        :flux_type: (string) integral or differential
        :west_detector: (string 1xp array) array indicating which GOES
            detector is facing westward for each time point
        :options: (string array) possible options to apply to data
            (GOES)
            :spacecraft: (string) primary or secondary if experiment is GOES_RT
        
        OUTPUTS:
        
        :energy_bins: (float nx2 array) appropriate energy bins for the
            experiment specified by the user
       
    """
    energy_bins = None
    energy_bin_centers = []
    
    #use corrected proton flux for GOES eps or epead; include hepad
    #-1 indicates infinity ([700, -1] means all particles above 700 MeV)
    if experiment == "SEPEM":
        energy_bins = [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                       [21.87,31.62],[31.62,45.73],[45.73,66.13],
                       [66.13,95.64],[95.64,138.3],[138.3,200.0],
                       [200.0,289.2]]
        energy_bin_centers = calculate_geometric_means(energy_bins)

    if experiment == "SEPEMv3":
        energy_bins = [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                       [21.87,31.62],[31.62,45.73],[45.73,66.13],
                       [66.13,95.64],[95.64,138.3],[138.3,200.0],
                       [200.0,289.2],[289.2,418.3],[418.3,604.9],
                       [604.9,874.7]]
        energy_bin_centers = calculate_geometric_means(energy_bins)

    if experiment == "CalGOES":
        energy_bins = [[10.0,10.0],[15.8,15.8],[25.1,25.1],[39.8,39.8],
                        [65.1,65.1],[100.0,100.0],[158.5,158.5],[251.2,251.2],
                        [398.1,398.1],[630.9,630.9],[1000.0,1000.0]]
        energy_bin_centers = [10.0, 15.8, 25.1, 39.8, 65.1, 100.0, 158.5, 251.2,
                        398.1,630.9,1000.0]
   
    if experiment == "ERNEf10":
        #The highest energy 3 channels tend to saturate and show incorrect
        #values during high intensity SEP events. For this reason, only the
        #>10 MeV integral fluxes should be used from ERNE data during the
        #most intense part of intense SEP events.
        #f10 format, from 2 December 1996
        #f10     2 Dec 1995     Original launch format
        energy_bins = [[1.5,1.8],[1.8,2.2],[2.2,2.7],[2.7,3.3],[3.3,4.1],
                       [4.1,5.1],[5.1,6.4],[6.4,8.1],[8.1,10],
                       [10.0,13.0],[14.0,17.0],[17.0,22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[41.0,51.0],
                       [51.0,67.0],[54.0,79.0],[79.0,114.0],[111.0,140.]]
        energy_bin_centers = [1.6, 2.0, 2.4, 3.0, 3.7, 4.6, 5.8, 7.2, 9.1, 11,
                        15.4, 18.9, 23.3, 29.0, 36.4, 45.6, 54.1, 67.5, 95.0, 116]

    if experiment == "ERNEf40":
        #f40 format, from 19 May 2000
        #f40    19 Apr 2000     Major update of the on-board program
        energy_bins = [[1.6,1.8],[1.8,2.2],[2.2,2.7],[2.7,3.3],[3.3,4.1],
                       [4.1,5.1],[5.1,6.4],[6.4,8.1],[8.1,10],
                       [10.0,13.0],[14.0,17.0],[17.0, 22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[40.0,51.0],[51.0,67.0],
                       [64.0,80.0],[80.0,101.0],[101.0,131.0]]
        energy_bin_centers = [1.7, 2.0, 2.4, 3.0, 3.7, 4.7, 5.7, 7.2, 9.1, 11,
                        15.4, 18.9, 23.3, 29.1, 36.4, 45.6, 57.4, 72.0, 90.5, 108]
                       
    if experiment == "ERNEf40brk":
        #f40brk format, from 21 Nov 2000
        #f40brk 21 Nov 2000     HED S1X H2 E-amplifier breakdown at 00:15:44.833
        #NOTE: all the f40brk HED data are left out from this data set. The
        #corresponding files are provided but are empty.
        energy_bins = [[1.6,1.8],[1.8,2.2],[2.2,2.7],[2.7,3.3],[3.3,4.1],
                       [4.1,5.1],[5.1,6.4],[6.4,8.1],[8.1,10],
                       [10.0,13.0],[14.0,17.0],[17.0, 22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[40.0,51.0],[51.0,67.0],
                       [64.0,80.0],[80.0,101.0],[101.0,131.0]]
        energy_bin_centers = [1.7, 2.0, 2.4, 3.0, 3.7, 4.7, 5.7, 7.2, 9.1, 11,
                        15.4, 18.9, 23.3, 29.1, 36.4, 45.6, 57.4, 72.0, 90.5, 108]
 
    if experiment == "ERNEf50":
        #f50 format, from 3 Jul 2001
        #f50     3 Jul 2001     On-board SW updated to (partially) fix the HED S1X
        #               E-amplifier breakdown.
        energy_bins = [[1.6,1.8],[1.8,2.2],[2.2,2.7],[2.7,3.3],[3.3,4.1],
                       [4.1,5.1],[5.1,6.4],[6.4,8.1],[8.1,10],
                       [10.0,13.0],[14.0,17.0],[17.0, 22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[40.0,51.0],[51.0,67.0],
                       [64.0,80.0],[80.0,101.0],[101.0,131.0]]
        energy_bin_centers = [1.7, 2.0, 2.4, 3.0, 3.7, 4.7, 5.7, 7.2, 9.1, 11,
                        15.4, 18.9, 23.3, 29.1, 36.4, 45.6, 57.4, 72.0, 90.5, 108]
                       
    if experiment == "EPHIN":
        #http://ulysses.physik.uni-kiel.de/costep/level3/l3i/
        #DOCUMENTATION-COSTEP-EPHIN-L3-20181002.pdf
        energy_bins = [[4.3,7.8],[7.8,25],[25,40.9],[40.9,53]]
        energy_bin_centers = calculate_geometric_means(energy_bins)

    if experiment == "EPHIN_REleASE":
        #This data may be downloaded by hand through a web interface at:
        #https://www.hesperia.astro.noa.gr/index.php/results/real-time-prediction-tools/data-retrieval-tool
        energy_bins = [[4.0,9.0],[9.0,15.8],[15.8,39.6],[20.0,35.5]]
        energy_bin_centers = calculate_geometric_means(energy_bins)

    ##### GOES SPACECRAFT ########
    #### GOES DOCUMENTATION USES THE GEOMETRIC MEAN AS THE BIN CENTER####
    #GOES-05 EPS starts in 1984-01-01 (magneto only previously)
    #Available files do not contain integral fluxes
    if experiment == "GOES-05":
        if flux_type == "differential":
            energy_bins = [[4.2,8.7],[8.7,14.5],[15.0,44.0],
                           [39.0,82.0],[84.0,200.0],[110.0,500.0]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
            if "S14" in options:
                    energy_bins = [[5.0,7.6],[8.7,13.1],[14.2,20.7],
                               [40.1,56.5],[96.3,133.0],[177.0,247.0],
                               [350.0,420.0],[420.0,510.0],[510.0,700.0],
                               [700.0,-1]]
                    energy_bin_centers = [6.3, 11.1, 17.9, 48.7, 114.0, 218.0,
                                math.sqrt(350.0*420.0), math.sqrt(420.0*510.0),
                                math.sqrt(510.0*700.0), 700]
        if (flux_type == "integral"):
             energy_bins = []
             energy_bin_centers = []




    if experiment == "GOES-06":
        #The highest energy bin is >685 MeV. For the
        #purposes of this code, it will be labeled as
        #>700 MeV and treated as consistent with the
        #other GOES spacecraft that specify >700 MeV
        print("define_energy_bins: Note that GOES-06 highest energy channel "
            "is >685 MeV. Setting to >700 MeV here to enable combining with other "
            "GOES spacecraft if selecting --Experiment GOES option.")
        if flux_type == "differential":
            energy_bins = [[4.2,8.7],[8.7,14.5],[15.0,44.0],
                           [39.0,82.0],[84.0,200.0],[110.0,500.0],
                           [375, 375],[465,465],[605,605], [700,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
        if (flux_type == "integral"):
            energy_bins = [[5.0,-1],[10.0,-1],[30.0,-1],
                            [50.0,-1],[60.0,-1],[100.0,-1], [700,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)


    if experiment == "GOES-07":
        print("define_energy_bins: Note that no HEPAD data is available at NCEI "
            "for GOES-07.")
        if flux_type == "differential":
            energy_bins = [[4.2,8.7],[8.7,14.5],[15.0,44.0],
                           [39.0,82.0],[84.0,200.0],[110.0,500.0]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
            if "S14" in options:
                    energy_bins = [[4.4,8.2],[7.8,14.6],[18.0,24.2],
                               [41.0,57.3],[89.5,139.0],[166.0,299.0]]
                    energy_bin_centers = [6.6, 11.2, 21.1, 50.5, 114.0, 243.0]
        if (flux_type == "integral"):
            energy_bins = [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)


    if (experiment == "GOES-08" or experiment == "GOES-09"
        or experiment == "GOES-10" or experiment == "GOES-11"
        or experiment == "GOES-12"):
        if (flux_type == "integral"):
            energy_bins = [[5.0,-1],[10.0,-1],[30.0,-1],
                            [50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
        if (flux_type == "differential"):
            #files named e.g. g08_eps_5m_yyyymmdd_yyyymmdd.csv
            energy_bins = [[4.0,9.0],[9.0,15.0],[15.0,44.0],
                           [40.0,80.0],[80.0,165.0],[165.0,500.0],
                           [350.0,420.0],[420.0,510.0],[510.0,700.0],
                           [700.0,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
            if "S14" in options:
                if experiment == "GOES-08":
                    energy_bins = [[4.0,7.9],[7.4,15.0],[13.3,21.3],
                               [37.0,53.6],[91.5,113.0],[119.0,179.0],
                               [350.0,420.0],[420.0,510.0],[510.0,700.0],
                               [700.0,-1]]
                    energy_bin_centers = [6.05, 10.6, 19.0, 47.8, 107.0, 153.0,
                                math.sqrt(350.0*420.0), math.sqrt(420.0*510.0),
                                math.sqrt(510.0*700.0), 700]
                if experiment == "GOES-11":
                    energy_bins = [[5.0,7.9],[9.4,15.9],[16.7,23.2],
                               [32.5,56.4],[89.8,114.0],[120.0,186.0],
                               [350.0,420.0],[420.0,510.0],[510.0,700.0],
                               [700.0,-1]]
                    energy_bin_centers = [6.4, 12.5, 20.8, 46.1, 104.0, 148.0,
                                math.sqrt(350.0*420.0), math.sqrt(420.0*510.0),
                                math.sqrt(510.0*700.0), 700]

    if (experiment == "GOES-13" or experiment == "GOES-14" or
        experiment == "GOES-15"):
        if (flux_type == "integral"):
            energy_bins = [[5.0,-1],[10.0,-1],[30.0,-1],
                            [50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
        if (flux_type == "differential"):
            energy_bins = [[4.2,8.7],[8.7,14.5],[15.0,40.0],
                            [38.0,82.0],[84.0,200.0],[110.0,900.0],
                            [330.0,420.0],[420.0,510.0],[510.0,700.0],
                            [700.0,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
            if "S14" in options:
                #S14 is not specifically calibrated to these GOES, but
                #apply the GOES-11 EPS energy bins for P2 - P7
                energy_bins = [[5.0,7.9],[9.4,15.9],[16.7,23.2],
                           [32.5,56.4],[89.8,114.0],[120.0,186.0],
                           [330.0,420.0],[420.0,510.0],[510.0,700.0],
                           [700.0,-1]]
                energy_bin_centers = [6.4, 12.5, 20.8, 46.1, 104.0, 148.0,
                            math.sqrt(350.0*420.0), math.sqrt(420.0*510.0),
                            math.sqrt(510.0*700.0), 700]
            if "Bruno2017" in options:
                #EPEAD CHANNELS P6 and P7
                if "uncorrected" in options:
                    if west_detector.count("A") >= west_detector.count("B"):
                        print("Choosing Bruno2017 energy bins for uncorrected A detector.")
                        #A detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [93.3,129.0]
                            energy_bin_centers[4] = 114.4
                            energy_bins[5] = [145.2,203.9]
                            energy_bin_centers[5] = 179.6
                        if experiment == "GOES-15":
                            energy_bins[4] = [97.5,134.8]
                            energy_bin_centers[4] = 119.5
                            energy_bins[5] = [142.2,199.0]
                            energy_bin_centers[5] = 175.5
                    if west_detector.count("B") > west_detector.count("A"):
                        print("Choosing Bruno2017 energy bins for uncorrected B detector.")
                        #B detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [92.3,127.5]
                            energy_bin_centers[4] = 113.1
                            energy_bins[5] = [143.5,200.5]
                            energy_bin_centers[5] = 177.0
                        if experiment == "GOES-15":
                            energy_bins[4] = [95.9,132.3]
                            energy_bin_centers = 117.4
                            energy_bins[5] = [144.6,202.3]
                            energy_bin_centers[5] = 178.5
                if "corrected" in options or "uncorrected" not in options: #Z89 applied
                    if west_detector.count("A") >= west_detector.count("B"):
                        print("Choosing Bruno2017 energy bins for corrected A detector.")
                        #A detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [93.3,129.1]
                            energy_bin_centers[4] = 114.5
                            energy_bins[5] = [146.7,205.6]
                            energy_bin_centers[5] = 181.2
                        if experiment == "GOES-15":
                            energy_bins[4] = [97.9,135.3]
                            energy_bin_centers[4] = 120.0
                            energy_bins[5] = [145.0,202.3]
                            energy_bin_centers[5] = 178.6
                    if west_detector.count("B") > west_detector.count("A"):
                        print("Choosing Bruno2017 energy bins for corrected B detector.")
                        #B detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [92.4,127.8]
                            energy_bin_centers[4] = 113.3
                            energy_bins[5] = [144.6,202.5]
                            energy_bin_centers[5] = 178.5
                        if experiment == "GOES-15":
                            energy_bins[4] = [96.1,132.8]
                            energy_bin_centers[4] = 117.7
                            energy_bins[5] = [147.3,205.7]
                            energy_bin_centers[5] = 181.5
                if experiment == "GOES-13": #Should use with bg-subtracted flux
                    energy_bins[6] = [273.9,387.5] #P8
                    energy_bin_centers[6] = 337.3
                    energy_bins[7] = [330.0,458.0] #P9
                    energy_bin_centers[7] = 407.1
                    energy_bins[8] = [418.7,566.0] #P10
                    energy_bin_centers[8] = 508.8
                    energy_bins[9] = [852.6,1081.2] #P11
                    energy_bin_centers[9] = 1002.4
                if experiment == "GOES-15":
                    energy_bins[6] = [240.4,335.6] #P8
                    energy_bin_centers[6] = 297.3
                    energy_bins[7] = [325.3,464.6] #P9
                    energy_bin_centers[7] = 407.4
                    energy_bins[8] = [420.4,573.1] #P10
                    energy_bin_centers[8] = 516.0
                    energy_bins[9] = [878.6,1230.0] #P11
                    energy_bin_centers[9] = 1094.7


    if experiment in goes_R:
        if flux_type == "differential":
            energy_bins = [[1.02,1.86],[1.9,2.3],[2.31,3.34],
                           [3.4,6.48],[5.84,11.0],[11.64,23.27],
                           [24.9,38.1],[40.3,73.4],[83.7,98.5],
                           [99.9,118.0],[115.0,143.0],[160.0,242.0],
                           [276.0,404.0],[500.0,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)
        if flux_type == "integral":
            energy_bins = [[1,-1],[5,-1],[10,-1],[30,-1],[50,-1],[100,-1],
                            [60,-1],[500,-1]]
            energy_bin_centers = calculate_geometric_means(energy_bins)


    if experiment == "GOES_RT":
        if flux_type == "integral":
            if spacecraft == "primary":
                energy_bins = [[1,-1],[5,-1],[10,-1],[30,-1],[50,-1],[100,-1],
                                [60,-1],[500,-1]]
            if spacecraft == "secondary":
                energy_bins = [[1,-1],[5,-1],[10,-1],[30,-1],[50,-1],[100,-1],
                                [60,-1],[500,-1]]
                 #GOES_RT when spacecraft is prior to GOES-R will have upper
                 #energy bins of >700 MeV!!!!!
            energy_bin_centers = calculate_geometric_means(energy_bins)

    
    if experiment == "STEREO-A" or experiment == "STEREO-B":
        #Uses the SUMMED LET bins with the HET bins
        energy_bins = [[1.8,3.6],[4.0,6.0],[6.0,10.0],[13.6,15.1],
                        [14.9,17.1],[17.0,19.3],[20.8,23.8],
                        [23.8,26.4],[26.3,29.7],[29.5,33.4],[33.4,35.8],
                        [35.5,40.5],[40.0,60.0],[60.0,100.0]]


    if experiment == "ACE_SIS":
        #https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/
        energy_bins = [[30,-1],[60,-1]]
        energy_bin_centers = calculate_geometric_means(energy_bins)

    if experiment == "ACE_EPAM_electrons":
        #https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/
        energy_bins = [[0.175,0.315]]
        energy_bin_centers = calculate_geometric_means(energy_bins)

    if experiment == "IMP8_CPME":
        #http://sd-www.jhuapl.edu/IMP/data/imp8/cpme/cpme_330s/protons/
        energy_bins = [[4.60, 15.0], [15.0, 25.0], [25.0, 48.0], [48.0, 96.0],
                        [96.0, 145.0], [145.0, 440.0]]
        energy_bin_centers = calculate_geometric_means(energy_bins)


    if experiment == "user":
        #modify to match your energy bins or integral channels
        #use -1 in the second edge of the bin for integral channel (infinity)
        if user_bins:
            energy_bins = user_bins #input into subroutine
        else:
            energy_bins = cfg.user_energy_bins #global from cfg
        energy_bin_centers = calculate_geometric_means(energy_bins)



    return energy_bins, energy_bin_centers
