from ..utils import config as cfg
import datetime
import argparse
import os
import sys


__version__ = "0.1"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"



def read_detector_list(detector_list):
    """ Read in the detector list with columns:
        year,month,goes detector
        
        comma separated
        
    """
    detector_dates = []
    detector = []
    with open(detector_list, 'r') as detfile:
        for line in detfile:
            if "#" in line: continue
            if line == "": continue
            line = line.strip().split(",")
            year = int(line[0])
            month = int(line[1])
            det = line[2]
            
            date = datetime.datetime(year=year,month=month,day=1)
            detector_dates.append(date)
            detector.append(det)
    
    return detector_dates, detector



def read_sep_times_file(septimes_file):
    """ Read in SEP times file with columns
        Start Time End Time
        
        white space separated
        
    """
    sep_sttimes = []
    sep_endtimes = []
    with open(septimes_file,'r') as timesfile:
        for line in timesfile:
            if "#" in line: continue
            if line == "": continue
            
            line = line.strip().split()
            stdt = line[0] + " " + line[1]
            enddt= line[2] + " " + line[3]
            
            stdate = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M:%S')
            enddate = datetime.datetime.strptime(enddt, '%Y-%m-%d %H:%M:%S')

            sep_sttimes.append(stdate)
            sep_endtimes.append(enddate)
            
    return sep_sttimes, sep_endtimes
    

def revise_sep_times(startdate, enddate, sep_sttimes, sep_endtimes):
    """ Adjust the start and end times of each SEP event (enhancement) to
        include more time before and after the event.
        
        Will check to make sure that events did not happen in quick
        succession.
        
        This subroutine is useful if SEPs will be run through OpSEP which
        benefits from some lead time prior to the SEP enhancement.
        
        INPUT:
        
            :sep_sttimes: (1xn datetime array) n SEP start times
            :sep_endtimes: (1xn datetime array) n SEP end times
            
        OUTPUT:
        
            :rev_sttimes: (1xn datetime array) n revised SEP start times
            :rev_endtimes: (1xn datetime array) n revised SEP end times

    """
    rev_sttimes = sep_sttimes
    rev_endtimes = sep_endtimes
    
    dt = datetime.timedelta(hours=6)
    dt2 = dt + dt
    
    #First SEP start and last SEP end
    diff = sep_sttimes[0] - startdate
    if diff > dt:
        rev_sttimes[0] = sep_sttimes[0] - dt
    diff = enddate - sep_endtimes[-1]
    if diff > dt:
        rev_endtimes[-1] = sep_endtimes[-1] + dt

    #Remaining SEP times
    for i in range(len(sep_sttimes)-1,0,-1):
        diff = sep_sttimes[i] - sep_endtimes[i-1]
        if diff > dt:
            rev_sttimes[i] = sep_sttimes[i] - dt
    
        if diff > dt2:
            rev_endtimes[i-1] = sep_endtimes[i-1] + dt

    return rev_sttimes, rev_endtimes


def get_detector(date, detector_dates, detector):
    """ For a given date, find the right GOES detector

    """
    year = date.year
    month = date.month
    checkdate = datetime.datetime(year=year,month=month,day=1)
    
    try:
        idx = detector_dates.index(checkdate)
    except:
        print("get_detector: Requested date not in list.")
        return None
    
    return detector[idx]


def make_event_list(str_startdate, str_enddate, septimes_file,
            detector_list, experiment, flux_type, options, outfile,
            revise=False):
    """ Read in SEP times and detector list, if applicable,
        and create a file in the appropriate format to run
        operational_sep_quantities.py in batch mode.
        
        Output file will create a file that has the
        SEP event start and end times, but will also
        list the time between the end of one SEP event to
        the start of the next.
        
        In the file, Flags and Options are separated by a
        semi-colon. Flags includes TwoPeaks, SubtractBG, and
        DetectPreviousEvent.Options includes S14, Bruno2017,
        corrected, uncorrected.
        
        Set revise = True to adjust the SEP start and end times
        to include some padding at the front of the event. This is
        useful if submitting times to OpSEP for processing.
        
    """
    try:
        startdate = datetime.datetime.strptime(str_startdate, "%Y-%m-%d %H:%M:%S")
    except:
        startdate = datetime.datetime.strptime(str_startdate, "%Y-%m-%d")
    try:
        enddate = datetime.datetime.strptime(str_enddate, "%Y-%m-%d %H:%M:%S")
    except:
        enddate = datetime.datetime.strptime(str_enddate, "%Y-%m-%d")


    header = "#Start Time,End Time,Experiment,Flux Type,Flags," \
            + "User Experiment Name,User Filename,Options,"\
            +"BGStart,BGEnd\n"
    
    out = open(cfg.outpath + "/idsep/" + outfile,'w')
    out.write(header)
    
    if detector_list != '':
        detector_dates, detector = read_detector_list(detector_list)
    
    sep_sttimes, sep_endtimes = read_sep_times_file(septimes_file)
    if revise:
        sep_sttimes, sep_endtimes = revise_sep_times(startdate, enddate,
                sep_sttimes, sep_endtimes)
    
    ####Start of full time range to first SEP time
    line = ''
    if detector_list != '':
        det = get_detector(startdate, detector_dates, detector)
        if det != None:
            line = str_startdate + "," + str(sep_sttimes[0])\
                + "," + det + "," + flux_type + ",,,," + options + ",,\n"
    else:
        line = str_startdate + "," + str(sep_sttimes[0]) + "," \
            + experiment + "," + flux_type + ",,,," + options + ",,\n"
    
    out.write(line)
    
    ###SEP Times
    for i in range(len(sep_sttimes)):
        #Print out line for SEP event
        if detector_list != '':
            det = get_detector(sep_sttimes[i], detector_dates, detector)
            if det == None: continue
            line = str(sep_sttimes[i]) + "," + str(sep_endtimes[i]) + ","\
                + det + "," + flux_type + ",,,," + options + ",,\n"
        else:
            line = str(sep_sttimes[i]) + "," + str(sep_endtimes[i]) + "," \
                + experiment + "," + flux_type + ",,,," + options + ",,\n"
            
        out.write(line)
        
        #Print out line for time period between SEP events
        if i == len(sep_sttimes) - 1: continue
        
        if detector_list != '':
            det = get_detector(sep_endtimes[i], detector_dates,
                    detector)
            if det == None: continue
            line = str(sep_endtimes[i]) + "," + str(sep_sttimes[i+1])\
                + "," + det + "," + flux_type + ",,,," + options + ",,\n"
        else:
            line = str(sep_endtimes[i]) + "," + str(sep_sttimes[i+1]) + ","\
                + experiment + "," + flux_type + ",,,," + options + ",,\n"
        
        out.write(line)


    ####End of full time range from last SEP time
    line = ''
    if detector_list != '':
        det = get_detector(enddate, detector_dates, detector)
        if det != None:
            line = str(sep_endtimes[-1]) + "," + str_enddate\
                + "," + det + "," + flux_type + ",,,," + options + ",,\n"
    else:
        line = str(sep_endtimes[-1]) + "," + str_enddate + "," \
            + experiment + "," + flux_type + ",,,," + options + ",,\n"
    
    out.write(line)
        
    out.close()
            
    
