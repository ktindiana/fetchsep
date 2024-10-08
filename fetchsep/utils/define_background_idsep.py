from ..utils import config as cfg
from . import date_handler as dateh
import re
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#2021-10-14, v0.1: This subroutine will attempt to
#   automatically identify background levels
#2021-12-07, v0.2: Using this version to generate
#   SEP and high point event lists for BON Extension
#   project.


def remove_none_one(vals):
    ''' Remove None values.
        Only values that are real are kept.
        
        INPUTS:

        :vals: (1xn float array)
        
        OUTPUTS:
        
        :clean_vals: (1xm float array) None values removed
        
    '''

    bad_index = [bad for bad, value in enumerate(vals) if value == None]
    vals_clean = list(vals)
    for bad in sorted(bad_index, reverse=True):
        del vals_clean[bad]

    return vals_clean, bad_index
    
    
def remove_nan_one(vals):
    ''' Remove nan values.
        Only values that are real are kept.
        
        INPUTS:

        :vals: (1xn float array)
        
        OUTPUTS:
        
        :clean_vals: (1xm float array) None values removed
        
    '''

    bad_index = [bad for bad, value in enumerate(vals) if math.isnan(value)]
    vals_clean = list(vals)
    for bad in sorted(bad_index, reverse=True):
        del vals_clean[bad]

    return vals_clean, bad_index
    
   
   
   
def remove_zero_one(vals):
    ''' Remove None values.
        Only values that are real are kept.
        
        INPUTS:

        :vals: (1xn float array)
        
        OUTPUTS:
        
        :clean_vals: (1xm float array) None values removed
        
    '''

    bad_index = [bad for bad, value in enumerate(vals) if value == 0]
    vals_clean = list(vals)
    for bad in sorted(bad_index, reverse=True):
        del vals_clean[bad]

    return vals_clean, bad_index



def remove_none_sigma(vals, sigma_low, sigma_high):
    ''' Remove None values.
        Only values that are real are kept.
        
        INPUTS:

        :vals: (1xn float array)
        
        OUTPUTS:
        
        :clean_vals: (1xm float array) None values removed
        
    '''

    bad_index = [bad for bad, value in enumerate(vals) if value == None]
    vals_clean = list(vals)
    sigma_low_clean = list(sigma_low)
    sigma_high_clean = list(sigma_high)
    for bad in sorted(bad_index, reverse=True):
        del vals_clean[bad]
        del sigma_low_clean[bad]
        del sigma_high_clean[bad]

    return vals_clean, sigma_low_clean, sigma_high_clean
    
def remove_zero_sigma(vals,sigma_low, sigma_high):
    ''' Remove None values.
        Only values that are real are kept.
        
        INPUTS:

        :vals: (1xn float array)
        
        OUTPUTS:
        
        :clean_vals: (1xm float array) None values removed
        
    '''

    bad_index = [bad for bad, value in enumerate(vals) if value == 0]
    vals_clean = list(vals)
    sigma_low_clean = list(sigma_low)
    sigma_high_clean = list(sigma_high)
    for bad in sorted(bad_index, reverse=True):
        del vals_clean[bad]
        del sigma_low_clean[bad]
        del sigma_high_clean[bad]

    return vals_clean, sigma_low_clean, sigma_high_clean


def remove_none(dates, vals):
    ''' Remove None values from corresponding dates and vals array.
        Only values that are real are kept.
        
        INPUTS:
        
        :dates: (1xn datetime array)
        :vals: (1xn float array) values associated with the dates
        
        OUTPUTS:
        
        :clean_dates: (1xm datetime array) dates with None values removed
        :clean_vals: (1xm float array) None values removed
        
    '''
    #Error checking
    if len(dates) != len(vals):
        sys.exit('remove_none: Both input arrays must be the same length! '
                'Exiting.')
    #Clean None values from observations and remove correponding entries in
    #the model
    bad_index = [bad for bad, value in enumerate(vals) if value == None]
    dates_clean = list(dates)
    vals_clean = list(vals)
    for bad in sorted(bad_index, reverse=True):
        del dates_clean[bad]
        del vals_clean[bad]

    return dates_clean, vals_clean



def remove_above_value(value, dates, vals):
    ''' Remove fluxes above "value" from corresponding dates and
        vals array.
  
        
        INPUTS:
        
        :value: (float) flux value threshold
        :dates: (1xn datetime array)
        :vals: (1xn float array) values associated with the dates
        
        OUTPUTS:
        
        :clean_dates: (1xm datetime array) dates with None values removed
        :clean_vals: (1xm float array) None values removed
        
    '''
    #Error checking
    if len(dates) != len(vals):
        sys.exit('remove_above_value: Both input arrays must be the same length! '
                'Exiting.')
    #Clean values > value from observations and remove correponding entries in
    #the model
    bad_index = [bad for bad, val in enumerate(vals) if val > value]
    dates_clean = list(dates)
    vals_clean = list(vals)
    for bad in sorted(bad_index, reverse=True):
        del dates_clean[bad]
        del vals_clean[bad]

    return dates_clean, vals_clean


def remove_above_value_one(value, vals):
    ''' Remove fluxes above "value" from corresponding dates and
        vals array.
  
        
        INPUTS:
        
        :value: (float) flux value threshold
        :vals: (1xn float array) values associated with the dates
        
        OUTPUTS:
        
        :clean_vals: (1xm float array) None values removed
        
    '''

    #Clean values > value from observations and remove correponding entries in
    #the model
    n = len(vals)
    vals_clean = list(vals)
    bad_index = []
    for i in range(n-1, -1, -1):
        val = vals[i]
        if val == None: continue
        if val > value:
            del vals_clean[i]
            if bad_index == []:
                bad_index = [i]
            else:
                bad_index.append(i)
    
    
    
#    bad_index = [bad for bad, val in enumerate(vals) if val > value]
#    vals_clean = list(vals)
#    for bad in sorted(bad_index, reverse=True):
#        print("deleted: " + str(vals_clean[bad]))
#        del vals_clean[bad]

    return vals_clean, bad_index


def separate_with_threshold(threshold, dates, vals):
    ''' Separate fluxes above and below a threshold value
  
        
        INPUTS:
        
        :threshold: (1xn float array) flux value threshold
        :dates: (1xn datetime array)
        :vals: (1xn float array) values associated with the dates
        
        OUTPUTS:
        
        :flux_below: (1xn float array) fluxes below threshold, with
            time periods where fluxes above threshold set to zero
        :flux_above: (1xn float array) fluxes above threshold, with
            time periods where fluxes below threshold set to zero
        
    '''
    #Error checking
    if len(dates) != len(vals):
        sys.exit('separate_with_threshold: Dates and values must be the same length! Exiting.')
    if len(dates) != len(threshold):
        sys.exit('separate_with_threshold: Threshold values must correspond to the dates array and both arrays must be the same length! Exiting.' +
            'dates: ' + str(len(dates)) + ', threshold: ' + str(len(threshold)))
    
    #Arrays above and below threshold
    npts = len(dates)
    flux_above = [0]*npts
    flux_below = [0]*npts
    
    for i in range(npts):
        val = vals[i]
        if val == None or math.isnan(val): continue
        if val <= threshold[i]:
            flux_below[i] = val
        else:
            flux_above[i] = val

    return flux_below, flux_above



def separate_with_dates(dates, vals, starttimes, endtimes, padding):
    ''' Separate fluxes above and below with a set of dates.
        Assume that the start and end dates indicate enhanced fluxes,
        e.g. SEP start and end times.
  
        
        INPUTS:
        
        :dates: (1xn datetime array)
        :vals: (1xn float array) values associated with the dates
        :starttimes: (1xm datetime array) list of start times for removing
            time periods from vals
        :endtimes: (1xm datetime array) list of end times for removing
            time periods from vals
        :padding: (int) number of days on either side of SEP start and end
            to exclude
        
        OUTPUTS:
        
        :flux_below: (1xn float array) fluxes below threshold, with
            time periods where fluxes above threshold set to zero
        :flux_above: (1xn float array) fluxes above threshold, with
            time periods where fluxes below threshold set to zero
    '''
    #Error checking
    if len(dates) != len(vals):
        sys.exit('separate_with_dates: Dates and values must be the same length! Exiting.')

    if len(starttimes) != len(endtimes):
        sys.exit('separate_with_dates: Start and end times must be the same '
            'lengths. Start: ' + str(starttimes) + ', END: ' + str(endtimes) +
            '. Exiting!!')
  
    #Arrays above and below threshold
    npts = len(dates)
    flux_above = [0]*npts #default to 0
    flux_below = np.array(vals)
    
    for j in range(len(starttimes)):
        indices = np.argwhere((np.array(dates) >= (starttimes[j]-datetime.timedelta(days=padding))) &
                (np.array(dates) < (endtimes[j]+datetime.timedelta(days=padding))))
        for ix in indices:
            flux_above[ix[0]] = vals[ix[0]]
            flux_below[ix[0]] = 0.

    return flux_below, flux_above



def monthly_average(dates, fluxes):
    """ Make a monthly average of flux.
    
        INPUTS:
        
        :dates: (1xn datetime array)
        :fluxes: (pxn float array) fluxes for p energy channels
        
        OUTPUTS:
        
        :monthly_dates: (1xm datetime array) dates at the midpoint of each
            time period
        :monthly_fluxes: (pxm float array) average flux for the month
        :monthly_sigma: (px2xm float array) sigma of flux for the month
        
    """
    #If entered array of fluxes for multiple energy channels
    nchan = len(fluxes)
    
    #If entered a 1D array of fluxes
    if len(dates) == len(fluxes): nchan = 1
    
    monthly_dates = []
    monthly_fluxes = [[]]*nchan
    monthly_sigma = [[]]*nchan
    
    year = dates[0].year
    month = dates[0].month
    day = 1
    
    #Average the flux between the 1st of each month
    firstdate = datetime.datetime(year=year,month=month,day=day)
    nextdate = dateh.get_next_month(firstdate)
    
    #Align the average flux with the middle of the averaged time period
    savedate = datetime.datetime(year=year,month=month,day=15)
    
    stidx = 0
    endidx = 0
    for i in range(len(dates)):
        if dates[i]<nextdate:
            endidx = i
        else:
            monthly_dates.append(savedate)
            for j in range(nchan):
                #Calculate and save needed values
                flux_arr = np.array(fluxes[j][stidx:endidx+1])
                flux_arr, bad_index = tools.remove_none(flux_arr)
                flux_arr, bad_index = remove_nan_one(flux_arr)
                flux_arr, bad_index = remove_zero_one(flux_arr)
                if flux_arr == []:
                    mean = 0
                    sigma = 0
                    sigma_low = 0
                    sigma_high = 0
                else:
                    #Convert to log space
                    flux_arr = [math.log10(x) for x in flux_arr]
                    mean = np.mean(flux_arr)
                    sigma = np.std(flux_arr)
                    
                    #transform to error bars that make sense
                    #in linear space
                    val_low = mean - sigma
                    val_high = mean + sigma
                    
                    sigma_low = 10**mean - 10**val_low
                    sigma_high = 10**val_high - 10**mean
                    mean = 10**mean
                
                if not monthly_fluxes[j]:
                    monthly_fluxes[j] = [mean]
                    monthly_sigma[j] = [[sigma_low],[sigma_high]]
                else:
                    monthly_fluxes[j].append(mean)
                    monthly_sigma[j][0].append(sigma_low)
                    monthly_sigma[j][1].append(sigma_high)
            
            #Advance to next month
            stidx = i
            firstdate = dates[i]
            nextdate = dateh.get_next_month(firstdate)
            savedate = datetime.datetime(year=firstdate.year,\
                                month=firstdate.month,day=15)
    
    #If entered a 1D array of fluxes, return a 1D array of fluxes
    if nchan == 1:
        monthly_fluxes = monthly_fluxes[0]
        monthly_sigma = monthly_sigma[0]
        
    return monthly_dates, monthly_fluxes, monthly_sigma


def find_last_good(idx, arr):
    """ Find the closest previous value that was non-zero.
    
        INPUTS:
        
        :idx: (int) index of location in array
        :arr: (1xn float array) array containing values
        
        OUTPUTS:
        
        :i: (int) location of the closest previous index
            that corresponds to a non-zero value in the array
    """
    
    for i in range(idx-1,-1,-1):
        if arr[i] != 0 and arr[i] != None:
            return i

    return None



def ndays_average(N, dates, fluxes, nsigma, remove_above):
    """ Average flux over N days.
    
        INPUTS:
        
        :N: (integer) number of days over which to average
        :dates: (1xn datetime array)
        :fluxes: (pxn float array) fluxes for p energy channels
        :remove_above: (float) - remove all flux values above this value
            (assumed that these are bad or incorrect fluxes)
        
        OUTPUTS:
        
        :ave_dates: (1xm datetime array) dates at the midpoint of each
            time period
        :ave_fluxes: (pxm float array) average flux for the N days
        :ave_sigma: (px2xm float array) sigma of flux for the N days
        
    """
    #If entered array of fluxes for multiple energy channels
    nchan = len(fluxes)
    
    #If entered a 1D array of fluxes
    if len(dates) == len(fluxes): nchan = 1
    
    ave_dates = []
    ave_fluxes = [[]]*nchan
    ave_sigma = [[]]*nchan
    threshold = [[]]*nchan
    threshold_dates = []
    
    year = dates[0].year
    month = dates[0].month
    day = dates[0].day
    
    #Average the flux between the 1st of each month
    firstdate = datetime.datetime(year=year,month=month,day=day)
    td = datetime.timedelta(days=N)
    nextdate = firstdate + td
    
    #Align the average flux with the middle of the averaged time period
    td2 = datetime.timedelta(days=int(N/2.))
    savedate = firstdate + td2
    
    stidx = 0
    endidx = 0
    for i in range(len(dates)):
        if dates[i]<nextdate:
            endidx = i
        else:
            if nextdate == dates[-1]: endidx = len(dates)-1
            ave_dates.append(savedate)
            for j in range(nchan):
                #Calculate and save needed values
                flux_arr = np.array(fluxes[j][stidx:endidx+1])
                flux_arr, bad_index = remove_none_one(flux_arr)
                flux_arr, bad_index = remove_nan_one(flux_arr)
                flux_arr, bad_index = remove_zero_one(flux_arr)
                flux_arr, bad_index = remove_above_value_one(remove_above,flux_arr)
                if flux_arr == []:
                    if stidx > 0: #Use previous good value
                        mean = ave_fluxes[j][-1]
                        sigma_low = ave_sigma[j][0][-1]
                        sigma_high = ave_sigma[j][1][-1]
                    else:
                        mean = 0
                        sigma_low = 0
                        sigma_high = 0
                else:
                    mean = np.mean(flux_arr)
                    sigma = np.std(flux_arr)
                    sigma_low = sigma
                    sigma_high = sigma
                
                if not ave_fluxes[j]:
                    ave_fluxes[j] = [mean]
                    ave_sigma[j] = [[sigma_low],[sigma_high]]
                else:
                    ave_fluxes[j].append(mean)
                    ave_sigma[j][0].append(sigma_low)
                    ave_sigma[j][1].append(sigma_high)
            
            #Fill in threshold - one value for each of the input dates
            #mean + n*sigma
#            print("stidx " + str(stidx) + ", endidx " + str(endidx) +
#                ",stdate " + str(dates[stidx]) + ", enddate " + str(dates[endidx])
#                +", ave flux: " + str(ave_fluxes[j][-1]) + ", sigma " +
#                str(ave_sigma[j][1][-1]))
            for k in range(stidx, endidx+1):
                threshold_dates.append(dates[k])
                for j in range(nchan):
                    if not threshold[j]:
                        threshold[j] = [ave_fluxes[j][-1] + nsigma*ave_sigma[j][1][-1]]
                    else:
                        threshold[j].append(ave_fluxes[j][-1] + nsigma*ave_sigma[j][1][-1])
                        
            
            
            #Advance to next month
            stidx = i
            firstdate = dates[i]
            td = datetime.timedelta(days=N)
            nextdate = firstdate + td
            if nextdate > dates[-1]: nextdate = dates[-1]
            td2 = datetime.timedelta(days=int(N/2.))
            savedate = firstdate + td2
            if savedate > dates[-1]: savedate = dates[-1]
        
    return ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold




def backward_window_background(N, dates, fluxes, nsigma):
    """ Average over a backward sliding window of N days.
        Estimate the value of the mean background (GCR) flux,
        sigma, and a threshold to separate GCR from SEP for
        every single time step of the data set.
        
        Expect that fluxes contain mostly background values
        and have already had some portion of SEPs excluded.
        This is meant to be a second iteration that further
        cleans the first attempt to remove SEPs.
    
        INPUTS:
        
        :N: (integer) number of days to smooth over
        :dates: (1xn datetime array)
        :fluxes: (pxn float array) fluxes for p energy channels
        :nsigma: (float) number of sigma to calculate threshold
        
        OUTPUTS:
        
        :smooth_dates: (1xm datetime array) currently unchanged
        :mean_background: (pxm float array) average flux for the month
        :ave_sigma: (px2xm float array) sigma of flux for the month
        
    """
    print("backward_window_background: Calculating the background using "
        "a backward smoothing window of " + str(N) + " days.")
    #If entered array of fluxes for multiple energy channels
    nchan = len(fluxes)
    fluxes_cp = [[]]*len(fluxes)
    fluxes_bad = [[]]*len(fluxes)
    for i in range(len(fluxes)):
        fluxes_cp[i] = fluxes[i][:]
        fluxes_bad[i] = [0]*len(fluxes[i])
        
    
    #Figure out how many data points are in the window defined by
    #N days
    time_res = dates[1] - dates[0]
    print("Time resolution of the data set is: "
            + str(time_res.total_seconds()) + " seconds.")
    time_res_sec = time_res.total_seconds()
    
    window = datetime.timedelta(days=N)
    window_sec = window.total_seconds()
    
    if time_res_sec > window_sec:
        sys.exit("backward_window_backgroud: The specified time window "
            "for smoothing is shorter than the time resolution between "
            "data points! Exiting. Window: " + str(window) +
            ", Data sets resolution: " + str(time_res))
    
    nwin_pts = int(window_sec/time_res_sec)
    print("backward_window_background: There are " + str(nwin_pts)
        + " data points in the " + str(N) + " days time window.")
    
    #If entered a 1D array of fluxes
    if len(dates) == len(fluxes): nchan = 1
    
    mean_background = [[]]*nchan
    ave_sigma = [[]]*nchan
    threshold = [[]]*nchan
    diff_fluxes = [[]]*nchan
    
    year = dates[0].year
    month = dates[0].month
    day = dates[0].day
    
    #Average the flux in a sliding window of N days
    firstdate = datetime.datetime(year=year,month=month,day=day)
    td = datetime.timedelta(days=N)
    nextdate = firstdate + td
    
    
    #######BEGINNING OF DATA SET#####
    stidx = 0
    endidx = 0
    #Apply a window to the beginning of the data set
    for i in range(len(dates)):
        if dates[i] < nextdate:
            endidx = i
        else:
            break
    
    #calculate values using first N days of data
    #iterate one time to have the opportunity to remove particularly
    #high values
    for j in range(nchan):
        for jj in range(2): #iterate to remove 3sigma+ values
            #Calculate and save needed values
            flux_arr = np.array(fluxes_cp[j][stidx:endidx+1])
            if jj == 1:
                thresh = mean + nsigma*sigma_high
                if thresh != 0 and thresh != None:
                    flux_arr, bad_index = remove_above_value_one(thresh,flux_arr)
                    for bad in bad_index:
                        fluxes_cp[j][stidx+bad] = None #remove
                        fluxes_bad[j][stidx+bad] = fluxes[j][stidx+bad]
            
            flux_arr, bad_index = remove_none_one(flux_arr)
            flux_arr, bad_index = remove_nan_one(flux_arr)
            flux_arr, bad_index = remove_zero_one(flux_arr)
            
            
            if len(flux_arr) < 2:
                mean = 0
                sigma = 0
                sigma_low = 0
                sigma_high = 0
                
            else:
                mean = np.mean(flux_arr)
                sigma = np.std(flux_arr)
                sigma_low = sigma
                sigma_high = sigma
            
            #Fill in first N days of data with the same threshold value
            if jj == 1:
                for i in range(endidx+1):
                    if not mean_background[j]:
                        mean_background[j] = [mean]
                        ave_sigma[j] = [[sigma_low],[sigma_high]]
                        threshold[j] = [mean_background[j][-1] + nsigma*ave_sigma[j][1][-1]]
                        diff_fluxes[j] = [fluxes[j][i] - mean]
                    else:
                        mean_background[j].append(mean)
                        ave_sigma[j][0].append(sigma_low)
                        ave_sigma[j][1].append(sigma_high)
                        threshold[j].append(mean_background[j][-1] + nsigma*ave_sigma[j][1][-1])
                        diff_fluxes[j].append(fluxes[j][i] - mean)

    #######REST OF DATA SET#######
    #Remaining dates
    #Get background and sigma for each day
    reftime = nextdate #from above
    for i in range(endidx+1,len(dates)):
        
        #Only calculate the background for each day
        #If there are multiple data points within 24
        #hours, assign them all the same mean and sigma
        checkdate = datetime.datetime(year=dates[i].year,
                        month=dates[i].month, day=dates[i].day)
        if checkdate == reftime:
            for j in range(nchan):
                mean_background[j].append(mean_background[j][-1])
                ave_sigma[j][0].append(ave_sigma[j][0][-1])
                ave_sigma[j][1].append(ave_sigma[j][1][-1])
                threshold[j].append(mean_background[j][-1] + nsigma*ave_sigma[j][1][-1])
                diff_fluxes[j].append(fluxes[j][i] - mean_background[j][-1])
            continue #go to next time step
        
        #If it's a new day, calculate a new sigma and bg
        if checkdate != reftime:
            reftime = datetime.datetime(year=dates[i].year,
                        month=dates[i].month, day=dates[i].day)
        firstdate = reftime - td
        for k in range(i):
            if dates[k] <= firstdate:
                stidx = k
            else:
                break
        
        for j in range(nchan):
            #Calculate and save needed values
            flux_arr = np.array(fluxes_cp[j][stidx:i+1])
            
            test_arr, bad_index_test = remove_none_one(fluxes[j][stidx:i+1])
            test_arr, bad_index_test = remove_nan_one(test_arr)
            test_arr, bad_index_test = remove_zero_one(test_arr)
            
            if len(test_arr) == []:
                mean = 0
                sigma_low = 0
                sigma_high = 0
            else:
                #If have points in the time frame, remove any that
                #were above the previous threshold value
                if threshold[j][i-1] != 0:
                    flux_arr, bad_index = remove_above_value_one(threshold[j][i-1],flux_arr)
                    #Remove the values that are too high from the
                    #data set
                    for bad in bad_index:
                        fluxes_cp[j][stidx+bad] = None #remove
                        fluxes_bad[j][stidx+bad] = fluxes[j][stidx+bad]
                
                if flux_arr != []:
                    flux_arr, bad_index = remove_none_one(flux_arr)
                    flux_arr, bad_index = remove_nan_one(flux_arr)
                    flux_arr, bad_index = remove_zero_one(flux_arr)

                #Require that at least some number of points in the time
                #window are good for calculating the background
                #otherwise use previous background level
                if flux_arr == []:
                    idx = find_last_good(i,mean_background[j])
                    if idx != None:
                        mean = mean_background[j][idx]
                        sigma_low = ave_sigma[j][0][idx]
                        sigma_high = ave_sigma[j][1][idx]
                    else:
                        mean = 0
                        sigma_low = 0
                        sigma_high = 0
                elif len(flux_arr) < cfg.percent_points*nwin_pts: #Use previous good value
                                #< 0.15*nwin_pts: Pioneer, Ulysses, STEREOA,B
                                #       IMP-8/CRNC,GME
                                #< 4: works for Voyager 1,2
                    idx = find_last_good(i, mean_background[j])
                    if idx != None:
                        mean = np.mean(flux_arr)
                        sigma_low = ave_sigma[j][0][idx]
                        sigma_high = ave_sigma[j][1][idx]
                    else:
                        mean = 0
                        sigma_low = 0
                        sigma_high = 0
                else:
                    mean = np.mean(flux_arr)
                    sigma = np.std(flux_arr)
                    sigma_low = sigma
                    sigma_high = sigma
                
            mean_background[j].append(mean)
            ave_sigma[j][0].append(sigma_low)
            ave_sigma[j][1].append(sigma_high)
            threshold[j].append(mean_background[j][-1] + nsigma*ave_sigma[j][1][-1])
            diff_fluxes[j].append(fluxes[j][i] - mean)
    
    #TESTING
#    for i in range(len(fluxes_bad)):
#        figname = "Bad fluxes " + str(i)
#        fig = plt.figure(figname,figsize=(12,4))
#        maskfluxes = np.ma.masked_less_equal(fluxes_bad[i], 0)
#        plt.plot(dates,maskfluxes,'.-')
#        plt.yscale("log")
        
    return mean_background, ave_sigma, threshold, diff_fluxes







#def ndays_average(N, dates, fluxes, nsigma, remove_above):
#    """ Average flux over N days.
#    
#        INPUTS:
#        
#        :N: (integer) number of days over which to average
#        :dates: (1xn datetime array)
#        :fluxes: (pxn float array) fluxes for p energy channels
#        :remove_above: (float) - remove all flux values above this value
#            (assumed that these are enhanced or undesirable fluxes)
#        
#        OUTPUTS:
#        
#        :ave_dates: (1xm datetime array) dates at the midpoint of each
#            time period
#        :ave_fluxes: (pxm float array) average flux for the N days
#        :ave_sigma: (pxm float array) sigma of flux for the N days
#        
#    """
#
#    #Average the flux between N days starting at midnight
#    firstdate = datetime.datetime(year=dates[0].year, month=dates[0].month, day=dates[0].day)
#    td = datetime.timedelta(days=N)
#    Ntd = int((dates[-1] - firstdate)/td) + 1
#    nextdate = firstdate + td
#    
#    
#    #Put dates and fluxes into a dataframe
#    dict = {'dates': dates}
#    cols = []
#    for ii in range(len(fluxes)):
#        dict.update({'fluxes'+str(ii): fluxes[ii]})
#        cols.append('fluxes'+str(ii))
#    df = pd.DataFrame(dict)
#    
#    means = []
#    sigmas = []
#    df_thresholds = pd.DataFrame()
#    for i in range(Ntd):
#        starttime = firstdate + i*td
#        endtime = firstdate + (i+1)*td
#        
#        sub = df.loc[(df['dates'] >= starttime) & (df['dates'] < endtime)]
#        #Number of dates in this range - need for threshold later
#        ndates = len(sub)
#        selected_dates = sub['dates'].to_list()
#        
#        #Replace all zero values
#        sub = sub.replace(0,np.nan)
#        #Replace all values above remove_above
#        for col in cols:
#            sub.loc[(sub[col] > remove_above)] = np.nan
#
#        
#        if sub.empty:
#            if not means: #No good data encountered yet
#                zeroes = [starttime + td/2.]
#                zeroes.extend([0.]*len(cols))
#                colnames = ['dates']
#                colnames.extend(cols)
#                mean = pd.Series(data=zeroes,index=colnames)
#                sigma = pd.Series(data=zeroes,index=colnames)
#                threshold = pd.Series(data=zeroes,index=colnames)
#            else: #Use last good data
#                mean = means[-1]
#                sigma = sigmas[-1]
#                threshold = thresholds[-1]
#        else:
#            #Take the mean and standard deviation
#            mean = sub.mean()
#            sigma = sub.std()
#            threshold = mean[cols] + sigma[cols]*nsigma
#
#
#        #One mean and sigma per averaged time period
#        means.append(mean)
#        sigmas.append(sigma)
#        
#        #Threshold for every date in the original data set
#        #Make small dataframes and concatenate for each averaged timeframe
#        df_thresh = pd.DataFrame([threshold]*ndates)
#        df_thresh.insert(0,'dates',selected_dates)
#        df_thresholds = pd.concat([df_thresholds,df_thresh],ignore_index=True)
#
#
#    df_means = pd.DataFrame(means)
#    df_sigmas = pd.DataFrame(sigmas)
#
#    ave_dates = df_means['dates'].to_list()
#    ave_fluxes = df_means[cols].T.to_numpy()
#    ave_sigma = df_sigmas[cols].T.to_numpy()
#
#    threshold_dates = df_thresholds['dates'].to_list()
#    threshold = df_thresholds[cols].T.to_numpy()
#    
#    return ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold
#
#
#
#
#def backward_window_background(N, dates, fluxes, nsigma):
#    """ Average over a backward sliding window of N days.
#        Estimate the value of the mean background (GCR) flux,
#        sigma, and a threshold to separate GCR from SEP for
#        every single time step of the data set.
#        
#        Expect that fluxes contain mostly background values
#        and have already had some portion of SEPs excluded.
#        This is meant to be a second iteration that further
#        cleans the first attempt to remove SEPs.
#    
#        INPUTS:
#        
#        :N: (integer) number of days to smooth over
#        :dates: (1xn datetime array)
#        :fluxes: (pxn float array) fluxes for p energy channels
#        :nsigma: (float) number of sigma to calculate threshold
#        
#        OUTPUTS:
#        
#        :smooth_dates: (1xm datetime array) currently unchanged
#        :mean_background: (pxm float array) average flux for the month
#        :ave_sigma: (px2xm float array) sigma of flux for the month
#        
#    """
#    print("backward_window_background: Calculating the background using "
#        "a backward smoothing window of " + str(N) + " days.")
##    #If entered array of fluxes for multiple energy channels
##    nchan = len(fluxes)
##    fluxes_cp = [[]]*len(fluxes)
##    fluxes_bad = [[]]*len(fluxes)
##    for i in range(len(fluxes)):
##        fluxes_cp[i] = fluxes[i][:]
##        fluxes_bad[i] = [0]*len(fluxes[i])
#        
#    
#    #Figure out how many data points are in the window defined by
#    #N days
#    time_res = dates[1] - dates[0]
#    print("Time resolution of the data set is: "
#            + str(time_res.total_seconds()) + " seconds.")
#    time_res_sec = time_res.total_seconds()
#    
#    window = datetime.timedelta(days=N)
#    window_sec = window.total_seconds()
#    
#    if time_res_sec > window_sec:
#        sys.exit("backward_window_backgroud: The specified time window "
#            "for smoothing is shorter than the time resolution between "
#            "data points! Exiting. Window: " + str(window) +
#            ", Data sets resolution: " + str(time_res))
#    
#    nwin_pts = int(window_sec/time_res_sec)
#    print("backward_window_background: There are " + str(nwin_pts)
#        + " data points in the " + str(N) + " days time window.")
#    
#
#    #Average the flux between N days
#    td = datetime.timedelta(days=N)
#
#    
#    
#    #Put dates and fluxes into a dataframe
#    dict = {'dates': dates}
#    cols = []
#    for ii in range(len(fluxes)):
#        dict.update({'fluxes'+str(ii): fluxes[ii]})
#        cols.append('fluxes'+str(ii))
#    df = pd.DataFrame(dict)
#    
#    #Set up the starting dates and time steps
#    #For data resolution less than 1 day
#    firstdate = dates[0]
#    lastdate = dates[-1]
#    td_step = datetime.timedelta(days=1)
#    Nsteps = int((lastdate-firstdate)/td_step) + 1
#    
#    #For data resolution greater than 1 day
#    #One mean per time step
#    if time_res_sec > 86400:
#        td_step = datetime.timedelta(seconds=time_res_sec)
#        firstdate = dates[0]
#        lastdate = dates[-1]
#        Nsteps = len(dates)
#
#    #Number of steps inside of a averaging window
#    Nstart = int(td/td_step)
#
#    df_means = pd.DataFrame()
#    df_sigmas = pd.DataFrame()
#    df_thresholds = pd.DataFrame()
#    df_diff_fluxes = pd.DataFrame()
#    for i in range(Nstart,Nsteps+1,1):
#        #Start N days into the calculation so can use the
#        #Specify a backwards window from Ndays earlier up to current date
#        endtime = firstdate + i*td_step
#        starttime = endtime-td
#        
#        sub = df.loc[(df['dates'] >= starttime) & (df['dates'] < endtime)]
#        selected_dates = sub['dates'].to_list()
#        #All dates in the current time step; all get assigned the same mean
#        #and threshold values
#        current_dates = df['dates'].loc[(df['dates'] >= endtime-td_step) & (df['dates'] < endtime)].to_list()
#        
#        #Replace all zero values
#        sub = sub.replace(0,np.nan)
#        
#        insert_dates = []
#        if i == Nstart: #All dates from the start
#            insert_dates = selected_dates
#        else:
#            insert_dates = current_dates
#        
#        #Check that there are enough points in the selected data to calculate
#        #reliable background and sigma values
#        use_sub = True
#        if len(sub) < cfg.percent_points*nwin_pts:
#            use_sub=False
#            #< 0.15*nwin_pts: Pioneer, Ulysses, STEREOA,B
#            #       IMP-8/CRNC,GME
#            #< 4: works for Voyager 1,2
#            
#            
#        if sub.empty or not use_sub:
#            if df_means.empty: #No good data encountered yet
#                zeroes = [0.]*len(cols)
#                mean = pd.Series(data=zeroes,index=cols)
#                sigma = pd.Series(data=zeroes,index=cols)
#                threshold = pd.Series(data=zeroes,index=cols)
#                
#            else: #Use last good data
#                mean = df_means[cols].iloc[-1]
#                sigma = df_sigmas[cols].iloc[-1]
#                threshold = df_thresholds[cols].iloc[-1]
#        else:
#            #Take the mean and standard deviation
#            mean = sub.mean(numeric_only=True)
#            sigma = sub.std(numeric_only=True)
#            threshold = (mean[cols] + sigma[cols]*nsigma).to_list()
#
#
#        df_mean = pd.DataFrame([mean]*len(insert_dates))
#        df_mean.insert(0,'dates',insert_dates)
#        df_means = pd.concat([df_means,df_mean],ignore_index=True)
#
#        df_sigma = pd.DataFrame([sigma]*len(insert_dates))
#        df_sigma.insert(0,'dates',insert_dates)
#        df_sigmas = pd.concat([df_sigmas,df_sigma],ignore_index=True)
#
#        df_thresh = pd.DataFrame([threshold]*len(insert_dates))
#        df_thresh.insert(0,'dates',insert_dates)
#        df_thresholds = pd.concat([df_thresholds,df_thresh],ignore_index=True)
#
#
#    df_diffs = df[cols] - df_means[cols]
#
#    ave_dates = df_means['dates'].to_list()
#    mean_background = df_means[cols].T.to_numpy()
#    ave_sigma = df_sigmas[cols].T.to_numpy()
#    threshold = df_thresholds[cols].T.to_numpy()
#    diff_fluxes = df_diffs[cols].T.to_numpy()
#
#    return mean_background, ave_sigma, threshold, diff_fluxes
#
#
#
##
##
##    #If entered a 1D array of fluxes
##    if len(dates) == len(fluxes): nchan = 1
##    
##    mean_background = [[]]*nchan
##    ave_sigma = [[]]*nchan
##    threshold = [[]]*nchan
##    diff_fluxes = [[]]*nchan
##    
##    year = dates[0].year
##    month = dates[0].month
##    day = dates[0].day
##    
##    #Average the flux in a sliding window of N days
##    firstdate = datetime.datetime(year=year,month=month,day=day)
##    td = datetime.timedelta(days=N)
##    nextdate = firstdate + td
##    
##    
##    #######BEGINNING OF DATA SET#####
##    stidx = 0
##    endidx = 0
##    #Apply a window to the beginning of the data set
##    for i in range(len(dates)):
##        if dates[i] < nextdate:
##            endidx = i
##        else:
##            break
##    
##    #calculate values using first N days of data
##    #iterate one time to have the opportunity to remove particularly
##    #high values
##    for j in range(nchan):
##        for jj in range(2): #iterate to remove 3sigma+ values
##            #Calculate and save needed values
##            flux_arr = np.array(fluxes_cp[j][stidx:endidx+1])
##            if jj == 1:
##                thresh = mean + nsigma*sigma_high
##                if thresh != 0 and thresh != None:
##                    flux_arr, bad_index = remove_above_value_one(thresh,flux_arr)
##                    for bad in bad_index:
##                        fluxes_cp[j][stidx+bad] = None #remove
##                        fluxes_bad[j][stidx+bad] = fluxes[j][stidx+bad]
##            
##            flux_arr, bad_index = remove_none_one(flux_arr)
##            flux_arr, bad_index = remove_nan_one(flux_arr)
##            flux_arr, bad_index = remove_zero_one(flux_arr)
##            
##            
##            if len(flux_arr) < 2:
##                mean = 0
##                sigma = 0
##                sigma_low = 0
##                sigma_high = 0
##                
##            else:
##                mean = np.mean(flux_arr)
##                sigma = np.std(flux_arr)
##                sigma_low = sigma
##                sigma_high = sigma
##            
##            #Fill in first N days of data with the same threshold value
##            if jj == 1:
##                for i in range(endidx+1):
##                    if not mean_background[j]:
##                        mean_background[j] = [mean]
##                        ave_sigma[j] = [[sigma_low],[sigma_high]]
##                        threshold[j] = [mean_background[j][-1] + nsigma*ave_sigma[j][1][-1]]
##                        diff_fluxes[j] = [fluxes[j][i] - mean]
##                    else:
##                        mean_background[j].append(mean)
##                        ave_sigma[j][0].append(sigma_low)
##                        ave_sigma[j][1].append(sigma_high)
##                        threshold[j].append(mean_background[j][-1] + nsigma*ave_sigma[j][1][-1])
##                        diff_fluxes[j].append(fluxes[j][i] - mean)
##
##    #######REST OF DATA SET#######
##    #Remaining dates
##    #Get background and sigma for each day
##    reftime = nextdate #from above
##    for i in range(endidx+1,len(dates)):
##        
##        #Only calculate the background for each day
##        #If there are multiple data points within 24
##        #hours, assign them all the same mean and sigma
##        checkdate = datetime.datetime(year=dates[i].year,
##                        month=dates[i].month, day=dates[i].day)
##        if checkdate == reftime:
##            for j in range(nchan):
##                mean_background[j].append(mean_background[j][-1])
##                ave_sigma[j][0].append(ave_sigma[j][0][-1])
##                ave_sigma[j][1].append(ave_sigma[j][1][-1])
##                threshold[j].append(mean_background[j][-1] + nsigma*ave_sigma[j][1][-1])
##                diff_fluxes[j].append(fluxes[j][i] - mean_background[j][-1])
##            continue #go to next time step
##        
##        #If it's a new day, calculate a new sigma and bg
##        if checkdate != reftime:
##            reftime = datetime.datetime(year=dates[i].year,
##                        month=dates[i].month, day=dates[i].day)
##        firstdate = reftime - td
##        for k in range(i):
##            if dates[k] <= firstdate:
##                stidx = k
##            else:
##                break
##        
##        for j in range(nchan):
##            #Calculate and save needed values
##            flux_arr = np.array(fluxes_cp[j][stidx:i+1])
##            
##            test_arr, bad_index_test = remove_none_one(fluxes[j][stidx:i+1])
##            test_arr, bad_index_test = remove_nan_one(test_arr)
##            test_arr, bad_index_test = remove_zero_one(test_arr)
##            
##            if len(test_arr) == []:
##                mean = 0
##                sigma_low = 0
##                sigma_high = 0
##            else:
##                #If have points in the time frame, remove any that
##                #were above the previous threshold value
##                if threshold[j][i-1] != 0:
##                    flux_arr, bad_index = remove_above_value_one(threshold[j][i-1],flux_arr)
##                    #Remove the values that are too high from the
##                    #data set
##                    for bad in bad_index:
##                        fluxes_cp[j][stidx+bad] = None #remove
##                        fluxes_bad[j][stidx+bad] = fluxes[j][stidx+bad]
##                
##                if flux_arr != []:
##                    flux_arr, bad_index = remove_none_one(flux_arr)
##                    flux_arr, bad_index = remove_nan_one(flux_arr)
##                    flux_arr, bad_index = remove_zero_one(flux_arr)
##
##                #Require that at least some number of points in the time
##                #window are good for calculating the background
##                #otherwise use previous background level
##                if flux_arr == []:
##                    idx = find_last_good(i,mean_background[j])
##                    if idx != None:
##                        mean = mean_background[j][idx]
##                        sigma_low = ave_sigma[j][0][idx]
##                        sigma_high = ave_sigma[j][1][idx]
##                    else:
##                        mean = 0
##                        sigma_low = 0
##                        sigma_high = 0
##                elif len(flux_arr) < cfg.percent_points*nwin_pts: #Use previous good value
##                                #< 0.15*nwin_pts: Pioneer, Ulysses, STEREOA,B
##                                #       IMP-8/CRNC,GME
##                                #< 4: works for Voyager 1,2
##                    idx = find_last_good(i, mean_background[j])
##                    if idx != None:
##                        mean = np.mean(flux_arr)
##                        sigma_low = ave_sigma[j][0][idx]
##                        sigma_high = ave_sigma[j][1][idx]
##                    else:
##                        mean = 0
##                        sigma_low = 0
##                        sigma_high = 0
##                else:
##                    mean = np.mean(flux_arr)
##                    sigma = np.std(flux_arr)
##                    sigma_low = sigma
##                    sigma_high = sigma
##                
##            mean_background[j].append(mean)
##            ave_sigma[j][0].append(sigma_low)
##            ave_sigma[j][1].append(sigma_high)
##            threshold[j].append(mean_background[j][-1] + nsigma*ave_sigma[j][1][-1])
##            diff_fluxes[j].append(fluxes[j][i] - mean)
##    
##    #TESTING
###    for i in range(len(fluxes_bad)):
###        figname = "Bad fluxes " + str(i)
###        fig = plt.figure(figname,figsize=(12,4))
###        maskfluxes = np.ma.masked_less_equal(fluxes_bad[i], 0)
###        plt.plot(dates,maskfluxes,'.-')
###        plt.yscale("log")
##        
##    return mean_background, ave_sigma, threshold, diff_fluxes
##
##
##
