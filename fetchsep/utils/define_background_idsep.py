from ..utils import config as cfg
from . import date_handler as dateh
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
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


######BACKGROUND AND ENHANCEMENT IDENTIFICATION##############
#Note that the original versions of those code are ndays_average()
#and backward_window_background().
#The algorithms were rewritten to take advantage of dataframes to
#attempt to speed up the computation. These subroutines are labeled
#_optimized().

def write_df(df, name, verbose=True):
    """Writes a pandas dataframe to the standard location in multiple formats
    """
    dataformats = (('pkl' , getattr(df, 'to_pickle'), {}),
                   ('csv',  getattr(df, 'to_csv'), dict(index=False)))
    for ext, write_func, kwargs in dataformats:
        filepath = os.path.join(cfg.outpath, 'idsep', ext, name + '.' + ext)
        write_func(filepath, **kwargs)
        if verbose:
            print('Wrote ' + filepath)



#####################################
#OPTIMIZED ALGORITHM
#####################################
def ndays_average_optimized(N, dates, fluxes, nsigma, remove_above):
    """ Average flux over N days.

        INPUTS:

        :N: (integer) number of days over which to average
        :dates: (1xn datetime array)
        :fluxes: (pxn float array) fluxes for p energy channels
        :remove_above: (float) - remove all flux values above this value
            (assumed that these are enhanced or undesirable fluxes)

        OUTPUTS:

        :ave_dates: (1xm datetime array) dates at the midpoint of each
            time period
        :ave_fluxes: (pxm float array) average flux for the N days
        :ave_sigma: (pxm float array) sigma of flux for the N days

    """

    #Average the flux between N days starting at midnight
    firstdate = datetime.datetime(year=dates[0].year, month=dates[0].month, day=dates[0].day)
    td = datetime.timedelta(days=N)
    Ntd = int((dates[-1] - firstdate)/td) + 1
    nextdate = firstdate + td


    #Put dates and fluxes into a dataframe
    dict = {'dates': dates}
    cols = []
    for ii in range(len(fluxes)):
        dict.update({'fluxes'+str(ii): fluxes[ii]})
        cols.append('fluxes'+str(ii))
    df = pd.DataFrame(dict)

    means = []
    sigmas = []
    ave_dates = []
    df_thresholds = pd.DataFrame()
    for i in range(Ntd):
        starttime = firstdate + i*td
        endtime = firstdate + (i+1)*td

        sub = df.loc[(df['dates'] >= starttime) & (df['dates'] < endtime)]
        #Number of dates in this range - need for threshold later
        ndates = len(sub)
        selected_dates = sub['dates'].to_list()
        ave_date = sub['dates'].mean

        #Replace all zero values
        sub = sub.replace(0,np.nan)
        #Replace all values above remove_above
        for col in cols:
            sub.loc[(sub[col] > remove_above)] = np.nan


        if sub.empty:
            if not means: #No good data encountered yet
                zeroes = [starttime + td/2.]
                zeroes.extend([0.]*len(cols))
                colnames = ['dates']
                colnames.extend(cols)
                mean = pd.Series(data=zeroes,index=colnames)
                sigma = pd.Series(data=zeroes,index=colnames)
                threshold = pd.Series(data=zeroes,index=colnames)
            else: #Use last good data
                mean = means[-1]
                sigma = sigmas[-1]
                threshold = thresholds[-1]
        else:
            #Take the mean and standard deviation
            mean = sub.mean()
            sigma = sub.std()
            threshold = mean[cols] + sigma[cols]*nsigma


        #One mean and sigma per averaged time period
        means.append(mean)
        sigmas.append(sigma)
        ave_dates.append(ave_date)

        #Threshold for every date in the original data set
        #Make small dataframes and concatenate for each averaged timeframe
        df_thresh = pd.DataFrame([threshold]*ndates)
        df_thresh.insert(0,'dates',selected_dates)
        df_thresholds = pd.concat([df_thresholds,df_thresh],ignore_index=True)


    df_means = pd.DataFrame(means)
    df_sigmas = pd.DataFrame(sigmas)

    #ave_dates = df_means['dates'].to_list()
    ave_fluxes = df_means[cols].T.to_numpy()
    ave_sigma = df_sigmas[cols].T.to_numpy()

    threshold_dates = df_thresholds['dates'].to_list()
    threshold = df_thresholds[cols].T.to_numpy()

    #Write fluxes to file for testing and use
    write_df(df_means,'mean_background_fluxes_ndays_optimized')
    write_df(df_sigmas,'background_sigma_ndays_optimized')
    write_df(df_thresholds,'threshold_ndays_optimized')

    return ave_dates, ave_fluxes, ave_sigma, threshold_dates, threshold




#####################################
##OPTIMIZED AGLORITHM
#####################################
def backward_window_background_optimized(N, dates, fluxes, nsigma,iteration=0):
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
    
    #Figure out how many data points are in the window defined by
    #N days
    time_res = dates[1] - dates[0]
    print("Time resolution of the data set is: "
            + str(time_res.total_seconds()) + " seconds.")
    time_res_sec = time_res.total_seconds()
    
    window_sec = datetime.timedelta(days=N).total_seconds()
    
    if time_res_sec > window_sec:
        sys.exit("backward_window_background_optimized: The specified time window "
            "for smoothing is shorter than the time resolution between "
            "data points! Exiting. Window: " + str(window) +
            ", Data sets resolution: " + str(time_res))
    
    nwin_pts = int(window_sec/time_res_sec)
    print("backward_window_background_optimized: There are " + str(nwin_pts)
        + " data points in the " + str(N) + " days time window.")
    print("backward_window_background_optimized: Require " + str(cfg.percent_points*nwin_pts)
        + " point to calculate background.")

    #Average the flux between N days
    td_win = datetime.timedelta(days=N)

    
    #Put dates and fluxes into a dataframe
    dict = {'dates': dates}
    cols = []
    for ii in range(len(fluxes)):
        dict.update({'fluxes'+str(ii): fluxes[ii]})
        cols.append('fluxes'+str(ii))
    df = pd.DataFrame(dict)
    
    #Set up the starting dates and time steps
    #For data resolution less than 1 day
    firstdate = dates[0]
    lastdate = dates[-1]
    td_step = datetime.timedelta(days=1)
    Nsteps = int((lastdate-firstdate)/td_step) + 1
    
    #For data resolution greater than 1 day
    #One mean per time step
    if time_res_sec > 86400:
        td_step = datetime.timedelta(seconds=time_res_sec)
        Nsteps = len(dates)

    #Number of steps inside of an averaging window
    Nstart = int(td_win/td_step)

    df_means = pd.DataFrame()
    df_sigmas = pd.DataFrame()
    df_thresholds = pd.DataFrame()
    df_diff_fluxes = pd.DataFrame()
    for i in range(Nstart,Nsteps+1,1):
        #Start N days into the calculation so can use the
        #Specify a backwards window from Ndays earlier up to current date
        endtime = firstdate + i*td_step
        starttime = endtime - td_win
        
        #Get all timesteps before removing NaN values so can generate a
        #background for all timesteps in the original data set
        sub = df.loc[(df['dates'] >= starttime) & (df['dates'] < endtime)]
        selected_dates = sub['dates'].to_list()
        #All dates in the current time step (e.g. 1 day) get assigned the same mean
        #and threshold values
        #Getting background value for the last day in the 27 day window*****
        current_dates = df['dates'].loc[(df['dates'] >= endtime-td_step) & (df['dates'] < endtime)].to_list()
        
        insert_dates = []
        if i == Nstart: #All dates from the start
            insert_dates = selected_dates
        else:
            insert_dates = current_dates

        #Replace all zero values with nan
        #nan values are ignored by pd.mean and pd.sigma
        sub = sub.replace(0,np.nan)
        #print(f"Start Time: {starttime}, End Time: {endtime}, All points: {len(sub)}, Required: {cfg.percent_points*nwin_pts}")
        #For each column of fluxes, calculate the mean and sigma.
        #Check that there are enough points in the selected data to calculate
        #reliable background and sigma values
        means = []
        sigmas = []
        thresholds = []
        for col in cols:
            #Set points above the previous threshold to nan
            if not df_thresholds.empty:
                prev_thresh = df_thresholds[col].iloc[-1]
                #print(f"Column: {col}, Start Time: {starttime}, End Time: {endtime}, Previous Threshold: {prev_thresh}")
                if not pd.isnull(prev_thresh) and prev_thresh != 0:
                    sub.loc[(sub[col] > prev_thresh),col] = np.nan

            ngood = len(sub[col].dropna())
            #print(f"Start Time: {starttime}, End Time: {endtime}, Number of good points: {ngood}, Required: {cfg.percent_points*nwin_pts}")
            if ngood < cfg.percent_points*nwin_pts:
                #If no good points yet, then set to zero
                if df_means.empty:
                    mean = 0
                    sigma = 0
                    threshold = 0
                else:
                    mean = df_means[col].iloc[-1]
                    sigma = df_sigmas[col].iloc[-1]
                    threshold = df_thresholds[col].iloc[-1]

            else:
                mean = sub[col].mean()
                sigma = sub[col].std(ddof=0) #1/N
                threshold = mean + sigma*nsigma

            
            means.append(mean)
            sigmas.append(sigma)
            thresholds.append(threshold)

        smean = pd.Series(means,index=cols)
        ssigma = pd.Series(sigmas,index=cols)
        sthreshold = pd.Series(thresholds,index=cols)


        df_mean = pd.DataFrame([smean]*len(insert_dates))
        df_mean.insert(0,'dates',insert_dates)
        df_means = pd.concat([df_means,df_mean],ignore_index=True)

        df_sigma = pd.DataFrame([ssigma]*len(insert_dates))
        df_sigma.insert(0,'dates',insert_dates)
        df_sigmas = pd.concat([df_sigmas,df_sigma],ignore_index=True)

        df_thresh = pd.DataFrame([sthreshold]*len(insert_dates))
        df_thresh.insert(0,'dates',insert_dates)
        df_thresholds = pd.concat([df_thresholds,df_thresh],ignore_index=True)

    ave_dates = df_means['dates'].to_list()
    mean_background = df_means[cols].T.to_numpy()
    ave_sigma = df_sigmas[cols].T.to_numpy()
    threshold = df_thresholds[cols].T.to_numpy()

    #Write fluxes to file for testing and use
    write_df(df_means,'mean_background_fluxes_optimized_it'+str(iteration))
    write_df(df_sigmas,'background_sigma_optimized_it'+str(iteration))
    write_df(df_thresholds,'threshold_optimized_it'+str(iteration))

    return mean_background, ave_sigma, threshold



