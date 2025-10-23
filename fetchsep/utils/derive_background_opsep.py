from . import read_datasets as datasets
from ..utils import config as cfg
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
#import urllib2
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
import scipy.integrate
from numpy import exp
import array as arr
import pandas as pd
import scipy

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

#2021-09-25, changes in 0.2: print out means and sigmas in derive_background


def about_derive_background_opsep():
    """ About derive_background_opsep.py
        
        Perform background-subtraction and return only SEP flux values.

        1. Read in a user-specified time period for the background flux
        2. Estimate the mean flux and level of variation (sigma)
        3. Define all flux above mean+nsigma*sigma (nsigma is defined in global_vars)as SEP flux and all flux below as background flux
        4. Separate the SEP fluxes and subtract the mean background value
        5. Return the total flux, the background flux, and the background-subtracted SEP fluxes
    """
    
def remove_none(flux):
    """ Takes 1D array of flux and removes None values.
        None values in a list are converted to NaN values in a numpy array.
        So check for NaN and remove.
        
        INPUTS:
        
        :flux: (float 1xn array)
        
        OUTPUTS:
        
        :clean_flux: (float 1xm array) with None or NaN values removed
        
    """
    bad_index = np.argwhere(np.isnan(flux))
    clean_flux = np.delete(flux, bad_index)
    return clean_flux, bad_index


def remove_zero(flux):
    """Takes 1D array of flux and removes zero values.
    
        INPUTS:
        
        :flux: (float 1xn array)
        
        OUTPUTS:
        
        :clean_flux: (float 1xm array) with zero values removed
        
    """
    bad_index = np.argwhere(flux == 0)
    clean_flux = np.delete(flux, bad_index)
    return clean_flux, bad_index



def remove_above(flux, val):
    """ Remove any flux above a specific value, val.
        
        INPUTS:
        
        :flux: (float 1xn array)
        :val: (float) lower threshold value
        
        OUTPUTS:
        
        :strip_flux: (float 1xm array) with flux > val removed
    
    """
    indices = np.nonzero(flux > val)
    strip_flux = np.delete(flux, indices)
    return strip_flux


def remove_below(flux, val):
    """ Remove any flux below a specific value, val.
    
        INPUTS:
        
        :flux: (float 1xn array)
        :val: (float) upper threshold value
        
        OUTPUTS:
        
        :strip_flux: (float 1xm array) with flux < val removed
        
    """
    indices = np.nonzero(flux < val)
    strip_flux = np.delete(flux, indices)
    return strip_flux


def separate_sep_and_background(fluxes, means, sigmas,
    nsigma=cfg.opsep_nsigma, doBGSub=True):
    """ Take the input fluxes, separate them into arrays containing
        the background flux and SEP flux. Values above mean + Nsigma*sigma is
        considered SEP flux while values below are considered the background.
        Perform a background subtraction on the SEP flux by subtracting the
        mean background value.
        The input flux array is a numpy array, but the output will be a list
        for flexibility.
        Nsigma is specified in config/config_opsep.py.
        
        INPUTS:
        
        :fluxes: (float nxm array) fluxes for n energy channels and m time points
        :dates: (datetime 1xm array) time points for flux time profile
        :means: (float 1xn array) mean background flux for n energy channels
        :sigmas: (float 1xn array) expected variability sigma for n energy channels
        
        OUTPUTS:
        
        :bgfluxes: (float nxm array) background fluxes for n energy channels and
            m time points
        :sepfluxes: (float nxm array) SEP fluxes for n energy channels and
            m time points
        
    """
    nflx = len(fluxes)
    sepfluxes = []
    bgfluxes = []

    for i in range(nflx):
        bgflux = [-1]
        sepflux = [-1]
        for j in range(len(fluxes[i])):
            if fluxes[i][j] <= means[i] + nsigma*sigmas[i]:
                bgflux.append(fluxes[i][j])
                sepflux.append(0)
            if fluxes[i][j] > means[i] + nsigma*sigmas[i]:
                bgflux.append(0)
                bgsubflux = fluxes[i][j] - means[i]
                if bgsubflux < 0: bgsubflux = 0
                if doBGSub:
                    sepflux.append(bgsubflux)
                else:
                    sepflux.append(fluxes[i][j])
            if np.isnan(fluxes[i][j]):
                bgflux.append(np.nan)
                sepflux.append(np.nan)

        bgflux.pop(0)
        sepflux.pop(0)

        if not bgfluxes:
            bgfluxes = [bgflux]
        else:
            bgfluxes.append(bgflux)

        if not sepfluxes:
            sepfluxes = [sepflux]
        else:
            sepfluxes.append(sepflux)

    return bgfluxes, sepfluxes


def separate_sep_and_background_idsep(fluxes, means, sigmas,
    nsigma=cfg.opsep_nsigma, doBGSub=True):
    """ Take the input fluxes, separate them into arrays containing
        the background flux and SEP flux. Values above mean + Nsigma*sigma is
        considered SEP flux while values below are considered the background.
        Perform a background subtraction on the SEP flux by subtracting the
        mean background value.
        The input flux array is a numpy array, but the output will be a list
        for flexibility.
        Nsigma is specified in config/config_opsep.py.
        
        INPUTS:
        
        :fluxes: (float nxm array) fluxes for n energy channels and m time points
        :dates: (datetime 1xm array) time points for flux time profile
        :means: (float nxm array) mean background flux for n energy channels and m time points
        :sigmas: (float nxm array) expected variability sigma for n energy channels and m time points
        
        OUTPUTS:
        
        :bgfluxes: (float nxm array) background fluxes for n energy channels and
            m time points
        :sepfluxes: (float nxm array) SEP fluxes for n energy channels and
            m time points
        
    """
    nflx = len(fluxes)
    sepfluxes = []
    bgfluxes = []

    for i in range(nflx):
        bgflux = [-1]
        sepflux = [-1]
        for j in range(len(fluxes[i])):
            if np.isnan(fluxes[i][j]):
                bgflux.append(np.nan)
                sepflux.append(np.nan)
            elif pd.isnull(means[i][j]): #idsep did not get a good background
                bgflux.append(np.nan)
                sepflux.append(np.nan)
            elif fluxes[i][j] <= means[i][j] + nsigma*sigmas[i][j]:
                bgflux.append(fluxes[i][j])
                sepflux.append(0)
            elif fluxes[i][j] > means[i][j] + nsigma*sigmas[i][j]:
                bgflux.append(0)
                bgsubflux = fluxes[i][j] - means[i][j]
                if bgsubflux < 0: bgsubflux = 0
                if doBGSub:
                    sepflux.append(bgsubflux)
                else:
                    sepflux.append(fluxes[i][j])
        bgflux.pop(0)
        sepflux.pop(0)

        if not bgfluxes:
            bgfluxes = [bgflux]
        else:
            bgfluxes.append(bgflux)

        if not sepfluxes:
            sepfluxes = [sepflux]
        else:
            sepfluxes.append(sepflux)

    return bgfluxes, sepfluxes



def define_hist_bins(flux):
    """Takes a 1D numpy array of flux and defines a set of histogram bins
        between the min and max flux values equally spaced in log space.
        
        INPUTS:
        
        :flux: (float 1xm array) flux time profile for m time steps
        
        OUTPUTS:
        
        :hist_bins: (float 1x20 (nbins) array) bins for a histogram equally
            spaced in log space over 20 betweens between the minimum and
            maximum flux values in flux
        
    """
    if flux.any(0):
        flux, bad_index = remove_zero(flux)
    max_val = flux.max()
    min_val = flux.min()

    nbins = 20;
    bins = []
    #Create nbins bins in log space between min_val and max_val
    logmax = math.log10(max_val)
    logmin = math.log10(min_val)
    log_bins = np.linspace(start=logmin, stop=logmax, num=nbins + 1,\
                        endpoint=True)

    for bin in log_bins:
        if not bins:
            bins = [10**bin]
        else:
            bins.append(10**bin)

    hist_bins = np.array(bins)
    return hist_bins


def create_histogram(flux, energy_bin, iteration):
    """Take a list of flux with time and generate a histogram of the values.
        NaN values are removed prior to creating histogram.
        The histogram is created with bins extending from the min flux to the
        max flux and equally spaced in log space.
        Estimate the mean by averaging the bin centers weighted by the
        frequency. Calculate the variance and take the square root to estimate
        sigma. Return the mean and sigma.
        
        INPUTS:
        
        :flux: (float 1xm array) flux time profile for m time points
        :energy_bin: (float 2x1 array) energy bin that defines the energy
            channel for the flux array (only used for plotting)
        :iteration: (integer) indicates how many times the flux has gone
            through this subroutine (only used for plotting)
            
        OUTPUTS:
        
        :hist_mean: (float) weighted mean of the histogram
        :sigma: (float) standard deviation (sqrt(variance)) of
            the histogram
        
    """
    #Bad values in the data were set to None
    #Remove None values to calculate flux distribution in a histogram
    clean_flux, bad_index = remove_none(flux) #remove any None values in numpy array

    #Define a set of logarithmic bins that span the flux values
    hist_bins = define_hist_bins(clean_flux)
    hist, _ = np.histogram(clean_flux, bins=hist_bins)

    #Check if the lowest histogram bin represents a data floor (as in GOES).
    #This is indicated by the lowest energy bin containing the most values
    if hist[0] == hist.max():
        if hist.max() <= 0.7*np.sum(hist):
            #remove lowest energy bin and recalculate histogram
            hist_bins = np.delete(hist_bins, 0)
            hist, _ = np.histogram(clean_flux, bins=hist_bins)

    centers = []
    for i in range(len(hist_bins)-1):
        center = math.sqrt(hist_bins[i]*hist_bins[i+1])
        if not centers:
            centers = [center]
        else:
            centers.append(center)
    bin_centers = np.array(centers)
    hist_mean = np.average(bin_centers, weights=hist)
    variance = np.average((bin_centers - hist_mean)**2, weights=hist)
    sigma = np.sqrt(variance)

    # Plot histogram
#    figname = 'FluxHistogram_'+ str(energy_bin[0]) + '_' \
#            + str(energy_bin[1]) + '_it' +str(iteration)
#    fig = plt.figure(figname,figsize=(8,5))
#    ax = plt.subplot(111)
#    n, bins, patches = plt.hist(x=clean_flux, bins=hist_bins, color='#0504aa', alpha=0.7, rwidth=0.85)
#    plt.grid(axis='y', alpha=0.75)
#    plt.xlabel('Flux')
#    plt.ylabel('Frequency')
#    plt.title('Distribution of Flux Values for '
#        + str(energy_bin[0]) + ' to ' + str(energy_bin[1]) + ' MeV')
#    plt.text(0.8, 0.5, (r'$\mu$= %.4g, sigma= %.4g'%(hist_mean,
#        sigma)), horizontalalignment='center',
#        verticalalignment='center', transform=ax.transAxes)
#    maxfreq = n.max()
#     # Set a clean upper y-axis limit.
#    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
#    plt.show()

    return hist_mean, sigma


def calc_mean_sigma(flux):
    """ Calculate the mean and the sigma of the flux using simple
        mean and sigma definitions.
    """
    #Bad values in the data were set to None
    #Remove None values to calculate flux distribution in a histogram
    clean_flux, bad_index = remove_none(flux) #remove any None values in numpy array
    clean_flux, bad_index = remove_zero(flux) #remove any None values in numpy array
    mean  = sum(clean_flux)/len(clean_flux)
    sigmasq = sum([(x-mean)**2 for x in clean_flux])/len(clean_flux)
    sigma = math.sqrt(sigmasq)

    return mean, sigma
    

def iterate_background(fluxes, energy_bins):
    """Bin fluxes into histograms to calculate the background mean and sigma.
        Exclude fluxes above and below mean +- 3sigma and recalculate mean
        and sigma. Use thes values as the final estimates of the background flux
        and the expected level of variability in the background.
        
        INPUTS:
        
        :fluxes: (float nxm array) flux time profiles for n energy channels and
            m time points
        :energy_bins: (float nx2 array) energy bins for n energy channels
        
        OUTPUTS:
        
        :means: (float 1xn array) mean values of histogram for n energy channels
        :sigmas: (float 1xn array) sigma of histograms for n energy channels
        
    """
    means = []
    sigmas = []
    for i in range(len(fluxes)):
        strip_flux = fluxes[i]
        for it in range(1): #number of iterations
            #First iteration: Use all fluxes in the background time period to
            #estimate mean and sigma.
            mean, sigma = create_histogram(strip_flux, energy_bins[i],it)
            #mean, sigma = calc_mean_sigma(strip_flux)
            #print("mean: " + str(mean) + ", sigma: "+ str(sigma))
            #exclude values above mean + 3sigma
            highval = mean + 3*sigma
            #print("highval: " + str(highval))
            strip_flux = remove_above(strip_flux,highval)
            #exclude values below mean - 3sigma
            lowval = mean - 3*sigma
            #print("lowval: "  + str(lowval))
            strip_flux = remove_below(strip_flux,lowval)

        if not means:
            means = [mean]
        else:
            means.append(mean)
        if not sigmas:
            sigmas = [sigma]
        else:
            sigmas.append(sigma)

    return means, sigmas



def derive_background(experiment, flux_type, options,
    bgstartdate, bgenddate, dates, fluxes, energy_bins,
    showplot, saveplot, nsigma=cfg.opsep_nsigma, doBGSub=True):
    """ Derive the background using fluxes in the time period between
        background start and end dates specified by the user. Derive the
        mean background value along with an expected level of variation (sigma)
        in the background. Make two separate arrays containing 1) background
        flux, 2) background-subtracted SEP fluxes.
        The fluxes will be separate by selecting fluxes above and below
        mean + Nsigma*sigma. The value of Nsigma is specified in
        config/config_opsep.py.
        Return the background and background-subtracted SEP flux arrays along
        with a date array. The fluxes and dates will extend from BGStartdate to
        SEPEndDate. The fluxes will be numpy arrays and the dates are a list.
        
        If doBGSub is False, then don't subtract the background, just return 
        original fluxes but with the background values set to zero and only
        the enhanced fluxes as nonzero.
        
        INPUTS:
        
        :bgstartdate: (datetime) starting date of background time period
        :bgenddate: (datetime) ending date of background time period 
        :dates: (datetime 1xm array) 
        :fluxes: (numpy float nxm array) - fluxes for n energy channels and m
            time steps
        :energy_bins: (array nx2 for n thresholds)
        :showplot: (bool) True to print plot to screen
        :saveplot: (bool) True to save plot automatically
        
        OUTPUTS:
        
        :bgfluxes: (float nxm array) background fluxes for n energy channels
            and m time points
        :sepfluxes: (float nxm array) background-subtracted SEP fluxes for n
            energy channels and m time points
        :dates: (datetime 1xm array) m time points extending from
            str_bgstartdate to str_enddate containing background flux time period
            and SEP flux time period
        
    """

    sepem_end_date = datetime.datetime(2015,12,31,23,55,00)
    if(experiment == "SEPEM" and (enddate > sepem_end_date)):
        sys.exit('The SEPEM (RSDv2) data set only extends to '
                  + str(sepem_end_date) +
            '. Please change your requested dates. Exiting.')

#    if experiment[0:4] == "GOES" and flux_type == "integral":
#        sys.exit("Do not perform background subtraction on GOES integral "
#                "fluxes. Integral fluxes have already been derived by "
#                "applying corrections for cross-contamination and removing "
#                "the instrument background levels.")
#
#    if experiment[0:4] == "GOES" and "uncorrected" not in options:
#        print("Warning: GOES corrected fluxes have already been derived by "
#                "applying corrections for cross-contamination and removing "
#                "the instrument and GCR background levels up to channel P6. "
#                "Please be sure it makes sense to perform a background "
#                "subtraction of this data (e.g. for HEPAD energies). "
#                "Otherwise, please add --options uncorrected to perform "
#                "background subtracion on GOES uncorrected fluxes. Continuing.")
#
#    if experiment[0:4] == "GOES" and "uncorrected" in options:
#        print("Note: Background-subtraction of uncorrected GOES fluxes "
#                "does not remove the effects of spurious increases in the low "
#                "energy channels due to cross-talk from the high energy "
#                "channels, particularly at the onset of well-connected, "
#                "intense SEP events. It also does not remove any contamination "
#                "due to particles entering the GOES detectors from the sides.")


    #Pull out the fluxes in the time period to be used for calculating
    #background
    bg_dates, bg_fluxes = datasets.extract_date_range(bgstartdate,
            bgenddate, dates, fluxes)
    print(f"Calculating background with data from {bgstartdate} to {bgenddate}.")

    means, sigmas = iterate_background(bg_fluxes, energy_bins)
    bgfluxes, sepfluxes = separate_sep_and_background(fluxes, means, sigmas,
                        doBGSub=doBGSub)
                     
    print("=====BACKROUND IDENITIFCATION=====")
    for k in range(len(means)):
        print(f"Mean: {means[k]} +- {cfg.opsep_nsigma} * {sigmas[k]}")
    
    #To get into a consistent format as used in opsep.py
    bgfluxes = np.array(bgfluxes)
    sepfluxes = np.array(sepfluxes)

    return bgfluxes, sepfluxes, means, sigmas

