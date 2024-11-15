from . import plotting_tools as plt_tools
from . import config as cfg
import pandas as pd
import datetime
import os
import math
import numpy as np
from statistics import mode

def sort_bin_order(all_fluxes, energy_bins):
    """Check the order of the energy bins. Usually, bins go from
        low to high energies, but some modelers or users may
        go in reverse order. Usually expect:
        [[10,20],[20,30],[30,40]]
        But user may instead input files with fluxes in order of:
        [[30,40],[20,30],[10,20]]
        
        This subroutine will reorder the fluxes and energy_bins
        to go in increasing order. If differential fluxes were input,
        this reordering will ensure that integral fluxes are
        estimated properly.
        
        INPUTS:
        
        :all_fluxes: (float nxm array) - fluxes for n energy channels
            and m time points
        :energy_bins: (float 2xn array) - energy bins for each of the
            energy channels
            
        OUTPUTS:
        
        :sort_fluxes: (float nxm array) - same as above, but sorted so
            that the lowest energy channel is first and highest is last
        :sort_bins: (float 2xn array) - same as above, but sorted
        
    """

    nbins = len(energy_bins)
    #Rank energy bins in order of lowest to highest effective
    #energies
    eff_en = []
    for i in range(nbins):
        if energy_bins[i][1] == -1:
            eff_en.append(energy_bins[i][0])
        else:
            midpt = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
            eff_en.append(midpt)
            
    eff_en_np = np.array(eff_en)
    sort_index = np.argsort(eff_en_np) #indices in sorted order
    
    sort_fluxes = np.array(all_fluxes)
    sort_bins = []
    for i in range(nbins):
        sort_fluxes[i] = all_fluxes[sort_index[i]]
        sort_bins.append(energy_bins[sort_index[i]])
    
    return sort_fluxes, sort_bins


def determine_time_resolution(dates):
    """ The time resolution is found by taking the difference between
        every consecutive data point. The most common difference is
        taken as the time resolution. Even if the data set has gaps,
        if there are enough consecutive time points in the observational
        or model output, the correct time resolution should be identified.
        
        INPUTS:
        
        :dates: (datetime 1xm array) - dates associated with flux time profile
        
        OUTPUTS:
        
        :time_resolution: (time delta object)
        
    """
    ndates = len(dates)
    time_diff = [a - b for a,b in zip(dates[1:ndates],dates[0:ndates-1])]
    if not time_diff:
        sys.exit("determine_time_resolution: Require more than 1 data point "
                "to determine time resolution. Please extend your "
                "requested time range and try again. Exiting.")
    time_resolution = mode(time_diff)
    return time_resolution


def write_fluxes(experiment, flux_type, options, energy_bins, dates, fluxes, module):
    """ Write dates, fluxes to a standard format csv file with datetime in the 
        first column and fluxes in the remaining column. The energy bins will 
        be indicated in the file comments.
        
        INPUTS:

        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        :options: (string array) options that may be applied to GOES data
        :dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :fluxes: (float nxm array) fluxes for n energy channels and m
            time points
        :module: (string) "idsep" or "opsep"
        
        OUTPUTS:
            
            csv file written to outpath
                 
    """
    modifier, title_modifier = plt_tools.setup_modifiers(options,False)
    stdate = dates[0].strftime("%Y%m%d")
    enddate = dates[-1].strftime("%Y%m%d")
    fname = (f"fluxes_{experiment}_{flux_type}{modifier}_{stdate}_{enddate}.csv")
    fname = os.path.join(cfg.outpath,module,fname)
        
    keys = []
    for bin in energy_bins:
        keys.append((f"{bin[0]}-{bin[1]}"))
        
    dict = {"dates":dates}
    for i in range(len(fluxes)):
        dict.update({keys[i]:fluxes[i]})
        
    df = pd.DataFrame(dict)
    df.to_csv(fname, index=False)
    print("Wrote " + fname + " to file.")
