from . import config as cfg
from . import date_handler as dh
from . import names
import pandas as pd
import datetime
import os
import sys
import math
import numpy as np

def sort_bin_order(all_fluxes, energy_bins, energy_bin_centers=[]):
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
        :energy_bin_centers: (float 1xn array) - corresponding energy bin centers
            
        OUTPUTS:
        
        :sort_fluxes: (float nxm array) - same as above, but sorted so
            that the lowest energy channel is first and highest is last
        :sort_bins: (float 2xn array) - same as above, but sorted
        
    """

    nbins = len(energy_bins)
    #Rank energy bins in order of lowest to highest effective
    #energies
    if len(energy_bin_centers) == 0:
        for i in range(nbins):
            if energy_bins[i][1] == -1:
                energy_bin_centers.append(energy_bins[i][0])
            else:
                midpt = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
                energy_bin_centers.append(midpt)
            
    eff_energies = np.array(energy_bin_centers)
    sort_index = np.argsort(eff_energies) #indices in sorted order

    sort_fluxes = np.array(all_fluxes)
    sort_bins = []
    sort_centers = []
    for i in range(nbins):
        sort_fluxes[i] = all_fluxes[sort_index[i]]
        sort_bins.append(energy_bins[sort_index[i]])
        sort_centers.append(energy_bin_centers[sort_index[i]])
    
    return sort_fluxes, sort_bins, sort_centers



def write_fluxes(experiment, flux_type, exp_name, energy_bins, dates, fluxes,
    module='', subdir='', modifier='', savepath='', suffix=''):
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
    name = ''
    stdate = dates[0].strftime("%Y%m%d")
    enddate = dates[-1].strftime("%Y%m%d")
    if module == 'idsep' or module == 'download':
        name = idsep_naming_scheme(experiment, flux_type, exp_name, options,
                    spacecraft=spacecraft)
        fname = (f"fluxes_{name}_{stdate}_{enddate}.csv")
        fname = os.path.join(savepath, fname)
    elif module == 'opsep':
        #name like the json schema
        tzulu = dh.time_to_zulu(dates[0])
        tzulu = tzulu.replace(":","")
        fname = f"{subdir}.{tzulu}.csv"
        if suffix != '':
            fname = f"{subdir}.{tzulu}_{suffix}.csv"
        fname = os.path.join(savepath, fname)
        
    keys = []
    for bin in energy_bins:
        keys.append(names.energy_bin_key(bin))
    
    dict = {"dates":dates}
    for i in range(len(fluxes)):
        dict.update({keys[i]:fluxes[i]})
        
    df = pd.DataFrame(dict)
    df.to_csv(fname, index=False)
    print("Wrote " + fname + " to file.")

    return fname


def make_lists_array(lists):
    """ If lists is a string, make an array of the filenames inside. """

    if not os.path.isfile(lists):
        sys.exit(f"make_lists_array: File does not exist. {lists}")
    
    arr = []
    with open(lists, 'r') as file:
        for list in file:
            list = list.strip()
            arr.append(list)
            
    return arr


