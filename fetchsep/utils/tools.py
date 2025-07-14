from . import plotting_tools as plt_tools
from . import config as cfg
import pandas as pd
import datetime
import os
import math
import numpy as np
from statistics import mode
import scipy
from numpy import exp

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


def write_fluxes(experiment, flux_type, options, energy_bins, dates, fluxes, module,
    spacecraft=""):
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
    modifier, title_modifier = plt_tools.setup_modifiers(options,False,spacecraft=spacecraft)
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
    
    
def from_differential_to_integral_flux(experiment, min_energy, energy_bins,
    fluxes, bruno2017=False):
    """ If user selected differential fluxes, convert to integral fluxes to
        caluculate operational threshold crossings (>10 MeV protons exceed 10
        pfu, >100 MeV protons exceed 1 pfu).
        Assume that the measured fluxes correspond to the center of the energy
        bin and use power law interpolation to extrapolate integral fluxes
        above user input min_energy.
        The intent is to calculate >10 MeV and >100 MeV fluxes, but leaving
        flexibility for user to define the minimum energy for the integral flux.
        An integral flux will be provided for each timestamp (e.g. every 5 mins).
        
        Note that if bad values were not interpolated in previous steps,
        they will have been set to None in check_for_bad_data, which translated
        to NaN in numpy arrays.
       
       If GOES-13, -14, -15, then the 110 - 900 MeV bin is not included when
       estimating the integral fluxes (unless Bruno2017 is selected).
       This bin's energy appears to be unreliable and the other energy
       bins cover this range.
       
        INPUTS:
        
        :experiment: (string)
        :min_energy: (float) - bottom energy for integral flux calculation
        :energy_bins: (float nx2 array) - bins for each energy channel
            [[Emin1, Emax1], [Emin2, Emax2], [], ...]
        :fluxes: (float nxm array) - fluxes with time for each energy channel
        :bruno2017: (bool) - apply Bruno 2017 correction to GOES energy bins?
        
        OUTPUTS:
        
        :integral_fluxes: (float 1xm array) - estimate integral flux for >min_energy
            (Returns all zero values if no energy bins above min_energy)
            
    """
    print(f"Converting differential flux to integral flux for >{min_energy} MeV.")
    nbins = len(energy_bins)
    nflux = len(fluxes[0])
    #Check requested min_energy inside data energy range
    if min_energy < energy_bins[0][0] or min_energy >= energy_bins[nbins-1][0]:
        print(f"The selected minimum energy {min_energy} to create "
               + "integral fluxes is outside of the range of the data: "
                + f"{energy_bins[0][0]}  - {max(energy_bins[nbins-1][0],energy_bins[nbins-1][1])}")
        print(f"Setting all >{min_energy} fluxes to zero.")
        integral_fluxes = [0]*nflux
        return integral_fluxes

    #Calculate bin center in log space for each energy bin
    bin_center = []
    for i in range(nbins):
        if energy_bins[i][1] != -1:
            centerE = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
        else:
            centerE = -1
        bin_center.append(centerE)

    #The highest energy EPEAD bin overlaps with all of the HEPAD bins
    #For this reason, excluding the highest energy EPEAD bin in
    #integral flux estimation, 110 - 900 MeV
    if (experiment == "GOES-13" or experiment == "GOES-14" or
        experiment == "GOES-15") and not bruno2017:
        remove_bin = -1
        for i in range(nbins):
            if energy_bins[i][0] == 110 and energy_bins[i][1] == 900:
                remove_bin = i
        if remove_bin == -1:
            sys.exit("Attempting to remove 110 - 900 MeV bin for "
                    + experiment + '. Cannot locate bin. Please check '
                    'define_energy_bins to see if GOES-13 to GOES-15 bins '
                    'include 110 - 900 MeV. if not please comment out this '
                    'section in tools.py from_differential_to_integral_flux.')
        fluxes = np.delete(fluxes,remove_bin,0)
        energy_bins = np.delete(energy_bins,remove_bin,0)
        bin_center = np.delete(bin_center,remove_bin,0)
        nbins = nbins-1

    #integrate by power law interpolation; assume spectrum the shape of a
    #power law from one bin center to the next. This accounts for the case
    #where the minimum energy falls inside of a bin or there are overlapping
    #energy bins or gaps between bins.
    #An energy bin value of -1 (e.g. [700,-1]) indicates infinity - already an
    #integral channel. This happens for HEPAD. If the integral channel is the
    #last channel, then the flux will be added. If it is a middle bin, it will
    #be skipped.
    integral_fluxes = []
    integral_fluxes_check = []
    for j in range(nflux):  #flux at each time
        sum_flux = 0
        ninc = 0 #number of energy bins included in integral flux estimate
        for i in range(nbins-1):
            if bin_center[i+1] < min_energy:
                continue
            else:
                if energy_bins[i][1] == -1 or energy_bins[i+1][1] == -1:
                    #bin is already an integral flux, e.g. last HEPAD bin
                    continue
                if pd.isnull(fluxes[i,j]) or pd.isnull(fluxes[i+1,j]): #data gap
                    continue
                
                if fluxes[i,j] < 0 or fluxes[i+1,j] < 0: #bad data
                    sys.exit("from_differential_to_integral_flux: "
                            + f"Bad flux data value of {fluxes[i,j]} and "
                            + f"{fluxes[i+1,j]} found for bin {i}, {j}. "
                            + "This should not happen. Did you call check_for_bad_data() first?")


                if fluxes[i,j] == 0 or fluxes[i+1,j] == 0: #add 0 flux
                    ninc = ninc + 1
                    continue

                F1 = fluxes[i,j]
                F2 = fluxes[i+1,j]
                if F1 == 0 or pd.isnull(F1):
                    sys.exit("from_differential_to_integral_flux: found bin flux of "
                            + f"{fluxes[i,j]}. Should not happen here, "
                            + f"bin [i,j] [{i},{j}]." )

                if F2 == 0 or pd.isnull(F2):
                    sys.exit("from_differential_to_integral_flux: found bin flux of "
                            + f"{fluxes[i+1,j]}. Should not happen here, "
                            + f"bin [i+1,j] [{i+1},{j}]." )
                
                logF1 = np.log(F1)
                logF2 = np.log(F2)
                logE1 = np.log(bin_center[i])
                logE2 = np.log(bin_center[i+1])
                endE = bin_center[i+1]
                if i+1 == nbins-1:
                    endE = energy_bins[nbins-1][1] #extend to edge of last bin

                f = lambda x:exp(logF1
                            + (np.log(x)-logE1)*(logF2-logF1)/(logE2-logE1))
                startE = max(bin_center[i],min_energy)
                fint = scipy.integrate.quad(f,startE,endE)
                if math.isnan(fint[0]):
                    print("from_differential_to_integral_flux: flux integral"
                        "across bins is NaN. Setting to zero. Bin values are "
                        + str(F1) + ' and ' + str(F2))
                    fint = [0]
                if fint[0] < 1e-10:
                    fint = [0]
                sum_flux = sum_flux + fint[0]
                ninc = ninc + 1

        #if last bin is integral, add (HEPAD)
        if energy_bins[nbins-1][1] == -1 and fluxes[nbins-1,j] >= 0:
            intflx = fluxes[nbins-1,j]
            sum_flux = sum_flux + intflx
            ninc = ninc + 1

        if ninc == 0:
            sum_flux = -1
        integral_fluxes.append(sum_flux)

    return integral_fluxes



######### Functions for calculating the onset peak ########
def residual(fit, data):
    """ Calculate difference between fit and data.
    
    """
    
    resid = []
    for i in range(len(fit)):
        resid.append(data[i] - fit[i])
    
    return resid


def modified_weibull(times, Ip, a, b):
    """ Create a Weibull for times in seconds since first
        date with params a and b.
        
    """
    
    weibull = []
    for t in times:
        W = Ip*(-a/b)*(t/b)**(a-1)*math.exp(-(t/b)**a)
        weibull.append(W)
        
    return weibull


def lognormal(times, Ip, a, b):
    """ Create a lognormal fit for times in seconds since first
        date with params a and b. (a=sigma and b=mu)
        
    """
    
    func = []
    for t in times:
        LN = Ip/(t*a*(2.*math.pi)**0.5) * math.exp(-((math.log(t)-b)**2)/(2*a**2))
        func.append(LN)
        
    return func


def func_residual(params, *args):
    """ Caluate the residual of the Weibull fit
        compared to data.
        
    """
    
    pars = params.valuesdict()
    a = pars['alpha']
    b = pars['beta']
    Ip = pars['peak_intensity']
    
    times = args[0]
    data = args[1]
    
    fit = modified_weibull(times, Ip, a, b)
#    fit = lognormal(times, Ip, a, b)
    
    resid = residual(fit, data)

    return resid
    


def find_max_curvature(x, y):
    """ Calculate the curvature along a curve
        and find the maximum curvature location.
        
        https://undergroundmathematics.org/glossary/curvature
        
        INPUTS:
        
        :y: (float 1xn array) weibull fit points
    
    """
    xarr = np.array(x)
    yarr = np.array(y)
    yderiv = yarr[1:] - yarr[:-1]
    yderiv2 = yderiv[1:] - yderiv[:-1]
    
    k_x = yderiv2 #/((1 + yderiv[1:]**2))**(3./2.)
    
    max_k_idx= np.argmin(k_x)
    
    #rescale the curvature to overplot
    max_y = np.max(yarr)
    max_y_idx = np.argmax(yarr)
    k_x = (np.max(yarr)/np.max(k_x))*k_x
            
    return max_k_idx+2
