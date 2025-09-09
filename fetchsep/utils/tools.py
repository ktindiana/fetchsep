from . import config as cfg
from ..json import ccmc_json_handler as ccmc_json
import pandas as pd
import datetime
import os
import sys
import math
import numpy as np
from statistics import mode
import scipy
from numpy import exp

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


##### NAMING SCHEMA #####
def setup_modifiers(options, spacecraft="", doBGSubOPSEP=False, doBGSubIDSEP=False,
        OPSEPEnhancement=False, IDSEPEnhancement=False):
    """ Add modifier strings according to options.
    
    """
    modifier = '' #for appending to filenames
    title_mod = '' #for appending to plot titles

    doBGSub = False
    if doBGSubOPSEP == True or doBGSubIDSEP == True:
        doBGSub = True

    module = ''
    if doBGSubOPSEP or OPSEPEnhancement:
        module = 'opsep'
    if doBGSubIDSEP or IDSEPEnhancement:
        module = 'idsep'

    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '
    if spacecraft:
        modifier = modifier + '_' + spacecraft
        title_mod = title_mod + spacecraft + ' '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if IDSEPEnhancement or OPSEPEnhancement:
        modifier = modifier + '_enhancement'

    if module != '':
        modifier = modifier + '_' + module
        title_mod = title_mod + ' (' + module + ')'

    return modifier, title_mod


def setup_energy_bin_label(energy_bin):
    """ Label for a single energy bin.
    
    """
    label = ""
    if energy_bin[1] != -1:
        label = (f"{energy_bin[0]}-{energy_bin[1]} {cfg.energy_units}")
    else:
        label = (f">{energy_bin[0]} {cfg.energy_units}")

    return label


def energy_bin_key(bin):
    """ Create key for dataframe or columns header in dataframe and
        csv files.
        
    """
    return f"{bin[0]}-{bin[1]}"


def idsep_naming_scheme(experiment, flux_type, exp_name, options, spacecraft=''):
    """ Create naming scheme for subfolders in output/idsep """

    modifier, title_mod = setup_modifiers(options, spacecraft=spacecraft)
    name = (f"{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        name = (f"{exp_name}_{flux_type}{modifier}")
    
    return name
    

def opsep_subdir(experiment, flux_type, exp_name, options,
    spacecraft='', doBGSubOPSEP=False, doBGSubIDSEP=False, OPSEPEnhancement=False,
    IDSEPEnhancement=False):

    modifier, title_mod = setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    #Return a corresponding directory name
    dir = (f"{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        name = (f"{exp_name}_{flux_type}{modifier}")
 
    return dir


def opsep_naming_scheme(date, suffix, experiment, flux_type, exp_name, options,
    spacecraft='', doBGSubOPSEP=False, doBGSubIDSEP=False, OPSEPEnhancement=False,
    IDSEPEnhancement=False):
    """ Create naming scheme for subfolders in output/idsep """

    tzulu = ccmc_json.make_ccmc_zulu_time(date)
    tzulu = tzulu.replace(":","")

    modifier, title_mod = setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    name = (f"{tzulu}_{experiment}_{flux_type}{modifier}_{suffix}")
    if experiment == 'user' and exp_name != '':
        name = (f"{tzulu}_{exp_name}_{flux_type}{modifier}_{suffix}")

    #Return a corresponding directory name
    dir = opsep_subdir(experiment, flux_type, exp_name, options, spacecraft=spacecraft,
        doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    return name, dir


def write_fluxes(experiment, flux_type, exp_name, options, energy_bins, dates, fluxes,
    module, spacecraft=""):
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
    name = idsep_naming_scheme(experiment, flux_type, exp_name, options,
                    spacecraft=spacecraft)
    stdate = dates[0].strftime("%Y%m%d")
    enddate = dates[-1].strftime("%Y%m%d")
    fname = (f"fluxes_{name}_{stdate}_{enddate}.csv")
    if module == 'idsep':
        fname = os.path.join(cfg.outpath, module, name, fname)
    else:
        fname = os.path.join(cfg.outpath,module,fname)
        
    keys = []
    for bin in energy_bins:
        keys.append(energy_bin_key(bin))
    
    dict = {"dates":dates}
    for i in range(len(fluxes)):
        dict.update({keys[i]:fluxes[i]})
        
    df = pd.DataFrame(dict)
    df.to_csv(fname, index=False)
    print("Wrote " + fname + " to file.")
    
    
def from_differential_to_integral_flux(experiment, min_energy, energy_bins,
    fluxes, bruno2017=False, energy_bin_centers=[]):
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
    if len(energy_bin_centers) != 0:
        bin_center = energy_bin_centers
    else:
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



######### SEP event identification ################
def identify_sep_above_background_one(dates, fluxes):
    """ Identify which increases above backgrounds
        are SEP events.
        
        Used in IDSEP and OpSEP.
        
        INPUTS:
        
        :dates: (1xn datetime array) dates for each flux point
        :fluxes: (1xn float array) m energy channels and n time points
        
        OUTPUTS:
        
        :dates: (1xn datetime array) same as in
        :fluxes_sep: (1xn float array) all points set to zero except
            those identified as SEPs
            
    """
    time_res = determine_time_resolution(dates)
    print("Time resolution of the data set is: "
            + str(time_res.total_seconds()) + " seconds.")
    time_res_sec = time_res.total_seconds()
    
    #DEPENDS ON TIME RESOLUTION
    #CAN BE DIFFICULTIES IN IDENTIFYING SEP EVENTS IN VERY
    #GAPPY DATA
    time_increase = 86400/4 #86400 #Require an increase above threshold for 6 hrs
    if time_res_sec <= 60*60:
        time_increase = 3*60*60 #3 hr increase for hi-res data
    global nconsec
    nconsec = max(math.ceil(time_increase/time_res_sec) + 1,3) #num points
        #3 consecutive points for Voyager
    global allow_miss
    #In this stricter SEP event identification used in OpSEP, allow less
    #missed points for the identification of SEP onset
    #20 minutes of points for 5 minute data
    allow_miss = max(math.ceil(nconsec/10),1)
    if nconsec == 2 or nconsec == 3: allow_miss = 0

    #This algorithm is being used to identify enhancements above background using
    #points that are n*sigma above background. There will be sporadic points
    #that are above this level just due to statistical fluctuations. The dwell
    #time can't be too long, or those random fluctuations will extend events
    #artificially. For use in OpSEP, which aims to more precisely identify
    #SEP events, use a shorter dwell time.
    global dwell_pts
    dwell_time = 1*60.*60. #hr
    dwell_pts = int(dwell_time/time_res_sec)
    if time_res_sec > dwell_time:
        dwell_pts = 2

    print("Requiring " + str(nconsec) + " points (" + str(nconsec*time_res_sec/(60.*60.)) + " hours) to define an onset.")
    print("Allowing " + str(allow_miss) + " points to be missed in onset definition.")
    print("Event ends after " + str(dwell_pts) + " points are below threshold (dwell time).")
                
    npts = len(dates)

    IsSPE = False
    SPEflag = False
    
    SPEstart = pd.NaT
    SPEend = pd.NaT
    SPEfluxes = [0]*npts
    stidx = 0
    endidx = 0

    if npts < (nconsec+dwell_pts+allow_miss):
        print(f"identify_sep_above_background_one: Time series is too short {npts} to "
            f"identify a SEP event with the necessary requirements ({nconsec+dwell_pts+allow_miss}).")
        return SPEstart, SPEend, SPEfluxes

    for i in range(npts-nconsec):
        #Condition to identify start of SEP
        if fluxes[i] > 0 and not IsSPE:
            #Check that the increase continues
            #Flux below threshold or nan will be counted as a miss
            nhit = 0
            for k in range(i, i+nconsec):
                chk_flux = fluxes[k]
                if chk_flux > 0:
                    nhit = nhit + 1
                elif pd.isnull(chk_flux):
                    continue
                elif chk_flux == cfg.badval:
                    continue
            
            #too many zero points, not an SPE
            if nhit >= (nconsec - allow_miss): IsSPE = True
        
        #Identify the start of an SPE
        if IsSPE and not SPEflag:
            SPEflag = True
            stidx = i
            SPEstart = dates[i]
            i = i+nconsec #jump to end of required consecutive points
 
 
        #ONGOING SPE with allowed gap
        if IsSPE and SPEflag:
            if fluxes[i] == 0: #don't consider nan or badval; bg set to zero
                IsSPE = False
                end_dwell = min(i+dwell_pts,npts-1)
                for ii in range(i,end_dwell):
                    chk_flux = fluxes[ii]
                    if chk_flux > 0: IsSPE = True

            if not IsSPE or i == npts-1:
                SPEflag = False
                endidx = i
                SPEend = dates[i]
                if i== npts-1:
                    print("WARNING!! identify_sep_above_background_one: SEP event ended "
                        "at end of file. Consider extending the timeframe to get a "
                        "more accurate end time.")
            
                #Fill in the SEP flux array with the SEP points
                for kk in range(stidx, endidx+1):
                    SPEfluxes[kk] = fluxes[kk]

                #Return the first SEP found. If not returned here, code
                #continue looking for next SEP.
                print("identify_sep_above_background_one: SEP event found from "
                    f"{SPEstart} and {SPEend}.")
                return SPEstart, SPEend, SPEfluxes


    return SPEstart, SPEend, SPEfluxes



def identify_sep_above_background(dates, fluxes):
    """ Identify which increases above backgrounds
        are SEP events.
        
        Used in IDSEP and OpSEP.
        
        INPUTS:
        
        :dates: (1xn datetime array) dates for each flux point
        :fluxes: (mxn float array) m energy channels and n time points
        
        OUTPUTS:
        
        :dates: (1xn datetime array) same as in
        :fluxes_sep: (mxn float array) all points set to zero except
            those identified as SEPs
            
    """
    time_res = determine_time_resolution(dates)
    print("Time resolution of the data set is: "
            + str(time_res.total_seconds()) + " seconds.")
    time_res_sec = time_res.total_seconds()
    
    #DEPENDS ON TIME RESOLUTION
    #CAN BE DIFFICULTIES IN IDENTIFYING SEP EVENTS IN VERY
    #GAPPY DATA
    time_increase = 86400/4 #86400 #Require an increase above threshold for duration
    if time_res_sec <= 60*60:
        time_increase = 3*60*60 #6*60*60 #6 hr increase for hi-res data
    global nconsec
    nconsec = max(math.ceil(time_increase/time_res_sec) + 1,3) #num points
        #3 consecutive points for Voyager
    global allow_miss
    allow_miss = max(math.ceil(nconsec/4),1)
    if nconsec == 2 or nconsec == 3: allow_miss = 0
    global dwell_pts
    dwell_pts = max(math.ceil(nconsec/2),2)
    #Rosetta? needs --> dwell_pts = max(math.ceil(nconsec),2)

    print("Requiring " + str(nconsec) + " points (" + str(nconsec*time_res_sec/(60.*60.)) + " hours) to define an onset.")
    print("Allowing " + str(allow_miss) + " points to be missed in onset definition.")
    print("Event ends after " + str(dwell_pts) + " points are below threshold (dwell time).")
                
    npts = len(dates)

    IsSPE = False
    SPEflag = False
    
    SPEstart = [[]]*len(fluxes)
    SPEend = [[]]*len(fluxes)
    SPEfluxes = [[]]*len(fluxes)
    stidx = 0
    endidx = 0
    for j in range(len(fluxes)):
        SPEfluxes[j] = [0]*npts
        for i in range(npts-nconsec):
            if fluxes[j][i] > 0 and not IsSPE:
                IsSPE = True
                
                #Check that the increase continues
                #Flux below threshold will be set to zero
                nmiss = 0
                for k in range(i, i+nconsec):
                    chk_flux = fluxes[j][k]
                    if chk_flux <= 0 or pd.isnull(chk_flux):
                        nmiss = nmiss + 1
                
                #too many zero points, not an SPE
                if nmiss > allow_miss: IsSPE = False
            
            #Identify the start of an SPE
            if IsSPE and not SPEflag:
                SPEflag = True
                stidx = i
                if not SPEstart[j]:
                    SPEstart[j] = [dates[i]]
                else:
                    SPEstart[j].append(dates[i])
                i = min(i+nconsec,npts-nconsec-1) #jump to end of required consecutive points
              
            #ONGOING SPE with allowed gap
            if IsSPE and SPEflag:
                if fluxes[j][i] <= 0:
                    IsSPE = False
                    end_dwell = min(i+dwell_pts,npts-nconsec-1)
                    for ii in range(i,end_dwell):
                        chk_flux = fluxes[j][ii]
                        if chk_flux > 0: IsSPE = True
                        
                if not IsSPE or i == npts-nconsec-1:
                    SPEflag = False
                    endidx = i
                    if not SPEend[j]:
                        SPEend[j] = [dates[i]]
                    else:
                        SPEend[j].append(dates[i])
                
                    #Fill in the SEP flux array with the SEP points
                    for kk in range(stidx, endidx+1):
                        SPEfluxes[j][kk] = fluxes[j][kk]
                    
    return SPEstart, SPEend, SPEfluxes



def identify_sep_noaa(dates, fluxes, threshold):
    """ Follow SWPC approach to identifying event start and end 
        above threshold.
        
        Used in OpSEP.
    
    """
    threshold_crossed = False
    event_ended = False
    ndates = len(dates)
    sep_start_time = pd.NaT
    sep_end_time = pd.NaT

    end_threshold = cfg.endfac*threshold
            #endfac = 1.0 to get SWPC definition of event end for 5 min data
            #endfac = 0.85 used by SRAG operators in an alarm code
            #Included for flexibility, but 1.0 is used for operational values

    time_res = determine_time_resolution(dates)
    print("Time resolution of the data set is: "
            + str(time_res.total_seconds()) + " seconds.")
    time_res_sec = time_res.total_seconds()
    
    npoints = 3 #require 3 points above threshold as employed by SWPC
    if time_res_sec/60. > 15:
        npoints = 1 #time resolution >15 mins, require one point above threshold

    for i in range(ndates):
        if not threshold_crossed:
            if(fluxes[i] >= threshold):
                start_counter = 0
                if i+(npoints-1) < ndates:
                    for ii in range(npoints):
                        if fluxes[i+ii] >= threshold:
                            start_counter = start_counter + 1
                if start_counter == npoints:
                    sep_start_time = dates[i]
                    threshold_crossed = True
        if threshold_crossed and not event_ended:
            if (fluxes[i] >= end_threshold):
                end_counter = 0  #reset if go back above threshold
                end_tm0 = dates[i] #will catch the last date above threshold
            if (fluxes[i] <= end_threshold): #flux drops below endfac*threshold
                end_counter = end_counter + 1
                elapse = (dates[i]  - end_tm0).total_seconds()
                #When looking for an end of an event below an operational threshold
                #or threshold that isn't the background, apply a dwell time to ensure
                #the fluxes don't fluctuate above threshold again after a little while.
                if elapse > cfg.dwell_time: #N consecutive points longer than dwell time
                    event_ended = True
                    sep_end_time = dates[i-(end_counter-1)] #correct back time steps
                    #Double checked some calculated event end times with SWPC and
                    #this logic gave the correct end times. 2023-04-10 KW

    return sep_start_time, sep_end_time


######### Functions for calculating the onset peak ########
def residual(fit, data):
    """ Calculate difference between fit and data.
    
    """
    
    resid = []
    for i in range(len(fit)):
        resid.append(data[i] - fit[i])
    
    return resid


def ratio(fit, data):
    """ Calculate difference between fit and data.
    
    """
    
    resid = []
    for i in range(len(fit)):
        if data[i] != 0:
            err = abs(data[i] - fit[i])/data[i]
            resid.append(err)
 
    resid = sum(resid)/len(resid)
 
    return resid


def normchisq(fit, data, sigma=np.nan):
    """ Calculate difference between fit and data.
    
    """
    
    resid = []
    for i in range(len(fit)):
        if not pd.isnull(sigma):
            err = ((data[i] - fit[i])**2)/sigma
            resid.append(err)
    
        elif data[i] != 0:
            err = ((data[i] - fit[i])**2)/data[i]
            resid.append(err)
 
    resid = sum(resid)/(len(resid)-1)
 
    return resid


def determination(fit, data):
    """ Calculate difference between fit and data.
    
    """
    
    fiterr = []
    dataerr = []
    meandata = [x for x in data if (x > 0 and not pd.isnull(x))]
    meandata = sum(meandata)/len(meandata)
    for i in range(len(fit)):
        if data[i] != 0 and not pd.isnull(data[i]):
            fiterr.append((data[i] - fit[i])**2)
            dataerr.append((data[i] - meandata)**2)
    
    resid = 1 -  sum(fiterr)/sum(dataerr)
 
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
#    resid = normchisq(fit, data)

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
