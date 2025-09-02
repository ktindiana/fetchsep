from . import config as cfg
from . import error_check
from . import read_datasets as datasets
from . import tools
from . import derive_background_opsep as bgsub
from . import plotting_tools as plt_tools
from ..json import ccmc_json_handler as ccmc_json
import pandas as pd
import numpy as np
import os
import sys
from lmfit import minimize, Parameters
import datetime
import math
import pickle

#########################################################
################# General-use Classes ###################
#########################################################
class EnergyBin:
    def __init__(self, min, max, units, bin_center=np.nan):
        """ Energy Bin 
        
            INPUT:
            
                :min: (float) low edge of energy bin
                :max: (float) high edge of energy bin
                    set to -1 for integral channel
                    
            OUTPUT:
            
                EnergyBins object
        
        """
    
        self.min = float(min)
        self.max = float(max)
        self.units = units
        self.bin_center = bin_center
        if pd.isnull(self.bin_center):
            if self.max != -1 and (not pd.isnull(self.min) and
                not pd.isnull(self.max)):
                self.bin_center = math.sqrt(min*max) #geometric mean
            
        self.energy_channel = {'min': float(min), 'max': float(max), 'units': units}
        
        return


class Threshold:
    def __init__(self, threshold, units):
        """ Float flux threshold value and units. """
        
        self.threshold = float(threshold)
        self.threshold_units = units
        self.threshold_dict = {'threshold':float(threshold), 'threshold_units': units}
        
        return



################################################################
############## Data Class: Flux Data for OpSEP #################
################################################################
class Data:
    def __init__(self):
        """Read in input data. Contains:
            - all fluxes read into OpSEP
            - all SEP event definitions
            - any estimated integral fluxes (if differential input)
            - all fluxes in the energy bins relevant for the event 
                definitions, i.e. evaluated fluxes
           
           INPUT:
                
                :experiment: (str) name of spacecraft or "user"
                :flux_type: (str) integral or differential
                
            OUTPUT:
            
                a Data object
        
        """

        self.label = None
        self.experiment = None
        self.flux_type = None
        self.spacecraft = '' #GOES only; primary or secondary
        self.startdate = pd.NaT
        self.enddate = pd.NaT
        self.min_energy = np.nan #If set, only use bins above this value
        self.max_energy = np.nan #If set, only use bins below this value

        self.location = None #earth, mars, etc
        self.species = None #protons, electrons

        #Paths for data source, plots and output files
        self.datapath = cfg.datapath
        self.outpath = cfg.outpath
        self.plotpath = cfg.plotpath
        
        #Variables for user-input files
        self.user_delim = cfg.user_delim #delimeter between columns
        self.user_col = cfg.user_col #columns in input data file containing fluxes
        #self.err_col = cfg.err_col #columns containing error bars
        self.user_energy_bins = cfg.user_energy_bins #energy bins for each column

        #User-input data file
        self.user = False #user file?
        self.user_name = None
        self.user_filename = None #If user-input file
        
        #Original fluxes before any background subtraction or interpolation
        #Bad values are set to None
        self.original_dates = []
        self.original_fluxes = []
        
        #Do background subtraction on the fluxes?
        #With OPSEP
        self.doBGSubOPSEP = False #True to do background-subtraction
        self.bgstartdate = pd.NaT #defaults to idsep output if not set
        self.bgenddate = pd.NaT #defaults to idsep output if not set
        #With IDSEP
        self.doBGSubIDSEP = False
        self.idsep_path = ''
        
        #BACKGROUND
        self.nsigma = cfg.opsep_nsigma #N * sigma to use in background subtraction
        self.bgmeans = [] #Mean fluxes for each energy channel
        self.bgsigmas = [] #Sigmas for each energy channel
        self.bgdates = []
        self.bgfluxes = []

        #ENHANCEMENT ABOVE BACKGROUND
        self.OPSEPenhancement = False
        self.IDSEPenhancement = False #Get threshold from idsep output files
        #Use threshold to set background to zero and leave only enhancement

        #Apply linear interpolation in time to data gaps when
        #time steps exist that have bad data values (negative, NaN)
        self.do_interpolation = True
        
        #Two peaks may extend an event if it crosses threshold,
        #temporarily drops below, then increases above threshold again
        self.two_peaks = False
        
        #Plotting settings
        self.showplot = False
        self.saveplot = False
        
        #GOES-specific options
        self.options = []
        self.goes_datatype = None #corrected or uncorrected
        self.goes_S14 = None #True to apply Sandberg et al. 2014 effective energies
        self.goes_Bruno2017 = None #True to apply Bruno et al. 2017 effective energies
        
        #Energy channels and thresholds used for SEP event definitions
        #Dictionaries with energy channel and associated threshold
        #{'energy_channel': energy_bin_obj, 'threshold': threshold_obj}
        self.event_definitions = []

        #Filenames of background fluxes and thresholds from IDSEP
        self.idsep_background = None
        self.idsep_threshold = None

        #The flux timeseries after interpolation and any background subtraction
        self.energy_bin_objects = [] #include bin centers
        self.energy_bins = []
        self.energy_bin_centers = []
        self.dates = []
        self.fluxes = []
        self.time_resolution = np.nan #seconds
        
        #Flux timeseries to be evaluated from event definitions
        self.evaluated_energy_bins = []
        self.evaluated_dates = []
        self.evaluated_fluxes = []

        #Collect the results as individual Analyze objects for
        #each event definition
        self.results = []
        #If a SEP event happens in any channel, day is stored here to return
        self.sep_year = np.nan
        self.sep_month = np.nan
        self.sep_day = np.nan

        return


    #Allow the user to change various values that are in the config file
    def set_datapath(self, datapath): #location to measurement data (e.g. GOES)
        """ Set the path containing the data downloaded and read by FetchSEP """
        self.datapath = datapath
        return


    def set_outpath(self, outpath):
        """ Set the path to output files """
        self.outpath = outpath
        return


    def set_plotpath(self, plotpath):
        """ Set the path to the output plots """
        self.plotpath = plotpath
        return


    def set_user_delim(self, user_path):
        """ Set the path to the output plots """
        self.user_delim = user_delim
        return


    def set_user_col(self, user_col):
        """ Set the path to the output plots """
        self.user_col = user_col
        return


    def str_to_datetime(self, date):
        """ Convert string date to datetime. """
        
        if date == '':
            return pd.NaT
        
        if len(date) == 10: #only YYYY-MM-DD
            date = date  + ' 00:00:00'
        dt = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
        return dt


    def set_dates(self, startdate, enddate):
        """ Set the start and end dates.
        
            INPUT:
            
                :startdate: (str) YYYY-MM-DD or YYYY-MM-DD HH:MM:SS
                :enddate: (str) YYYY-MM-DD or YYYY-MM-DD HH:MM:SS
                
            OUTPUT:
            
                Set self.startdate and self.enddate as datetime
        
        """
        self.startdate = self.str_to_datetime(startdate)
        self.enddate = self.str_to_datetime(enddate)
        return


    def set_options(self, options):
        """ Set options (arr) if any specified """

        #First, reset
        self.options = []
        self.goes_datatype = None
        self.goes_S14 = None
        self.goes_Bruno2017 = None
    
        options = options.split(";")
        if options[0] != "": self.options = options
        
        if 'GOES' in self.experiment:
            if 'uncorrected' in options:
                self.goes_datatype = 'uncorrected'
            else:
                self.goes_datatype = 'corrected'
                
            if 'S14' in options:
                self.goes_S14 = True
            else:
                self.goes_S14 = False
                
            if 'Bruno2017' in options:
                self.goes_Bruno2017 = True
            else:
                self.goes_Bruno2017 = False

        return


    def set_opsep_background_info(self, doBGSubOPSEP, OPSEPEnhancement,
        bgstartdate, bgenddate):
        """ Indicate whether to perform background-subtraction.
            If start and end dates aren't set, then code
            will look for idsep output to use for mean background.
        
            INPUT:
            
                :doBGSubOPSEP: (bool) bg subtraction if True
                :OPSEPEnhancement: (bool) will use background mean + n*sigma
                    to separate background and SEP, with or without background subtraction
                :bg_startdate: (str) start of time period to use for 
                    background calculation
                :bg_enddate: (str) end of time period to use for
                    background calculation
                    
            OUTPUT:
            
                Set background attributes in Data object
        
        """
        self.doBGSubOPSEP = doBGSubOPSEP
        self.bgstartdate = self.str_to_datetime(bgstartdate)
        self.bgenddate = self.str_to_datetime(bgenddate)
        self.OPSEPEnhancement = OPSEPEnhancement
        #IF choose to do background subtraction, then automatically choose
        #to calculate enhancement above background
        if self.doBGSubOPSEP: self.OPSEPEnhancement = True
        
        if doBGSubOPSEP or OPSEPEnhancement:
            if pd.isnull(self.bgstartdate) or pd.isnull(self.bgenddate):
                sys.exit("WARNING!!! User selected to perform background-subtraction, but did not provide dates. Please provide dates of a quiet background period to use this feature.")
        
        return


    def set_idsep_background_info(self, doBGSubIDSEP, idsep_path, IDSEPEnhancement):
        """ Specify whether to use background calculated by idsep """

        #If want to use IDSEP files, but no path specified, try the default
        if (IDSEPEnhancement or doBGSubIDSEP) and idsep_path == '':
            name = tools.idsep_naming_scheme(self.experiment, self.flux_type, self.user_name,
                    self.options, spacecraft=self.spacecraft)
            idsep_path = os.path.join(cfg.outpath, 'idsep', name, 'csv')

        self.doBGSubIDSEP = doBGSubIDSEP
        self.idsep_path = idsep_path
        self.IDSEPEnhancement = IDSEPEnhancement
        #IF choose to do background subtraction, then automatically choose
        #to calculate enhancement above background
        if self.doBGSubIDSEP: self.IDSEPEnhancement = True

        return
    

    def create_event_definition(self, emin, emax, eunits, thresh, thresh_units):
        """ Easily create a single event definition for flexible
            use of methods from the command line. This subroutine
            is not used in opsep. It is meant to be a helper when 
            interfacing with the object manually.
            
            For units, refers to units specified in fetchsep.cfg.
           
            INPUT:
            
                :energy_bin: (arr) [Emin, Emax]
                :threshold: (float) flux threshold
            
        """

        energy_bin_obj = EnergyBin(emin,emax,eunits)
        threshold_obj = Threshold(thresh, thresh_units)
        
        event_definition = {'energy_channel': energy_bin_obj, 'threshold': threshold_obj}

        return event_definition


    def set_event_definitions(self, definitions, reset=True):
        """ Set the energy channel and applied threshold
            for each event definition requested by the user.
            FORMAT: "30,1;50,1;4.5-9.2,0.1"
            
            user-input thresholds in the format "30,1" for >30 MeV exceeds 1 pfu, 
            "4-7,0.01" for 4-7 MeV differential channel exceeds 0.01.  
            "30,1;4-7,0.01" multiple thresholds separated by semi-colon.
            
            Always apply default operational definitions:
                >10 MeV exceeds 10 pfu
                >100 MeV exceed 1 pfu 
            
            Units will be pulled from the config file for the
            appropriate flux type.
            
            A threshold set to -1 in the event definition indicates to use
            the mean and sigma calculated by idsep as a threshold. The 
            threshold can vary with time as a value is calculated for each
            timestep. Therefore a constant value cannot be applied.
            
            INPUT:
            
                :definitions: (str) user-input definitions in opsep
                    format e.g. "10,1;4.5-9.2,0.1"
                    Definitions separated by semi-colons and
                    energy channel and threshold separated by
                    commas
                :reset: (bool) True means all existing event definitions
                    will be reset and filled in again. False means 
                    provided event definitions will be appended to
                    existing.
                
            Output:
            
                Fill self.event_definitions Data attribute
                
        """
        #First, reset
        if reset:
            self.event_definitions = []
            #Default operational thresholds
            bins = [[10.0,-1], [100.0,-1]]
            thresholds = [10.0, 1.0]
        else:
            bins = []
            thresholds = []
        
        definitions = definitions.strip().split(";")
        if definitions[0] != "":
            for evdef in definitions:
                evdef = evdef.strip().split(",")
                if "-" in evdef[0]:
                    bin_edge = evdef[0].strip().split("-")
                    bins.append([float(bin_edge[0]), float(bin_edge[1])])
                else:
                    bins.append([float(evdef[0]), -1])
                    
                thresholds.append(float(evdef[1]))
        

        for index, (bin, thresh) in enumerate(zip(bins, thresholds)):
            energy_units = cfg.energy_units
            flux_units = ''
            if -1 in bin:
                flux_units = cfg.flux_units_integral
            else:
                flux_units = cfg.flux_units_differential
                
            energy_bin_obj = EnergyBin(bin[0],bin[1],energy_units)
            threshold_obj = Threshold(thresh, flux_units)
            
            self.event_definitions.append({'energy_channel': energy_bin_obj, 'threshold': threshold_obj})

            #If use IDSEP thresholds, append event definitions to all
            #previously selected energy channels and use -1 as threshold
            if self.IDSEPEnhancement or self.OPSEPEnhancement:
                idsep_threshold_obj = Threshold(-1, flux_units)
                self.event_definitions.append({'energy_channel': energy_bin_obj, 'threshold': idsep_threshold_obj})

        return

        
    def error_check(self):
        """ Error check the inputs and options. """
            
        error_check.error_check_options(self.experiment, self.flux_type, self.options, spacecraft=self.spacecraft)
        error_check.error_check_inputs(self.startdate, self.enddate, self.experiment, self.flux_type)
        error_check.error_check_background(self.experiment, self.flux_type, self.doBGSubOPSEP,
            self.doBGSubIDSEP, self.OPSEPEnhancement, self.IDSEPEnhancement)

        return


    def load_info(self, startdate, enddate, experiment, flux_type,
        user_name='', user_file='', spase_id='', showplot=False, saveplot=False,
        two_peaks=False, definitions='', options='',
        doBGSubOPSEP=False, OPSEPEnhancement=False, bgstartdate='', bgenddate='',
        doBGSubIDSEP=False, IDSEPEnhancement=False, idsep_path='output/idsep/csv/',
        nointerp=False, spacecraft='', location='earth', species='proton'):
        """ Create new Data object and load with all values.
        
            INPUT:
                :startdate: (string) - user input start date 
                    "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS"
                :enddate: (string) - user input end date "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS"
                :experiment: (string) - "GOES-05" up to "GOES-19", "SEPEM", "SEPEMv3","EPHIN", "EPHIN_REleASE", or "user"
                :flux_type: (string) - "integral" or "differential" 
                    indicates the type of flux to read in
                :user_name: (string) - If model is "user", set 
                    user_name to describe your model or data set (e.g. 
                    MyModel), otherwise set to ''.
                :user_file: (string) - Default is ''. If "user" is 
                    selected for experiment, specify name of flux file.
                :showplot: (bool) - True to show plots
                :saveplot: (bool) - True to save plots 
                :two_peaks: (bool) - option for extending event length
                :definitions: (string) - user-input thresholds in the 
                    format "30,1;4-7,0.01" multiple thresholds
                    separated by semi-colon. Same as user_thresholds in opsep
                :nointerp: (boolean) - True to fill in negative fluxes 
                    with None instead of linear interpolation in time
                :spacecraft: (string) primary or secondary 
                
            OUTPUT:
            
                :input_data: (Data Object)
        
        """
        
        if experiment == "user":
            self.user = True
            self.label = f"{user_name} {flux_type}"
        else:
            self.label = f"{experiment} {flux_type}"

        self.experiment = experiment
        self.flux_type = flux_type
        self.user_name = user_name #user input experiment (model or obs name)
        self.user_filename = user_file
        self.set_dates(startdate, enddate)
        self.user_filename = user_file #default tmp.txt
        self.spacecraft = spacecraft
        self.location = location
        self.species = species
        self.set_options(options)
        self.set_opsep_background_info(doBGSubOPSEP, OPSEPEnhancement, bgstartdate, bgenddate)
        self.set_idsep_background_info(doBGSubIDSEP, idsep_path, IDSEPEnhancement)
        self.set_event_definitions(definitions)
        self.two_peaks = two_peaks
        self.showplot = showplot
        self.saveplot = saveplot
        self.do_interpolation = not(nointerp)
        self.error_check()
 
        #Create subdirectory to hold values
        subdir = tools.opsep_subdir(self.experiment, self.flux_type,
            self.user_name, self.options, spacecraft=self.spacecraft,
            doBGSubOPSEP=self.doBGSubOPSEP, doBGSubIDSEP=self.doBGSubIDSEP,
            OPSEPEnhancement=self.OPSEPEnhancement,
            IDSEPEnhancement=self.IDSEPEnhancement)
        if not os.path.exists(os.path.join(cfg.outpath,'opsep', subdir)):
            os.mkdir(os.path.join(cfg.outpath,'opsep', subdir))
        if not os.path.exists(os.path.join(cfg.plotpath,'opsep', subdir)):
            os.mkdir(os.path.join(cfg.plotpath,'opsep', subdir))

        
 
        return


    def plot_background_subtraction(self, showplot=False):
        """ Make plots of background-subtracted fluxes """

        plt_tools.opsep_plot_bgfluxes("Total_Fluxes",self.experiment, self.flux_type, self.options,
                self.user_name, self.original_fluxes, self.original_dates,
                self.energy_bins, self.bgmeans, self.bgsigmas, self.saveplot,
                spacecraft = self.spacecraft, doBGSubOPSEP=self.doBGSubOPSEP,
                doBGSubIDSEP=self.doBGSubIDSEP,
                OPSEPEnhancement=self.OPSEPEnhancement,
                IDSEPEnhancement=self.IDSEPEnhancement)
        plt_tools.opsep_plot_bgfluxes("Background_Fluxes", self.experiment, self.flux_type,
                self.options, self.user_name, self.bgfluxes, self.bgdates,
                self.energy_bins, self.bgmeans, self.bgsigmas, self.saveplot,
                spacecraft = self.spacecraft, doBGSubOPSEP=self.doBGSubOPSEP,
                doBGSubIDSEP=self.doBGSubIDSEP,
                OPSEPEnhancement=self.OPSEPEnhancement,
                IDSEPEnhancement=self.IDSEPEnhancement)
        plt_tools.opsep_plot_bgfluxes("SEP_Fluxes", self.experiment,
                self.flux_type, self.options, self.user_name, self.fluxes, self.dates,
                self.energy_bins, self.bgmeans, self.bgsigmas, self.saveplot,
                spacecraft = self.spacecraft, doBGSubOPSEP=self.doBGSubOPSEP,
                doBGSubIDSEP=self.doBGSubIDSEP,
                OPSEPEnhancement=self.OPSEPEnhancement,
                IDSEPEnhancement=self.IDSEPEnhancement)
        
        if showplot:
            plt.show()
        
        return


    def opsep_background_and_sep_separation(self, all_dates, all_fluxes):
        """ Perform background separation using OpSEP for a
            small time period specified by the user. 
            If OpSEPdoBGSub == True, then background will be subtracted
            from SEP fluxes and background values will be set to zero.
            Otherwise, background values will be set to zero and SEP 
            fluxes will remain the same.
            
        """
        #sepfluxes are background subtracted fluxes
        #Previous version of subroutine had read in the data
        #again and did not apply linear interpolation on bad points
        nointerp_fluxes = datasets.check_for_bad_data(all_dates,all_fluxes,
            self.energy_bins, dointerp=False)
        bgfluxes, sepfluxes, means, sigmas = bgsub.derive_background(self.experiment,
            self.flux_type, self.options, self.bgstartdate, self.bgenddate,
            all_dates, nointerp_fluxes, self.energy_bins, self.showplot,
            self.saveplot, nsigma=self.nsigma, doBGSub=self.doBGSubOPSEP)
        self.bgmeans = means
        self.bgsigmas = sigmas
        self.bgfluxes = bgfluxes
        self.bgdates = all_dates
        #Extract the date range specified by the user for the
        #background-subtracted fluxes
        dates, fluxes = datasets.extract_date_range(self.startdate, self.enddate,
                            all_dates, sepfluxes)

        return dates, fluxes


    def read_idsep_files(self):
        """ Read in the idsep files needed for background subtraction and
            for applying the thresholds calculated by idsep as an event
            definition.
            
        """
        bgfilename = os.path.join(self.idsep_path,'background_mean_fluxes_optimized_FINAL.csv')
        sigmafilename = os.path.join(self.idsep_path, 'background_sigma_optimized_FINAL.csv')
        threshfilename = os.path.join(self.idsep_path, 'background_threshold_optimized_FINAL.csv')

        df_mean = pd.read_csv(bgfilename)
        df_mean['dates'] =pd.to_datetime(df_mean['dates'])
        df_sigma = pd.read_csv(sigmafilename)
        df_sigma['dates'] =  pd.to_datetime(df_sigma['dates'])
        df_thresh = pd.read_csv(threshfilename)
        df_thresh['dates'] = pd.to_datetime(df_thresh['dates'])

        df_mean = df_mean.replace(1e6,np.nan)
        df_thresh = df_thresh.replace(1e6,np.nan)
        
        #Trim to date range
        df_mean = df_mean.loc[(df_mean['dates'] >= self.startdate) & (df_mean['dates'] <= self.enddate)]
        df_sigma = df_sigma.loc[(df_sigma['dates'] >= self.startdate) & (df_sigma['dates'] <= self.enddate)]
        df_thresh = df_thresh.loc[(df_thresh['dates'] >= self.startdate) & (df_thresh['dates'] <= self.enddate)]
        if df_mean.empty:
            sys.exit("The idsep file containing the mean background does not cover the "
                    f"dates required. {bgfilename}")
                    
        mean_dates = df_mean['dates'].to_list()
        if mean_dates != self.original_dates:
            sys.exit("The idsep file containing the mean background doesn't have the same "
                    f"dates length. {bgfilename}")
 
        df_mean = df_mean.drop('dates',axis=1) #only flux columns
        
        #Check that the energy bins for the columns match the energy bins
        #of the current data
        columns = df_mean.columns
        for bin in self.energy_bins:
            key = tools.energy_bin_key(bin)
            if key not in columns:
                sys.exit("read_idsep_files: IDSEP files don't contain energy bins that "
                    f"match the data. IDSEP columns: {columns}, Data energy bins: {self.energy_bins}")
        
        means = []
        sigmas = []
        thresholds = []
        for i, col in enumerate(columns):
            bg_flux = df_mean[col].to_list()
            bg_sigma = df_sigma[col].to_list()
            bg_thresh = df_thresh[col].to_list()
            if i==0:
                means = [bg_flux]
                sigmas = [bg_sigma]
                thresholds = [bg_thresh]
            else:
                means.append(bg_flux)
                sigmas.append(bg_sigma)
                thresholds.append(bg_thresh)

        return means, sigmas, thresholds



    def idsep_background_and_sep_separation(self):
        """ Read in idsep files containing mean background with time,
            e.g. output/idsep/
            
            You must have run idsep for exactly the same experiment and 
            flux_type for date ranges that contain your period of interest.
            In this way, the energy bins and data cadence will match.
            
            If IDSEP background subtraction is selected, will 
            return background-subtracted fluxes with background
            set to zero.
            
            If background subtraction not selected, will return
            original fluxes, but with background set to zero.
            
        """

        means, sigmas, thresholds = self.read_idsep_files()
        
        bgfluxes, sepfluxes = bgsub.separate_sep_and_background_idsep(self.original_fluxes,
                            means, sigmas, nsigma=self.nsigma,
                            doBGSub=self.doBGSubIDSEP)
        
        bgfluxes = np.array(bgfluxes)
        sepfluxes = np.array(sepfluxes)
 
        self.bgmeans = means
        self.bgsigmas = sigmas
        self.bgfluxes = bgfluxes
        self.bgdates = self.original_dates
 
        return self.original_dates, sepfluxes
  

    def read_in_flux(self):
        """ Read in the appropriate data or user files. Performs
            background subtraction or background and SEP separation, 
            if requested. Trims to dates between start time and end time. 
            Interpolates bad points with linear interpolation in time.
            
            Loads into self.fluxes
            
            OUTPUTS:
            
            :dates: (datetime 1xm array) - times in flux time profile trimmed
                between startdate and enddate
            :fluxes: (numpy float nxm array) - fluxes for n energy channels and m
                time steps; these are background subtracted fluxes if background
                subtraction was selected.
            :energy_bins: (array nx2 for n thresholds)
            
        """

        detector= []
        west_detector = []
        startdate = self.startdate
        enddate = self.enddate
        #If want to do background subtraction with a specified date range,
        #make sure to read in data for the full date range required
        if self.doBGSubOPSEP:
            if not pd.isnull(self.bgstartdate) and not pd.isnull(self.bgenddate):
                startdate = min(self.startdate,self.bgstartdate)
                enddate = max(self.enddate, self.bgenddate)
        
        
        if self.experiment == "GOES": #Extra output
            filenames1, filenames2, filenames_orien, detector = \
                datasets.check_data(startdate, enddate, self.experiment, self.flux_type, self.user_filename, spacecraft=self.spacecraft)
        else:
            filenames1, filenames2, filenames_orien = datasets.check_data(startdate,
                    enddate, self.experiment, self.flux_type, self.user_filename, spacecraft=self.spacecraft)

                                        
        #read in flux files
        if not self.user:
            if self.experiment == "GOES":
                all_dates, all_fluxes, west_detector, energy_bins, energy_bin_centers = \
                    datasets.read_in_files(self.experiment, self.flux_type,
                            filenames1, filenames2, filenames_orien, self.options,
                            detector=detector, spacecraft=self.spacecraft)
            else:
                all_dates, all_fluxes, west_detector = \
                    datasets.read_in_files(self.experiment, self.flux_type,
                        filenames1, filenames2, filenames_orien, self.options,
                        detector=detector, spacecraft=self.spacecraft)
 
        else:
            all_dates, all_fluxes = datasets.read_in_user_files(filenames1,
                        delim=self.user_delim, flux_col=self.user_col)


        #Define energy bins
        if self.experiment == "ERNE":
            version = datasets.which_erne(startdate, enddate)
            energy_bins, energy_bin_centers = datasets.define_energy_bins(version,
                self.flux_type, west_detector, self.options)
        elif self.experiment != "GOES":
            energy_bins, energy_bin_centers = datasets.define_energy_bins(self.experiment, self.flux_type,
                        west_detector, self.options, spacecraft=self.spacecraft,
                        user_bins=self.user_energy_bins)


        if len(all_dates) <= 1:
            sys.exit(f"read_in_flux: The specified start and end dates ({startdate} to {enddate}) were not present in the specified input file or were too restrictive. Exiting.")

        #Full flux and date range for specified input files, not yet trimmed in date
        all_fluxes, energy_bins, energy_bin_centers = tools.sort_bin_order(all_fluxes, energy_bins, energy_bin_centers)
        
        self.energy_bins = energy_bins
        self.energy_bin_centers = energy_bin_centers
        
        ####Save original fluxes with bad points set to None
        #Extract date range that covers any background-subtraction periods
        print("Reading in original fluxes including any background subtraction periods with no interpolation.")
        orig_dates, orig_fluxes = datasets.extract_date_range(startdate, enddate,
                                        all_dates, all_fluxes)
        orig_fluxes = datasets.check_for_bad_data(orig_dates,orig_fluxes,energy_bins,dointerp=False)
        self.original_dates = orig_dates
        self.original_fluxes = orig_fluxes

        #IF BACKGROUND SUBTRACTION
        if self.doBGSubOPSEP or self.OPSEPEnhancement:
            #Background-subtracted fluxes with date range specified by the user
            dates, fluxes = self.opsep_background_and_sep_separation(all_dates, all_fluxes)
        elif self.doBGSubIDSEP or self.IDSEPEnhancement:
             dates,fluxes = self.idsep_background_and_sep_separation()
        #NO BACKGROUND SUBTRACTION OR USE OF IDSEP BACKGROUND IDENTIFICATION
        else:
            #Extract the date range specified by the user
            dates, fluxes = datasets.extract_date_range(self.startdate,
                    self.enddate, all_dates, all_fluxes)
        
        #Handle bad data points
        print("Removing bad points from final fluxes after any background "
              "subtraction and trimming. "
              f"Performing interpolation? {self.do_interpolation}")
        fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins, dointerp=self.do_interpolation)
         
        if len(dates) <= 1:
            print("read_in_flux: The specified start and end dates were not "
                f"present in the specified input file. Exiting. {startdate} to {enddate}")
            sys.exit()
        
        self.fluxes = fluxes
        self.dates = dates
        
        time_res = tools.determine_time_resolution(dates)
        self.time_resolution = time_res.total_seconds()

        #Plot background and SEP separation, may or may not include background
        #subtraction
        if self.doBGSubOPSEP or self.doBGSubIDSEP or self.OPSEPEnhancement\
        or self.IDSEPEnhancement:
            if self.showplot or self.saveplot:
                self.plot_background_subtraction()

        return


    def estimate_integral_fluxes(self):
        """ If input data is differential, then estimate 
            integral fluxes by doing a linear interpolation
            in log space from energy bin center to bin center.
            Append these estimated fluxes to self.fluxes.
            
        """

        if self.flux_type == "integral":
            return
            
        for evdef in self.event_definitions:
            #If integral channel, estimate integral fluxes
            if evdef['energy_channel'].max == -1:
                energy_threshold = evdef['energy_channel'].min
                integral_flux = tools.from_differential_to_integral_flux(self.experiment,
                            energy_threshold, self.energy_bins, self.fluxes,
                            bruno2017=self.goes_Bruno2017,
                            energy_bin_centers=self.energy_bin_centers)
                
                #Add estimated integral fluxes and associated energy bin to self
                self.evaluated_fluxes.append(integral_flux)
                bin = [evdef['energy_channel'].min, evdef['energy_channel'].max]
                self.evaluated_energy_bins.append(bin)
            
        return



    def extract_fluxes_to_evaluate(self):
        """ Pull out the flux timeseries for the requested event definitions. 
            
            If the input flux is differential, then estimate integral fluxes
            for integral channels specified in event definitions by calling
            estimate_integral_fluxes().
            
            Store the fluxes to be analyzed in:
            self.evaluated_energy_bins
            self.evaluated_fluxes
            self.evaluated_dates
        
        """
        
        #First, reset
        self.evaluated_energy_bins = []
        self.evaluated_dates = self.dates
        self.evaluated_fluxes = []
        to_remove = [] #event definitions with energy bins not in the data
        
        if self.flux_type == 'differential':
            self.estimate_integral_fluxes()
 
        for evdef in self.event_definitions:
            bin = [evdef['energy_channel'].min, evdef['energy_channel'].max]

            #Users may apply multiple thresholds to the same energy channel.
            #Only need to extract one copy of the flux and energy bins in
            #that channel, so skip if already included.
            if bin in self.evaluated_energy_bins:
                continue
            
            try:
                idx = self.energy_bins.index(bin)
            except:
                print("extract_fluxes_to_evaluate: Energy bin for requested event "
                    f"definition is not present in the data, {bin}. Skipping.")
                #Remove from event definitions
                to_remove.append(evdef)
                continue
            
            self.evaluated_energy_bins.append(bin)
            self.evaluated_fluxes.append(self.fluxes[idx])
 
        #Clean up by removing any event definitions that weren't in the data
        if len(to_remove) > 0:
            for evdef in to_remove:
                if evdef in self.event_definitions:
                    self.event_definitions.remove(evdef)
 
 
        return


    def add_results(self, analyze):
        """ Append Analyze object to Data object """
        self.results.append(analyze)
        return


    def get_sep_date(self):
        """ Get the year, month, day of SEP event if one is
            recorded for any of the event definitions.
            
        """
        sep_date = pd.NaT
        for analyze in self.results:
            if not pd.isnull(analyze.sep_start_time):
                self.sep_year = analyze.sep_start_time.year
                self.sep_month = analyze.sep_start_time.month
                self.sep_day = analyze.sep_start_time.day
                return
                
        return
        



################################################################
##### Analyze Class: Analysis of Flux Data for OpSEP ###########
################################################################
class Analyze:
    def __init__(self, data, event_definition):
        """ The Analyze class contains the functions that work on the
            fluxes in the Data class. These functions calculate the 
            SEP event characteristics produced by OpSEP.
            
            There is one Analyze class object per event definition.
            
            The values calculated in multiple Analyze objects (i.e. event
            definitions) can be combined with the source Data object
            to fill an Observation or Forecast object, which can be used
            directly by SPHINX, and produce a JSON file in the 
            CCMC SEP Scoreboard format.
            
            INPUT:

                :data: (object) a filled Data object
                :event_definition: (dict) dict of EnergyBin and Threshold objects
                    {'energy_channel': EnergyBin, 'threshold': Threshold}
            
        """
        self.event_definition = event_definition
        
        
        #Specific dates and fluxes for this event definition
        #Flux in a single energy channel
        self.dates = []
        self.flux = []
        
        #Derived values
        self.sep_start_time = pd.NaT
        self.sep_end_time = pd.NaT
        self.onset_peak = np.nan
        self.onset_peak_time = pd.NaT
        self.onset_rise_time = np.nan
        self.max_flux = np.nan
        self.max_flux_time = pd.NaT
        self.max_flux_rise_time = np.nan
        self.duration = np.nan
        self.fluence = np.nan
        self.fluence_spectrum = []

        self.energy_units = event_definition['energy_channel'].units
        self.flux_units = None
        self.rise_time_units = 'minutes'
        self.duration_units = 'hours'
        self.fluence_units = None
        self.fluence_spectrum_units = None

        self.sep_profile = None #name of output file continaining
            #datetime column and flux column for this event definition

        self.isgood = self.check_event_definition(data)

 
    def check_event_definition(self, data):
        """ Check if the energy channels in the data correspond to  
            the requested event definition. Exit if not.
            
        """
        energy_bin = self.make_energy_bin()
        try:
            data.evaluated_energy_bins.index(energy_bin)
        except:
            print(f"Analyze init: Requested energy bin in the event definition {energy_bin} "
                f"is not present in the data: {data.evaluated_energy_bins}. Exiting.")
            return False
        
        print(f"Analyze init: Applying event definition: "
            f"[{self.event_definition['energy_channel'].min}, "
            f"{self.event_definition['energy_channel'].max}] exceeds "
            f"{self.event_definition['threshold'].threshold} {self.event_definition['threshold'].threshold_units}")


        if self.event_definition['threshold'].threshold == -1:
            if not data.OPSEPEnhancement and not data.IDSEPEnhancement \
            and not data.doBGSubOPSEP and not data.doBGSubIDSEP:
                sys.exit("You are requesting to identify enhancement "
                    "above background, but background and SEP separation "
                    "has not been performed. Choose to do background subtraction "
                    "--OPSEPSubtractBG, variable=doBGSubOPSEP or --IDSEPSubtractBG, "
                    "variable=doBGSubIDSEP) "
                    "or set flags to identify enhancements above background "
                    "(--OPSEPEnhancement, variable=OPSEPEnhancement or "
                    "--IDSEPEnhancement, variable=IDSEPEnhancement)")
                    
            print("-1 threshold indicates you have selected to identify "
                "enhancement above background. Setting threshold to "
                f"opsep_min_threshold in fetchsep.cfg: {cfg.opsep_min_threshold}.")
            self.event_definition['threshold'].threshold = cfg.opsep_min_threshold
            
        return True

 
 
    def make_energy_channel_dict(self):
        """ {'min': 10, 'max': -1, 'units': 'MeV'} from event_definition 
            
            This form of the energy channel is used in SPHINX and the 
            Observation and Forecast objects.
        
        """
        energy_channel = {'min':self.event_definition['energy_channel'].min,
                        'max':self.event_definition['energy_channel'].max,
                        'units': self.event_definition['energy_channel'].units}
        return energy_channel



    def make_threshold_dict(self):
        """ {'threshold': 10, 'threshold_units': 'pfu'} from event_definition """
        threshold = {'threshold': self.event_definition['threshold'].threshold,
                    'threshold_units': self.event_definition['threshold'].threshold_units}
        return threshold



    def make_energy_bin(self):
        """ [Emin, Emax] from event_definition  """
        energy_bin = [self.event_definition['energy_channel'].min,
                    self.event_definition['energy_channel'].max]
        return energy_bin


    def select_fluxes(self, data, event_definition):
        """ Pull out the specific fluxes in the energy 
            channel being analyzed. 
        
        """
        energy_bin = [event_definition['energy_channel'].min,
                        event_definition['energy_channel'].max]
        idx = data.evaluated_energy_bins.index(energy_bin)
        fluxes = data.evaluated_fluxes[idx]
        dates = data.evaluated_dates
        
        self.dates = dates
        self.flux = fluxes

        #check this energy bin to determine units
        if energy_bin[1] == -1:
            self.flux_units = cfg.flux_units_integral
            self.fluence_units = cfg.fluence_units_integral
        else:
            self.flux_units = cfg.flux_units_differential
            self.fluence_units = cfg.fluence_units_differential

        #Check original input data to determine fluence spectrum units
        if data.flux_type == "integral":
            self.fluence_spectrum_units = cfg.fluence_units_integral
        if data.flux_type == "differential":
            self.fluence_spectrum_units = cfg.fluence_units_differential

        return


#FOR TESTING
#    def plot_sep_separation(self, data, dates, original_fluxes, sepfluxes,
#        bgfluxes, energy_bins, means, sigmas, showplot=False):
#        """ Make plots of background-subtracted fluxes """
#
#        plt_tools.opsep_plot_bgfluxes(f"Total_{data.experiment}", data.flux_type,
#                data.options, data.user_name, original_fluxes, dates, energy_bins,
#                means, sigmas, data.saveplot, spacecraft = data.spacecraft)
#        plt_tools.opsep_plot_bgfluxes(f"BackgroundFluxIDSEP_{data.experiment}",
#                data.flux_type, data.options, data.user_name, bgfluxes, dates,
#                energy_bins, means, sigmas, data.saveplot, spacecraft = data.spacecraft)
#        plt_tools.opsep_plot_bgfluxes(f"BGSubSEPFluxIDSEP_{data.experiment}",
#                data.flux_type, data.options, data.user_name, sepfluxes, dates,
#                energy_bins, means, sigmas, data.saveplot,
#                spacecraft = data.spacecraft)
#        
#        if showplot:
#            plt.show()
#        
#        return
#
#
#    def get_flux_above_idsep_threshold(self, data):
#        """ If threshold == -1, this indicates to apply
#            the IDSEP mean + nsigma. n is defined here and
#            optimized to identify SEP onsets above background.
#            The nsigma used here may not be the same that is 
#            used for background and SEP separation in IDSEP or
#            background subtraction step in OpSEP.
#            
#        """
#        thresh_nsigma = cfg.opsep_nsigma #allow for flexibility
#
#        fluxes = self.flux
#        dates = self.dates
#        energy_bin = self.make_energy_bin()
#        idx = data.energy_bins.index(energy_bin) #Mean and sigma index
#
#        means, sigmas, thresholds = data.read_idsep_files()
#        means = means[idx]
#        sigmas = sigmas[idx]
#        zeroes = [0]*len(means)
#
#        #If IDSEP mean background and sigma have already been
#        #applied to separate fluxes from background, then
#        #self.flux is already background subtracted.
#        #Want to do a background, SEP separation by applying
#        #thresh_nsigma * sigmas
#        if data.doBGSubIDSEP:
#            bgfluxes, sepfluxes = bgsub.separate_sep_and_background_idsep([fluxes],
#                            [zeroes], [sigmas], nsigma=thresh_nsigma,
#                            doBGSub=False)
#            if data.showplot or data.saveplot:
#                self.plot_sep_separation(data, dates, [fluxes], sepfluxes,
#                    bgfluxes, [energy_bin], [zeroes], [sigmas])
#
#        else:
#            bgfluxes, sepfluxes = bgsub.separate_sep_and_background_idsep([fluxes],
#                            [means], [sigmas], nsigma=thresh_nsigma,
#                            doBGSub=False)
#            if data.showplot or data.saveplot:
#                self.plot_sep_separation(data, dates, [fluxes], sepfluxes,
#                    bgfluxes, [energy_bin], [means], [sigmas])
#
#        #all background fluxes have been set to zero and non-zero fluxes
#        #remain only over the idsep threshold values with time
#        self.flux = sepfluxes[0]
#
#        return



    def calculate_threshold_crossing(self, data, event_definition):
        """ Calculate the threshold crossing times for a given energy bin
            and flux threshold. 

            An SEP event is considered to start if 3 consecutive points (depending on 
            dataset time resolution) are above threshold. The start time is set to the 
            first point that crossed threshold.
            
            If IDSEPEnhancement, OPSEPEnhancement or doBGSubIDSEP, doBGSubOPSEP are 
            selected, then require 5 consecutive points to attempt to avoid false 
            identifications of background time periods.
            
            Start time will be calculated with respect to the threshold. The end
            time can be calculated for a different threshold by applying a factor
            specified in the fetchsep.cfg file called endfac. The end threshold 
            will be set to cfg.endfac*threshold. 
            
            The flexibility of adding a factor for end time was inspired by 
            a SRAG alarm code that uses 0.85*threshold as the ending condition
            to avoid fluxes bouncing above and below the operational threshold.
            
            The end time is specified as the last point above threshold, even if
            there have been some drops below threshold. This is done by applying
            a dwell time, also specified in fetchsep.cfg. The dwell time is set 
            to best match NOAA SWPC's published end times.
        
            INPUT:
            
                :data: (Data object) contains flux information
                :event_definition: (dict) dict of EnergyBin and Threshold objects
                
            OUTPUT:
            
                :tc_start_time: (datetime) threshold crossing start time
        
        """
        energy_bin = [event_definition['energy_channel'].min,
                        event_definition['energy_channel'].max]
        threshold = event_definition['threshold'].threshold
        
        npoints = 3 #require 3 points above threshold as employed by SWPC
        if (data.IDSEPEnhancement or data.OPSEPEnhancement \
        or data.doBGSubIDSEP or data.doBGSubOPSEP)\
        and threshold == cfg.opsep_min_threshold:
            npoints = 8
        if data.time_resolution/60. > 15:
            npoints = 1 #time resolution >15 mins, require one point above threshold
            #If identifying enhancement above background, use more points because
            #all points above mean+3sigma will be present
            if (data.IDSEPEnhancement or data.OPSEPEnhancement \
            or data.doBGSubIDSEP or data.doBGSubOPSEP) \
            and threshold == cfg.opsep_min_threshold:
                npoints = 3
            
        if energy_bin not in data.evaluated_energy_bins:
            print(f"calculated_threshold_crossing: Requested energy bin {energy_bin} not "
                "specified in event definitions or not present in data. Skipping.")
            return

        fluxes = self.flux
        dates = self.dates
        
        threshold_crossed = False
        event_ended = False
        ndates = len(dates)
        sep_start_time = pd.NaT
        sep_end_time = pd.NaT
        
        end_threshold = cfg.endfac*threshold
                #endfac = 1.0 to get SWPC definition of event end for 5 min data
                #endfac = 0.85 used by SRAG operators in an alarm code
                #Included for flexibility, but 1.0 is used for operational values

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
                    #If identifying the end of an enhancement above background,
                    #the background fluxes are set to zero. The remaining fluxes
                    #are 3sigma above the mean background. The dwell time may not
                    #appropriately identify the event end in this case, so end the
                    #event after npoints are below the background threshold.
                    if( data.IDSEPEnhancement or data.OPSEPEnhancement \
                    or data.doBGSubIDSEP or data.doBGSubOPSEP) \
                    and threshold == cfg.opsep_min_threshold:
                        if end_counter > npoints:
                            event_ended = True
                            sep_end_time = dates[i-(end_counter-1)]
                    #When looking for an end of an event below an operational threshold
                    #or threshold that isn't the background, apply a dwell time to ensure
                    #the fluxes don't fluctuate above threshold again after a little while.
                    elif elapse > cfg.dwell_time: #N consecutive points longer than dwell time
                        event_ended = True
                        sep_end_time = dates[i-(end_counter-1)] #correct back time steps
                        #Double checked some calculated event end times with SWPC and
                        #this logic gave the correct end times. 2023-04-10 KW


        #In case that date range ended before fell before threshold,
        #use the last time in the file
        if not pd.isnull(sep_start_time) and pd.isnull(sep_end_time):
            sep_end_time = dates[ndates-1]
            print(f"WARNING !!!!File ended before SEP event ended for [{energy_bin[0]},{energy_bin[1]}], "
                 f"{threshold} {event_definition['threshold'].threshold_units}! "
                "Using the last time in the date range as the event end time. "
                "Extend your date range to get an improved estimate of the event "
                "end time and duration.")
        
        return sep_start_time, sep_end_time


    def trim_to_date_range(self, startdate, enddate, dates, array):
        """ Trim array to between startdate and enddate.
            dates corresponds to the time steps in array.
            
        """
        
        indices = [i for i in range(len(dates)) if (dates[i] >= startdate and dates[i] <= enddate)]
        
        if len(indices) == 0:
            return []
        
        nst = indices[0]
        nend = indices[-1] + 1
        
        return array[nst:nend]



    def calculate_max_flux(self, data):
        """ Identify maximum flux value and time during SEP event. """
 
        max_flux = np.nan
        max_flux_time = pd.NaT
        
        energy_bin = self.make_energy_bin()
 
        if energy_bin not in data.evaluated_energy_bins:
            print(f"calculated_threshold_crossing: Requested energy bin {energy_bin} not "
                "specified in event definitions or not present in data. Skipping.")
            return

        fluxes = self.flux
        dates = self.dates
        
        #If there is a SEP EVENT, extract between start and end times.
        #If NO EVENT, don't trim and take the maximum of the full timeseries
        if not pd.isnull(self.sep_start_time) and not pd.isnull(self.sep_end_time):
            fluxes = self.trim_to_date_range(self.sep_start_time, self.sep_end_time,
                                    dates, fluxes)
            dates = self.trim_to_date_range(self.sep_start_time, self.sep_end_time,
                                    dates, dates)

        max_flux = np.nanmax(fluxes)
        ix = np.nanargmax(fluxes) #First instance, if multiple
        try:
            max_flux_time = dates[ix]
        except:
            pass
        
        self.max_flux = max_flux
        self.max_flux_time = max_flux_time
        
        return max_flux, max_flux_time


    def calculate_onset_peak_from_fit(self, data):
        """Calculate the peak associated with the initial SEP onset. This subroutine
            searches for the rollover that typically occurs after the SEP onset.
            The peak value will be specified as the flux value at the rollover
            location.
            
            The onset peak may provide a more physically appropriate comparison
            with models.
            
            If code cannot identify onset peak, it will return values:
                onset_peak = np.nan
                onset_peak_time = pd.NaT
            
            The onset peak is found by fitting a Weibull function to the SEP
            time profile. The fit is performed using fluxes up to 6 hours prior to
            the threshold crossing up to 24 hours after the threshold crossing.
            
            A peak time is estimated using the second derivative
            of the fitted Weibull. The maximum measured flux value within 1 hour
            of the estimated onset peak location is taken to be the measured
            onset peak.
            
            A least 6 hours of flux measurements are required to fit the Weibull.
            If the duration of the time profile is shorter, then the onset peak
            will not be calcualated.
            
            INPUT:
            
                :event_definition: (dict) dict of EnergyBin and Threshold objects
                :sep_start_time: (datetime) SEP event start time
                :sep_end_time: (datetime) SEP event end time
                :max_flux: (float) maximum flux during SEP event, used to 
                    check for magnetitude of event to guide fitting
            
            OUTPUT:
            
                :onset_peak: (float) value of the onset peak
                :onset_peak_time: (datetime) time of onset peak
            
        """
        energy_bin = self.make_energy_bin()
        energy_units = self.event_definition['energy_channel'].units
        threshold = self.event_definition['threshold'].threshold
        threshold_units = self.event_definition['threshold'].threshold_units

        onset_peak = np.nan
        onset_peak_time = pd.NaT

        if energy_bin not in data.evaluated_energy_bins:
            print(f"calculated_onset_peak_from_fit: Requested energy bin {energy_bin} not "
                "specified in event definitions or not present in data. Skipping.")
            self.onset_peak = onset_peak
            self.onset_peak_time = onset_peak_time
            return onset_peak, onset_peak_time

        #NO SEP EVENT
        if pd.isnull(self.sep_start_time) or pd.isnull(self.sep_end_time):
            self.onset_peak = onset_peak
            self.onset_peak_time = onset_peak_time
            return onset_peak, onset_peak_time

        fluxes = self.flux
        dates = self.dates

        #Do a fit of the Weibull function for each time profile
        params_fit = Parameters()
        params_fit.add('alpha', value = -3, min = -20, max = -0.1) #-3, -5, -0.1
        params_fit.add('beta', value = 10, min = 1, max =500) #10, 1, 100
        params_fit.add('peak_intensity', value = 100, min = 1e-3, max =1e6) #100, 1e-3, 1e6

        flxratio = self.max_flux/threshold
        
        #For mid to strong SEP events, time according to prompt onset and the
        #possibility of a CME arriving around 24 hours later (set timing
        #to approximately exclude CME arrival)
        start_fit_time = max(dates[0], self.sep_start_time - datetime.timedelta(hours=3))
        end_fit_time = min(self.sep_end_time, self.sep_start_time + datetime.timedelta(hours=24))

        
        #For lower intensity events, extend back to initial elevation to
        #capture more of the event and also capture more time
        #after threshold crossing as likely a slower CME or CME that
        #will not arrive at Earth (for very gradual rises)
        if flxratio < 10:
            #Try to find the background level and evaluate the event
            #from the time it first deviates from background
            ratio = [0.1, 0.2, 0.3, 0.4, 0.5]
            for rat in ratio:
                low_thresh = rat*threshold
                low_evdef = data.create_event_definition(energy_bin[0],energy_bin[1],
                    energy_units, low_thresh, threshold_units)
                low_start_time, low_end_time = self.calculate_threshold_crossing(data, low_evdef)
                if low_start_time > dates[0] and low_start_time < start_fit_time:
                    start_fit_time = low_start_time
                    break

        print(f"calculate_onset_peak_from_fit: FITTING BETWEEN {start_fit_time} to {end_fit_time}")
        
        trim_fluxes = self.trim_to_date_range(start_fit_time, end_fit_time,
                                    dates, fluxes)
        trim_dates = self.trim_to_date_range(start_fit_time, end_fit_time,
                                    dates, dates)
 
        #Convert dates into a series of times in hours for fitting
        trim_times = [((t - dates[0]).total_seconds() + 60)/(60*60) for t in trim_dates]
 
        minimize_func = minimize(tools.func_residual, params_fit,
                    args = [trim_times, trim_fluxes],
                    nan_policy= 'propagate', max_nfev=np.inf)
                
        #Get Weibull fit parameters
        best_pars = minimize_func.params.valuesdict()
        best_a = best_pars['alpha']
        best_b = best_pars['beta']
        best_Ip = best_pars['peak_intensity']
        best_fit = tools.modified_weibull(trim_times, best_Ip, best_a, best_b)
        err = tools.ratio(best_fit, trim_fluxes)
    
        print(f"calculate_onset_peak_from_fit ==== {energy_bin} MeV =====")
        print(f"Best fit Weibull for onset peak Ip: {best_Ip}, a: {best_a}, b: {best_b}")
        print(f"Error in fit: {err}")

        if pd.isnull(best_Ip) or pd.isnull(best_a) or pd.isnull(best_b):
            print("calculate_onset_peak_from_fit: Fit failed for "
                f"{energy_bin}, {threshold}. Returning null values.")
            self.onset_peak = onset_peak
            self.onset_peak_time = onset_peak_time
            return onset_peak, onset_peak_time

        #IF WEIBULL FIT SUCCESSFUL
        #Find maximum curvature in the fit using the second derivative
        max_curve_idx = tools.find_max_curvature(trim_times, best_fit)

        max_curve_model_time = trim_times[max_curve_idx]
        max_curve_model_date = trim_dates[max_curve_idx]
        max_curve_model_peak = best_fit[max_curve_idx] #fit value
        
        #Get the maximum measured value near the location of the fit peak
        #Pull out max measured value around the maximum of curvature
        #Search +- 1 hour from the fit max time
        max_curve_meas_peak = 0
        max_curve_meas_time = pd.NaT
        max_curve_meas_date = pd.NaT
        dt = datetime.timedelta(hours=1)
        for k in range(len(trim_dates)):
            if trim_dates[k] >= max_curve_model_date - dt \
                and trim_dates[k] <= max_curve_model_date + dt:
                if trim_fluxes[k] > max_curve_meas_peak:
                    max_curve_meas_peak = trim_fluxes[k]
                    max_curve_meas_time = trim_times[k]
                    max_curve_meas_date = trim_dates[k]

        onset_peak = max_curve_meas_peak
        onset_peak_time = max_curve_meas_date

        ### VALUES ONLY USED IN PLOTS
        ####FIND PEAK BY JUST TAKING MAXIMUM OF WEIBULL
        max_val = np.max(best_fit)
        max_idx = np.where(best_fit == max_val)
        max_time = trim_times[max_idx[0][0]]
        
        #Pull out max measured value around this identified maximum in the fit
        model_max_date = datetime.timedelta(seconds=(max_time*60*60 - 60)) + dates[0]
        max_meas = 0
        max_meas_time = 0
        max_date = 0
        dt = datetime.timedelta(hours=1)
        for k in range(len(trim_dates)):
            if trim_dates[k] >= model_max_date - dt and trim_dates[k] <= model_max_date + dt:
                if trim_fluxes[k] > max_meas:
                    max_meas = trim_fluxes[k]
                    max_meas_time = trim_times[k]
                    max_date = trim_dates[k]

        
        #PLOT
        if data.saveplot or data.showplot:
            plt_tools.plot_weibull_fit(energy_bin, threshold, data.experiment,
                data.flux_type, data.user_name, data.options,
                self.sep_start_time, trim_times, trim_fluxes, best_pars, best_fit, max_time,
                max_val, max_meas_time, max_meas, max_curve_model_time, max_curve_model_peak,
                max_curve_meas_time, max_curve_meas_peak,
                data.saveplot, data.showplot, spacecraft=data.spacecraft,
                doBGSubOPSEP=data.doBGSubOPSEP, doBGSubIDSEP=data.doBGSubIDSEP,
                OPSEPEnhancement=data.OPSEPEnhancement,
                IDSEPEnhancement=data.IDSEPEnhancement)

        self.onset_peak = onset_peak
        self.onset_peak_time = onset_peak_time
        
        return onset_peak, onset_peak_time


    def derived_timing_values(self):
        """ Calculate duration and rise time in minutes. """
        
        onset_rise_time = pd.NaT
        max_rise_time = pd.NaT
        duration = pd.NaT
        
        if not pd.isnull(self.sep_start_time):
            if not pd.isnull(self.onset_peak_time):
                onset_rise_time = (self.onset_peak_time - self.sep_start_time).total_seconds()/60.
            if not pd.isnull(self.max_flux_time):
                max_rise_time = (self.max_flux_time - self.sep_start_time).total_seconds()/60.
            if not pd.isnull(self.sep_end_time):
                duration = (self.sep_end_time - self.sep_start_time).total_seconds()/(60.*60.)

        self.onset_rise_time = onset_rise_time #minutes
        self.max_rise_time = max_rise_time #minutes
        self.duration = duration #hours
        
        return onset_rise_time, max_rise_time, duration


    def calculate_fluence(self,fluxes, time_resolution):
        """ Calculate fluence for one energy bin of fluxes already trimmed
            to the time period of interest. Sum all the fluxes in the 
            array in time and return a single fluence value.
            
        """
        clean_flux = [fx for fx in fluxes if not pd.isnull(fx) and fx >=0]
        
        if len(clean_flux) == 0:
            return np.nan
        
        fluence = sum(clean_flux)*time_resolution*4.0*math.pi #multiply 4pi steradians
    
        return fluence
 
 
    def calculate_channel_fluence(self, data):
        """  Calculate the fluence for the specified event definition
        """
 
        flux = self.flux
        dates = self.dates

        if pd.isnull(self.sep_start_time) or pd.isnull(self.sep_end_time):
            print("calculated_channel_fluence: The SEP start or end time is null. No event. Returning NaN for fluence.")
            return np.nan

        trim_flux = self.trim_to_date_range(self.sep_start_time, self.sep_end_time, dates, flux)
        fluence = self.calculate_fluence(trim_flux, data.time_resolution)
        self.fluence = fluence
        
        return fluence



    def calculate_fluence_spectrum(self, data):
        """ Calculate the spectrum between the start and end times 
            for a specific event definition.
            
            The spectrum is calculated by summing the flux in all
            of the ORIGINAL energy bins between the start and end
            times. If the input flux was in differential energy
            channels, then the spectrum will be in differential 
            channels, even if the, e.g., >10 meV, 10 pfu event
            definition is used. If the input fluxes were integral,
            then the energy spectrum will be in integral flux bins.
            
        """
        
        fluxes = data.fluxes #all fluxes for all energy bins in input data
        dates = data.dates
        
        #Trim to the SEP start and end times
        sep_dates, sep_fluxes = datasets.extract_date_range(self.sep_start_time,
                                self.sep_end_time, dates, fluxes)

        fluence_spectrum = []
        for flux in sep_fluxes:
            fluence = self.calculate_fluence(flux, data.time_resolution)
            fluence_spectrum.append(fluence)
            
        self.fluence_spectrum = fluence_spectrum
        
        return fluence_spectrum



    def calculate_event_info(self, data):
        """ Calculate SEP event characteristics for a single event
            definition. Calculate from the fluxes and energy bins 
            that have been extracted for evaluation for the specified
            event definitions.
            
            INPUT:
            
                :event_definition: (dict) dict of EnergyBin and Threshold objects
                
            OUTPUT:
            
                :event_info: (dict) SEP event values
                
        """
        #calculate event values and fill in a dictionary that will
        #save info needed for Observation or Forecast objects
        self.select_fluxes(data, self.event_definition) #Load fluxes to obj
        sep_start_time, sep_end_time = self.calculate_threshold_crossing(data, self.event_definition)
        self.sep_start_time = sep_start_time
        self.sep_end_time = sep_end_time
        self.calculate_max_flux(data)
        self.calculate_onset_peak_from_fit(data)
        self.derived_timing_values()

        #If the onset peak time is AFTER the max flux time, set the onset peak
        #to the max flux value and time.
        if self.onset_peak_time > self.max_flux_time:
            self.onset_peak = self.max_flux
            self.onset_peak_time = self.max_flux_time

        #Fluence
        self.calculate_channel_fluence(data)
        self.calculate_fluence_spectrum(data)

        energy_bin = self.make_energy_bin()
        energy_units = self.event_definition['energy_channel'].units
        threshold = self.event_definition['threshold'].threshold
        threshold_units = self.event_definition['threshold'].threshold_units

        energy_label = f"{energy_bin[0]} - {energy_bin[1]} {energy_units}"
        if energy_bin[1] == -1:
            energy_label = f">{energy_bin[0]} {energy_units}"
        threshold_label = f"{threshold} {threshold_units}"
        print()
        print(f"====SEP Event Characteristics for {energy_label}, {threshold_label}====")
        print(f"SEP Start Time: {self.sep_start_time}")
        print(f"SEP End Time: {self.sep_end_time}")
        print(f"Onset Peak: {self.onset_peak} {self.flux_units} at {self.onset_peak_time}")
        print(f"Max Flux: {self.max_flux} {self.flux_units} at {self.max_flux_time}")
        print(f"Rise time to Onset: {self.onset_rise_time} {self.rise_time_units}")
        print(f"Rise time to Max: {self.max_rise_time} {self.rise_time_units}")
        print(f"Duration: {self.duration} {self.duration_units}")
        print(f"Channel Fluence: {self.fluence} {self.fluence_units}")
        print(f"Fluence Spectrum: {self.fluence_spectrum} {self.fluence_spectrum_units}")
        print(f"Fluence Energy Bins: {data.energy_bins}")
        print()
        
        return
        


class Output:
    def __init__(self, data, json_type, spase_id=None):
        """Output the data generated by OpSEP and write out the data
            files. Multiple Analyze objects with individual event 
            definitions are combined to create summary text files,
            plots, CCMC SEP Scoreboard and SPHINX jsons, and SPHINX
            Observation or Forecast objects.
           
           INPUT:
                
                :data: (Data object) Data object loaded up with 
                    Analyze objects
                
            OUTPUT:
            
                various output files including CCMC JSONs and 
                Forecast or Observation objects used by SPHINX.
        
        """

        self.data = data
        
        #Check if any Analysis objects were created
        if len(data.results) == 0:
            sys.exit("Output init: No event definitions were applied to the data. "
                "There are no results to report. Available energy bins are "
                f"{self.data.energy_bins} Exiting.")
    
        self.spase_id = spase_id
        self.json_type = json_type #observation or forecast
        self.mode = None #modes for CCMC json when forecasts: allowed values: forecast, historical, nowcast,
                         #simulated_realtime_forecast, simulated_realtime_nowcast
        self.json_dict = {} #json dictionary from template
        self.json_filename = None #output path and filename
        self.issue_time = pd.NaT
        
        # Subdirectory with unique string to hold data
        self.subdir = tools.opsep_subdir(self.data.experiment, self.data.flux_type,
            self.data.user_name, self.data.options, spacecraft=self.data.spacecraft,
            doBGSubOPSEP=self.data.doBGSubOPSEP, doBGSubIDSEP=self.data.doBGSubIDSEP,
            OPSEPEnhancement=self.data.OPSEPEnhancement,
            IDSEPEnhancement=self.data.IDSEPEnhancement)

        if not os.path.exists(os.path.join(cfg.outpath,'opsep', self.subdir)):
            os.mkdir(os.path.join(cfg.outpath,'opsep', self.subdir))
        if not os.path.exists(os.path.join(cfg.plotpath,'opsep', self.subdir)):
            os.mkdir(os.path.join(cfg.plotpath,'opsep', subdir))


    def set_json_type(self, json_type):
        self.json_type = json_type
        return
        
    
    def set_spase_id(self, spase_id):
        self.spase_id = spase_id
        return


    def set_json_filename(self):
        """ Filename in CCMC SEP Scoreboard format. 
            Set filenames and set issue time.
            
        """

        modifier, title_mod = tools.setup_modifiers(self.data.options, spacecraft=self.data.spacecraft, doBGSubOPSEP=self.data.doBGSubOPSEP,
            doBGSubIDSEP=self.data.doBGSubIDSEP, OPSEPEnhancement=self.data.OPSEPEnhancement,
            IDSEPEnhancement=self.data.IDSEPEnhancement)

        #Get issue time of forecast (now)
        now = datetime.datetime.now()
        self.issue_time = now
        
        issue_time = ccmc_json.make_ccmc_zulu_time(now)
        issue_time = issue_time.replace(":","")
        zstdate = ccmc_json.make_ccmc_zulu_time(self.data.startdate)
        zstdate = zstdate.replace(":","")

        #Filenames for observations don't include issue time
        fnameprefix = ""
        if self.json_type == "observations":
            fnameprefix = f"{self.data.experiment}_{self.data.flux_type}{modifier}.{zstdate}"
            if not pd.isnull(self.data.user_name) and self.data.user_name != "":
                fnameprefix = f"{self.data.user_name}_{self.data.flux_type}{modifier}.{zstdate}"

        #Filenames for model output do include issue time
        if self.json_type == "model":
            fnameprefix = f"{self.data.experiment}_{self.data.flux_type}{modifier}.{zstdate}.{issue_time}"
            if not pd.isnull(self.data.user_name) and self.data.user_name != "":
                fnameprefix = f"{self.data.user_name}_{self.data.flux_type}{modifier}.{zstdate}.{issue_time}"

        ####JSON FILE
        self.json_filename = fnameprefix + ".json"
    
        return

    
    def set_sep_profile_filename(self, analyze):

        fnameprefix = self.json_filename.strip().split(".json")
        fnameprefix = fnameprefix[0]
        
        ####TIME PROFILE
        energy_bin = analyze.make_energy_bin()
        if energy_bin[1] == -1: #integral
           profname = f"{fnameprefix}.{energy_bin[0]}.{analyze.energy_units}.txt"
        else:
           profname = f"{fnameprefix}.{energy_bin[0]}-{energy_bin[1]}.{analyze.energy_units}.txt"
        analyze.sep_profile = profname
        
        return analyze


    def write_zulu_time_profile(self, analyze):
        """ Write out the time profile with the date in the
            first column as the ISO standard and flux in the
            second column as:
            
            YYYY-MM-DDTHH:MM:SSZ    Float
            
            INPUTS:
            
            :Filename: (string) - name of file to write
            :date: (datetime 1xn array) - list of dates
            :fluxes: (float 1xn array) - corresponding fluxes
            
            OUTPUTS:
            
            None but writes output file with filename
            
        """

        fname = os.path.join(cfg.outpath,'opsep', self.subdir, analyze.sep_profile)
        outfile = open(fname, "w")
        for i in range(len(analyze.dates)):
            zdate = ccmc_json.make_ccmc_zulu_time(analyze.dates[i])
            outfile.write(zdate + "    " + str(analyze.flux[i]) + "\n")
            
        outfile.close()

        print("write_zulu_time_profile: Wrote file --> " + fname)

        return


    def fill_event_info_dict(self, analyze):
        """ Initialize dictionary that contains SEP event info
            saved in a single Analyze object combined with Data. 
            
            This dictionary contains the derived values and 
            supporting contextual information for a single
            event definition.
            
        """
        
        dict = {'experiment': self.data.experiment, #Experiment or model name; GOES-13
                'flux_type': self.data.flux_type, #ORIGINAL input data - integral or differential
                'startdate': self.data.startdate, #Start of analyzed time period
                'enddate': self.data.enddate, #End of analyzed time period
                'background_subtraction': self.data.doBGSubOPSEP, #bool doBGSubOPSEP
                'options': self.data.options, #options applied to data
                'original_energy_bins': self.data.energy_bins, #All original energy bins for input data
               # 'event_definition': None, #Dictionary of Energy Channel and Threshold obj
                'energy_channel': analyze.make_energy_channel_dict(), #{'min': 10, 'max': -1, 'units': 'MeV'}
                'energy_bin': analyze.make_energy_bin(), #[Emin, Emax]
                'threshold_dict': analyze.make_threshold_dict(), #{'threshold': 10, 'threshold_units': 'pfu'}
                'threshold': analyze.event_definition['threshold'].threshold, #float
                'sep_start_time': analyze.sep_start_time, #SEP start time
                'sep_end_time': analyze.sep_end_time, #SEP end time
                'onset_peak': analyze.onset_peak, #Onset peak
                'onset_peak_time': analyze.onset_peak_time, #Time of onset peak
                'onset_rise_time': analyze.onset_rise_time, #Time from sep_start_time to sep_onset_peak_time
                'max_flux': analyze.max_flux, #Maximum flux during SEP event
                'max_flux_time': analyze.max_flux_time, #Time of maximum flux
                'max_rise_time': analyze.max_rise_time, #Time from sep_start_time to sep_max_flux_time
                'duration': analyze.duration, #sep_end_time - sep_start_time
                'fluence': analyze.fluence, #fluence in single energy channel summed between sep_start_time and sep_end_time
                'fluence_spectrum': analyze.fluence_spectrum, #fluence in all_energy_bins summed between sep_start_time and sep_end_time
                'flux_units': analyze.flux_units, #str
                'fluence_units': analyze.fluence_units, #str
                'fluence_spectrum_units': analyze.fluence_spectrum_units, #str
                'rise_time_units': analyze.rise_time_units, #str
                'duration_units': analyze.duration_units, #str
                'sep_profile': analyze.sep_profile #str
            
        }

        return dict


    def event_info_dict_for_csv(self, analyze):
        """ Create a flat dictionary with all event info with the
            dictionary keys labeled according to energy channel
            and threshold information.
            
            Useful to ultimately export to csv.
            
        """
        energy_bin = analyze.make_energy_bin()
        energy_units = analyze.event_definition['energy_channel'].units
        threshold = analyze.event_definition['threshold'].threshold
        threshold_units = analyze.event_definition['threshold'].threshold_units

        if energy_bin[1] == -1:
            channel_label = f">{energy_bin[0]} {energy_units}"
        else:
            channel_label = f"{energy_bin[0]}-{energy_bin[1]} {energy_units}"
    
        threshold_label = f"{threshold} {threshold_units}"

        fluence_spectrum_str = str(analyze.fluence_spectrum)
        fluence_spectrum_str = fluence_spectrum_str.replace(",", ";")
        fluence_spectrum_energy_bins = str(self.data.energy_bins)
        fluence_spectrum_energy_bins = fluence_spectrum_energy_bins.replace(",", ";")
        fluence_spectrum_energy_bin_centers = str(self.data.energy_bin_centers)
        fluence_spectrum_energy_bin_centers = fluence_spectrum_energy_bin_centers.replace(",", ";")

        if not pd.isnull(analyze.sep_start_time):
            sttime = analyze.sep_start_time.strftime("%Y-%m-%d %H:%M:%S")
        else:
            sttime = None

        if not pd.isnull(analyze.sep_end_time):
            endtime = analyze.sep_end_time.strftime("%Y-%m-%d %H:%M:%S")
        else:
            endtime = None

        if not pd.isnull(analyze.onset_peak_time):
            optime = analyze.onset_peak_time.strftime("%Y-%m-%d %H:%M:%S")
        else:
            optime = None

        if not pd.isnull(analyze.max_flux_time):
            mftime = analyze.max_flux_time.strftime("%Y-%m-%d %H:%M:%S")
        else:
            mftime = None

        dict = {f"{channel_label} {threshold_label} SEP Start Time": sttime,
                f"{channel_label} {threshold_label} SEP End Time": endtime,
                f"{channel_label} {threshold_label} SEP Duration ({analyze.duration_units})": analyze.duration,
                f"{channel_label} {threshold_label} Onset Peak ({analyze.flux_units})": analyze.onset_peak,
                f"{channel_label} {threshold_label} Onset Peak Time": optime,
                f"{channel_label} {threshold_label} Rise Time to Onset ({analyze.rise_time_units})": analyze.onset_rise_time,
                f"{channel_label} {threshold_label} Max Flux ({analyze.flux_units})": analyze.max_flux,
                f"{channel_label} {threshold_label} Max Flux Time": mftime,
                f"{channel_label} {threshold_label} Rise Time to Max ({analyze.rise_time_units})": analyze.max_rise_time,
                f"{channel_label} {threshold_label} Fluence ({analyze.fluence_units})": analyze.fluence,
                f"{channel_label} {threshold_label} Fluence Spectrum ({analyze.fluence_spectrum_units})": fluence_spectrum_str,
                f"{channel_label} {threshold_label} Fluence Spectrum Energy Bins ({energy_units})": fluence_spectrum_energy_bins,
                f"{channel_label} {threshold_label} Fluence Spectrum Energy Bin Centers ({energy_units})": fluence_spectrum_energy_bin_centers
            }

        return dict


    def event_info_dict_for_pkl(self, analyze):
        """ Create a flat dictionary with all event info with the
            dictionary keys labeled according to energy channel
            and threshold information.
            
            Useful to ultimately export to pkl.
            
        """
        energy_bin = analyze.make_energy_bin()
        energy_units = analyze.event_definition['energy_channel'].units
        threshold = analyze.event_definition['threshold'].threshold
        threshold_units = analyze.event_definition['threshold'].threshold_units

        if energy_bin[1] == -1:
            channel_label = f">{energy_bin[0]} {energy_units}"
        else:
            channel_label = f"{energy_bin[0]}-{energy_bin[1]} {energy_units}"
    
        threshold_label = f"{threshold} {threshold_units}"


        dict = {f"{channel_label} {threshold_label} SEP Start Time": analyze.sep_start_time,
                f"{channel_label} {threshold_label} SEP End Time": analyze.sep_end_time,
                f"{channel_label} {threshold_label} SEP Duration ({analyze.duration_units})": analyze.duration,
                f"{channel_label} {threshold_label} Onset Peak ({analyze.flux_units})": analyze.onset_peak,
                f"{channel_label} {threshold_label} Onset Peak Time": analyze.onset_peak_time,
                f"{channel_label} {threshold_label} Rise Time to Onset ({analyze.rise_time_units})": analyze.onset_rise_time,
                f"{channel_label} {threshold_label} Max Flux ({analyze.flux_units})": analyze.max_flux,
                f"{channel_label} {threshold_label} Max Flux Time": analyze.max_flux_time,
                f"{channel_label} {threshold_label} Rise Time to Max ({analyze.rise_time_units})": analyze.max_rise_time,
                f"{channel_label} {threshold_label} Fluence ({analyze.fluence_units})": analyze.fluence,
                f"{channel_label} {threshold_label} Fluence Spectrum ({analyze.fluence_spectrum_units})": analyze.fluence_spectrum,
                f"{channel_label} {threshold_label} Fluence Spectrum Energy Bins ({energy_units})": self.data.energy_bins,
                f"{channel_label} {threshold_label} Fluence Spectrum Energy Bin Centers ({energy_units})": self.data.energy_bin_centers
            }

        return dict



    def fill_event_info_list(self, analyze):
        """ Return a list that contains SEP event info
            saved in a single Analyze object combined with Data. 
            
            This list contains the derived values and 
            supporting contextual information for a single
            event definition.
            
            The list is in the same order as the event info dict.
            
        """
        event_info_list = [self.data.experiment,
                            self.data.flux_type,
                            self.data.startdate,
                            self.data.enddate,
                            self.data.doBGSubOPSEP,
                            self.data.options,
                            self.data.energy_bins,
                            self.data.energy_bin_centers,
                            analyze.make_energy_channel_dict(),
                            analyze.make_energy_bin(),
                            analyze.make_threshold_dict(),
                            analyze.event_definition['threshold'].threshold,
                            analyze.sep_start_time,
                            analyze.sep_end_time,
                            analyze.onset_peak,
                            analyze.onset_peak_time,
                            analyze.onset_rise_time,
                            analyze.max_flux,
                            analyze.max_flux_time,
                            analyze.max_rise_time,
                            analyze.duration,
                            analyze.fluence,
                            analyze.fluence_spectrum,
                            analyze.flux_units,
                            analyze.fluence_units,
                            analyze.fluence_spectrum_units,
                            analyze.rise_time_units,
                            analyze.duration_units,
                            analyze.sep_profile
                            ]
 
        return event_info_list


    def fill_json_header(self):
        """ Fill the json dictionary with the extracted SEP values. """
        self.json_dict = ccmc_json.fill_json_header(self.json_type,
            self.issue_time, self.data.experiment, self.data.flux_type,
            self.data.options, self.spase_id,
            user_name=self.data.user_name, spacecraft=self.data.spacecraft)

        return


    def fill_json_block(self, analyze):
        """ Store SEP event values into the json block associated with
            the event definition in the Analyze object.
            
        """

        energy_channel_dict = analyze.make_energy_channel_dict()
        threshold_dict = analyze.make_threshold_dict()
        
        self.json_dict = ccmc_json.fill_json_block(self.json_dict,
                                    self.json_type,
                                    energy_channel_dict,
                                    threshold_dict,
                                    self.data.startdate,
                                    self.data.enddate,
                                    analyze.sep_start_time,
                                    analyze.sep_end_time,
                                    analyze.onset_peak,
                                    analyze.onset_peak_time,
                                    analyze.max_flux,
                                    analyze.max_flux_time,
                                    analyze.flux_units,
                                    analyze.fluence,
                                    analyze.fluence_units,
                                    analyze.fluence_spectrum,
                                    analyze.fluence_spectrum_units,
                                    self.data.energy_bins,
                                    analyze.sep_profile,
                                    self.data.location,
                                    self.data.species)
        
        return


    def clean_json(self):
        """ Remove empty fields or fields with bad values """
        
        self.json_dict = ccmc_json.clean_json(self.json_dict, self.data.experiment,
                self.json_type)
        return


    def write_json(self):
        filename = os.path.join(cfg.outpath, 'opsep', self.subdir, self.json_filename)
        is_good = ccmc_json.write_json(self.json_dict, filename)
        return filename


    def write_ccmc_json(self):
        """ Write all event definitions out to CCMC json file 
            https://ccmc.gsfc.nasa.gov/publicData/sepsb/files/sepscoreboard_visual_schema.pdf
        
        """
        self.set_json_filename()
        self.fill_json_header()
        
        #Cycle through all Analyze objects for the various event definitions
        for i, analyze in enumerate(self.data.results):
            #Set SEP profile filename
            analyze = self.set_sep_profile_filename(analyze)
            #Write out SEP profile in CCMC format and save to Analyze obj
            self.write_zulu_time_profile(analyze)
            self.fill_json_block(analyze)
            self.data.results[i] = analyze
        
        self.clean_json()
        filename = self.write_json()

        return filename


    def create_csv_dict(self):
        """ A flat dictionary of values for all event definitions. """

        exp_name = self.data.experiment
        if not pd.isnull(self.data.user_name) and self.data.user_name != "":
            exp_name = self.data.user_name

        bgsub = 'None'
        if self.data.doBGSubOPSEP: bgsub = 'OPSEP'
        if self.data.doBGSubIDSEP: bgsub = 'IDSEP'

        dict = {"Experiment": exp_name,
                "Flux Type": self.data.flux_type,
                "Options": str(self.data.options).replace(",",";"),
                "Background Subtraction": bgsub,
                "Time Period Start": self.data.startdate.strftime("%Y-%m-%d %H:%M:%S"),
                "Time Period End": self.data.enddate.strftime("%Y-%m-%d %H:%M:%S")
                }
        
        for analyze in self.data.results:
            analyze_dict = self.event_info_dict_for_csv(analyze)
            dict.update(analyze_dict)
            
        header = ''
        row = ''
        for key in dict.keys():
            header += key + ","
            row += str(dict[key]) + ","
        
        header = header[:-1] + "\n"
        row = row[:-1] + "\n"
        
        filename = self.json_filename
        filename = filename.replace(".json",".csv")
        filename = os.path.join(cfg.outpath,"opsep",self.subdir, filename)
        fout = open(filename,"w+")
        fout.write(header)
        fout.write(row)
        fout.close()
        print(f"create_csv_dict: Wrote {filename}")

        return dict


    def create_pkl_dict(self):
        """ A flat dictionary of values for all event definitions. """

        exp_name = self.data.experiment
        if not pd.isnull(self.data.user_name) and self.data.user_name != "":
            exp_name = self.data.user_name

        bgsub = 'None'
        if self.data.doBGSubOPSEP: bgsub = 'OPSEP'
        if self.data.doBGSubIDSEP: bgsub = 'IDSEP'

        dict = {"Experiment": exp_name,
                "Flux Type": self.data.flux_type,
                "Options": self.data.options,
                "Background Subtraction": bgsub,
                "Time Period Start": self.data.startdate,
                "Time Period End": self.data.enddate
                }
        
        for analyze in self.data.results:
            analyze_dict = self.event_info_dict_for_pkl(analyze)
            dict.update(analyze_dict)
            
        header = ''
        row = ''
        for key in dict.keys():
            header += key + ","
            row += str(dict[key]) + ","
        
        header = header[:-1] + "\n"
        row = row[:-1] + "\n"
        
        filename = self.json_filename
        filename = filename.replace(".json",".pkl")
        filename = os.path.join(cfg.outpath,"opsep",self.subdir, filename)
        with open(filename, 'wb') as file:
            pickle.dump(dict, file)
            print(f"create_pkl_dict: Wrote {filename}")
        
        return dict


    def extract_analyze_lists(self):
        """ Pull out the SEP start and end times, onset peaks,
            max fluxes, fluences, and fluence spectra for plotting.
        
        """
        event_definitions = []
        fluxes = []
        sep_start_times = []
        sep_end_times = []
        onset_peaks = []
        onset_peak_times = []
        max_fluxes = []
        max_flux_times = []
        fluences = []
        fluence_spectra = []
        fluence_spectra_units = []

        for analyze in self.data.results:
            event_definitions.append(analyze.event_definition)
            fluxes.append(analyze.flux)
            sep_start_times.append(analyze.sep_start_time)
            sep_end_times.append(analyze.sep_end_time)
            onset_peaks.append(analyze.onset_peak)
            onset_peak_times.append(analyze.onset_peak_time)
            max_fluxes.append(analyze.max_flux)
            max_flux_times.append(analyze.max_flux_time)
            fluences.append(analyze.fluence)
            fluence_spectra.append(analyze.fluence_spectrum)
            fluence_spectra_units.append(analyze.fluence_spectrum_units)

        return event_definitions, fluxes, sep_start_times, sep_end_times, onset_peaks,\
            onset_peak_times, max_fluxes, max_flux_times, fluences, fluence_spectra,\
            fluence_spectra_units


    def plot_event_definitions(self):
        """ Plot the fluxes used for event definitions with threshold,
            start and end times, onset peak and max flux.
        
        """
        #Collect calculated values from Analyze objects
        event_definitions, analyzed_fluxes, sep_start_times, sep_end_times,\
        onset_peaks, onset_peak_times, max_fluxes, max_flux_times, fluences,\
        fluence_spectra, fluence_spectra_units = self.extract_analyze_lists()
        
        plt_tools.opsep_plot_event_definitions(self.data.experiment,
            self.data.flux_type, self.data.user_name, self.data.options,
            self.data.evaluated_dates, analyzed_fluxes,
            self.data.evaluated_energy_bins, event_definitions,
            sep_start_times, sep_end_times, onset_peaks, onset_peak_times,
            max_fluxes, max_flux_times, self.data.showplot, self.data.saveplot,
            spacecraft=self.data.spacecraft, doBGSubOPSEP=self.data.doBGSubOPSEP,
            doBGSubIDSEP=self.data.doBGSubIDSEP, OPSEPEnhancement=self.data.OPSEPEnhancement,
            IDSEPEnhancement=self.data.IDSEPEnhancement)


    def plot_all_fluxes(self):
        """ Plot threshold crossings on top of all fluxes """
 
        #Collect calculated values from Analyze objects
        event_definitions, analyzed_fluxes, sep_start_times, sep_end_times,\
        onset_peaks, onset_peak_times, max_fluxes, max_flux_times, fluences,\
        fluence_spectra, fluence_spectra_units = self.extract_analyze_lists()

        plt_tools.opsep_plot_all_bins(self.data.experiment, self.data.flux_type,
            self.data.user_name, self.data.options,
            self.data.dates, self.data.fluxes, self.data.energy_bins,
            self.data.event_definitions, sep_start_times, sep_end_times,
            self.data.showplot, self.data.saveplot, spacecraft=self.data.spacecraft,
            doBGSubOPSEP=self.data.doBGSubOPSEP, doBGSubIDSEP=self.data.doBGSubIDSEP,
            OPSEPEnhancement=self.data.OPSEPEnhancement,
            IDSEPEnhancement=self.data.IDSEPEnhancement)


    def plot_fluence_spectra(self):
        """ Plots the fluence spectra from all event definitions """

        #Collect calculated values from Analyze objects
        event_definitions, analyzed_fluxes, sep_start_times, sep_end_times,\
        onset_peaks, onset_peak_times, max_fluxes, max_flux_times, fluences,\
        fluence_spectra, fluence_spectra_units = self.extract_analyze_lists()
        
        #Check that at least one SEP crossed threshold, otherwise
        #don't need to generate fluence plot.
        plot_fluence = False
        for time in sep_start_times:
            if not pd.isnull(time):
                plot_fluence = True

        if not plot_fluence: return

        plt_tools.opsep_plot_fluence_spectrum(self.data.experiment, self.data.flux_type,
            self.data.user_name, self.data.options,
            self.data.event_definitions, self.data.evaluated_dates,
            self.data.energy_bin_centers,
            fluence_spectra, fluence_spectra_units, self.data.showplot,
            self.data.saveplot, spacecraft=self.data.spacecraft,
            doBGSubOPSEP=self.data.doBGSubOPSEP, doBGSubIDSEP=self.data.doBGSubIDSEP,
            OPSEPEnhancement=self.data.OPSEPEnhancement,
            IDSEPEnhancement=self.data.IDSEPEnhancement)
        
        

###############################################
######### CLASSES FOR SEP VALUES ##############
###############################################
### Correspond to json entries in CCMC json ###
### format used by the SEP Scoreboards and ####
### the SPHINX validation framework. ##########
### Observation or Forecast classes. ##########
###############################################

###---OBSERVATION CLASS-------
#If the input data is observational data
#measured by a spacecraft.
class Observation():
    def __init__(self, energy_channel):
        
        self.label = 'observation'
        self.energy_channel = energy_channel #dict
        self.short_name = None
        self.issue_time = pd.NaT

        
        #General info
        self.species = None
        self.location = None
        self.observation_window_start = pd.NaT
        self.observation_window_end = pd.NaT
        
        
        #Observed Values
        self.all_clear = All_Clear(None, np.nan, None, np.nan) #All_Clear object
        self.peak_intensity = Flux_Intensity('peak_intensity', np.nan, None, np.nan, np.nan, np.nan, pd.NaT) #Flux_Intensity object
        self.peak_intensity_max = Flux_Intensity('peak_intensity_max', np.nan, None, np.nan, np.nan, np.nan, pd.NaT)#Flux_Intensity object
        self.event_lengths = []
        self.fluences = []
        self.fluence_spectra = []
        self.threshold_crossings = []
        self.sep_profile = None

        #Triggers, if known
        #Triggers
        self.cmes = []
        self.flares = []

        return



##-----FORECAST CLASS------
class Forecast():
    def __init__(self, energy_channel):
        
        self.label = 'forecast'
        self.energy_channel = energy_channel #dict
        #all thresholds applied in forecasted quantities for this
        #energy channel. Will be a list if present.
        self.all_thresholds = []
        self.short_name = None
        self.original_short_name = None
        self.issue_time = pd.NaT
        self.valid = None #indicates whether prediction window starts
                          #at the same time or after triggers/inputs
        self.invalid_reason = '' #If invalid, save reason
        
        #General info
        self.species = None
        self.location = None
        self.prediction_window_start = pd.NaT
        self.prediction_window_end = pd.NaT
        
        
        #Triggers
        self.cmes = []
        self.cme_simulations = []
        self.flares = []
        self.particle_intensities = []
        
        #Inputs
        self.magnetic_connectivity = []
        self.magnetograms = []
        self.coronagraphs = []
        self.human_evaluations = []
        
        
        #Forecasts
        self.source = None #source from which forcasts ingested
                        #JSON filename or perhaps database in future
        self.path = None #Path to JSON file
        self.all_clear = All_Clear(None, np.nan, None, np.nan) #All_Clear object
        self.point_intensity = Flux_Intensity('point_intensity', np.nan, None, np.nan, np.nan, np.nan, pd.NaT)
        self.peak_intensity = Flux_Intensity('peak_intensity', np.nan, None, np.nan, np.nan, np.nan, pd.NaT) #Flux_Intensity object
        self.peak_intensity_max = Flux_Intensity('peak_intensity_max', np.nan, None, np.nan, np.nan, np.nan, pd.NaT) #Flux_Intensity object
        self.event_lengths = []
        self.fluences = []
        self.fluence_spectra = []
        self.threshold_crossings = []
        self.probabilities = []
        self.sep_profile = None

        return


########### SUB CLASSES ###############
#Classes for all the types of values in
#Observation and Forecast classes
class All_Clear:
    def __init__(self, all_clear, threshold, threshold_units,
                probability_threshold):
        """
        Input:
            :self: (object) All_Clear object
            :all_clear: (boolean) all clear value
            :threshold: (float) threshold applied to get all clear value
            :threshold_units: (astropy units)
            :probability_threshold: (float) threshold applied to derive
                all clear from probabilistic model
        
        Output: an All_Clear object
        
        """
        self.label = 'all_clear'
        self.all_clear_boolean = all_clear
        self.threshold = threshold
        self.threshold_units = threshold_units
        self.probability_threshold = probability_threshold
    
        return

    def SetValues(self, all_clear, threshold, threshold_units,
        probability_threshold):
        self.all_clear = all_clear
        self.threshold = threshold
        self.threshold_units = threshold_units
        self.probability_threshold = probability_threshold
        
        return


class Flux_Intensity:
    def __init__(self, label, intensity, units, uncertainty, uncertainty_low,
        uncertainty_high, time):
        """
        Input:
            :self: (object) Peak_Intensity object
            :intensity: (float) intensity value
            :units: (astropy) units
            :uncertainty: (float)
            :uncertainty_low: (float)
            :uncertainty_high: (float)
            :time: (datetime)
        
        Output: a Flux_Intensity object
        
        """
        self.label = label
        self.intensity = intensity
        self.units = units
        self.uncertainty = uncertainty
        self.uncertainty_low = uncertainty_low
        self.uncertainty_high = uncertainty_high
        self.time = time
        
        return



class Event_Length:
    def __init__(self, start_time, end_time, threshold, threshold_units):
        """
        Input:
            :self: (object) Event_Length object
            :start_time: (datetime) start of SEP event
            :end_time: (datetime) end of SEP event
            :threshold: (float) threshold to determine start of SEP event
            :threshold_units: (astropy) units
        
        Output: an Event_Length object
        
        """
        self.label = 'event_length'
        self.start_time = start_time
        self.end_time = end_time
        self.threshold = threshold
        self.threshold_units = threshold_units
        
        return
        

class Fluence:
    def __init__(self, id, fluence, units, threshold, threshold_units,
        uncertainty_low, uncertainty_high):
        """
        Input:
            :self: (object) Fluence object
            :id: (string) unique id
            :fluence: (float)
            :units: (astropy) units
            :uncertainty_low: (float)
            :uncertainty_high: (float)
        
        Output: a Fluence object
        
        """
        self.label = 'fluence'
        self.id = id
        self.fluence = fluence
        self.units = units
        self.threshold = threshold
        self.threshold_units = threshold_units
        self.uncertainty_low = uncertainty_low
        self.uncertainty_high = uncertainty_high
        
        return
        
class Fluence_Spectrum:
    def __init__(self, start_time, end_time, threshold_start, threshold_end, threshold_units, fluence_units, fluence_spectrum):
        """
        Input:
            :self: (object) Fluence_Spectrum object
            :id: (string) unique id
            :start_time: (datetime) start of SEP event
            :end_time: (datetime) end of SEP event
            :threshold_start: (float) threshold to determine start of SEP event
            :threshold_end: (float) threshold to determine end of
                SEP event; if not present, then assume
                threshold_end = threshold_start
            :threshold_units: (astropy units)
            :fluence_units: (astropy) units
            :fluence_spectrum: (array of dict)
            e.g. [{"energy_min": 5.0, "energy_max": -1, "fluence": 78527636.38502692}, {"energy_min": 10.0, "energy_max": -1, "fluence": 46371821.92788475}, {"energy_min": 30.0, "energy_max": -1, "fluence": 16355421.889077082}, {"energy_min": 50.0, "energy_max": -1, "fluence": 7673363.706302568}, {"energy_min": 60.0, "energy_max": -1, "fluence": 5425386.382761811}, {"energy_min": 100.0, "energy_max": -1, "fluence": 2085984.6018625232}, {"energy_min": 700.0, "energy_max": -1, "fluence": 187.6881309476662}]}]
 
        Output: a Fluence_Spectrum object
        
        """
        self.label = 'fluence_spectrum'
        self.start_time = start_time
        self.end_time = end_time
        self.threshold_start = threshold_start
        self.threshold_end = threshold_end
        self.threshold_units = threshold_units
        self.fluence_units = fluence_units
        self.fluence_spectrum = fluence_spectrum
        
        return


class Threshold_Crossing:
    def __init__(self, crossing_time, uncertainty, threshold, threshold_units):
        """
        Input:
            :self: (object) Threshold_Crossing object
            :crossing_time: (datetime) start of SEP event
            :threshold: (float) threshold to determine start of SEP event
            :threshold_units: (astropy units)

        Output: a Threshold_Crossing object
        
        """
        self.label = 'threshold_crossing'
        self.crossing_time = crossing_time
        self.uncertainty = uncertainty
        self.threshold = threshold
        self.threshold_units = threshold_units
        
        return


class Probability:
    def __init__(self, probability_value, uncertainty, threshold,
        threshold_units):
        """
        Input:
            :self: (object) Threshold_Crossing object
            :probability: (float)
            :uncertainty: (float)
            :threshold: (float) threshold to determine start of SEP event
            :threshold_units: (astropy units)

        Output: a Probability object
        
        """
        self.label = 'probability'
        self.probability_value = probability_value
        self.uncertainty = uncertainty
        self.threshold = threshold
        self.threshold_units = threshold_units
        
        return


#If user inputs flux data output by a model to generate
#jsons for a model forecast, can add these triggers and
#inputs associated with the model forecast.
#For observed SEP events with known flare and CME triggers,
#these classes can be populated in the Observation class.
#------ TRIGGERS ------
class CME:
    def __init__(self, cme_id, start_time, liftoff_time, lat, lon, pa,
                half_width, coordinates, speed, catalog, catalog_id):
        """
        Input:
            :self: (object) CME object
            :cme_id: (string) unique id for CME (useful if multiple CMEs
                used as triggers for a forecast)
            :start_time: (datetime)
            :liftoff_time: (datetime)
            :lat: (float) latitude of CME source eruption
            :lon: (float) longitude of CME source eruption
            :pa: (float) position angle
            :half_width: (float) half width of CME cone
            :speed: (float) CME speed
            :catalog: (string) CME catalog
            :catalog_id: (string) id of CME catalog
        
        Output: a CME object
        
        """
        self.label = 'cme'
        self.id = cme_id
        self.start_time = start_time
        self.liftoff_time = liftoff_time
        self.lat = lat
        self.lon = lon
        self.pa = pa
        self.half_width = half_width
        self.speed = speed
        self.coordinates = coordinates
        self.catalog = catalog
        self.catalog_id = catalog_id
        self.cme_allowed_tags = ['start_time', 'liftoff_time', 'lat',
                'lon', 'pa', 'half_width', 'speed', 'coordinates',
                'catalog', 'catalog_id']
 
        return
        
    def SetValues(self, cme_id, start_time, liftoff_time, lat, lon, pa,
                half_width, coordinates, speed, catalog, catalog_id):
        
        self.id = cme_id
        self.start_time = start_time
        self.liftoff_time = liftoff_time
        self.lat = lat
        self.lon = lon
        self.pa = pa
        self.half_width = half_width
        self.speed = speed
        self.coordinates = coordinates
        self.catalog = catalog
        self.catalog_id = catalog_id
        
        return


class CME_Simulation:
    def __init__(self, cmesim_id, model, sim_completion_time):
        """
        Input:
            :model: (string) model
            :sim_completion_time: (datetime) simulation completion time
            
        Output: a CME_Simulation object
        
        """
        self.label = 'cme_simulation'
        self.id = cmesim_id
        self.model = model
        self.completion_time = sim_completion_time
        
        return

    def SetValues(self, cmesim_id, model, sim_completion_time):
        self.id = cmesim_id
        self.model = model
        self.completion_time = sim_completion_time
        
        return


class Flare:
    def __init__(self, flare_id, last_data_time, start_time, peak_time, \
        end_time, location, lat, lon, intensity, integrated_intensity, \
        noaa_region):
        """
        Input:
            :last_data_time: (datetime)
            :start_time: (datetime)
            :peak_time: (datetime)
            :end_time: (datetime)
            :location: (string) location of flare source eruption N00W00
            :lat: (int) latitude
            :lon: (int) longitude
            :intensity: (float) X-ray intensity of flare at last_data_time
            :integrated_intensity: (float) X-ray intensity summed from start  
                to last
            :noaa_region: (string) identifier of NOAA active region
       
        Ouput: a Flare object
            
        """
        self.label = 'flare'
        self.id = flare_id
        self.last_data_time = last_data_time
        self.start_time = start_time
        self.peak_time = peak_time
        self.end_time = end_time
        self.location = location
        self.lat = lat
        self.lon = lon
        self.intensity = intensity
        self.integrated_intensity = integrated_intensity
        self.noaa_region = noaa_region
        
        return
        
    def SetValues(self, flare_id, last_data_time, start_time, peak_time, end_time, \
        location, lat, lon, intensity, integrated_intensity, noaa_region):
        self.id = flare_id
        self.last_data_time = last_data_time
        self.start_time = start_time
        self.peak_time = peak_time
        self.end_time = end_time
        self.location = location
        self.lat = lat
        self.lon = lon
        self.intensity = intensity
        self.integrated_intensity = integrated_intensity
        self.noaa_region = noaa_region
        
        return
            

class Particle_Intensity:
    def __init__(self, part_id, observatory, instrument, last_data_time,
        ongoing_events):
        """
        Inputs:
            :part_id: (string) unique identifier for this measurement
            :observatory: (string)
            :instrument: (string)
            :last_data_time: (datetime)
            :ongoing_events: (array)
            
        Ouput: a Particle_Intensity object
        
        """
        self.label = 'particle_intensity'
        self.id = part_id
        self.observatory = observatory
        self.instrument = instrument
        self.last_data_time = last_data_time
        self.ongoing_events = ongoing_events #array of dict
        
        return
        
    def SetValues(part_id, observatory, instrument, last_data_time,
        ongoing_events):
        self.id = part_id
        self.observatory = observatory
        self.instrument = instrument
        self.last_data_time = last_data_time
        self.ongoing_events = ongoing_events
    
        return

class HumanEvaluation:
    def __init__(self, human_evaluation_id, last_data_time):
        """
        Inputs:
            :human_evaluation_id: (string) unique identifier for human evaluation entry
            :last_data_time: (datetime)

        Output: A HumanEvaluation object
        """
        self.label = 'human_evaluation'
        self.id = human_evaluation_id
        self.last_data_time = last_data_time

        return


##INPUTS
class Magnetic_Connectivity:
    def __init__(self, magcon_id, method, lat, lon, connection_angle,
        solar_wind):
        """
        Input:
            :method: (string)
            :lat: (float)
            :lon: (float)
            :connection_angle: (dict)
            :solar_wind: (dict)
            
        Output: A Magnetic_Connectivity object
        
        """
        self.label = 'magnetic_connectivity'
        self.id = magcon_id
        self.method = method
        self.lat = lat
        self.lon = lon
        self.connection_angle = connection_angle
        self.solar_wind = solar_wind
        
        return
        
    def SetValues(self, magcon_id, method, lat, lon, connection_angle,
        solar_wind):
        self.id = magcon_id
        self.method = method
        self.lat = lat
        self.lon = lon
        self.connection_angle = connection_angle
        self.solar_wind = solar_wind
        
        return



class Magnetogram:
    def __init__(self, magneto_id, observatory, instrument, products):
        """
        Input:
            :magneto_id: (string) unique identifier for magnetogram entry
            :observatory: (string)
            :instrument: (string)
            :products: (array of dict)
        
        
        products has the format e.g.,
        [{"product":"hmi.M_45s_nrt","last_data_time":"2022-03-26T00:59Z"},
         {"product":"hmi.sharp_720s_nrt","last_data_time":"2022-03-26T00:58Z"}]
        
        Output: A Magnetogram object
        
        """
        self.label = 'magnetogram'
        self.id = magneto_id
        self.observatory = observatory
        self.instrument = instrument
        self.products = products
    
        return


class Coronagraph:
    def __init__(self, corona_id, observatory, instrument, products):
        """
        Input:
            :corona_id: (string) unique identifier for coronagraph entry
            :observatory: (string)
            :instrument: (string)
            :products: (array of dict)
        
        Output: A Coronagraph object
        
        """
        self.label = 'coronagraph'
        self.id = corona_id
        self.observatory = observatory
        self.instrument = instrument
        self.products = products
    
        return
