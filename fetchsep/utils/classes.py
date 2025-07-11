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

        #Determine Observation or Forecast class according
        #to type of input data
        self.datatype = 'observation' #or 'forecast' for model output

        #Paths for data source, plots and output files
        self.datapath = cfg.datapath
        self.outpath = cfg.outpath
        self.plotpath = cfg.plotpath
        
        #Variables for user-input files
        self.user_delim = cfg.user_delim #delimeter between columns
        self.use_col = cfg.user_col #columns in input data file containing fluxes
        #self.err_col = cfg.err_col #columns containing error bars
        self.user_energy_bins = cfg.user_energy_bins #energy bins for each column

        #User-input data file
        self.user = False #user file?
        self.model_name = None
        self.user_filename = None #If user-input file
        
        #Original fluxes before any background subtraction or interpolation
        #Bad values are set to None
        self.original_dates = []
        self.original_fluxes = []
        
        #Do background subtraction on the fluxes?
        self.doBGSub = False #True to do background-subtraction
        self.bgstartdate = pd.NaT #defaults to idsep output if not set
        self.bgenddate = pd.NaT #defaults to idsep output if not set
        self.bgmeans = [] #Mean fluxes for each energy channel
        self.bgsigmas = [] #Sigmas for each energy channel
        self.bgdates = []
        self.bgfluxes = []
        
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
        self.energy_bins = []
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


    def set_background_subtraction_info(self, doBGSub, bgstartdate, bgenddate):
        """ Indicate whether to perform background-subtraction.
            If start and end dates aren't set, then code
            will look for idsep output to use for mean background.
        
            INPUT:
            
                :doBGSub: (bool) bg subtraction if True
                :bg_startdate: (str) start of time period to use for 
                    background calculation
                :bg_enddate: (str) end of time period to use for
                    background calculation
                    
            OUTPUT:
            
                Set background attributes in Data object
        
        """
        self.doBGSub = doBGSub
        self.bgstartdate = self.str_to_datetime(bgstartdate)
        self.bgenddate = self.str_to_datetime(bgenddate)

        if doBGSub:
            if pd.isnull(self.bgstartdate) or pd.isnull(self.bgenddate):
                print("WARNING!!! User selected to perform background-subtraction, but did not provide dates. Will look for idsep output files to extract mean background.")
        
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
                    bins.append([float(bin_edge[0]), float(bins_edge[1])])
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
            
        return

        
    def error_check(self):
        """ Error check the inputs and options. """
            
        error_check.error_check_options(self.experiment, self.flux_type, self.options, self.doBGSub, spacecraft=self.spacecraft)
        error_check.error_check_inputs(self.startdate, self.enddate, self.experiment, self.flux_type)

        return


    def load_info(self, startdate, enddate, experiment, flux_type,
        model_name, user_file, showplot, saveplot, two_peaks,
        definitions, options, doBGSub, bgstartdate,
        bgenddate, nointerp, spacecraft):
        """ Create new Data object and load with all values.
        
            INPUT:
                :startdate: (string) - user input start date 
                    "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS"
                :enddate: (string) - user input end date "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS"
                :experiment: (string) - "GOES-05" up to "GOES-19", "SEPEM", "SEPEMv3","EPHIN", "EPHIN_REleASE", or "user"
                :flux_type: (string) - "integral" or "differential" 
                    indicates the type of flux to read in
                :model_name: (string) - If model is "user", set 
                    model_name to describe your model or data set (e.g. 
                    MyModel), otherwise set to ''.
                :user_file: (string) - Default is ''. If "user" is 
                    selected for experiment, specify name of flux file.
                :showplot: (bool) - True to show plots
                :saveplot: (bool) - True to save plots 
                :two_peaks: (bool) - option for extending event length
                :definitions: (string) - user-input thresholds in the 
                    format "30,1;4-7,0.01" multiple thresholds
                    separated by semi-colon. Same as str_thresh in opsep
                :nointerp: (boolean) - True to fill in negative fluxes 
                    with None instead of linear interpolation in time
                :spacecraft: (string) primary or secondary 
                
            OUTPUT:
            
                :input_data: (Data Object)
        
        """
        
        if experiment == "user":
            self.user = True
            self.label = f"{model_name} {flux_type}"
        else:
            self.label = f"{experiment} {flux_type}"

        self.experiment = experiment
        self.flux_type = flux_type
        self.model_name = model_name #user input experiment (model or obs name)
        self.user_filename = user_file
        self.set_dates(startdate, enddate)
        self.user_filename = user_file #default tmp.txt
        self.spacecraft = spacecraft
        self.set_options(options)
        self.set_background_subtraction_info(doBGSub, bgstartdate, bgenddate)
        self.set_event_definitions(definitions)
        self.two_peaks = two_peaks
        self.showplot = showplot
        self.saveplot = saveplot
        self.dointerpolation = not(nointerp)
        self.error_check()
        
        return


    def plot_background_subtraction(self, showplot=False):
        """ Make plots of background-subtracted fluxes """

        plt_tools.opsep_plot_bgfluxes(f"Total_{experiment}", self.flux_type, self.options,
                self.model_name, self.original_fluxes, self.original_dates,
                self.energy_bins, self.means, self.sigmas, self.saveplot,
                spacecraft = self.spacecraft)
        plt_tools.opsep_plot_bgfluxes(f"BackgroundFlux_{experiment}", self.flux_type,
                self.options, self.model_name, self.bgfluxes, self.bgdates,
                self.energy_bins, self.means, self.sigmas, self.saveplot,
                spacecraft = self.spacecraft)
        plt_tools.opsep_plot_bgfluxes(f"BGSubSEPFlux_{experiment}", self.flux_type,
                self.options, self.model_name, self.fluxes, self.dates,
                self.energy_bins, self.means, self.sigmas, self.saveplot,
                spacecraft = self.spacecraft)
        
        if showplot:
            plt.show()
        
        return



    def read_in_flux_files(self):
        """ Read in the appropriate data or user files. Performs
            background subtraction, if requested. Trims to dates
            between start time and end time. Interpolates bad
            points with linear interpolation in time.
            
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
        if self.doBGSub:
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
                all_dates, all_fluxes, west_detector, energy_bins = \
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
                        user_delim=self.user_delim, flux_col=self.user_col)


        #Define energy bins
        if self.experiment == "ERNE":
            version = datasets.which_erne(startdate, enddate)
            energy_bins = datasets.define_energy_bins(version, self.flux_type,
                                west_detector, self.options)
        elif self.experiment != "GOES":
            energy_bins = datasets.define_energy_bins(self.experiment, self.flux_type,
                        west_detector, self.options, spacecraft=self.spacecraft,
                        user_bins=self.user_energy_bins)


        if len(all_dates) <= 1:
            sys.exit(f"read_in_flux_files: The specified start and end dates ({startdate} to {enddate}) were not present in the specified input file or were too restrictive. Exiting.")

        #Full flux and date range for specified input files, not yet trimmed in date
        all_fluxes, energy_bins = tools.sort_bin_order(all_fluxes, energy_bins)
        
        self.energy_bins = energy_bins
        
        ####Save original fluxes with bad points set to None
        #Extract date range that covers any background-subtraction periods
        orig_dates, orig_fluxes = datasets.extract_date_range(startdate, enddate,
                                        all_dates, all_fluxes)
        orig_fluxes = datasets.check_for_bad_data(orig_dates,orig_fluxes,energy_bins,False)
        self.original_dates = orig_dates
        self.original_fluxes = orig_fluxes


        #IF BACKGROUND SUBTRACTION
        if self.doBGSub:
            #sepfluxes are background subtracted fluxes
            #Previous version of subroutine had read in the data
            #again and did not apply linear interpolation on bad points
            nointerp_fluxes = datasets.check_for_bad_data(all_dates,all_fluxes,
                energy_bins,False)
            bgfluxes, sepfluxes, means, sigmas = bgsub.derive_background(self.bgstartdate, self.bgenddate, all_dates, nointerp_fluxes, energy_bins, self.showplot,
                self.saveplot)
            self.bgmeans = means
            self.bgsigmas = sigmas
            self.bgfluxes = bgfluxes
            self.bgdates = all_dates
            #Extract the date range specified by the user for the
            #background-subtracted fluxes
            dates, fluxes = datasets.extract_date_range(startdate, enddate,
                                        all_dates, sepfluxes)

        #NO BACKGROUND SUBTRACTION
        if not self.doBGSub:
            #Extract the date range specified by the user
            dates, fluxes = datasets.extract_date_range(self.startdate, self.enddate,
                                    all_dates, all_fluxes)
        
        #Handle bad data points
        fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,self.do_interpolation)
         
        if len(dates) <= 1:
            print("read_in_flux_files: The specified start and end dates were not "
                f"present in the specified input file. Exiting. {startdate} to {enddate}")
            sys.exit()
        
        self.fluxes = fluxes
        self.dates = dates
        
        time_res = tools.determine_time_resolution(dates)
        self.time_resolution = time_res.total_seconds()

        if self.doBGSub and (self.showplot or self.saveplot):
            plot_background_subtraction()

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
                            bruno2017=self.goes_Bruno2017)
                
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
                continue
            
            self.evaluated_energy_bins.append(bin)
            self.evaluated_fluxes.append(self.fluxes[idx])
    
        return


    def add_results(self, analyze):
        """ Append Analyze object to Data object """
        self.results.append(analyze)
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

        self.flux_units = None
        self.rise_time_units = 'minutes'
        self.duration_units = 'hours'
        self.fluence_units = None
        self.fluence_spectrum_units = None

        self.check_event_definition(data)


 
    def check_event_definition(self, data):
        """ Check if the energy channels in the data correspond to  
            the requested event definition. Exit if not.
            
        """
        energy_bin = self.make_energy_bin()
        try:
            data.evaluated_energy_bins.index(energy_bin)
        except:
            sys.exit(f"Analyze init: Requested energy bin in the event definition {energy_bin} "
                f"is not present in the data: {data.evaluated_energy_bins}. Exiting.")
        
        print(f"Analyze init: Applying event definition: "
            f"[{self.event_definition['energy_channel'].min}, "
            f"{self.event_definition['energy_channel'].min} exceeds "
            f"{self.event_definition['threshold'].threshold} {self.event_definition['threshold'].threshold_units}")
 
 
 
    def make_energy_channel_dict(self):
        """ {'min': 10, 'max': -1, 'units': 'MeV'} from event_definition 
            
            This form of the energy channel is used in SPHINX and the 
            Observation and Forecast objects.
        
        """
        energy_channel = {'min':self.event_definition['energy_channel'].min,
                            'max':self.event_definition['energy_channel'].max,
                            'units': elf.event_definition['energy_channel'].units}
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



    def calculate_threshold_crossing(self, data, event_definition):
        """ Calculate the threshold crossing times for a given energy bin
            and flux threshold. 

            An SEP event is considered to start if 3 consecutive points (depending on 
            dataset time resolution) are above threshold. The start time is set to the 
            first point that crossed threshold.
            
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
        
        npoints = 3 #require 3 points above threshold
        if data.time_resolution/60. > 15:
            npoints = 1 #time resolution >15 mins, require one point above threshold
            
        if energy_bin not in data.evaluated_energy_bins:
            print(f"calculated_threshold_crossing: Requested energy bin {energy_bin} not "
                "specified in event definitions or not present in data. Skipping.")
            return

        idx = data.evaluated_energy_bins.index(energy_bin)
        fluxes = data.evaluated_fluxes[idx]
        dates = data.evaluated_dates
        
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
                    if elapse > cfg.dwell_time: #N consecutive points longer than dwell time
                        event_ended = True
                        sep_end_time = dates[i-(end_counter-1)] #correct back time steps
                        #Double checked some calculated event end times with SWPC and
                        #this logic gave the correct end times. 2023-04-10 KW

        self.sep_start_time = sep_start_time
        self.sep_end_time = sep_end_time
        
        return sep_start_time, sep_end_time


    def trim_to_date_range(self, startdate, enddate, dates, array):
        """ Trim array to between startdate and enddate.
            dates corresponds to the time steps in array.
            
        """
        
        indices = [i for i in range(len(dates)) if (dates[i] >= startdate and dates[i] <= enddate)]
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

        idx = data.evaluated_energy_bins.index(energy_bin)
        fluxes = data.evaluated_fluxes[idx]
        dates = data.evaluated_dates
        
        #If there is a SEP EVENT, extract between start and end times.
        #If NO EVENT, don't trim and take the maximum of the full timeseries
        if not pd.isnull(self.sep_start_time) and not pd.isnull(self.sep_end_time):
            fluxes = self.trim_to_date_range(self.sep_start_time, self.sep_end_time,
                                    dates, fluxes)
            dates = self.trim_to_date_range(self.sep_start_time, self.sep_end_time,
                                    dates, dates)


        max_flux = max(fluxes)
        ix = np.where(fluxes == max_flux) #First instance, if multiple
        try:
            max_flux_time = dates[ix[0][0]]
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

        idx = data.evaluated_energy_bins.index(energy_bin)
        fluxes = data.evaluated_fluxes[idx]
        dates = data.evaluated_dates

        #Do a fit of the Weibull function for each time profile
        params_fit = Parameters()
        params_fit.add('alpha', value = -3, min = -5, max = -0.1)
        params_fit.add('beta', value = 10, min = 1, max =100)
        params_fit.add('peak_intensity', value = 100, min = 1e-3, max =1e6)

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
                    energy_units, low_thresh, thresh_units)
                low_start_time, low_end_time = calculate_threshold_crossing(data, low_evdef)
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

        print(f"calculate_onset_peak_from_fit ==== {energy_bin} MeV =====")
        print(f"Best Fit Ip: {best_Ip}, a: {best_a}, b: {best_b}")
        
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
                self.sep_start_time,trim_times, trim_fluxes, best_pars, best_fit, max_time,
                max_val, max_meas_time, max_meas, max_curve_model_time, max_curve_model_peak,
                max_curve_meas_time, max_curve_meas_peak,
                data.saveplot, data.showplot, spacecraft=data.spacecraft)

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
        fluence = sum(clean_flux)*time_resolution*4.0*math.pi #multiply 4pi steradians
    
        return fluence
 
 
    def calculate_channel_fluence(self, data):
        """  Calculate the fluence for the specified event definition
        """
        energy_bin = self.make_energy_bin()
        idx = data.evaluated_energy_bins.index(energy_bin)

        flux = data.evaluated_fluxes[idx]
        dates = data.evaluated_dates

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
        energy_bin = self.make_energy_bin()
        energy_units = self.event_definition['energy_channel'].units
        threshold = self.event_definition['threshold'].threshold
        threshold_units = self.event_definition['threshold'].threshold_units
        
        #calculate event values and fill in a dictionary that will
        #save info needed for Observation or Forecast objects
        self.calculate_threshold_crossing(data, self.event_definition)
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

        energy_label = f"{energy_bin[0]} - {energy_bin[1]} {energy_units}"
        if energy_bin[1] == -1:
            energy_label = f">{energy_bin[0]} {energy_units}"
        threshold_label = f"{threshold} {threshold_units}"
        print(f"====SEP Event Characteristics for {energy_label}, {threshold_label}====")
        print(f"SEP Start Time: {self.sep_start_time}")
        print(f"SEP End Time: {self.sep_end_time}")
        print(f"Onset Peak: {self.onset_peak} {self.flux_units} at {self.onset_peak_time}")
        print(f"Max Flux: {self.max_flux} {self.flux_units} at {self.max_flux_time}")
        print(f"Rise time to Onset: {self.onset_rise_time} {self.rise_time_units}")
        print(f"Rise time to Max: {self.max_rise_time} {self.rise_time_units}")
        print(f"Duration: {self.duration} {self.duration_units}")
        print(f"Channel Fluence: {self.fluence} {self.fluence_units}")
        print(f"Fluence Spectrum: {self.fluence_spectrum} {self.fluence_units}")

        #Create a dictionary containing all of the calculated values
        dict = tools.fill_event_info(data.experiment, data.flux_type, self.event_definition,
            data.startdate, data.enddate, data.energy_bins, data.doBGSub, data.options,
            self.sep_start_time, self.sep_end_time,  self.onset_peak, self.onset_peak_time, self.onset_rise_time,
            self.max_flux, self.max_flux_time, self.max_rise_time, self.duration, self.fluence, self.fluence_spectrum)
        

        return dict




class Output:
    def __init__(self, data, json_type, spase_id=None, location="earth",
        species="proton"):
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
    
        self.spase_id = None
        self.json_type = None #observation or forecast
        self.json_dict = {} #json dictionary from template
        self.location = location #earth, mars, etc
        self.species = species #protons, electrons
        

    def set_json_type(self, json_type):
        self.json_type = json_type
        return
        
    
    def set_spase_id(self, spase_id):
        self.spase_id = spase_id
        return


    def intialize_json(self, type):
        """ Construct the main components of the json files """
        dict = ccmc_json.initialize_json(self.json_type)
        self.json_dict = dict
        
        return template


    def fill_json_header(self):
        """ Fill the json dictionary with the extracted SEP values. """
        issue_time = datetime.datetime.now()
        
        self.json_template = ccmc_json.fill_json_header(self.json_template, self.json_type,
            issue_time, self.data.experiment, self.data.flux_type, self.spase_id,
            model_name=self.data.model_name)

        return


    def fill_json_block(self, analyze):
        """ Store SEP event values into the json block associated with
            the event definition in the Analyze object.
            
        """

#        self.event_definition = event_definition
#        
#        #Derived values
#        self.sep_start_time = pd.NaT
#        self.sep_end_time = pd.NaT
#        self.onset_peak = np.nan
#        self.onset_peak_time = pd.NaT
#        self.onset_rise_time = np.nan
#        self.max_flux = np.nan
#        self.max_flux_time = pd.NaT
#        self.max_flux_rise_time = np.nan
#        self.duration = np.nan
#        self.fluence = np.nan
#        self.fluence_spectrum = []
#
#        self.flux_units = None
#        self.rise_time_units = 'minutes'
#        self.duration_units = 'hours'
#        self.fluence_units = None
#        self.fluence_spectrum_units = None


        



    def write_ccmc_json(self):
        """ Write all event definitions out to CCMC json file 
            https://ccmc.gsfc.nasa.gov/publicData/sepsb/files/sepscoreboard_visual_schema.pdf
        
        """
        self.read_in_template(self.json_type, filename=self.json_template_filename)
        self.fill_json_header()
        



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
