from ..utils import config as cfg
from ..utils import names
from ..utils import date_handler as dh
from ..utils import experiments as expts
from ..utils import directories as dirs
from ..utils import error_check
import os
import pandas as pd

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

""" Class to store the choices made by the user at run
    and any information derivative of those choices.
    
"""

class Parameters:
    def __init__(self, module, startdate, enddate, experiment):
        """ Parameters that may be set by the user. 
            Start with default values and change if user specifies a different value.
            
        """

        #### INPUT BY USER #####
        self.module = module #opsep, idsep, download
        self.experiment = experiment

        #Dates
        if startdate == "" or enddate == "" or startdate == None or enddate == None:
            sys.exit('You must enter start and end dates. Exiting.')
        self.startdate = dh.str_to_datetime(startdate)
        self.enddate = dh.str_to_datetime(enddate)
        print(f"Set analysis dates {self.startdate} to {self.enddate}")

        self.flux_type = None
        self.spacecraft = '' #GOES only; primary or secondary

        #User-input data file
        self.user = False #user file?
        self.user_name = ''
        self.user_filename = '' #If user-input file
        self.is_unixtime = False

        #Directory behavior
        self.directory_depth = 2
        self.use_absolute_datapath = False

        #GOES-specific options
        self.options = []
        self.goes_datatype = 'corrected' #corrected or uncorrected
        self.goes_S14 = False #True to apply Sandberg et al. 2014 effective energies
        self.goes_Bruno2017 = False #True to apply Bruno et al. 2017 effective energies

        #Plotting settings
        self.showplot = False
        self.saveplot = False
 
        #### Directory and plotting names ####
        self.modifier = '' #for plots and filenames
        self.title_modifier = '' #For plot titles
        self.idsep_subdir = '' #Need idsep paths in opsep and download too
        self.idsep_outpath = ''
        self.idsep_plotpath = ''

        self.idsep_path = '' #May be set by user as path to background_mean_fluxes_FINAL.csv, etc

        self.module_subdir = ''
        self.module_outpath = ''
        self.module_plotpath = ''

        ######## IDSEP-SPECIFIC PARAMETERS #######
        self.remove_above=999999
        self.for_inclusive=False
        self.idsep_nsigma=cfg.idsep_nsigma
        self.init_win=cfg.init_win
        self.sliding_win=cfg.sliding_win
        self.percent_points=cfg.percent_points
        self.plot_timeseries_only=False
        self.write_fluxes=True
        #When calculating the background in idsep, some time periods may have very
        #non-Gaussian flux distributions. This may cause problems with the
        #background calculation. A cut may be applied on kurtosis of the
        #histogrammed flux distributions to exclude the data from the background
        #calculation. This is la configurable parameter as different data sets may
        #result in different numerical values of kurtosis.
        #
        #Kurtosis values are saved with time in output/idsep/kurtosis*.csv
        #
        #A kurtosis value >= kurtosis_cut will be excluded from background calculation.
        #If no cut is desired, set this variable to a very large value.
        #
        #DEFAULT VALUE SET HIGH SO THAT KURTOSIS CUT IS NOT APPLIED.
        #Values in experiments.py will override this default.
        #A user-input value when running idsep will override all.
        self.kurtosis_cut=999

        ######## OPSEP-SPECIFIC PARAMETERS #######
        self.color_scheme=1
        self.no_goes_colors=False #Set to True to turn off SWPC colors
        
        self.location = None #earth, mars, etc
        self.species = None #protons, electrons
        
        self.json_type=''
        self.json_mode=''
        self.spase_id=''
 
        #Fill bad fluxes with nan values (default)
        self.do_interpolation = False
 
        self.user_thresholds = ''
 
        #Do background subtraction on the fluxes?
        #BACKGROUND
        self.opsep_nsigma = cfg.opsep_nsigma #N * sigma to use in background subtraction
        #With OPSEP
        self.doBGSubOPSEP = False #True to do background-subtraction
        self.bgstartdate = pd.NaT #defaults to idsep output if not set
        self.bgenddate = pd.NaT #defaults to idsep output if not set
        self.OPSEPEnhancement = False #enhancement above background
        #With IDSEP
        self.doBGSubIDSEP = False
        self.IDSEPEnhancement = False #Get threshold from idsep output files
        #Use threshold to set background to zero and leave only enhancement
        #Filenames of background fluxes and thresholds from IDSEP

        #Two peaks may extend an event if it crosses threshold,
        #temporarily drops below, then increases above threshold again
        self.two_peaks = False
        self.detect_prev_event = False
    

    def set_options(self, options):
        """ Set options (arr) if any specified """

        #First, reset
#        self.options = []
#        self.goes_datatype = None
#        self.goes_S14 = None
#        self.goes_Bruno2017 = None
    
        if not isinstance(options, str):
            return
    
        options = options.split(";")
        if options[0] != "":
            self.options = options
        
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


    def configure_idsep(self, remove_above=None, kurtosis_cut=None, idsep_nsigma=None,
        init_win=None, sliding_win=None, percent_points=None):
        if remove_above != None:
            self.remove_above = remove_above
            print(f"parameters: Setting remove_above to {remove_above}.")
        if kurtosis_cut != None:
            self.kurtosis_cut = kurtosis_cut
            print(f"parameters: Setting kurtosis_cut to {kurtosis_cut}.")
            print(f"parameters: Superceding any kurtosis_cut value set in experiments.py")
        if idsep_nsigma != None:
            self.idsep_nsigma = idsep_nsigma
            print(f"parameters: Setting idsep_nsigma to {idsep_nsigma}.")
        if init_win != None:
            self.idsep_init_win = init_win
            print(f"parameters: Setting idsep init_win to {init_win}.")
        if sliding_win != None:
            self.idsep_sliding_win = sliding_win
            print(f"parameters: Setting idsep sliding_win to {sliding_win}.")
        if percent_points != None:
            self.idsep_percent_points = percent_points
            print(f"parameters: Setting idsep percent_points to {percent_points}.")


    def set_opsep_background_info(self):
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
        #IF choose to do background subtraction, then automatically choose
        #to calculate enhancement above background
        if self.doBGSubOPSEP: self.OPSEPEnhancement = True
        
        if self.doBGSubOPSEP or self.OPSEPEnhancement:
            if pd.isnull(self.bgstartdate) or pd.isnull(self.bgenddate):
                sys.exit("WARNING!!! User selected to perform background-subtraction, but did not provide dates. Please provide dates of a quiet background period to use this feature.")
        
        return


    def set_idsep_background_info(self):
        """ Specify whether to use background calculated by idsep """

        #If want to use IDSEP files, but no path specified, try the default
        if (self.IDSEPEnhancement or self.doBGSubIDSEP) and self.idsep_path == '':
            self.idsep_path = self.idsep_outpath

        #IF choose to do background subtraction, then automatically choose
        #to calculate enhancement above background
        if self.doBGSubIDSEP: self.IDSEPEnhancement = True

        return


    def error_check(self):
        """ Error check the inputs and options. """
            
        error_check.error_check_options(self.experiment, self.flux_type, self.options, spacecraft=self.spacecraft)
        error_check.error_check_inputs(self.startdate, self.enddate, self.experiment, self.flux_type, module=self.module, init_win=self.init_win)
        error_check.error_check_background(self.experiment, self.flux_type, self.doBGSubOPSEP,
            self.doBGSubIDSEP, self.OPSEPEnhancement, self.IDSEPEnhancement)

        return


    def check_json_info(self):
        json_type = expts.get_json_type(self.experiment)
        
        #If json type isn't specified for the experiment in experiments.py
        if json_type == '' or json_type == None:
            if self.json_type != '':
                pass
            else:
                sys.exit(f"check_json_info: You must specify a json type for {self.data.experiment} experiment. Choose from observations or model. Exiting.")
        
        #User specified wrong json type; override
        elif json_type != '' and self.json_type != '':
            if json_type != self.json_type:
                print(f"check_json_info: You specified a json type of {json_type}, however the "
                    f"correct json_type is {json_type}. Replacing.")
                self.json_type = json_type
        
        #User didn't specify json type and want to get it from experiments.py
        elif json_type != '' and self.json_type == '':
            print(f"check_json_info: Automatically setting json type {json_type}.")
            self.json_type = json_type

 
        json_mode = expts.get_json_mode(self.experiment)
        #If json mode isn't specified for the experiment in experiments.py
        if json_mode == '':
            if self.json_mode != '':
                pass
            else:
                sys.exit(f"check_json_info: You must specify a json mode for {self.data.experiment} experiment. Choose from measurements for observations or forecast, historical, nowcast, etc for model. Recommend to follow the options for the CCMC SEP Scoreboard JSON schema. Exiting.")

        #User specified wrong json mode; override
        elif self.json_mode != '' and json_mode != '':
            if json_mode != self.json_mode:
                print(f"check_json_info: You specified a json mode of {self.json_mode}, however the "
                    f"correct json_mode is {json_mode}. Replacing.")
                self.json_mode = json_mode

        #User didn't specify json mode and want to get it from experiments.py
        elif json_mode != '' and self.json_mode == '':
            print(f"check_json_info: Automatically setting json mode {json_mode}.")
            self.json_mode = json_mode
 
        return
 

    def print_parameters(self):
        print("########## SET PARAMETERS ###########")
        for key, value in vars(self).items():
            print(f"{key}: {value}")
        print("########## END PARAMETERS ###########")

    def set_values(self, flux_type=None,
        spacecraft=None,
        user_name=None,
        user_file=None,
        is_unixtime=None,
        options=None,
        dointerp=None,
        showplot=None,
        saveplot=None,
        directory_depth=None,
        use_absolute_datapath=None,
        write_fluxes=None,
        for_inclusive=None,
        remove_above=None,
        kurtosis_cut=None,
        idsep_nsigma=None,
        init_win=None,
        sliding_win=None,
        percent_points=None,
        opsep_nsigma=None,
        color_scheme=None,
        no_goes_colors=None,
        json_type=None,
        json_mode=None,
        spase_id=None,
        detect_prev_event=None,
        two_peaks=None,
        user_thresholds=None,
        doBGSubOPSEP=None,
        OPSEPEnhancement=None,
        bgstartdate=None,
        bgenddate=None,
        doBGSubIDSEP=None,
        IDSEPEnhancement=None,
        idsep_path=None,
        location=None,
        species=None):
        """ Set all the values related to user choices that are needed across fetchsep. 
            Values that are directly specified by the user (i.e. those not None) will 
            overwrite the default values set when the Parameter object was created.
            
        """
        print(f"flux_type is {flux_type}")
        if flux_type == '' or flux_type == None:
            flux_type = expts.get_flux_type(self.experiment)
        if flux_type != None: self.flux_type = flux_type

        #If user specifies a spacecraft but isn't relevant to experiment,
        #overrides and sets spacecraft to ''
        if spacecraft == '' or spacecraft == None:
            spacecraft = expts.get_spacecraft(self.experiment, spacecraft)
        if spacecraft != None: self.spacecraft = spacecraft #GOES only; primary or secondary

        if self.experiment != 'user':
            exp_info = expts.experiment_info(self.experiment)
            location = exp_info['location']
            species = exp_info['species']
        if location != None: self.location = location #earth, mars, etc
        if species != None: self.species = species #protons, electrons

        if directory_depth != None: self.directory_depth = directory_depth
        if color_scheme != None: self.color_scheme = color_scheme
        if no_goes_colors != None: self.no_goes_colors = no_goes_colors

        if use_absolute_datapath != None:
            self.use_absolute_datapath = use_absolute_datapath

        if write_fluxes != None:
            self.write_fluxes = write_fluxes

        #Plotting settings
        if showplot != None: self.showplot = showplot
        if saveplot != None: self.saveplot = saveplot

        #User-input data file
        if self.experiment == 'user':
            self.user = True
        if user_name != None: self.user_name = user_name
        if user_file != None: self.user_filename = user_file
        if is_unixtime != None: self.is_unixtime = is_unixtime

        #IDSEP VALUES
        if for_inclusive != None:
            self.for_inclusive = for_inclusive
        if kurtosis_cut == None:
            kurtosis_cut = expts.set_kurtosis_cut(self.experiment, self.flux_type)
        self.configure_idsep(remove_above=remove_above, kurtosis_cut=kurtosis_cut,
            idsep_nsigma=idsep_nsigma, init_win=init_win,
            sliding_win=sliding_win, percent_points=percent_points)


        ### JSON FILE INFO
        if json_type == '' or json_type == None:
            json_type = expts.get_json_type(self.experiment)
        self.json_type=json_type

        if json_mode == '' or json_type == None:
            json_mode = expts.get_json_mode(self.experiment)
        if json_mode != None: self.json_mode=json_mode

        self.check_json_info()

        if spase_id != None: self.spase_id = spase_id
 
        #Fill bad fluxes with nan values (default)
        if dointerp != None: self.do_interpolation = dointerp
 
        if user_thresholds!= None: self.user_thresholds = user_thresholds
 
        #GOES-specific options
        if options != None:
            self.set_options(options)
 
        #### Directory and plotting names FOR IDSEP OR DOWNLOAD ####
        #Directory name, filename, and plot title modifiers based on selections
        #Need idsep subdir before setting idsep background info
        id_modifier, id_title_mod = names.setup_modifiers(options, spacecraft=spacecraft)
 
        #Need idsep subdir and outpath in opsep and download, too
        self.idsep_subdir = names.idsep_naming_scheme(self.experiment, self.flux_type, self.user_name, modifier=id_modifier)
        self.idsep_outpath = dirs.create_subdirectories(cfg.outpath, module=self.module,
            subdir=self.idsep_subdir, directory_depth=self.directory_depth)
        
        if self.module != 'opsep':
            self.modifier = id_modifier
            self.title_modifier = id_title_mod
            self.module_subdir = self.idsep_subdir


        #Do background subtraction on the fluxes?
        #BACKGROUND
        if opsep_nsigma != None: self.opsep_nsigma = opsep_nsigma
        #With OPSEP
        if bgstartdate != None: self.bgstartdate = dh.str_to_datetime(bgstartdate)
        if bgenddate != None: self.bgenddate = dh.str_to_datetime(bgenddate)
        if doBGSubOPSEP != None: self.doBGSubOPSEP = doBGSubOPSEP
        if OPSEPEnhancement != None: self.OPSEPEnhancement = OPSEPEnhancement
        self.set_opsep_background_info()
        #With IDSEP
        if idsep_path != None: self.idsep_path = idsep_path
        if doBGSubIDSEP != None: self.doBGSubIDSEP = doBGSubIDSEP
        self.set_idsep_background_info()

        #### Directory and plotting names FOR OPSEP ####
        #Directory name, filename, and plot title modifiers based on selections
        if self.module == 'opsep':
            modifier, title_mod = names.setup_modifiers(options, spacecraft=spacecraft,
                doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
                OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
            self.modifier = modifier
            self.title_modifier = title_mod
     
            #Create subdirectory to hold values
            self.module_subdir = names.opsep_subdir(self.experiment, self.flux_type,
                    self.user_name, modifier=self.modifier)

        self.module_outpath = dirs.create_subdirectories(cfg.outpath, module=self.module,
            subdir=self.module_subdir, directory_depth=self.directory_depth)
        self.module_plotpath = dirs.create_subdirectories(cfg.plotpath, module=self.module,
            subdir=self.module_subdir, directory_depth=self.directory_depth)


        #Two peaks may extend an event if it crosses threshold,
        #temporarily drops below, then increases above threshold again
        if two_peaks != None: self.two_peaks = two_peaks
        if detect_prev_event != None: self.detect_prev_event = detect_prev_event
    
        #Quality controls
        self.error_check()

        self.print_parameters()
