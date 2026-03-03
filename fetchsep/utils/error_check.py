from ..utils import config as cfg
from . import experiments as expts
import datetime
import os
import sys

""" Check inputs and options for errors and incompatabilities."""


def error_check_options(experiment, flux_type, options, subroutine=None,
    spacecraft=""):
    """ Make sure the selected options make sense for the experiment.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral or differential
        :options: (string array) - various options applied to GOES data
        :doBGSub: (boolean) - indicates if background subtraction to be
            performed
        
        OUTPUTS:
        
        no outputs by system exit if error found
        
    """
    if "S14" in options and experiment[0:4] != "GOES":
        sys.exit("Sandberg et al. (2014) effective energies (S14) may only "
            "be applied to GOES data.")
    if "S14" in options and "uncorrected" not in options:
        sys.exit("Sandberg et al. (2014) effective energies (S14) may only be "
            "applied to GOES uncorrected fluxes. Please add "
            "\"uncorrected\" to options.")
    if "uncorrected" in options and flux_type == "integral":
        sys.exit("The uncorrected option cannot be used with integral fluxes. "
                "Please remove this option and run again. Exiting.")
    if "uncorrected" in options and experiment[0:4] != "GOES":
        sys.exit("The uncorrected option may only be specified for GOES "
                "differential fluxes. Exiting.")
    if "Bruno2017" in options and (experiment != "GOES-13" and \
        experiment != "GOES-15"):
        sys.exit("Bruno2017 effective energies may only be appied to GOES-13 "
                "or GOES-15 fluxes. Exiting.")
    if "S14" in options and (experiment == "GOES-13" or \
        experiment == "GOES-15"):
        print("Sandberg et al. (2014) effective energies found for GOES-11 "
            "will be applied to channels P2-P7. Continuing.")
    if "S14" in options and "Bruno2017" in options:
        print("Sandberg et al. (2014) effective energies from GOES-11 will be "
            "applied to P2-P5. Bruno (2017) effective energies will be applied "
            "to P6-P11.")
    if ("uncorrected" in options or "S14" in options or "Bruno2017" in options)\
                and experiment[0:4] != "GOES":
        sys.exit("The options you have selected are only applicable to GOES "
                "data. Please remove these options and run again: "
                "uncorrected, S14, or Bruno2017.")


def error_check_inputs(startdate, enddate, experiment, flux_type,
    json_type=None, json_mode='', is_diff_thresh=[], subroutine=None):
    """ Check that all of the user inputs make sense and fall within bounds.
        
        INPUTS:
        
        :startdate: (datetime) - start of time period entered by user
        :enddate: (datetime) - end of time period entered by user
        :experiment: (string) - name of experiment specifed by user
        :flux_type: (string) - integral or differential
        :json_type: (string) - model or observations, only needed if "user" experiment
        :json_mode: (string) - measurement (for observations), or values like forecast,
            historical, simulated_realtime_forecast for model
        :is_diff_thresh: (bool 1xn array) - where n indicates the number of
            thresholds input by the user, e.g. "30,1;50,1" n=2
            Indicates if the user-input thresholds apply to integral or
            differential channels
            
        OUTPUTS:
        
        None, but system exit if error found
        
    """
    #CHECKS ON INPUTS
    if (enddate < startdate):
        sys.exit('End time before start time! Enter a valid date range. '
                'Exiting.')
 
    if subroutine == 'idsep':
        dt = enddate - startdate
        if (dt.days < cfg.init_win):
            print(f'Date range from {startdate.date()} to {enddate.date()} ({dt.days} days) '
                'is less than the '
                f'length of the background subtraction window, init_win={cfg.init_win} days. '
                'fetchsep is extending the time frame automatically to run and '
                'trimming results to requested time frame at the end. Continuing.')
 
 
    if experiment != "user":
        exp_info = expts.experiment_info(experiment)
        if len(exp_info['flux_type']) > 1 and flux_type == "":
            sys.exit(f"User must indicate whether input flux is {exp_info['flux_type']}. Exiting.")

        if flux_type != "" and flux_type not in exp_info['flux_type']:
             sys.exit(f"User must specify flux type for {experiment} from the choices: {exp_info['flux_type']}. Exiting.")

        if exp_info['last_date'] != None:
            date = exp_info['last_date']+datetime.timedelta(hours=24)
            exclusive_end_date = datetime.datetime(date.year, date.month, date.day)
            if startdate < exp_info['first_date'] or enddate > exclusive_end_date:
                sys.exit(f"The {experiment} data is available from {exp_info['first_date']} to {exp_info['last_date']}. Please change your requested dates. Exiting.")
        else:
            if startdate < exp_info['first_date']:
                sys.exit(f"The {experiment} data is available from {exp_info['first_date']} to present. Please change your requested dates. Exiting.")

    if experiment == "GOES-RT" and flux_type == "integral":
        print('Using GOES primary satellite real time fluxes as provided by SWPC in their 3-day jsons '
            'and archived by CCMC. Available starting 2010-04-14.')

    if experiment == "GOES-RT" and flux_type == "differential":
        print('Using GOES primary satellite real time fluxes as provided by SWPC in their 7-day jsons.')

    goes_R = expts.goes_R()
    goes16_integral_stdate = datetime.datetime(2020,3,8)
    if(experiment in goes_R) and startdate < goes16_integral_stdate:
        if startdate >= datetime.datetime(2017,9,1) and \
            startdate <= datetime.datetime(2017,9,20):
            print("error_check_inputs: Only special event data for September 2017 is available for GOES-16.")
        else:
            sys.exit('The GOES-R real time integral fluxes are only available '
                    + 'starting on '+ str(goes16_integral_stdate) +
                '. Please change your requested dates and use GOES-RT for the experiment. Exiting.')
    elif (experiment in goes_R) and flux_type == "integral":
        #UNTIL NOAA PROVIDES A SUPPORTED INTEGRAL PRODUCT
        sys.exit('Note: The GOES-R integral fluxes are real-time fluxes archived at CCMC. '
            'When NOAA\'s official L2 integral fluxes become available, they will be included in FetchSEP. '
            'Please specify GOES-RT for --Experiment to use GOES-R integral fluxes.')
  
  


def error_check_background(experiment, flux_type, doBGSubOPSEP, doBGSubIDSEP,
    OPSEPEnhancement, IDSEPEnhancement):
    """ Check background separatation and background subtraction options.
    
    """
    if doBGSubOPSEP and doBGSubIDSEP:
        sys.exit("You chose background subtraction applied by both OPSEP and "
                "IDSEP. Please select only one method to perform background "
                "subtraction by selecting only --doBGSubOPSEP or --doBGSubIDSEP.")
                
    if doBGSubOPSEP and IDSEPEnhancement:
        sys.exit("You chose to perform background subtraction with OPSEP "
            "and to identify enhancements above background using IDSEP. "
            "Please choose only one method to perform background "
            "subtraction and identify enhancements.")
            
    if doBGSubIDSEP and OPSEPEnhancement:
        sys.exit("You chose to perform background subtraction with IDSEP "
            "and to identify enhancements above background using OPSEP. "
            "Please choose only one method to perform background "
            "subtraction and identify enhancements.")

    if (doBGSubOPSEP or doBGSubIDSEP) and experiment[0:4] == "GOES" and flux_type == "integral":
        print("WARNING: Do you want to perform background subtraction? "
                "Do not perform background subtraction on GOES integral "
                "fluxes. Integral fluxes have already been derived by "
                "applying corrections for cross-contamination and removing "
                "the instrument background levels.")
