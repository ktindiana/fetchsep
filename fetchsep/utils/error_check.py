import datetime
import os
import sys

""" Check inputs and options for errors and incompatabilities."""


def error_check_options(experiment, flux_type, options, doBGSub, subroutine=None):
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
    if doBGSub and experiment[0:4] == "GOES" and flux_type == "integral":
        sys.exit("Do not perform background subtraction on GOES integral "
                "fluxes. Integral fluxes have already been derived by "
                "applying corrections for cross-contamination and removing "
                "the instrument background levels.")
    if ("uncorrected" in options or "S14" in options or "Bruno2017" in options)\
                and experiment[0:4] != "GOES":
        sys.exit("The options you have selected are only applicable to GOES "
                "data. Please remove these options and run again: "
                "uncorrected, S14, or Bruno2017.")



def error_check_inputs(startdate, enddate, experiment, flux_type,
    json_type=None, is_diff_thresh=None, subroutine=None):
    """ Check that all of the user inputs make sense and fall within bounds.
        
        INPUTS:
        
        :startdate: (datetime) - start of time period entered by user
        :enddate: (datetime) - end of time period entered by user
        :experiment: (string) - name of experiment specifed by user
        :flux_type: (string) - integral or differential
        :json_type: (string) - model or observations, only needed if "user" experiment
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
        if (dt.days < init_win):
            print(f'Date range from {startdate.date()} to {enddate.date()} ({dt.days} days) '
                'is less than the '
                f'length of the background subtraction window, init_win={init_win} days. '
                'fetchsep is extending the time frame automatically to run and '
                'trimming results to requested time frame at the end. Continuing.')
 
 
    if flux_type == "":
        sys.exit('User must indicate whether input flux is integral or '
                'differential. Exiting.')

    goes_R = ["GOES-16", "GOES-17", "GOES-18"]

    if ("SEPEM" in experiment and flux_type == "integral"):
        sys.exit('The SEPEM data sets only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    if (("EPHIN" in experiment) and flux_type == "integral"):
        sys.exit('The SOHO/EPHIN data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')
    
    #UNTIL NOAA PROVIDES A SUPPORTED INTEGRAL PRODUCT
    if ((experiment in goes_R) and flux_type == "integral"):
        sys.exit('Note: The GOES-R integral fluxes are only available as real time fluxes '
            'for the primary GOES spacecraft served by NOAA and archived at CCMC (as of 2024-09-06). '
            'When NOAA\'s official L2 fluxes, become available, they will be included in FetchSEP. '
            'Please specify GOES_RT for --Experiment to use GOES-R integral fluxes.')

    if experiment == "GOES_RT" and flux_type == "differential":
        sys.exit('GOES_RT real time fluxes are only available as integral fluxes. Change FluxType to '
            '\"integral\" or select a different GOES experiment.')

    if experiment == "GOES_RT":
        print('Using GOES primary satellite real time fluxes as provided by SWPC in their 3-day jsons '
            'and archived by CCMC. Available starting 2010-04-14.')

    if experiment == "user" and (json_type != "model" and json_type != "observations"):
        sys.exit('User experiments must specify a JSONType of \"model\" or '
            '\"observations\". Please specify your JSONType. Exiting.')

    for diff_thresh in is_diff_thresh:
        if diff_thresh and flux_type == "integral":
            sys.exit('The input flux type is specified as integral, but you '
                    'have requested a threshold in a differential energy bin. '
                    'Flux must be differential to impelement a threshold on a '
                    'differential energy bin. Exiting.')

    sepem_end_date = datetime.datetime(2015,12,31,23,55,00)
    if(experiment == "SEPEM" and (startdate > sepem_end_date or
                   enddate > sepem_end_date + datetime.timedelta(days=1))):
        sys.exit('The SEPEM (RSDv2) data set only extends to '
                  + str(sepem_end_date) +
            '. Please change your requested dates. Exiting.')

    sepemv3_end_date = datetime.datetime(2017,12,31,23,55,00)
    if(experiment == "SEPEMv3" and (startdate > sepemv3_end_date or
                   enddate > sepemv3_end_date + datetime.timedelta(days=1))):
        sys.exit('The SEPEM (RSDv3) data set only extends to '
                  + str(sepemv3_end_date) +
            '. Please change your requested dates. Exiting.')

    
    goes16_integral_stdate = datetime.datetime(2020,3,8)
    if(experiment in goes_R) and startdate < goes16_integral_stdate:
        if startdate >= datetime.datetime(2017,9,1) and \
            startdate <= datetime.datetime(2017,9,20):
            print("error_check_inputs: Only special event data for September 2017 is available for GOES-16.")
        else:
            sys.exit('The GOES-R real time integral fluxes are only available '
                    + 'starting on '+ str(goes16_integral_stdate) +
                '. Please change your requested dates and use GOES_RT for the experiment. Exiting.')
            
    stereoB_end_date = datetime.datetime(2014,9,27,16,26,00)
    if(experiment == "STEREO-B" and (startdate > stereoB_end_date or
                   enddate > stereoB_end_date + datetime.timedelta(days=1))):
        sys.exit('The STEREO-B data set only extends to '
                  + str(stereoB_end_date) +
            '. Please change your requested dates. Exiting.')


