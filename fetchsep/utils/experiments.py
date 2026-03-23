import sys
import os
import datetime
import math
from . import config as cfg

def valid_experiments():
    """ Return a list of experiments that can be processed by FetchSEP """

    valid_experiments = ['user', 'GOES', 'GOES-05','GOES-06', 'GOES-07', 'GOES-08',
        'GOES-10', 'GOES-11', 'GOES-12', 'GOES-13', 'GOES-14', 'GOES-15',
        'GOES-16', 'GOES-17','GOES-18','GOES-19','GOES-RT', 'GOES-SWPC',
        'SEPEM', 'SEPEMv3', 'EPHIN', 'EPHIN_REleASE','ERNE', 'CalGOES',
        'STEREO-A', 'STEREO-B', 'ACE_SIS', 'ACE_EPAM_electrons', 'IMP8_CPME']
        
    valid_nm = valid_neutron_monitors()
    valid_experiments = valid_experiments + valid_nm

    return valid_experiments


def valid_neutron_monitors():
    """ Return a list of neutron monitors that can be processed by FetchSEP """

    valid_nm = ['AATA', 'AATB', 'AHMD', 'APTY', 'ARNM', 'ATHN', 'BKSN', 'BUDA',
                'CALG', 'CALM', 'CHAC', 'CLMX', 'DJON', 'DOMB', 'DOMC', 'DRBS',
                'DRHM', 'ESOI', 'FSMT', 'HRMS', 'HUAN', 'ICRB', 'ICRO', 'INVK',
                'IRK2', 'IRK3', 'IRKT', 'JBGO', 'JUNG', 'JUNG1', 'KERG', 'KGSN',
                'KIEL', 'KIEL2', 'LMKS', 'MCMU', 'MCRL', 'MGDN', 'MOSC', 'MRNY',
                'MWSB', 'MWSN', 'MXCO', 'NAIN', 'NANM', 'NEU3', 'NEWK', 'NRLK',
                'NVBK', 'OULU', 'PSNM', 'PTFM', 'PWNK', 'ROME', 'SANB', 'SNAE',
                'SOPB', 'SOPO', 'TERA', 'THUL', 'TSMB', 'TXBY', 'UFSZ', 'YKTK',
                'ZUGS']

    return valid_nm


def goes_R():
    """ Return a list of GOES-R+ satellites """

    goes_R = ['GOES-16', 'GOES-17', 'GOES-18', 'GOES-19']

    return goes_R


def goes_sc():
    """ Return a list of middle GOES satellites """

    goes_sc = ['GOES-08', 'GOES-09','GOES-10','GOES-11',
                'GOES-12','GOES-13','GOES-14','GOES-15']
 
    return goes_sc


def old_goes_sc():
    """ Return a list of early GOES satellites """
    
    old_goes_sc = ['GOES-05', 'GOES-06', 'GOES-07']
    
    return old_goes_sc


def get_flux_type(experiment):
    """ If an experiment has only one possible flux type, set value
        automatically. If there is more than one option, the user must
        specify the flux type.
        
    """
    exp_info = experiment_info(experiment)

    if len(exp_info['flux_type']) == 1:
        return exp_info['flux_type'][0]
        
    else:
        sys.exit(f"You must specify a flux type from options {exp_info['flux_type']}. Exiting.")


def set_config_energy_units(experiment):
    """ Set global energy bin units """

    if experiment == "user":
        #Do nothing, units should already be set in config
        return
    else:
        exp_info = experiment_info(experiment)
        energy_units = exp_info['energy_units']
        cfg.set_energy_units(energy_units)


def set_config_flux_units(experiment):
    """ Set global flux and fluence units """

    if experiment == "user":
        #Do nothing, units should already be set in config
        return
    else:
        exp_info = experiment_info(experiment)
        flux_units_integral = exp_info['integral']['flux_units']
        fluence_units_integral = exp_info['integral']['fluence_units']
        flux_units_differential = exp_info['differential']['flux_units']
        fluence_units_differential = exp_info['differential']['fluence_units']
        cfg.set_flux_units(flux_units_integral, fluence_units_integral,
            flux_units_differential, fluence_units_differential)


def set_config_kurtosis_cut(experiment, flux_type):
    """ Set global kurtosis_cut values for idsep """
    
    exp_info = experiment_info(experiment)
    kurtosis_cut = exp_info[flux_type]['kurtosis_cut']
    cfg.set_kurtosis_cut(kurtosis_cut)
    

def get_spacecraft(experiment, spacecraft):
    """ Checks whether user-specified spacecraft is a valid for experiment.
        If experiment doesn't have a spacecraft option, set to 
        an empty string.

    """
 
    if experiment == "user":
        #allow user to specify whatever they want for this field
        return spacecraft
    else:
        exp_info = experiment_info(experiment)
        #no spacecraft field if none needed for experiment
        #overrides user-specified value
        if 'spacecraft' not in exp_info:
            return ''
 
        if 'spacecraft' in exp_info:
            sc = exp_info['spacecraft']
            if spacecraft in sc:
                return spacecraft
            else:
                sys.exit(f"{experiment} spacecraft must be selected from {sc}. Please correct and run again.")


def get_json_type(experiment):
    exp_info = experiment_info(experiment)
    json_type = exp_info['json_type']
    return json_type


def get_json_mode(experiment):
    exp_info = experiment_info(experiment)
    json_mode = exp_info['json_mode']
    return json_mode


#first_date and last_date indicate data availability.
#last_date = None indicates experiment continues to the present
#Energy bins listed for GOES are the nominal energy bins provided by NOAA

def calculate_geometric_means(energy_bins):
    """ Define the bin centers as the geometric means:
        center = sqrt(Elow*Ehigh)
        
        This approach is typically used to define the bin center of 
        proton experiments if not better known.
        
        This should only be applied to differential channels. 
        If a channel is an integral channel, will return the bottom
        energy of the channel.
        
    """
    centers = []
    for bin in energy_bins:
        if bin[1] == -1:
            centers.append(bin[0])
        else:
            centers.append(math.sqrt(bin[0]*bin[1]))
        
    return centers

########
#Fluxes can be converted back and forth between integral and differential by fetchsep.
#Corresponding integral and differential units are provided for all experiments
#except neutron monitors.
########
def experiment_info(experiment):
    """ Contains information about each experiment 'native' to fetchsep.
        Returns a dictionary with parameters for requested experiment.
    
    """

    experiments = {
            'user':{
                'first_date': None,
                'last_date': None,
                'flux_type': ['integral', 'differential'],
                'json_type': '',
                'json_mode': '',
                'species': '',
                'location': '',
                'differential': {
                    'kurtosis_cut': cfg.kurtosis_cut,
                },
                'integral': {
                    'kurtosis_cut': cfg.kurtosis_cut,
                },

            },

            'ACE_SIS':{
                'first_date': datetime.datetime(2001,8,7),#2001-08-07
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'kurtosis_cut': 999, #idsep for calculating mean background
                },
                'integral': {
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[30.0,-1],[60.0,-1]],
                    'energy_bin_centers': [30.0,60.0],
                    'url': 'https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/',
                    'kurtosis_cut': 999,
                }
            },

            'ACE_EPAM_electrons':{
                'first_date': datetime.datetime(2001,8,7),#2001-08-07
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'electron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[0.175,0.315]],
                    'energy_bin_centers': calculate_geometric_means([[0.175,0.315]]),
                    'url': 'https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },


            'CalGOES':{
                'first_date': datetime.datetime(1986,1,1),#1986-01-01
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'year',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[10.0,10.0],[15.8,15.8],[25.1,25.1],[39.8,39.8],
                        [65.1,65.1],[100.0,100.0],[158.5,158.5],[251.2,251.2],
                        [398.1,398.1],[630.9,630.9],[1000.0,1000.0]],
                    'energy_bin_centers': [10.0, 15.8, 25.1, 39.8, 65.1, 100.0, 158.5, 251.2,
                        398.1,630.9,1000.0],
                    'url': 'Request from Shaowen Hu of the NASA Space Radiation Analysis Group',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },


            'EPHIN':{
                'first_date': datetime.datetime(1995,12,8),#1995-12-08
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'year',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.3,7.8],[7.8,25.0],[25.0,40.9],[40.9,53.0]],
                    'energy_bin_centers': calculate_geometric_means([[4.3,7.8],[7.8,25.0],[25.0,40.9],[40.9,53.0]]),
                    'url': 'http://ulysses.physik.uni-kiel.de/costep/level3/l3i/',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },
            
            'EPHIN_HESPERIA':{
                'first_date': datetime.datetime(2015,10,4),#2015-10-04
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'variable',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.0,9.0],[9.0,15.8],[15.8,39.6],[20.0,35.5]],
                    'energy_bin_centers': calculate_geometric_means([[4.0,9.0],[9.0,15.8],[15.8,39.6],[20.0,35.5]]),
                    'url': 'https://hesperia.astro.noa.gr/data-retrieval-tool/', #download by hand
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            'EPHIN_REleASE':{
                'first_date': datetime.datetime(1995,12,8),#1995-12-08
                'last_date': datetime.datetime(2016,12,31),#2016-12-31
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'year',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[3.98, 4.47],[4.47, 5.01],[5.01, 5.62],[5.62, 6.31],
                        [6.31, 7.08], [7.08, 7.94], [7.94, 8.91], [8.91, 10.00], [10.00, 11.22],
                        [11.22, 12.59], [12.59, 14.13], [14.13, 15.85], [15.85, 17.78],
                        [17.78, 19.95], [19.95, 22.39], [22.39, 25.12], [25.12, 28.18],
                        [28.18, 31.62], [31.62, 35.48], [35.48, 39.81], [39.81, 44.67],
                        [44.67, 50.12]],
                    'energy_bin_centers': calculate_geometric_means([[3.98, 4.47],
                        [4.47, 5.01],[5.01, 5.62],[5.62, 6.31],
                        [6.31, 7.08], [7.08, 7.94], [7.94, 8.91], [8.91, 10.00],
                        [10.00, 11.22],
                        [11.22, 12.59], [12.59, 14.13], [14.13, 15.85], [15.85, 17.78],
                        [17.78, 19.95], [19.95, 22.39], [22.39, 25.12], [25.12, 28.18],
                        [28.18, 31.62], [31.62, 35.48], [35.48, 39.81], [39.81, 44.67],
                        [44.67, 50.12]]),
                    'url': 'https://zenodo.org/records/14191918', #download by hand
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            'ERNE':{
                'info': 'The highest energy 3 channels of ERNE tend to saturate and show incorrect values during high intensity SEP events. For this reason, only the >10 MeV integral fluxes should be used from ERNE data during themost intense part of intense SEP events.',
                'first_date': datetime.datetime(1996,5,7),#1996-05-07
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'variable',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [], #Changes with time, see ERNEf10, ERNEf40, ERNEf50
                    'url': 'https://export.srl.utu.fi',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            #SUPPORTING INFORMATION FOR ERNE
            #f10 format, from 2 December 1996
            #f10     2 Dec 1996     Original launch format
            'ERNEf10':{
                'info': 'original launch format and energy bin definitions',
                'first_date': datetime.datetime(1996,5,7),#1996-05-07
                'last_date': datetime.datetime(2000,4,19), #2000-04-19
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'variable',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[1.5,1.8],[1.8,2.2],[2.2,2.7],[2.7,3.3],[3.3,4.1],
                       [4.1,5.1],[5.1,6.4],[6.4,8.1],[8.1,10],
                       [10.0,13.0],[14.0,17.0],[17.0,22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[41.0,51.0],
                       [51.0,67.0],[54.0,79.0],[79.0,114.0],[111.0,140.0]],
                    'energy_bin_centers': [1.6, 2.0, 2.4, 3.0, 3.7, 4.6, 5.8, 7.2, 9.1, 11.0,
                        15.4, 18.9, 23.3, 29.0, 36.4, 45.6, 54.1, 67.5, 95.0, 116.0],
                    'url': 'https://export.srl.utu.fi',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            #SUPPORTING INFORMATION FOR ERNE
            #f40brk format, from 21 Nov 2000
            #f40brk 21 Nov 2000     HED S1X H2 E-amplifier breakdown at 00:15:44.833
            #NOTE: all the f40brk HED data are left out from this data set. The
            #corresponding files are provided but are empty.
            'ERNEf40':{
                'info': 'Major update of the on-board program. From 21 Nov 2000, HED S1X H2 E-amplifier breakdown at 00:15:44.833. All HED data are left out of this dataset.',
                'first_date': datetime.datetime(2000,4,20), #2000-04-20
                'last_date': datetime.datetime(2001,7,2), #2001-07-02
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'variable',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[1.6,1.8],[1.8,2.2],[2.2,2.7],[2.7,3.3],[3.3,4.1],
                       [4.1,5.1],[5.1,6.4],[6.4,8.1],[8.1,10],
                       [10.0,13.0],[14.0,17.0],[17.0, 22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[40.0,51.0],[51.0,67.0],
                       [64.0,80.0],[80.0,101.0],[101.0,131.0]],
                    'energy_bin_centers': [1.7, 2.0, 2.4, 3.0, 3.7, 4.7, 5.7, 7.2, 9.1, 11.0,
                        15.4, 18.9, 23.3, 29.1, 36.4, 45.6, 57.4, 72.0, 90.5, 108.0],
                    'url': 'https://export.srl.utu.fi',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            #SUPPORTING INFORMATION FOR ERNE
            #f50 format, from 3 Jul 2001
            #f50     3 Jul 2001     On-board SW updated to (partially) fix the HED S1X
            #               E-amplifier breakdown.
            'ERNEf50':{
                'info': 'On-board SW updated to (partially) fix the HED S1X E-amplifier breakdown. ',
                'first_date': datetime.datetime(2001,7,3), #2000-04-20
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'variable',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[1.6,1.8],[1.8,2.2],[2.2,2.7],[2.7,3.3],[3.3,4.1],
                       [4.1,5.1],[5.1,6.4],[6.4,8.1],[8.1,10],
                       [10.0,13.0],[14.0,17.0],[17.0, 22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[40.0,51.0],[51.0,67.0],
                       [64.0,80.0],[80.0,101.0],[101.0,131.0]],
                    'energy_bin_centers': [1.7, 2.0, 2.4, 3.0, 3.7, 4.7, 5.7, 7.2, 9.1, 11,
                        15.4, 18.9, 23.3, 29.1, 36.4, 45.6, 57.4, 72.0, 90.5, 108.0],
                    'url': 'https://export.srl.utu.fi',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            'GOES':{
                'first_date': datetime.datetime(1984,1,1),#1984-01-01 GOES-05
                'last_date': datetime.datetime(2020,3,4), #2020-03-04, GOES-15
                'flux_type': ['integral'],
                'spacecraft': ['primary', 'secondary'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'various',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'url': 'various',
                    'kurtosis_cut': 50,
                },
                'integral': {
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins':[[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'various',
                    'kurtosis_cut': 65,
                }
            },


            #GOES-05 has very irregular coverage of particle measurements
            'GOES-05':{
                'first_date': datetime.datetime(1984,1,1),#1984-01-01
                'last_date': datetime.datetime(1985,12,31), #1985-12-31
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.2,8.7],[8.7,14.5],[15.0,44.0],[39.0,82.0],[84.0,200.0],[110.0,500.0]],
                    'energy_bin_centers': calculate_geometric_means([[4.2,8.7],[8.7,14.5],[15.0,44.0],[39.0,82.0],[84.0,200.0],[110.0,500.0]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-06':{
                'first_date': datetime.datetime(1986,1,1),#1986-01-01
                'last_date': datetime.datetime(1994,11,30), #1994-11-30
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.2,8.7],[8.7,14.5],[15.0,44.0],
                        [39.0,82.0],[84.0,200.0],[110.0,500.0],
                        [375.0, 375.0],[465.0,465.0],[605.0,605.0], [685.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.2,8.7],[8.7,14.5],[15.0,44.0],
                        [39.0,82.0],[84.0,200.0],[110.0,500.0],
                        [375.0, 375.0],[465.0,465.0],[605.0,605.0], [685.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
                'integral': {
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1], [685,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,685.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                },
            },

            'GOES-07':{
                'first_date': datetime.datetime(1987,3,1),#1987-03-01
                'last_date': datetime.datetime(1996,8,31), #1996-08-31
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.2,8.7],[8.7,14.5],[15.0,44.0],
                        [39.0,82.0],[84.0,200.0],[110.0,500.0]],
                    'energy_bin_centers': calculate_geometric_means([[4.2,8.7],[8.7,14.5],[15.0,44.0],
                        [39.0,82.0],[84.0,200.0],[110.0,500.0]]),
                    'url':'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0],
                    'url':'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
            },


            'GOES-08':{
                'first_date': datetime.datetime(1995,1,1),#1995-01-01
                'last_date': datetime.datetime(2003,6,30), #2003-05-31
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.0,9.0],[9.0,15.0],[15.0,44.0],
                               [40.0,80.0],[80.0,165.0],[165.0,500.0],
                               [350.0,420.0],[420.0,510.0],[510.0,700.0],
                               [700.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.0,9.0],[9.0,15.0],[15.0,44.0],
                               [40.0,80.0],[80.0,165.0],[165.0,500.0],
                               [350.0,420.0],[420.0,510.0],[510.0,700.0],
                               [700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
                'integral': {
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins':[[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                }
            },


            'GOES-09':{
                'first_date': datetime.datetime(1996,4,1),#1996-04-01
                'last_date': datetime.datetime(1998,7,31), #1998-07-31
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0], [350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0], [350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url':'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-10':{
                'first_date': datetime.datetime(1998,7,1),#1998-07-01
                'last_date': datetime.datetime(2004,6,30), #2004-06-30
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0],[350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0],[350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-11':{
                'first_date': datetime.datetime(2003,6,1),#2003-06-01
                'last_date': datetime.datetime(2011,2,28), #2011-02-28
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0],[350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0],[350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-12':{
                'first_date': datetime.datetime(2003,1,1),#2003-01-01
                'last_date': datetime.datetime(2010,8,31), #2010-08-31
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0],[350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.0,9.0],[9.0,15.0],[15.0,44.0],
                        [40.0,80.0],[80.0,165.0],[165.0,500.0],[350.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 50,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                    }
            },

            'GOES-13':{
                'first_date': datetime.datetime(2010,5,1),#2010-05-01
                'last_date': datetime.datetime(2017,12,31), #2017-12-14
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.2,8.7],[8.7,14.5],[15.0,40.0],
                        [38.0,82.0],[84.0,200.0],[110.0,900.0],[330.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.2,8.7],[8.7,14.5],[15.0,40.0],
                        [38.0,82.0],[84.0,200.0],[110.0,900.0],[330.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-14':{
                'first_date': datetime.datetime(2009,7,1),#2009-07-01
                'last_date': datetime.datetime(2020,3,4), #2020-03-04
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.2,8.7],[8.7,14.5],[15.0,40.0],
                        [38.0,82.0],[84.0,200.0],[110.0,900.0],[330.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]],
                    'energy_bin_centers': calculate_geometric_means([[4.2,8.7],[8.7,14.5],[15.0,40.0],
                        [38.0,82.0],[84.0,200.0],[110.0,900.0],[330.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-15':{
                'first_date': datetime.datetime(2011,1,1),#2011-01-01
                'last_date': datetime.datetime(2020,3,4), #2020-03-04
                'flux_type': ['integral', 'differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'month',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.2,8.7],[8.7,14.5],[15.0,40.0],
                        [38.0,82.0],[84.0,200.0],[110.0,900.0],[330.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]],
                     'energy_bin_centers': calculate_geometric_means([[4.2,8.7],[8.7,14.5],[15.0,40.0],
                        [38.0,82.0],[84.0,200.0],[110.0,900.0],[330.0,420.0],
                        [420.0,510.0],[510.0,700.0],[700.0,-1]]),
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,60.0,100.0,700.0],
                    'url': 'https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-16':{
                'first_date': datetime.datetime(2020,11,1),#2020-11-01
                'last_date': datetime.datetime(2025,4,6), #2025-04-06, there is a file on the 7th, but has a problem
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[3.4,6.48],[5.84,11.0],[11.64,23.27],[24.9,38.1],
                        [40.3,73.4],[83.7,98.5],[99.9,118.0],[115.0,143.0],[160.0,242.0],
                        [276.0,404.0],[500.0,-1]], #[1.02,1.86],[1.9,2.3],[2.31,3.34],
                     'energy_bin_centers': calculate_geometric_means([[3.4,6.48],[5.84,11.0],[11.64,23.27],[24.9,38.1],
                        [40.3,73.4],[83.7,98.5],[99.9,118.0],[115.0,143.0],[160.0,242.0],
                        [276.0,404.0],[500.0,-1]]),
                    'url': 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/sgps-l2-avg5m/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-17':{
                'first_date': datetime.datetime(2020,12,1),#2020-12-01
                'last_date': datetime.datetime(2023,3,14), #2023-03-14
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[3.4,6.48],[5.84,11.0],[11.64,23.27],[23.9,32.6],
                        [40.7,68.2],[83.9,98.4],[99.7,118.0],[123.0,148.0],[156.0,237.0],
                        [267.0,390.0],[500.0,-1]], #[1.02,1.86],[1.9,2.3],[2.31,3.34],
                     'energy_bin_centers': calculate_geometric_means([[3.4,6.48],[5.84,11.0],[11.64,23.27],[23.9,32.6],
                        [40.7,68.2],[83.9,98.4],[99.7,118.0],[123.0,148.0],[156.0,237.0],
                        [267.0,390.0],[500.0,-1]]),
                    'url': 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes17/l2/data/sgps-l2-avg5m/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-18':{
                'first_date': datetime.datetime(2022,9,13),#2022-09-13
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[3.4,6.48],[5.84,11.0],[11.64,23.27],[25.5,38.4],
                        [41.0,77.0],[80.9,97.6],[96.3,118.4],[114.88,138.4],[153.3,229.3],
                        [267.0,390.0],[500.0,-1]], #[1.02,1.86],[1.9,2.3],[2.31,3.34],
                     'energy_bin_centers': calculate_geometric_means([[3.4,6.48],[5.84,11.0],[11.64,23.27],[25.5,38.4],
                        [41.0,77.0],[80.9,97.6],[96.3,118.4],[114.88,138.4],[153.3,229.3],
                        [267.0,390.0],[500.0,-1]]),
                    'url': 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes18/l2/data/sgps-l2-avg5m/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-19':{
                'first_date': datetime.datetime(2024,8,22),#2024-08-22
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[3.4,6.48],[5.84,11.0],[11.64,23.27],[25.9,35.2],
                        [41.0,74.0],[78.0,100.7],[97.9,120.6],[114.6,142.4],[150.7,231.5],
                        [267.0,390.0],[500.0,-1]], #[1.02,1.86],[1.9,2.3],[2.31,3.34],
                     'energy_bin_centers': calculate_geometric_means([[3.4,6.48],[5.84,11.0],[11.64,23.27],[25.9,35.2],
                        [41.0,74.0],[78.0,100.7],[97.9,120.6],[114.6,142.4],[150.7,231.5],
                        [267.0,390.0],[500.0,-1]]),
                    'url': 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes19/l2/data/sgps-l2-avg5m/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-RT':{
                'first_date': datetime.datetime(2010,4,14),#2010-04-14
                'last_date': None,
                'flux_type': ['integral'],
                'spacecraft': ['primary','secondary'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[100.0,-1],[60.0,-1],[500.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,100.0,60.0,500.0],
                    'url': 'https://iswa.gsfc.nasa.gov/IswaSystemWebApp/hapi/',
                    'kurtosis_cut': 65,
                }
            },

            'GOES-SWPC':{
                'first_date': datetime.datetime.now() - datetime.timedelta(days=7),#previous 7 days
                'last_date': None,
                'flux_type': ['integral', 'differential'],
                'spacecraft': ['primary','secondary'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential':{
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[3.4,6.48],[5.84,11.00],[11.64,23.27],[25.90, 38.10],[40.30,73.40],
                        [83.70,98.50],[99.9,118.0],[115.0,143.0],[160.0,242.0],[276.0,404.0]],
                     'energy_bin_centers': calculate_geometric_means([[3.4,6.48],[5.84,11.00],[11.64,23.27],
                        [25.90, 38.10],[40.30,73.40],
                        [83.70,98.50],[99.9,118.0],[115.0,143.0],[160.0,242.0],[276.0,404.0]]),
                    'url': 'https://services.swpc.noaa.gov/json/goes/',
                    'kurtosis_cut': 65,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'energy_bins': [[5.0,-1],[10.0,-1],[30.0,-1],[50.0,-1],[100.0,-1],[60.0,-1],[500.0,-1]],
                    'energy_bin_centers': [5.0,10.0,30.0,50.0,100.0,60.0,500.0],
                    'url': 'https://services.swpc.noaa.gov/json/goes/',
                    'kurtosis_cut': 65,
                }
            },


            'IMP8_CPME':{
                'first_date': datetime.datetime(1974,2,17),#1974-02-17
                'last_date': datetime.datetime(2001,11,7), #2001-11-07
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'variable',
                'resolution': datetime.timedelta(seconds=330),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[4.60, 15.0], [15.0, 25.0],[25.0, 48.0], [48.0, 96.0],[96.0, 145.0], [145.0, 440.0]],
                    'energy_bin_centers': calculate_geometric_means([[4.60, 15.0], [15.0, 25.0],[25.0, 48.0], [48.0, 96.0],[96.0, 145.0], [145.0, 440.0]]),
                    'url': 'http://sd-www.jhuapl.edu/IMP/data/imp8/cpme/cpme_330s/protons/',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            'SEPEM':{
                'first_date': datetime.datetime(1974,7,1),#1974-07-01
                'last_date': datetime.datetime(2015,12,31), #2015-12-31
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'year', #after processing by FetchSEP
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                           [21.87,31.62],[31.62,45.73],[45.73,66.13],
                           [66.13,95.64],[95.64,138.3],[138.3,200.0],
                           [200.0,289.2]],
                    'energy_bin_centers': calculate_geometric_means([[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                           [21.87,31.62],[31.62,45.73],[45.73,66.13],
                           [66.13,95.64],[95.64,138.3],[138.3,200.0],
                           [200.0,289.2]]),
                    'url': 'http://sepem.eu/help/SEPEM_RDS_v2-00.zip',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            'SEPEMv3':{
                'first_date': datetime.datetime(1974,7,1),#1974-07-01
                'last_date': datetime.datetime(2017,12,31), #2017-12-31
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'earth',
                'cadence': 'year', #after processeing by FetchSEP
                'resolution': datetime.timedelta(minutes=5),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                           [21.87,31.62],[31.62,45.73],[45.73,66.13],
                           [66.13,95.64],[95.64,138.3],[138.3,200.0],
                           [200.0,289.2],[289.2,418.3],[418.3,604.9],
                           [604.9,874.7]],
                    'energy_bin_centers': calculate_geometric_means([[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                           [21.87,31.62],[31.62,45.73],[45.73,66.13],
                           [66.13,95.64],[95.64,138.3],[138.3,200.0],
                           [200.0,289.2],[289.2,418.3],[418.3,604.9],
                           [604.9,874.7]]),
                    'url': 'http://www.sepem.eu/help/SEPEM_RDS_v3.zip',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },
 
            'STEREO-A':{
                'first_date': datetime.datetime(2006,12,1), #2006-12-01
                'last_date': None,
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'stereoa',
                'cadence': 'multiple', #LET and HET package files differently
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[1.8,3.6],[4.0,6.0],[6.0,10.0],[13.6,15.1],
                        [14.9,17.1],[17.0,19.3],[20.8,23.8],
                        [23.8,26.4],[26.3,29.7],[29.5,33.4],[33.4,35.8],
                        [35.5,40.5],[40.0,60.0],[60.0,100.0]],
                    'energy_bin_centers': calculate_geometric_means([[1.8,3.6],[4.0,6.0],[6.0,10.0],[13.6,15.1],
                        [14.9,17.1],[17.0,19.3],[20.8,23.8],
                        [23.8,26.4],[26.3,29.7],[29.5,33.4],[33.4,35.8],
                        [35.5,40.5],[40.0,60.0],[60.0,100.0]]),
                    'url': 'https://izw1.caltech.edu/STEREO/DATA/',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

            'STEREO-B':{
                'first_date': datetime.datetime(2006,12,1), #2006-12-01
                'last_date': datetime.datetime(2014,9,27,16,26,0), #2014-09-27 16:26:00
                'flux_type': ['differential'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'proton',
                'location': 'stereob',
                'cadence': 'multiple', #LET and HET package files differently
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'MeV',
                'differential': {
                    'flux_units': 'MeV^-1*cm^-2*s^-1*sr^-1',
                    'fluence_units': 'MeV^-1*cm^-2',
                    'energy_bins': [[1.8,3.6],[4.0,6.0],[6.0,10.0],[13.6,15.1],
                        [14.9,17.1],[17.0,19.3],[20.8,23.8],
                        [23.8,26.4],[26.3,29.7],[29.5,33.4],[33.4,35.8],
                        [35.5,40.5],[40.0,60.0],[60.0,100.0]],
                    'energy_bin_centers': calculate_geometric_means([[1.8,3.6],[4.0,6.0],[6.0,10.0],[13.6,15.1],
                        [14.9,17.1],[17.0,19.3],[20.8,23.8],
                        [23.8,26.4],[26.3,29.7],[29.5,33.4],[33.4,35.8],
                        [35.5,40.5],[40.0,60.0],[60.0,100.0]]),
                    'url': 'https://izw1.caltech.edu/STEREO/DATA/',
                    'kurtosis_cut': 999,
                },
                'integral':{
                    'flux_units': 'pfu',
                    'fluence_units': 'cm^-2',
                    'kurtosis_cut': 999,
                }
            },

 
            ####NEUTRON MONITORS
            'BKSN': {
                'info': 'Baksan (R=5.70, Alt=1700 m)',
                'first_date': datetime.datetime(2009,5,18),#2009-05-18
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day', #after processeing by FetchSEP
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential':{
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[5.70,-1]],
                    'energy_bin_centers': [5.70],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            
            },

            'OULU': {
                'info': 'Oulu (R=0.81, Alt=15 m)',
                'first_date': datetime.datetime(1964,4,1),#1964-04-01
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day', #after processeing by FetchSEP
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential':{
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.81,-1]],
                    'energy_bin_centers': [0.81],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'SOPO': {
                'info': 'South Pole (R=0.10, Alt=2820 m)',
                'first_date': datetime.datetime(1964,3,1),#1964-04-01
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day', #after processeing by FetchSEP
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential':{
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.1,-1]],
                    'energy_bin_centers': [0.1],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },
            
            ##########TO UPDATE WITH START AND END TIMES ##############
            'AATA': {
                'info': 'Alma-Ata A (R=5.90, Alt=897 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[5.9,-1]],
                    'energy_bin_centers': [5.9],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'AATB': {
                'info': 'Alma-Ata B (R=5.90, Alt=3340 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[5.9,-1]],
                    'energy_bin_centers': [5.9],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'AHMD': {
                'info': 'Ahmedabad (R=15.94, Alt=50 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[15.94,-1]],
                    'energy_bin_centers': [15.94],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'APTY': {
                'info': 'Apatity (R=0.65, Alt=181 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.65,-1]],
                    'energy_bin_centers': [0.65],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'ARNM': {
                'info': 'Aragats (R=7.10, Alt=3200 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[7.1,-1]],
                    'energy_bin_centers': [7.1],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'ATHN': {
                'info': 'Athens (R=8.53, Alt=260 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[8.53,-1]],
                    'energy_bin_centers': [8.53],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'BKSN': {
                'info': 'Baksan (R=5.70, Alt=1700 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[5.7,-1]],
                    'energy_bin_centers': [5.7],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'BUDA': {
                'info': 'BUDA (R=, Alt= m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[None,-1]],
                    'energy_bin_centers': [None],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'CALG': {
                'info': 'Calgary (R=1.08, Alt=1123 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[1.08,-1]],
                    'energy_bin_centers': [1.08],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'CALM': {
                'info': 'NM de Castilla la Mancha (R=6.95, Alt=708 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[6.95,-1]],
                    'energy_bin_centers': [6.95],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'CHAC': {
                'info': 'CHACALTAYA (R=11.80, Alt=5270 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[11.8,-1]],
                    'energy_bin_centers': [11.8],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'CLMX': {
                'info': 'Climax (R=3.00, Alt=3400 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[3.0,-1]],
                    'energy_bin_centers': [3.0],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'DJON': {
                'info': 'Daejeon (R=11.20, Alt=200 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[11.2,-1]],
                    'energy_bin_centers': [11.2],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'DOMB': {
                'info': 'Dome C mini NM (bare) (R=0.01, Alt=3233 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.01,-1]],
                    'energy_bin_centers': [0.01],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'DOMC': {
                'info': 'Dome C mini NM (R=0.01, Alt=3233 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.01,-1]],
                    'energy_bin_centers': [0.01],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'DRBS': {
                'info': 'Dourbes (R=3.18, Alt=225 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[3.18,-1]],
                    'energy_bin_centers': [3.18],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'DRHM': {
                'info': 'Durham (R=2.21, Alt=20 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.21,-1]],
                    'energy_bin_centers': [2.21],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'ESOI': {
                'info': 'Emilio Segre Obs. Israel (R=10.75, Alt=2055 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[10.75,-1]],
                    'energy_bin_centers': [10.75],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'FSMT': {
                'info': 'Fort Smith (R=0.30, Alt=180 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.3,-1]],
                    'energy_bin_centers': [0.3],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'HRMS': {
                'info': 'Hermanus (R=4.58, Alt=26 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[4.58,-1]],
                    'energy_bin_centers': [4.58],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'HUAN': {
                'info': 'Huancayo (R=12.92, Alt=3400 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[12.92,-1]],
                    'energy_bin_centers': [12.92],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'ICRB': {
                'info': 'Izana Cosmic Ray Observatory (bare) (R=11.50, Alt=2390 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[11.5,-1]],
                    'energy_bin_centers': [11.5],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'ICRO': {
                'info': 'Izana Cosmic Ray Observatory (R=11.50, Alt=2390 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[11.5,-1]],
                    'energy_bin_centers': [11.5],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'INVK': {
                'info': 'Inuvik (R=0.30, Alt=21 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.3,-1]],
                    'energy_bin_centers': [0.3],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'IRK2': {
                'info': 'Irkustk 2 (R=3.64, Alt=2000 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[3.64,-1]],
                    'energy_bin_centers': [3.64],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'IRK3': {
                'info': 'Irkutsk 3 (R=3.64, Alt=3000 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[3.64,-1]],
                    'energy_bin_centers': [3.64],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'IRKT': {
                'info': 'Irkustk (R=3.64, Alt=435 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[3.64,-1]],
                    'energy_bin_centers': [3.64],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'JBGO': {
                'info': 'JangBogo (R=0.30, Alt=29 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.3,-1]],
                    'energy_bin_centers': [0.3],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'JUNG': {
                'info': 'IGY Jungfraujoch (R=4.49, Alt=3570 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[4.49,-1]],
                    'energy_bin_centers': [4.49],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'JUNG1': {
                'info': 'NM64 Jungfraujoch (R=4.49, Alt=3475 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[4.49,-1]],
                    'energy_bin_centers': [4.49],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'KERG': {
                'info': 'Kerguelen (R=1.14, Alt=33 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[1.14,-1]],
                    'energy_bin_centers': [1.14],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'KGSN': {
                'info': 'Kingston (R=1.88, Alt=65 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[1.88,-1]],
                    'energy_bin_centers': [1.88],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'KIEL': {
                'info': 'Kiel (R=2.36, Alt=54 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.36,-1]],
                    'energy_bin_centers': [2.36],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'KIEL2': {
                'info': 'KielRT (R=2.36, Alt=54 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.36,-1]],
                    'energy_bin_centers': [2.36],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'LMKS': {
                'info': 'Lomnicky Stit (R=3.84, Alt=2634 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[3.84,-1]],
                    'energy_bin_centers': [3.84],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MCMU': {
                'info': 'Mc Murdo (R=0.30, Alt=48 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.3,-1]],
                    'energy_bin_centers': [0.3],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MCRL': {
                'info': 'Mobile Cosmic Ray Laboratory (R=2.46, Alt=200 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.46,-1]],
                    'energy_bin_centers': [2.46],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MGDN': {
                'info': 'Magadan (R=2.10, Alt=220 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.1,-1]],
                    'energy_bin_centers': [2.1],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MOSC': {
                'info': 'Moscow (R=2.43, Alt=200 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.43,-1]],
                    'energy_bin_centers': [2.43],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MRNY': {
                'info': 'Mirny (R=0.03, Alt=30 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.03,-1]],
                    'energy_bin_centers': [0.03],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MWSB': {
                'info': 'Mawson Bare (R=0.22, Alt=30 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.22,-1]],
                    'energy_bin_centers': [0.22],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MWSN': {
                'info': 'Mawson (R=0.22, Alt=30 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.22,-1]],
                    'energy_bin_centers': [0.22],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'MXCO': {
                'info': 'Mexico (R=8.28, Alt=2274 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[8.28,-1]],
                    'energy_bin_centers': [8.28],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'NAIN': {
                'info': 'Nain (R=0.30, Alt=46 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.3,-1]],
                    'energy_bin_centers': [0.3],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'NANM': {
                'info': 'Nor-Amberd (R=7.10, Alt=2000 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[7.1,-1]],
                    'energy_bin_centers': [7.1],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'NEU3': {
                'info': 'Neumayer III mini neutron monitor (R=0.10, Alt=40 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.1,-1]],
                    'energy_bin_centers': [0.1],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'NEWK': {
                'info': 'Newark (R=2.40, Alt=50 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.4,-1]],
                    'energy_bin_centers': [2.4],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'NRLK': {
                'info': 'Norilsk (R=0.63, Alt=0 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.63,-1]],
                    'energy_bin_centers': [0.63],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'NVBK': {
                'info': 'Novosibirsk (R=2.91, Alt=163 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[2.91,-1]],
                    'energy_bin_centers': [2.91],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'PSNM': {
                'info': 'Doi Inthanon (Princess Sirindhorn NM) (R=16.80, Alt=2565 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[16.8,-1]],
                    'energy_bin_centers': [16.8],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'PTFM': {
                'info': 'Potchefstroom (R=6.98, Alt=1351 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[6.98,-1]],
                    'energy_bin_centers': [6.98],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'PWNK': {
                'info': 'Peawanuck (R=0.30, Alt=53 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.3,-1]],
                    'energy_bin_centers': [0.3],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'ROME': {
                'info': 'Rome (R=6.27, Alt=0 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[6.27,-1]],
                    'energy_bin_centers': [6.27],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'SANB': {
                'info': 'Sanae D (R=0.73, Alt=52 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.73,-1]],
                    'energy_bin_centers': [0.73],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'SNAE': {
                'info': 'Sanae IV (R=0.73, Alt=856 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.73,-1]],
                    'energy_bin_centers': [0.73],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'SOPB': {
                'info': 'South Pole Bare (R=0.10, Alt=2820 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.1,-1]],
                    'energy_bin_centers': [0.1],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'TERA': {
                'info': 'Terre Adelie (R=0.01, Alt=32 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.01,-1]],
                    'energy_bin_centers': [0.01],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'THUL': {
                'info': 'Thule (R=0.30, Alt=26 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.3,-1]],
                    'energy_bin_centers': [0.3],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'TSMB': {
                'info': 'Tsumeb (R=9.15, Alt=1240 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[9.15,-1]],
                    'energy_bin_centers': [9.15],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'TXBY': {
                'info': 'Tixie Bay (R=0.48, Alt=0 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[0.48,-1]],
                    'energy_bin_centers': [0.48],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'UFSZ': {
                'info': 'Zugspitze (R=4.10, Alt=2650 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[4.1,-1]],
                    'energy_bin_centers': [4.1],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'YKTK': {
                'info': 'Yakutsk (R=1.65, Alt=105 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[1.65,-1]],
                    'energy_bin_centers': [1.65],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },

            'ZUGS': {
                'info': 'Zugspitze (R=4.24, Alt=2960 m)',
                'first_date': datetime.datetime(1940,1,1),
                'last_date': None,
                'flux_type': ['integral'],
                'json_type': 'observations',
                'json_mode': 'measurement',
                'species': 'neutron',
                'location': 'earth',
                'cadence': 'day',
                'resolution': datetime.timedelta(minutes=1),
                'energy_units': 'GV',
                'differential': {
                    'flux_units': None,
                    'fluence_units': None,
                    'kurtosis_cut': 999,
                },
                'integral': {
                    'flux_units': 'counts*s^-1',
                    'fluence_units': 'counts',
                    'energy_bins': [[4.24,-1]],
                    'energy_bin_centers': [4.24],
                    'url': 'https://www.nmdb.eu/nest/',
                    'kurtosis_cut': 999,
                }
            },
            

    } #end of experiments json

    if experiment not in experiments:
        print(f"{experiment} is not a valid experiment choice. Returning empty dictionary.")
        return {}

    return experiments[experiment]
