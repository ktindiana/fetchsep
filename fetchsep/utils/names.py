from . import config as cfg
from . import experiments as expts
from . import date_handler as dh
import datetime
import os
import sys



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
        modifier = modifier + '_uncor'
        title_mod = title_mod + 'uncorrected '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_B17'
        title_mod = title_mod + 'Bruno2017 '
    if spacecraft:
        modifier = modifier + '_' + spacecraft
        title_mod = title_mod + spacecraft + ' '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if IDSEPEnhancement or OPSEPEnhancement:
        modifier = modifier + '_enhance'

    if module != '':
        modifier = modifier + '_' + module
        title_mod = title_mod + ' (' + module + ')'

    return modifier, title_mod


#####Refer to these subroutines to set units everywhere
def get_energy_units():
    energy_units = cfg.energy_units
    return energy_units

def get_flux_units(flux_type):
    if flux_type == "integral": flux_units = cfg.flux_units_integral
    if flux_type == "differential": flux_units = cfg.flux_units_differential
    return flux_units

def get_flux_units_bin(energy_bin):
    if energy_bin[1] == -1:
        flux_units = cfg.flux_units_integral
    else:
        flux_units = cfg.flux_units_differential
        
    return flux_units

def get_fluence_units(flux_type):
    if flux_type == "integral": fluence_units = cfg.fluence_units_integral
    if flux_type == "differential": fluence_units = cfg.fluence_units_differential
    return fluence_units
    
def get_fluence_units_bin(energy_bin):
    if energy_bin[1] == -1:
        fluence_units = cfg.fluence_units_integral
    else:
        fluence_units = cfg.fluence_units_differential
        
    return fluence_units
#########################


def setup_energy_bin_label(energy_bin):
    """ Label for a single energy bin.
    
    """
    label = ""
    energy_units = get_energy_units()
    if energy_bin[1] != -1:
        label = (f"{energy_bin[0]}-{energy_bin[1]} {energy_units}")
    else:
        label = (f">{energy_bin[0]} {energy_units}")

    return label


def energy_bin_key(bin):
    """ Create key for dataframe or columns header in dataframe and
        csv files.
        
    """
    return f"{bin[0]}-{bin[1]}"


def idsep_naming_scheme(experiment, flux_type, exp_name, modifier=''):
    """ Create naming scheme for subfolders in output/idsep 
    
        Used in:
        SEPTimes files
        HighPoints files
        subdirectories in output
        
    """

    name = (f"{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        name = (f"{exp_name}_{flux_type}{modifier}")
    
    return name
    

def opsep_subdir(experiment, flux_type, exp_name, modifier=''):
    #Return a corresponding directory name
    dir = (f"{experiment}_{flux_type}{modifier}")
    if experiment == 'user' and exp_name != '':
        dir = (f"{exp_name}_{flux_type}{modifier}")
 
    return dir


def opsep_naming_scheme(date, suffix, experiment, flux_type, exp_name, modifier=''):
    """ Create naming scheme for subfolders in output/idsep """

    tzulu = dh.time_to_zulu(date)
    tzulu = tzulu.replace(":","")

    name = (f"{tzulu}_{experiment}_{flux_type}{modifier}_{suffix}")
    if experiment == 'user' and exp_name != '':
        name = (f"{tzulu}_{exp_name}_{flux_type}{modifier}_{suffix}")

    return name
