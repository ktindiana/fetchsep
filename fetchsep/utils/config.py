import configparser
import argparse
import os.path
import shutil

"""fetchsep configuration

fetchsep is configured using configparser, however configuration
values are accessed in the code as package variables.  This package
reads in a number of configuration files in reverse order of
precedance and evaluates the values found into this package
namespace.  As an example, in a configuration file fetchsep.cfg

    foo = "bar"

Can be accessed like this

    from fetchsep.utils import config
    print(config.foo) # prints "bar"

The order of precedance is:

 1. fetchsep.cfg in the current working directory
 2. .fetchsep in the user's home directory
 3. fetchsep.cfg in the directory of this package

When a given config key is defined in multiple files, the one with the
highest precedence is used.  The higher precedence files do not need
to define all variables; e.g. the package configuration file can be
considered the default and the working directory configuration file
can only contain a few items that the user wishes to be different than
the default.

Note that while section names are required (a configparser
requirement), they are not used.  All variables from all sections are
dumped into the same package namespace.
"""

pkg_configfile = os.path.join(os.path.dirname(__file__), 'fetchsep.cfg')
home_configfile = os.path.expanduser('~/.fetchsep')
cwd_configfile = './fetchsep.cfg'

config = configparser.ConfigParser()
config.read([pkg_configfile, home_configfile, cwd_configfile])

#Sections to "squash" to bring variables to top level
sections = ['fetchsep', 'user_tseries', 'idsep', 'prepare_obs', 'opsep']
pkg_globals = globals()
for section in sections:
    if section in config.sections():
        for k in config[section]:
            pkg_globals[k] = eval(config[section][k])

def export_config(overwrite=False):
    if not overwrite and os.path.isfile(cwd_configfile):
        print("Not overwriting existing config file:", cwd_configfile)
    else:
        print("Copying default configuration file to ", cwd_configfile)
        shutil.copy(pkg_configfile, cwd_configfile)

def prepare_dirs():
    for dir in [datapath, outpath, plotpath, listpath]:
        if not os.path.isdir(dir):
            print('Making directory:', dir)
            os.makedirs(dir)

        if dir != datapath:
            for sub in ['idsep', 'opsep']:
                subdir = os.path.join(dir, sub)
                if not os.path.isdir(subdir):
                    print('Making directory:', subdir)
                    os.makedirs(subdir)

#def configure_for(experiment):
#    section = 'experiment.' + experiment
#    if section in config.sections():
#        for k in config[section]:
#            pkg_globals[k] = eval(config[section][k])

def print_configured_values():
    print()
    print("====== CONFIG VALUES SET TO ========================")
    for section in sections:
        if section in config.sections():
            for k in config[section]:
                print(f"{k}: {pkg_globals[k]}")
    print("======= END CONFIG VALUES =========================")
    print()

##### PATHS ####
def set_datapath(datapath):
    #allow user to set the datapath on the fly across all modules
    pkg_globals['datapath'] = datapath
    print(f"config: Setting global datapath to {datapath}.")
    
def set_outpath(outpath):
    #allow user to set the outpath on the fly across all modules
    pkg_globals['outpath'] = outpath
    print(f"config: Setting global outpath to {outpath}.")
    
def set_plotpath(plotpath):
    #allow user to set the plotpath on the fly across all modules
    pkg_globals['plotpath'] = plotpath
    print(f"config: Setting global plotpath to {plotpath}.")
    
def set_listpath(listpath):
    #allow user to set the listpath on the fly across all modules
    pkg_globals['listpath'] = listpath
    print(f"config: Setting global listpath to {listpath}.")

def set_config_paths(path_to_data=None, path_to_output=None, path_to_plots=None,
    path_to_lists=None):
    #set paths specified by user
    if path_to_data != None and path_to_data != '':
        set_datapath(path_to_data)
    if path_to_output != None and path_to_output != '':
        set_outpath(path_to_output)
    if path_to_plots != None and path_to_plots != '':
        set_plotpath(path_to_plots)
    if path_to_lists != None and path_to_lists != '':
        set_listpath(path_to_lists)
    #Create idsep and opsep subdirectories
    prepare_dirs()


#### IDSEP ####
def set_kurtosis_cut(kurtosis_cut):
    pkg_globals['kurtosis_cut'] = kurtosis_cut
    print(f"config: Setting global kurtosis_cut to {kurtosis_cut}.")

def set_idsep_nsigma(idsep_nsigma):
    pkg_globals['idsep_nsigma'] = idsep_nsigma
    print(f"config: Setting global idsep_nsigma to {idsep_nsigma}.")

def set_idsep_init_win(init_win):
    pkg_globals['init_win'] = init_win
    print(f"config: Setting global idsep init_win to {init_win}.")

def set_idsep_sliding_win(sliding_win):
    pkg_globals['sliding_win'] = sliding_win
    print(f"config: Setting global idsep sliding_win to {sliding_win}.")

def set_idsep_percent_points(percent_points):
    pkg_globals['percent_points'] = percent_points
    print(f"config: Setting global idsep percent_points to {percent_points}.")

def configure_idsep(kurtosis_cut=None, idsep_nsigma=None,
    init_win=None, sliding_win=None, percent_points=None):
    if kurtosis_cut != None:
        set_kurtosis_cut(kurtosis_cut)
        print(f"configure_idsep: Superceding any kurtosis_cut value set in experiments.py")
    if idsep_nsigma != None:
        set_idsep_nsigma(idsep_nsigma)
    if init_win != None:
        set_idsep_init_win(init_win)
    if sliding_win != None:
        set_idsep_sliding_win(sliding_win)
    if percent_points != None:
        set_idsep_percent_points(percent_points)



#### OPSEP ####
def set_opsep_nsigma(opsep_nsigma):
    pkg_globals['opsep_nsigma'] = opsep_nsigma
    print(f"config: Setting global opsep_nsigma to {opsep_nsigma}.")

def configure_opsep(opsep_nsigma=None):
    if opsep_nsigma != None:
        set_opsep_nsigma(opsep_nsigma)


#### UNITS ####
def set_energy_units(energy_units):
    #allow user to set the datapath on the fly across all modules
    pkg_globals['energy_units'] = energy_units
    print(f"config: Setting global energy_units to {energy_units}.")

def set_flux_units(flux_units_integral, fluence_units_integral,
    flux_units_differential, fluence_units_differential):
    #allow user to set the datapath on the fly across all modules
    pkg_globals['flux_units_integral'] = flux_units_integral
    pkg_globals['fluence_units_integral'] = fluence_units_integral
    pkg_globals['flux_units_differential'] = flux_units_differential
    pkg_globals['fluence_units_differential'] = fluence_units_differential
    print(f"config: Setting global flux and fluence units.")


if __name__ == '__main__':
    """If called from the command line, copy the package config file to the current
    working directory to give the user an example to start with."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--datapath", type=str, default=None,
            help=("Specify path where satellite data should be downloaded. Will default to "
                "datapath listed in fetchsep.cfg if a value is not specified."))

    parser.add_argument("--outpath", type=str, default=None,
            help=("Specify path where output files should be saved. Will default to "
                "outpath listed in fetchsep.cfg if a value is not specified."))
                
    parser.add_argument("--plotpath", type=str, default=None,
            help=("Specify path where plots should be saved. Will default to "
                "plotpath listed in fetchsep.cfg if a value is not specified."))

    parser.add_argument("--listpath", type=str, default=None,
            help=("Specify path where lists should be saved. Will default to "
                "listpath listed in fetchsep.cfg if a value is not specified."))

    args = parser.parse_args()
    path_to_data = args.datapath
    path_to_output = args.outpath
    path_to_plots = args.plotpath
    path_to_lists = args.listpath

    export_config()
    set_config_paths(path_to_data=path_to_data, path_to_output=path_to_output,
        path_to_plots=path_to_plots, path_to_lists=path_to_lists)
