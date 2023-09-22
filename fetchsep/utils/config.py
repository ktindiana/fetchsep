import configparser
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

pkg_globals = globals()
for section in config.sections():
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

if __name__ == '__main__':
    """If called from the command line, copy the package config file to the current
    working directory to give the user an example to start with."""
    export_config()
    prepare_dirs()
