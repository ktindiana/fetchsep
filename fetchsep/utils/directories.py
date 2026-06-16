from ..utils import config as cfg
from ..utils import experiments as expts
import os

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

""" Directory structure for data downloaded from the original source by fetchsep. 
    Directories containing output and plots created by idsep, opsep, and download.

    Current directory structures:
        - FetchSEP default
        - SEP Knowledgebase

    Default base directories are set and created in config.py or
    through keywords passed to idsep, opsep, and download:
        cfg.datapath (default = data/)
        cfg.outpath (default = output/)
        cfg.plotpath (default = plots/)
        cfg.listpath (default = lists/)
    
    The directories specified here are the subdirectories below these top levels.

    FetchSEP has the default behavior to create unique subdirectories that include
    the experiment name, flux type, and other modifiers. This behavior may be
    suppressed here.
    
"""


def add_subdir(path, subdir):
    """ Add a subdir under a path and create if doesn't exist. """
    fullpath = os.path.join(path, subdir)
    if not os.path.isdir(fullpath):
        print("Making directory:", fullpath)
        os.mkdir(fullpath)

    return fullpath
    

def add_subdirs(path,subdirs):
    """ Add multiple subdirectories to the path
    
        path/subdir1/subdir2
        
        INPUTS:
        
            :path: (string)
            :subdirs: (arr of strings) [subdir1, subdir2, ...]
    
    """
    fullpath = path
    for subdir in subdirs:
        add_subdir(fullpath, subdir)
        fullpath = os.path.join(fullpath, subdir)

    return fullpath

#################### FETCHSEP PATHS #######################
def create_subdirectories(basedir, module='', subdir='', directory_depth=2):
    """ Default paths within output/ plots/ for specified fetchsep module: 
            idsep
            opsep
            download
            
        INPUTS:
        
            :basedir: (string) base directory, e.g. cfg.outpath (output/)
            :module: (string) idsep, opsep, download
            :subdir: (string) name of subdirectory below module
            :directory_depth: (int) default = 2
                0 - Directory only to cfg.oupath (output/), cfg.plotpath (plots/) level
                1 - Include subdirectory to module, cfg.outpath/module (output/opsep)
                2 - Include subdirectory named according to experiment and
                    options, e.g. cfg.outpath/module/subdir
                    (output/opsep/GOES-13_integral/
            
    """
    
    if directory_depth == 0:
        return basedir
        
    elif directory_depth == 1:
        fullpath = add_subdir(basedir, module)
        return fullpath
        
    elif directory_depth == 2:
        fullpath = add_subdirs(basedir, [module, subdir])
        return fullpath
        
    else:
        sys.exit(f"create_subdirectories: directory_depth valid values are (0, 1, 2). You entered {directory_depth}. Exiting.")
    
####### DEFAULT DIRECTORIES WHERE FETCHSEP DOWNLOADS OBSERVATIONAL DATA #######
def check_path(dir):
    """Check that the paths that hold the data and output exist. If not, create. """
    if not os.path.isdir(dir):
        print(f"check_paths: Creating directory for observational data, {dir}")
        os.makedirs(dir, exist_ok=True)


def get_directory(experiment, use_absolute_datapath=False):
    """ Get directory for location where experiment data is saved. 
        Experiment directory structures are saves in fetchsep.cfg.
        
    """
    if use_absolute_datapath:
        return cfg.datapath
    elif experiment == 'user': #not relevant
        return cfg.datapath
    else:
        #variable has to be inside the subroutine in order to get cfg.datapath
        #set when opsep, idsep, or download were executed.
        """ Return the directory where the observational data is downloaded. """
        if experiment in expts.valid_neutron_monitors():
            path = os.path.join(cfg.datapath, cfg.experiment_directories['NeutronMonitor'], experiment)
        else:
            path = os.path.join(cfg.datapath,cfg.experiment_directories[experiment])

    return path
