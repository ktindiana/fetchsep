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
def check_paths(experiment):
    """Check that the paths that hold the data and output exist. If not, create.
    """
    if 'STEREO' in experiment:
        for det in ['LET', 'HET']:
            dir = get_directory(experiment + '_' + det)
            if not os.path.isdir(dir):
                print(f"check_paths: Creating directory for {experiment} observational data, {dir}")
                os.makedirs(dir, exist_ok=True)
    else:
        dir = get_directory(experiment)
        if not os.path.isdir(dir):
            print(f"check_paths: Creating directory for {experiment} observational data, {dir}")
            os.makedirs(dir, exist_ok=True)


def get_directory(experiment):

    #variable has to be inside the subroutine in order to get cfg.datapath
    #set when opsep, idsep, or download were executed.
    experiment_directories = {
        'GOES-RT': os.path.join(cfg.datapath, 'GOES-RT'),
        'GOES-SWPC': os.path.join(cfg.datapath, 'GOES-SWPC'),
        'GOES': os.path.join(cfg.datapath, 'GOES'),
        'STEREO-A_HET': os.path.join(cfg.datapath, 'STEREO-A','HET'),
        'STEREO-A_LET': os.path.join(cfg.datapath, 'STEREO-A','LET'),
        'STEREO-B_HET': os.path.join(cfg.datapath, 'STEREO-B','HET'),
        'STEREO-B_LET': os.path.join(cfg.datapath, 'STEREO-B','LET'),
        'ACE_SIS': os.path.join(cfg.datapath, 'ACE','SIS'),
        'ACE_EPAM_electrons': os.path.join(cfg.datapath, 'ACE','EPAM'),
        'ERNE': os.path.join(cfg.datapath, 'ERNE'),
        'EPHIN': os.path.join(cfg.datapath, 'EPHIN'),
        'NeutronMonitor': os.path.join(cfg.datapath, 'NeutronMonitor'),
        'SEPEM': os.path.join(cfg.datapath, 'SEPEM'),
        'SEPEMv3': os.path.join(cfg.datapath, 'SEPEMv3'),
        'CalGOES': os.path.join(cfg.datapath, 'CalGOES'),
        'CalGOES': os.path.join(cfg.datapath, 'CalGOES'),
        'EPHIN_HESPERIA': os.path.join(cfg.datapath, 'EPHIN_HESPERIA'),
        'EPHIN_REleASE': os.path.join(cfg.datapath, 'EPHIN_REleASE'),
        'IMP8_CPME': os.path.join(cfg.datapath, 'IMP8','CPME'),
    }
    """ Return the directory where the observational data is downloaded. """
    if experiment in expts.valid_neutron_monitors():
        path = os.path.join(experiment_directories['NeutronMonitor'], experiment)
    else:
        path = experiment_directories[experiment]
    print(f"THIS IS THE PATH {path}")
    return path
