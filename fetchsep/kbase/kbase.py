import os
import sys
import json


""" Class to interact with the SEP Knowledgebase developed by Ricky Egeland.

"""

class KBase:
    def __init__(self, event_dir, experiment, species):
        """ Set up information defining SEP Knowledgebase directory structure
            and naming scheme.
            
            Naming scheme list:
            SEP Knowledgebase/prototype/200501_01/insitu/proton/ACE_SIS/txt
            
            INPUTS:
            
                :event_dir: (str) name of the directory containing the 
                    observations and measurements for a specific SEP event.
                    Specified by yearmonth_evno: 200310_01, 200310_02
                :experiment: (str) name used by FetchSEP to specify observations
                :species: (str) particle type, e.g. proton, electron

        """
        
        
        self.base_dir = "/Users/kwhitman/Library/CloudStorage/Box-Box"
        self.kbase_dir = "SEP Knowledgebase"
        self.data_dir = "prototype"
        self.event_lists = os.path.join("event_lists","prototype.csv")
        self.event_dir = None
        self.file_type = None
        self.experiment = experiment_name(experiment)
        self.species = species
        self.fetchsep_output_dir = "derived"
        self.manifest = self.read_manifest()

    def set_base_dir(self, dir):
        """ Set the directory on the user computer to the location of 
            SEP Knowledgebase 
            
        """
        self.set_base_dir = dir


    def set_kbase_dir(self, dir):
        """ Set the name of the SEP Knowledgebase top level directory """
        self.set_kbase_dir = dir


    def set_data_dir(self, dir):
        """ Set the location of the data directory """
        self.data_dir = dir


    def set_event_dir(self, dir):
        """ Set the directory to the SEP event in the data directory. 
            Specified by yearmonth_evno: 200310_01, 200310_02
            
        """
        self.event_dir = dir



    def read_manifest(self):
        """ Read manifest.json containing in each event directory
            with info specific to each event period.
            
            ~/SEP Knowledgebase/prototype/200501_01/manifest.json
            
        """
        dir = os.path.join(self.base_dir, self.kbase_dir, self.data_dir, self.event_dir)
        fname = os.path.join(dir,"manifest.json")
        with open(fname,'r') as file:
            manifest = json.load(file)

        return manifest
        
        

    def experiment_name(self, experiment):
        """ Maps experiment names used by FetchSEP to experiment names
            used by KBase.
            
        """
        name = experiment
        map = {"EPHIN": "SOHO_EPHIN",
               "EPHIN_REleASE": "SOHO_EPHIN_HI-RES",
               "ERNE": "SOHO_ERNE",
               "ACE_EPAM_electrons": "ACE_EPAM"
            }

        if experiment in map:
            name = map[experiment]
            
        return name


