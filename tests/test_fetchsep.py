#!/usr/bin/env python

"""Tests for `fetchsep` package."""


import unittest
import os
import json
import shutil
import contextlib
import math
import datetime
import pandas as pd
import sys
import glob

import fetchsep
import fetchsep.opsep.opsep as opsep
import fetchsep.utils.read_datasets as datasets
import fetchsep.utils.date_handler as dh
import fetchsep.json.ccmc_json_handler as ccmc_json
import fetchsep.utils.experiments as expts
import fetchsep.utils.config as cfg
import fetchsep.utils.download as fsdl
import fetchsep.utils.parameters as fsparam


def utility_get_verbosity():
    """
    Sets verbosity for unit test suite.
    """
    if ('-v' in sys.argv) or ('--verbose' in sys.argv):
        verbosity = 2
    elif ('--quiet' in sys.argv):
        verbosity = 0
    else:
        verbosity = 1
    return verbosity


datapath = os.path.join('tests','data')
outpath = os.path.join('tests','output')
plotpath = os.path.join('tests','plots')
listpath = os.path.join('tests','lists')

#Clean paths
for pth in [datapath, outpath, plotpath, listpath]:
    if os.path.isdir(pth):
        shutil.rmtree(pth)

with contextlib.redirect_stdout(None):
    cfg.set_config_paths(path_to_data=datapath, path_to_output=outpath,
        path_to_plots=plotpath, path_to_lists=listpath)

class TestFetchsep(unittest.TestCase):
    """Tests for `fetchsep` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_version(self):
        """Test that version exists."""
        self.assertIsInstance(fetchsep.__version__, str)



class TestDownload(unittest.TestCase):
    """ Test download functions for each native experiment """

    def setUp(self):
        self.verbosity = utility_get_verbosity()


    def test_download_data(self):
        """ Test download and read-in functions in read_datasets """
 
        #SOHO/ERNE .dates files
        erne_path = os.path.join(cfg.datapath,'SOHO','ERNE')
        os.mkdir(os.path.join(cfg.datapath,'SOHO'))
        os.mkdir(erne_path)
        os.mkdir(os.path.join(erne_path,'export.src'))
        for file in glob.glob('tests/files/data/SOHO/ERNE/*.dates'):
            shutil.copy(file, erne_path)
            
        for file in glob.glob('tests/files/data/SOHO/ERNE/export.src/*'):
            shutil.copy(file, os.path.join(erne_path, 'export.src'))

        experiments = expts.valid_experiments()

        #Experiments to skip.
        #Don't test every neutron monitor. OULU is sufficient.
        all_nm = expts.valid_neutron_monitors()
        all_nm.remove('OULU')
        skip = ['user', 'GOES', 'EPHIN_HESPERIA', 'CalGOES', 'SEPEM', 'SEPEMv3'] + all_nm
        
        for experiment in experiments:
            if experiment in skip:
                continue
            
            exp_info = expts.experiment_info(experiment)

            #Dates to test data download
            start_date = exp_info['first_date']
            if experiment == "GOES-RT" or experiment == "GOES-16":
                start_date = datetime.datetime(2021,9,3) #Data gaps after end of GOES-15 on 2020-03-04
            elif experiment == "GOES-SWPC":
                start_date = datetime.datetime.now() - datetime.timedelta(hours = 73)
            end_date = start_date + datetime.timedelta(hours = 72)

            #Available flux_types for experiment
            flux_types = exp_info['flux_type']
            
            if 'spacecraft' in exp_info.keys():
                spacecrafts = exp_info['spacecraft']
            else:
                spacecrafts = [None]
            
            
            for flux_type in flux_types:
                for spacecraft in spacecrafts:
                    if self.verbosity == 2:
                        print(f"Downloading {experiment} {flux_type} spacecraft={spacecraft} for {start_date} to {end_date}")

                    with contextlib.redirect_stdout(None):
                        params = fsdl.load_parameters(start_date, end_date, experiment,
                            flux_type=flux_type, spacecraft=spacecraft)
                        
                    #Check that parameters were set appropriately
                    self.assertEqual(params.json_type, exp_info['json_type'])
                    self.assertEqual(params.json_mode, exp_info['json_mode'])
                    self.assertEqual(params.species, exp_info['species'])
                    self.assertEqual(params.location, exp_info['location'])
                    self.assertEqual(params.kurtosis_cut, exp_info[flux_type]['kurtosis_cut'])

                    #Try to download data
                    #with contextlib.redirect_stdout(None):
                    module_outpath, module_plotpath, dates, fluxes, energy_bins, energy_bin_centers = fsdl.get_data(params, saveplot=True)
            
                    
                    self.assertNotEqual(len(dates), 0)
                    if experiment != 'SOHO_ERNE':
                        self.assertEqual(sorted(exp_info[flux_type]['energy_bins']), energy_bins)
                        self.assertEqual(sorted(exp_info[flux_type]['energy_bin_centers']), energy_bin_centers)


    def test_local_data(self):
        """ Test read-in functions in read_datasets for experimental data that
            must already be present on the user's computer - CalGOES, SEPEM, SEPEMv3.
            This unit test contains a small sample of each dataset, but the
            user must download the full datasets to use them with fetchsep.
            
        """
        #Copy subset of data files that user has to download manually for use in FetchSEP
        #SEPEMv2
        sepemv2_path = os.path.join(cfg.datapath,'SEPEMv2')
        os.mkdir(sepemv2_path)
        shutil.copy('tests/files/data/SEPEMv2/SEPEM_H_reference_1974.csv', sepemv2_path)

        #SEPEMv3
        sepemv3_path = os.path.join(cfg.datapath,'SEPEMv3')
        os.mkdir(sepemv3_path)
        shutil.copy('tests/files/data/SEPEMv3/SEPEM_RDS_v3_H_1974.csv', sepemv3_path)

        #CalGOES
        cg_path = os.path.join(cfg.datapath,'CalGOES')
        os.mkdir(cg_path)
        shutil.copy('tests/files/data/CalGOES/srag12_1986.dat', cg_path)

        experiments = ['CalGOES', 'SEPEM', 'SEPEMv3']

        for experiment in experiments:
            exp_info = expts.experiment_info(experiment)

            #Dates to test data download
            start_date = exp_info['first_date']
            end_date = start_date + datetime.timedelta(hours = 72)

            #Available flux_types for experiment
            flux_types = exp_info['flux_type']
            spacecrafts = [None]
            
            for flux_type in flux_types:
                for spacecraft in spacecrafts:
                    if self.verbosity == 2:
                        print(f"Reading in {experiment} {flux_type} spacecraft={spacecraft} for {start_date} to {end_date}.")

                    with contextlib.redirect_stdout(None):
                        params = fsdl.load_parameters(start_date, end_date, experiment,
                            flux_type=flux_type, spacecraft=spacecraft)
                        
                    #Check that parameters were set appropriately
                    self.assertEqual(params.json_type, exp_info['json_type'])
                    self.assertEqual(params.json_mode, exp_info['json_mode'])
                    self.assertEqual(params.species, exp_info['species'])
                    self.assertEqual(params.location, exp_info['location'])
                    self.assertEqual(params.kurtosis_cut, exp_info[flux_type]['kurtosis_cut'])

                    #Try to download data
                    #with contextlib.redirect_stdout(None):
                    module_outpath, module_plotpath, dates, fluxes, energy_bins, energy_bin_centers = fsdl.get_data(params, saveplot=True)
            
                    self.assertNotEqual(len(dates), 0)
                    self.assertEqual(sorted(exp_info[flux_type]['energy_bins']), energy_bins)
                    self.assertEqual(sorted(exp_info[flux_type]['energy_bin_centers']), energy_bin_centers)




########### USER INPUT FILE ##############
def load_flux_timeseries(filename):
    """ Read in .txt output by opsep in format YYYY-MM-DDTHH:MM:SSZ FLUX """
    dates = []
    fluxes = []
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip().split()
            date = dh.str_to_datetime(line[0])
            flux = float(line[1])
            dates.append(date)
            fluxes.append(flux)
            
        file.close()
    return dates, fluxes


class TestExperimentOpsep(unittest.TestCase):
    """ Test native experiments in OpSEP and overall functionality. """
    @classmethod
    def setUpClass(cls):
        startdate = "2012-05-17 00:10:00"
        enddate = "2012-05-22"
        experiment = "GOES-13"
        flux_type = "integral"
        showplot = False
        saveplot = True
        user_thresholds = "30,1;50,1"
        associations = True

        #####INTEGRAL FLUXES
        with contextlib.redirect_stdout(None):
            cls.opsep_outputs_integral = opsep.run_opsep(startdate, enddate,
                experiment, flux_type=flux_type, showplot=showplot, saveplot=saveplot,
                user_thresholds=user_thresholds)
        
        if utility_get_verbosity() == 2:
            print(f"\n[setUpClass] Creating opsep outputs for native experiment ({experiment} {flux_type}) tests.")

        flux_type = "differential"
        user_thresholds = "30,1;50,1;38.0-82.0,0.1"

        #####DIFFERENTIAL FLUXES
        with contextlib.redirect_stdout(None):
            cls.opsep_outputs_differential = opsep.run_opsep(startdate, enddate,
                experiment, flux_type=flux_type, showplot=showplot, saveplot=saveplot,
                user_thresholds=user_thresholds)
        
        if utility_get_verbosity() == 2:
            print(f"\n[setUpClass] Creating opsep outputs for native experiment ({experiment} {flux_type}) tests.")



    @classmethod
    def tearDownClass(cls):
        if utility_get_verbosity() == 2:
            print("\n[tearDownClass] Cleaning up and deleting opsep native experiment outputs.")
        del cls.opsep_outputs_integral
        del cls.opsep_outputs_differential


#    def setUp(self):
#        self.verbosity = utility_get_verbosity()
#        
#        ref_file = 'tests/files/output/opsep/GOES-13_integral/'
#        with open(ref_file,"r") as f:
#            self.ref_json = json.load(f)
#            f.close()
#        
#        test_file = self.opsep_outputs["jsonfname"]
#        with open(test_file,"r") as f:
#            self.test_json = json.load(f)
#            f.close()
#
#        self.ref_pathnm = 'tests/files/user/'
#
# 
 


class TestUserOpsep(unittest.TestCase):
    """ Test user-input flux timeseries in OpSEP """
    @classmethod
    def setUpClass(cls):
        filename = 'tests/files/user/ZEUS+iPATH_CME_20260605_090900_20260605_132007_mars_differential-flux.csv'
        startdate = "2026-06-05 11:20:00"
        enddate = "2026-06-08 10:20:00"
        experiment = "user"
        user_name = "iPATH_test"
        flux_type = "differential"
        showplot = False
        saveplot = True
        location = "mars"
        json_type = "model"
        json_mode = "unit test"
        dointerp = False
        user_thresholds = "30,1;50,1"

        with contextlib.redirect_stdout(None):
            #Set iPATH user file values
            cfg.set_user_delimeter(',')
            cfg.set_user_columns([1,2,3,4,5,6,7,8,9,10,11,
                                  12,13,14,15,16,17,18,19,
                                  20,21,22,23,24,25])
            cfg.set_user_energy_bins(
                [[0.1,0.1],[0.14677993,0.14677993],[0.21544347,0.21544347],
                [0.31622777,0.31622777],[0.46415888,0.46415888],
                [0.68129207,0.68129207],[1,1],[1.4677993,1.4677993],
                [2.1544347,2.1544347],[3.1622777,3.1622777],[4.6415888,4.6415888],
                [6.8129207,6.8129207],[10,10],[14.677993,14.677993],
                [21.544347,21.544347],[31.622777,31.622777],[46.415888,46.415888],
                [68.129207,68.129207],[100,100],[146.77993,146.77993],
                [215.44347,215.44347],[316.22777,316.22777],[464.15888,464.15888],
                [681.29207,681.29207],[1000,1000]])
        
        cme_start_time = "2026-06-05T09:09:00Z"
        cme_half_width = 46.0
        cme_speed = 999.0
        cme_lat = -23.0
        cme_lon = -161.0
        cme_height = 21.5
        cme_time_at_height_time = "2026-06-05T12:42Z"
        cme_time_at_height_height = 21.5
        cme_coordinates = "HEEQ"
        cme_catalog = "DONKI"
        cme_catalog_id = "2026-06-05T09:09:00-CME-001"
        cme_urls = ["https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CMEAnalysis/46649/-1"]
        cme_derivation_process = "manual"
        cme_derivation_method = "SWPC_CAT"
        cme_measurement_type = "LE"

        with contextlib.redirect_stdout(None):
            cls.opsep_outputs = opsep.run_opsep(startdate, enddate,
                experiment, flux_type=flux_type, user_name=user_name,
                user_file=filename, json_type=json_type, json_mode=json_mode,
                dointerp=dointerp, showplot=showplot, saveplot=saveplot,
                user_thresholds=user_thresholds, location=location,
                cme_start_time=cme_start_time,
                cme_half_width=cme_half_width,
                cme_speed=cme_speed,
                cme_lat=cme_lat, cme_lon=cme_lon,
                cme_height=cme_height,
                cme_time_at_height_time=cme_time_at_height_time,
                cme_time_at_height_height=cme_time_at_height_height,
                cme_coordinates=cme_coordinates,
                cme_catalog=cme_catalog,
                cme_catalog_id=cme_catalog_id,
                cme_urls=cme_urls,
                cme_derivation_process=cme_derivation_process,
                cme_derivation_method=cme_derivation_method,
                cme_measurement_type=cme_measurement_type)
        
        if utility_get_verbosity() == 2:
            print("\n[setUpClass] Creating opsep outputs for user input file tests.")


    @classmethod
    def tearDownClass(cls):
        if utility_get_verbosity() == 2:
            print("\n[tearDownClass] Cleaning up and deleting opsep outputs.")
        del cls.opsep_outputs

    def setUp(self):
        self.verbosity = utility_get_verbosity()
        
        ref_file = 'tests/files/user/ZEUS+iPATH_CME.Mars.2026-06-05T090900Z.2026-06-05T132007Z.json'
        with open(ref_file,"r") as f:
            self.ref_json = json.load(f)
            f.close()
        
        test_file = self.opsep_outputs["jsonfname"]
        with open(test_file,"r") as f:
            self.test_json = json.load(f)
            f.close()

        self.ref_pathnm = 'tests/files/user/'


    def test_user_event_start_end(self):
        for i in range(len(self.ref_json["sep_forecast_submission"]["forecasts"][0]["event_lengths"])):
            ref = self.ref_json["sep_forecast_submission"]["forecasts"][0]["event_lengths"][i]["start_time"]
            test = self.test_json["sep_forecast_submission"]["forecasts"][0]["event_lengths"][i]["start_time"]
            self.assertEqual(ref, test)
            if self.verbosity == 2:
                print(f"------------ test_user_event_start_end START TIME: ref {ref}, test {test}")


            ref = self.ref_json["sep_forecast_submission"]["forecasts"][0]["event_lengths"][i]["end_time"]
            test = self.test_json["sep_forecast_submission"]["forecasts"][0]["event_lengths"][i]["end_time"]
            self.assertEqual(ref, test)
            if self.verbosity == 2:
                print(f"------------ test_user_event_start_end END TIME: ref {ref}, test {test}")


    def test_user_fluence(self):
        for i in range(len(self.ref_json["sep_forecast_submission"]["forecasts"][0]["fluences"])):
            ref = self.ref_json["sep_forecast_submission"]["forecasts"][0]["fluences"][i]
            test = self.test_json["sep_forecast_submission"]["forecasts"][0]["fluences"][i]
            self.assertEqual(ref, test)
            if self.verbosity == 2:
                print(f"------------ test_user_fluence: ref {ref}, test {test}")


    def test_user_fluence_spectrum(self):
        ref = self.ref_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["start_time"]
        test = self.test_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["start_time"]
        self.assertEqual(ref, test)
        if self.verbosity == 2:
            print(f"------------ test_user_fluence_spectrum: ref {ref}, test {test}")

        ref = self.ref_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["end_time"]
        test = self.test_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["end_time"]
        self.assertEqual(ref, test)
        if self.verbosity == 2:
            print(f"------------ test_user_fluence_spectrum: ref {ref}, test {test}")

        ref = self.ref_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["fluence_units"]
        test = self.test_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["fluence_units"]
        self.assertEqual(ref, test)
        if self.verbosity == 2:
            print(f"------------ test_user_fluence_spectrum: ref {ref}, test {test}")

        for i in range(len(self.ref_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["fluence_spectrum"])):
            ref = self.ref_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["fluence_spectrum"][i]["fluence"]
            test = self.test_json["sep_forecast_submission"]["forecasts"][0]["fluence_spectra"][0]["fluence_spectrum"][i]["fluence"]
            self.assertAlmostEqual(ref, test)
            if self.verbosity == 2:
                print(f"------------ test_user_fluence_spectrum: ref {ref}, test {test}")
            

    def test_user_max_flux(self):
        for i in range(len(self.ref_json["sep_forecast_submission"]["forecasts"])):
            ref = self.ref_json["sep_forecast_submission"]["forecasts"][i]["peak_intensity_max"]["intensity"]
            test = self.test_json["sep_forecast_submission"]["forecasts"][i]["peak_intensity_max"]["intensity"]
            self.assertAlmostEqual(ref, test)
            if self.verbosity == 2:
                print(f"------------ test_user_max_flux: ref {ref}, test {test}")
        
            ref = self.ref_json["sep_forecast_submission"]["forecasts"][i]["peak_intensity_max"]["units"]
            test = self.test_json["sep_forecast_submission"]["forecasts"][i]["peak_intensity_max"]["units"]
            self.assertAlmostEqual(ref, test)
            if self.verbosity == 2:
                print(f"------------ test_user_max_flux: ref {ref}, test {test}")

            ref = self.ref_json["sep_forecast_submission"]["forecasts"][i]["peak_intensity_max"]["time"]
            test = self.test_json["sep_forecast_submission"]["forecasts"][i]["peak_intensity_max"]["time"]
            self.assertAlmostEqual(ref, test)
            if self.verbosity == 2:
                print(f"------------ test_user_max_flux: ref {ref}, test {test}")


    def test_user_time_profile(self):
        test_pathnm = os.path.dirname(self.opsep_outputs["jsonfname"])
        
        for i in range(len(self.ref_json["sep_forecast_submission"]["forecasts"])):
            ref_fname = os.path.join(self.ref_pathnm, self.ref_json["sep_forecast_submission"]["forecasts"][i]["sep_profile"])
            ref_dates, ref_fluxes = load_flux_timeseries(ref_fname)
 
            test_fname = os.path.join(test_pathnm, self.test_json["sep_forecast_submission"]["forecasts"][i]["sep_profile"])
            test_dates, test_fluxes = load_flux_timeseries(test_fname)

            for i in range(len(ref_dates)):
                self.assertEqual(ref_dates[i], test_dates[i])
                #print(f"------------ test_user_time_profile: ref {ref_dates[i]}, test {test_dates[i]}")

                self.assertAlmostEqual(ref_fluxes[i], test_fluxes[i], delta=0.0001)
                #print(f"------------ test_user_time_profile: ref {ref_fluxes[i]}, test {test_fluxes[i]}")


    def test_user_trigger(self):
        keys = self.ref_json["sep_forecast_submission"]["triggers"][0]["cme"].keys()
        for key in keys:
            if isinstance(self.ref_json["sep_forecast_submission"]["triggers"][0]["cme"][key], dict):
                keys2 = self.ref_json["sep_forecast_submission"]["triggers"][0]["cme"][key].keys()
                for key2 in keys2:
                    ref = self.ref_json["sep_forecast_submission"]["triggers"][0]["cme"][key][key2]
                    test = self.test_json["sep_forecast_submission"]["triggers"][0]["cme"][key][key2]
                    if 'time' in key2:
                        ref = dh.str_to_datetime(ref)
                        test = dh.str_to_datetime(test)
                    self.assertEqual(ref, test)
                    if self.verbosity == 2:
                        print(f"------------ test_user_trigger: ref {ref}, test {test}")
            else:
                ref = self.ref_json["sep_forecast_submission"]["triggers"][0]["cme"][key]
                test = self.test_json["sep_forecast_submission"]["triggers"][0]["cme"][key]
                if 'time' in key:
                    ref = dh.str_to_datetime(ref)
                    test = dh.str_to_datetime(test)
                self.assertEqual(ref, test)
                if self.verbosity == 2:
                    print(f"------------ test_user_trigger: ref {ref}, test {test}")


class TestParameters(unittest.TestCase):
    """ Test setting and modifying fetchsep parameters. """

    def setUp(self):
        self.verbosity = utility_get_verbosity()

        ####EXPECTED DEFAULTS####
        self.expected_defaults = {
                    'experiment': 'user',
                    'user': False,
                    'user_name': '',
                    'user_filename': '',
                    'is_unixtime': False,
                    'directory_depth': 2,
                    'use_absolute_datapath': False,
                    'options': [],
                    'goes_datatype': 'corrected',
                    'goes_S14': False,
                    'goes_Bruno2017': False,
                    'showplot': False,
                    'saveplot': False,
                    'modifier': '',
                    'title_modifier': '',
                    'idsep_subdir': '',
                    'idsep_outpath': '',
                    'idsep_plotpath': '',
                    'idsep_path': '',
                    'module_subdir': '',
                    'module_outpath': '',
                    'module_plotpath': '',
                    'remove_above': 999999,
                    'for_inclusive': False,
                    'idsep_nsigma': cfg.idsep_nsigma,
                    'init_win': cfg.init_win,
                    'sliding_win': cfg.sliding_win,
                    'percent_points': cfg.percent_points,
                    'write_fluxes': True,
                    'kurtosis_cut': 999,
                    'location': 'earth',
                    'species': 'proton',
                    'json_type': '',
                    'json_mode': '',
                    'spase_id': '',
                    'do_interpolation': False,
                    'user_thresholds': '',
                    'opsep_nsigma': cfg.opsep_nsigma,
                    'doBGSubOPSEP': False,
                    'bgstartdate': pd.NaT,
                    'bgenddate': pd.NaT,
                    'OPSEPEnhancement': False,
                    'doBGSubIDSEP': False,
                    'IDSEPEnhancement': False,
                    'two_peaks': False,
                    'detect_prev_event': False
                }

        #### CHANGE PARAMETER VALUES ###
        self.expected_references = {
            'experiment': 'GOES-13',
            'flux_type': 'differential',
            'user': False,
            'user_name': 'Parameter_Unit_Test',
            'user_filename': 'test_file.txt',
            'is_unixtime': True,
            'directory_depth': 0,
            'use_absolute_datapath': True,
            'options': ['S14', 'Bruno2017', 'uncorrected'],
            'goes_datatype': 'uncorrected',
            'goes_S14': True,
            'goes_Bruno2017': True,
            'showplot': True,
            'saveplot': True,
            'modifier': '_uncor_S14_B17_bgsub_enhance_opsep',
            'title_modifier': 'uncorrected S14 Bruno2017 BG-subtracted  (opsep)',
            'idsep_subdir': 'GOES-13_differential_uncor_S14_B17',
            'module_subdir': 'GOES-13_differential_uncor_S14_B17_bgsub_enhance_opsep',
            'remove_above': 10.0,
            'for_inclusive': True,
            'idsep_nsigma': 4,
            'init_win': 25,
            'sliding_win': 10,
            'percent_points': 0.1,
            'write_fluxes': False,
            'kurtosis_cut': 60,
            'location': 'earth',
            'species': 'proton',
            'json_type': 'observations',
            'json_mode': 'measurement',
            'spase_id': 'test_spase_id',
            'do_interpolation': True,
            'user_thresholds': '30,1;50,1',
            'opsep_nsigma': 2.1,
            'doBGSubOPSEP': True,
            'bgstartdate': datetime.datetime(2012, 1, 1),
            'bgenddate': datetime.datetime(2012, 1, 5),
            'OPSEPEnhancement': True,
            'doBGSubIDSEP': False,
            'IDSEPEnhancement': False,
            'two_peaks': True,
            'detect_prev_event': True
        }


        self.expected_references_user = {
            'experiment': 'user',
            'flux_type': 'differential',
            'user': True,
            'user_name': 'Parameter_Unit_Test',
            'user_filename': 'test_file.txt',
            'directory_depth': 1,
            'modifier': '_bgsub_enhance_idsep',
            'title_modifier': 'BG-subtracted  (idsep)',
            'idsep_subdir': 'Parameter_Unit_Test_differential',
            'module_subdir': 'Parameter_Unit_Test_differential_bgsub_enhance_idsep',
            'location': 'mars',
            'species': 'electron',
            'json_type': 'observations',
            'json_mode': 'measurement',
            'doBGSubIDSEP': True,
            'IDSEPEnhancement': True,
        }

        #self.ref_idsep_outpath = ''
        #self.ref_idsep_plotpath = ''
        #self.ref_idsep_path = ''
        #self.ref_module_outpath = ''
        #self.ref_module_plotpath = ''

    def test_default_parameters(self):
        """ Test default Parameters class """
        with contextlib.redirect_stdout(None):
            test_param = fsparam.Parameters('unittest', '2012-01-05', '2012-01-10', 'user')
        for param, ref in self.expected_defaults.items():
            test = getattr(test_param, param)
            if self.verbosity == 2:
                print(f"DEFAULT Parameters {param} ref: {ref} test: {test}")
            if isinstance(ref, list):
                self.assertEqual(ref, test)
            elif pd.isna(ref):
                self.assertTrue(pd.isna(test))
            else:
                self.assertEqual(ref, test)
 

    def test_set_parameters(self):
        """ Test setting parameters in Parameter class with native experiment """
        with contextlib.redirect_stdout(None):
            test_param = fsparam.Parameters('opsep', '2012-01-05', '2012-01-10', 'GOES-13')
            test_param.set_values(
                flux_type=self.expected_references['flux_type'],
                spacecraft=None,  # Not present in self.expected_references
                user_name=self.expected_references['user_name'],
                user_file=self.expected_references['user_filename'],
                is_unixtime=self.expected_references['is_unixtime'],
                options='S14;Bruno2017;uncorrected',
                dointerp=self.expected_references['do_interpolation'],
                showplot=self.expected_references['showplot'],
                saveplot=self.expected_references['saveplot'],
                directory_depth=self.expected_references['directory_depth'],
                use_absolute_datapath=self.expected_references['use_absolute_datapath'],
                write_fluxes=self.expected_references['write_fluxes'],
                for_inclusive=self.expected_references['for_inclusive'],
                remove_above=self.expected_references['remove_above'],
                kurtosis_cut=self.expected_references['kurtosis_cut'],
                idsep_nsigma=self.expected_references['idsep_nsigma'],
                init_win=self.expected_references['init_win'],
                sliding_win=self.expected_references['sliding_win'],
                percent_points=self.expected_references['percent_points'],
                opsep_nsigma=self.expected_references['opsep_nsigma'],
                color_scheme=None,  # Not present in self.expected_references
                no_goes_colors=None,  # Not present in self.expected_references
                json_type=self.expected_references['json_type'],
                json_mode=self.expected_references['json_mode'],
                spase_id=self.expected_references['spase_id'],
                detect_prev_event=self.expected_references['detect_prev_event'],
                two_peaks=self.expected_references['two_peaks'],
                user_thresholds=self.expected_references['user_thresholds'],
                doBGSubOPSEP=self.expected_references['doBGSubOPSEP'],
                OPSEPEnhancement=self.expected_references['OPSEPEnhancement'],
                bgstartdate=self.expected_references['bgstartdate'],
                bgenddate=self.expected_references['bgenddate'],
                doBGSubIDSEP=self.expected_references['doBGSubIDSEP'],
                IDSEPEnhancement=self.expected_references['IDSEPEnhancement'],
                idsep_path=None,  # Not present in self.expected_references
                location=self.expected_references['location'],
                species=self.expected_references['species']
            )

        for param, ref in self.expected_references.items():
            test = getattr(test_param, param)
            if self.verbosity == 2:
                print(f"SET Parameters {param} ref: {ref} test: {test}")
            if isinstance(ref, list):
                for i in range(len(ref)):
                    self.assertEqual(ref[i], test[i])
            elif pd.isna(ref):
                self.assertTrue(pd.isna(test))
            else:
                self.assertEqual(ref, test)



    def test_set_parameters_user(self):
        """ Test setting additional parameters in Parameter class with 'user' experiment """
        with contextlib.redirect_stdout(None):
            test_param = fsparam.Parameters('opsep', '2012-01-05', '2012-01-10', 'user')
            test_param.set_values(
                flux_type=self.expected_references_user['flux_type'],
                user_name=self.expected_references_user['user_name'],
                user_file=self.expected_references_user['user_filename'],
                directory_depth=self.expected_references_user['directory_depth'],
                json_type=self.expected_references_user['json_type'],
                json_mode=self.expected_references_user['json_mode'],
                doBGSubIDSEP=self.expected_references_user['doBGSubIDSEP'],
                IDSEPEnhancement=self.expected_references_user['IDSEPEnhancement'],
                location=self.expected_references_user['location'],
                species=self.expected_references_user['species']
            )

        for param, ref in self.expected_references_user.items():
            test = getattr(test_param, param)
            if self.verbosity == 2:
                print(f"SET Parameters USER {param} ref: {ref} test: {test}")
            if isinstance(ref, list):
                for i in range(len(ref)):
                    self.assertEqual(ref[i], test[i])
            elif pd.isna(ref):
                self.assertTrue(pd.isna(test))
            else:
                self.assertEqual(ref, test)




