from ..opsep import opsep
from . import config as cfg
from ..json import keys
from ..json import ccmc_json_handler as ccmc_json
import datetime
import sys
import os
import pandas as pd
import numpy as np
import shutil

__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

""" Moves output produced by FetchSEP to a target directory.
    Lists containing SEP event dates and 'observation windows'
    are produced.
    
    If these lists containing the SEP event dates and observation
    windows already exist in the target directory, then there
    is the option to add only new files to the directory and
    append the lists. This allows the user to run fetchsep and
    incrementally add observation files to the target directory
    without overwriting.
    
    Structure of target directory:
        target_dir/observation_windows.txt
        target_dir/approved_SEPs.txt
        target_dir/data/*.json and *.txt - files containing observations
        target_dir/data/plots/*.png - plots of prepared time periods

"""

def read_batch_list():
    """ Read in outpath/idsep/batch_event_list.txt which contains
        the observation windows that were used to process time
        periods with OpSEP.
        
    """
    obs_st = []
    obs_end = []
    
    fname = os.path.join(cfg.outpath, "idsep","batch_event_list.txt")
    if not os.path.isfile(fname):
        sys.exit("read_batch_list: Cannot find " + fname + ". Exiting.")
    
    with open(fname,"r") as file:
        for line in file:
            line.strip()
            if '#' in line: continue
            if line == '': continue
    
            line = line.split(',')
            stdate = line[0].strip()
            enddate = line[1].strip()
            
            if len(stdate) == 10:
                stdate = stdate + " 00:00:00"
            if len(enddate) == 10:
                enddate = enddate + " 00:00:00"
                
            obs_st.append(datetime.datetime.strptime(stdate, "%Y-%m-%d %H:%M:%S"))
            obs_end.append(datetime.datetime.strptime(enddate, "%Y-%m-%d %H:%M:%S"))

    return obs_st, obs_end


def make_observation_window_list(path):
    """ Look in the output/opsep directory for observation files.
        Create a list of filenames and compile the associated 
        observation window start and end times.
        
    """
    print("make_observation_window_list: Identifying observation windows for files in " + path)
    
    #List all json files in the opsep output directory
    allfiles = os.listdir(path)
    jsonfiles = [os.path.join(path,f) for f in allfiles if '.json' in f]

    win_st = []
    win_end = []
    ek_id = keys.id_energy_channel
    id_pred_win = keys.id_prediction_window
    for fname in jsonfiles:
        data = ccmc_json.read_in_json(fname)
        if data == None:
            continue
        nblocks = ccmc_json.return_nforecasts(data)
        #Grab the observation window from the first observations block
        #The observation window is the same for all blocks (unlike forecasts,
        #where the prediction window may differ for different energy channels)
        if nblocks > 0:
            #Energy Channel
            energy_channel = ccmc_json.return_json_value_by_index(data,ek_id,0)
            if energy_channel == cfg.errval: continue
            energy_key = ccmc_json.energy_channel_to_key(energy_channel)

            #Observation Window
            obs_win = ccmc_json.return_json_value_by_index(data,id_pred_win,0)
            str_obs_win_st = obs_win["start_time"]
            obs_st = ccmc_json.zulu_to_time(str_obs_win_st)
            
            str_obs_win_end = obs_win["end_time"]
            obs_end = ccmc_json.zulu_to_time(str_obs_win_end)

            win_st.append(obs_st)
            win_end.append(obs_end)

    #SORT THE OBSERVATION WINDOWS IN TIME ORDER
    sort_idx = [sorted(win_st).index(x) for x in win_st]
    sort_obs_st = [win_st[sort_idx[ix]] for ix in sort_idx]
    sort_obs_end = [win_end[sort_idx[ix]] for ix in sort_idx]

    return sort_obs_st, sort_obs_end

 
 
#def read_target_win_list(target_dir):
#    """ Read in the list in the target directory that contains
#        the observation window start and end times for observations
#        that have already been processed and are already inside the
#        target directory.
#        
#        The file is:
#            target_dir/observation_windows.txt
#        Format:
#            YYYY-MM-DD HH:MM:SS, YYYY-MM-DD HH:MM:SS
#        
#    """
#    target_st = []
#    target_end = []
#
#    fname = target_dir + "/observation_windows.csv"
#    with open(fname,"r") as file:
#        for line in file:
#            line.strip()
#            if '#' in line: continue
#            if line == '': continue
#            line = line.split(',')
#            stdate = line[0].strip()
#            enddate = line[1].strip()
#            
#            target_st.append(datetime.datetime.strptime(stdate, "%Y-%m-%d %H:%M:%S"))
#            target_end.append(datetime.datetime.strptime(enddate, "%Y-%m-%d %H:%M:%S"))
# 
#    return target_st, target_end



def identify_new_obs(target_dir, enforce_new=True):
    """ Identify which of the time periods are new and should
        be added to the target_dir. Do not overwrite observations
        that are already in the directory.
    """
    new_st = []
    new_end = []
    
#    obs_st, obs_end = read_batch_list()
    sort_obs_st, sort_obs_end = make_observation_window_list(os.path.join(cfg.outpath,"opsep"))
    
    if not sort_obs_st:
        sys.exit("identify_new_obs: No new observations are available. "
            "Exiting.")


    if not enforce_new:
        return sort_obs_st, sort_obs_end

    target_st, target_end = make_observation_window_list(target_dir)
    
    #Target file is empty, so all observations are considered new
    if not target_st:
        return sort_obs_st, sort_obs_end

    #end of last observation window in the target directory
    last_time = max(target_end)

    for i in range(len(sort_obs_st)):
        if sort_obs_st[i] >= last_time:
            new_st.append(sort_obs_st[i])
            new_end.append(sort_obs_end[i])
            
    return new_st, new_end




def check_for_sep(path):
    """ Read in json files output by OpSEP and check if SEP events
        are present.
        
        json files are in outpath/opsep.
        
        INPUT:
            
            :stdates: (datetime array) start of observation window and
            time that is in the filename of the jsons
       
       OUTPUT:
        
            :df_sep: (pandas dataframe) list of SEP dates for associated
                energy channels and thresholds
            
    """
    ek_id = keys.id_energy_channel
    thresh_id = keys.id_threshold_crossings
    id_pred_win = keys.id_prediction_window
    
    dict = {"Energy Channel": [],
            "Threshold": [],
            "Threshold Crossing Time": [],
            "Observation Window Start": [],
            "Observation Window End": [],
            "SEP Occurred": []}
    

    onlyfiles = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    jsonfiles = [f for f in onlyfiles if '.json' in f]

    for fname in jsonfiles:
        data = ccmc_json.read_in_json(fname)
        if data == None:
            continue
        nblocks = ccmc_json.return_nforecasts(data)
        for i in range(nblocks):
            #Energy Channel
            energy_channel = ccmc_json.return_json_value_by_index(data,ek_id,i)
            if energy_channel == cfg.errval: continue
            energy_key = ccmc_json.energy_channel_to_key(energy_channel)

            #Observation Window
            obs_win = ccmc_json.return_json_value_by_index(data,id_pred_win,i)
            str_obs_win_st = obs_win["start_time"]
            obs_st = ccmc_json.zulu_to_time(str_obs_win_st)
            
            str_obs_win_end = obs_win["end_time"]
            obs_end = ccmc_json.zulu_to_time(str_obs_win_end)

            #Threshold Crossings
            thresh_crossings = ccmc_json.return_json_value_by_index(data,thresh_id,i)
            if thresh_crossings == cfg.errval:
                dict["Energy Channel"].append(energy_key)
                dict["Threshold"].append(None)
                dict["Threshold Crossing Time"].append(pd.NaT)
                dict["Observation Window Start"].append(obs_st)
                dict["Observation Window End"].append(obs_end)
                dict["SEP Occurred"].append(False)

            else:
                for tc in thresh_crossings:
                    threshold = {"threshold": tc["threshold"],
                                "threshold_units": tc["threshold_units"]}
                    thresh_key = ccmc_json.threshold_to_key(threshold)
                    
                    str_cross_time = tc["crossing_time"]
                    cross_time = ccmc_json.zulu_to_time(str_cross_time)
                    if pd.isnull(cross_time) or cross_time == 0\
                        or cross_time == '':
                        continue
                    
                    dict["Energy Channel"].append(energy_key)
                    dict["Threshold"].append(thresh_key)
                    dict["Threshold Crossing Time"].append(cross_time)
                    dict["Observation Window Start"].append(obs_st)
                    dict["Observation Window End"].append(obs_end)
                    dict["SEP Occurred"].append(True)


    df_sep = pd.DataFrame(dict)
    df_sep = df_sep.sort_values(by=["Observation Window Start"], ascending=[True])
    
    return df_sep




def read_approved_sep(target_dir):
    """ Read in the (manually) approved SEP events  in the
        list in the target directory.
        
        target_dir/approved_SEP.csv (csv of pandas DataFrame)
        
        DataFrame columns:
        Energy Channel, Threshold, Threshold Crossing Time, Observation Window Start, Observation Window End
    """
    fname = os.path.join(target_dir, "approved_SEP.csv")
    df = pd.read_csv(fname)
    df['Threshold Crossing Time'] = pd.to_datetime(df['Threshold Crossing Time'])
    df['Observation Window Start'] = pd.to_datetime(df['Observation Window Start'])
    df['Observation Window End'] = pd.to_datetime(df['Observation Window End'])

    return df



def check_target(target_dir):
    """ Check target_dir exists and other files.
    
    """
    if not os.path.isdir(target_dir):
        sys.exit("check_target: Specified target directory does not exist. "
            + target_dir + " Please create or check path and rerun. Exiting.")

    plot_dir = os.path.join(target_dir, "plots")
    if not os.path.isdir(plot_dir):
        print("check_target: plots subdirectory not found. Creating. "
                + plot_dir)
        os.makedirs(plot_dir)
       

    fname = os.path.join(target_dir, "observation_windows.csv")
    if not os.path.isfile(fname):
        print("check_target: Cannot find " + fname + ". Creating.")
        file = open(fname,"w+")
        file.write("#Observation Window Start, Observation Window End\n")
        file.close()

    fname = os.path.join(target_dir, "approved_SEP.csv")
    if not os.path.isfile(fname):
        print("check_target: Cannot find " + fname + ". Creating and populating with any SEP events already in directory.")
        
        df = check_for_sep(target_dir)
        df = df[["Energy Channel","Threshold","Threshold Crossing Time", "Observation Window Start", "Observation Window End"]]
        df = df.dropna()
        df.to_csv(fname,index=False)



def check_approved_sep(target_dir, df_sep, df_approved, obs_st, enforce_sep_stop=True):
    """ Check if there are SEP events in the new observations and 
        determine if they can be moved.
        
    """
    if df_approved.empty:
        if enforce_sep_stop:
            sys.exit("move_obs: The observation with observation "
                "window starting at " + str(obs_st) + " contains an "
                "SEP event with threshold crossing time\n "
                + str(df) + ". "
                "There are NO EVENTS in the approved SEP file "
                + target_dir + "/approved_SEP.csv. Outputs cannot "
                "be moved until the SEP event has been approved. "
                "Exiting.")


    #Check if there is an SEP event in this observation
    df = df_sep.loc[(df_sep["Observation Window Start"] == obs_st)]
    
    sep_time = pd.NaT
    sep = pd.NaT
    if not df.empty:
        sep_time = df["Threshold Crossing Time"]
      
        for j in range(len(df)):
            sep = df["Threshold Crossing Time"].iloc[j]
            if pd.isnull(sep):
                continue
            
            df_check = df_approved.loc[(df_approved["Threshold Crossing Time"] == str(sep))]
            if df_check.empty:
                if enforce_sep_stop:
                    sys.exit("move_obs: The observation with observation "
                        "window starting at " + str(obs_st) + " contains "
                        "an SEP event with threshold crossing time\n "
                        + str(df) + "\n "
                        "This SEP " +str(sep) +" is not in the approved SEP file "
                        + target_dir + "/approved_SEP.csv. Outputs "
                        "cannot be moved until the SEP event has been "
                        "approved. Exiting.")
            else:
                print("SEP Event is approved or no SEP present: " + str(sep))



def move_files(target_dir, obspath, pltpath, allfiles, allplots, df_sep, obs_st, obs_end):
    """ Move the approved files to the target directory.
    
    """
    #Move the quiet and approved observations
    zulu_st = ccmc_json.make_ccmc_zulu_time(obs_st)
    zulu_st = zulu_st.replace(':','')
    
    selectfiles = [f for f in allfiles if zulu_st in f]
    selectplots = [f for f in allplots if zulu_st in f]
    
    
    df = df_sep.loc[(df_sep["Observation Window Start"] == obs_st)]
    if not df.empty:
        print(df)
        for sep in df["Threshold Crossing Time"]:
            if not pd.isnull(sep):
                zulu_sep = ccmc_json.make_ccmc_zulu_time(sep)
                zulu_sep = zulu_sep.replace(':','')
                selectplots2 = [f for f in allplots if zulu_sep in f]
                if selectplots2:
                    selectplots.extend(selectplots2)
            
    success = False
    for file in selectfiles:
        fname = os.path.join(obspath, file)
        targetname = os.path.join(target_dir, file)
        try:
            shutil.move(fname, targetname)
            success = True
            print("move_output: Moved file " + fname)
        except:
            print("move_output: Could not move " + fname
                + ". Continuing.")
            

    for plot in selectplots:
        fname = os.path.join(pltpath, plot)
        targetname = os.path.join(target_dir, "plots", plot)
        try:
            shutil.move(fname, targetname)
            print("move_output: Moved file " + fname)
        except:
            print("move_output: Could not move " + fname
                + ". Continuing.")
            success = False

    #Append observation windows file in target directory with newly
    #added observations windows
    if success:
        target_win_file = open(os.path.join(target_dir, "observation_windows.csv"), "a")
        target_win_file.write(str(obs_st) + "," + str(obs_end) + "\n")
        target_win_file.close()


def cleanup_csv(obspath, allfiles):
    """ Remove unneccessary csv files produced by opsep from obspath.
    
    """
    selectfiles = [f for f in allfiles if '.csv' in f]
    for file in selectfiles:
        fname = os.path.join(obspath, file)
        if os.path.isfile(fname):
            os.remove(fname)



def print_target_info(target_dir):
    #Approved SEP in target directory
    df_approved = read_approved_sep(target_dir)
    print("Approved SEP in " + target_dir + "/approved_SEP.csv")
    print(df_approved)

    target_st, target_end = make_observation_window_list(target_dir)
    
    print("Observation Windows in " + target_dir)
    for i in range(len(target_st)):
        print(f"Start: {target_st[i]}, End: {target_end[i]}")
    


def move_output(target_dir, enforce_new=True, enforce_sep_stop=True):
    """ Move the json and supporting files created by OpSEP to
        the target directory. Move the associated plots into
        a plots/ dir inside the target directory.
        
        if enforce_new = False, then all files will be moved
            to the target directory and write over any files with
            the same name that are already present. (not recommended)
            
        if enforce_sep_stop = False, SEP event files will be allowed
            to be moved into the target directory without requiring
            them to be in the approved list. (not recommended)
            
    """
    
    check_target(target_dir)

    #Check which observations are not already present in the target dir
    new_obs_st, new_obs_end = identify_new_obs(target_dir, enforce_new)
    
    #SEP in fetchsep observations
    obspath = os.path.join(cfg.outpath, "opsep")
    df_sep = check_for_sep(obspath)
    
    #Approved SEP in target directory
    df_approved = read_approved_sep(target_dir)

    #Identify which files to move and move them
    allfiles = [f for f in os.listdir(obspath) if os.path.isfile(os.path.join(obspath, f))]
    
    pltpath = os.path.join(cfg.plotpath, "opsep")
    allplots = [f for f in os.listdir(pltpath) if os.path.isfile(os.path.join(pltpath, f))]

    for i in range(len(new_obs_st)):
        obs_st = new_obs_st[i]
        obs_end = new_obs_end[i]

        print("move_output: Checking files with observations windows "
            + str(obs_st) + " to " + str(obs_end))

        ###SEP CHECKS
        check_approved_sep(target_dir, df_sep, df_approved, obs_st, enforce_sep_stop)

        ###MOVE FILES
        move_files(target_dir, obspath, pltpath, allfiles, allplots, df_sep,
            obs_st, obs_end)

        ###REMOVE LEFTOVER csv FILES
        cleanup_csv(obspath, allfiles)



def update_observations(target_dir, start_date, end_date, experiment,
    flux_type, threshold, nointerp=False, two_peaks=False):
    """ Run opsep for a time period that starts where the last set of observations
        in target_dir end.
        
        The options and settings below are optimized for SEP Scoreboard validation.
    
        INPUTS:
        
            :target_dir: (string) after the opsep files are made, they will be
                transferred to the target_dir
            :end_date: (string) YYYY-MM-DD or "YYYY-MM-DD HH:MM:SS" in quotes, end date
                of the observation time period
            :experiment: (string) spacecraft, eg. GOES_RT, GOES-13, STEREO-A, etc
            :flux_type: (string) integral or differential
            :threshold: (string) additional thresholds to apply to fluxes, see opsep for
                more guidance
            :ninterp: (bool) False means that will not apply interpolation in data gaps
                
        OUTPUTS:
        
            None. Files will be moved if any SEP events in the time frame are approved.
    
    """
    if start_date == '':
        #Get the date of the most recent data in the target directory
        #Use the end of the last prediction window as the start time
        target_st, target_end = make_observation_window_list(target_dir)
        start_date = str(max(target_end))
    
    #Optimal settings for SEP Scoreboard
    showplot = False
    saveplot = True
    
    model_name = ''
    user_file = ''
    json_type = 'observations' #not relevant unless input a user file
    spase_id = ''
    detect_prev_event = False
    umasep = False
    option = ''
    doBGSub = False
    bgstartdate = ''
    bgenddate = ''
    
    
    sep_year, sep_month, \
    sep_day, jsonfname = opsep.run_all(start_date, end_date,
        experiment, flux_type, model_name, user_file, json_type,
        spase_id, showplot, saveplot, detect_prev_event,
        two_peaks, umasep, threshold, option, doBGSub, bgstartdate,
        bgenddate, nointerp)

