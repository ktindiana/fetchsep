from . import config as cfg
from ..json import keys
from ..json import ccmc_json_handler as ccmc_json
import datetime
import sys
import os
import pandas as pd
import shutil

__version__ = "1.0"
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

datapath = cfg.datapath
outpath = cfg.outpath
plotpath = cfg.plotpath

def read_batch_list():
    """ Read in outpath/idsep/batch_event_list.txt which contains
        the observation windows that were used to process time
        periods with OpSEP.
        
    """
    obs_st = []
    obs_end = []
    
    fname = outpath + "/idsep/batch_event_list.txt"
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
 
 
def read_target_win_list(target_dir):
    """ Read in the list in the target directory that contains
        the observation window start and end times for observations
        that have already been processed and are already inside the
        target directory.
        
        The file is:
            target_dir/observation_windows.txt
        Format:
            YYYY-MM-DD HH:MM:SS, YYYY-MM-DD HH:MM:SS
        
    """
    target_st = []
    target_end = []

    fname = target_dir + "/observation_windows.csv"
    with open(fname,"r") as file:
        for line in file:
            line.strip()
            if '#' in line: continue
            if line == '': continue
            line = line.split(',')
            stdate = line[0].strip()
            enddate = line[1].strip()
            
            target_st.append(datetime.datetime.strptime(stdate, "%Y-%m-%d %H:%M:%S"))
            target_end.append(datetime.datetime.strptime(enddate, "%Y-%m-%d %H:%M:%S"))
 
    return target_st, target_end



def identify_new_obs(target_dir, enforce_new=True):
    """ Identify which of the time periods are new and should
        be added to the target_dir. Do not overwrite observations
        that are already in the directory.
    """
    new_st = []
    new_end = []
    
    obs_st, obs_end = read_batch_list()
    if obs_st == []:
        sys.exit("identify_new_obs: No new observations are available. "
            "Exiting.")

    #SORT THE OBSERVATION WINDOWS IN TIME ORDER
    sort_idx = [sorted(obs_st).index(x) for x in obs_st]
    sort_obs_st = [obs_st[sort_idx[ix]] for ix in sort_idx]
    sort_obs_end = [obs_end[sort_idx[ix]] for ix in sort_idx]

    if not enforce_new:
        return sort_obs_st, sort_obs_end

    target_st, target_end = read_target_win_list(target_dir)
    
    #Target file is empty, so all observations are considered new
    if target_st == []:
        return sort_obs_st, sort_obs_end

    #end of last observation window in the target directory
    last_time = max(target_end)

    for i in range(len(sort_obs_st)):
        if sort_obs_st[i] >= last_time:
            new_st.append(sort_obs_st[i])
            new_end.append(sort_obs_end[i])
            
    return new_st, new_end



def check_for_sep():
    """ Read in json files output by OpSEP and check if SEP events
        are present.
        
        json files are in outpath/opsep.
        
        INPUT:
            
            :stdates: (datetime array) start of observation window and
            time that is in the filename of the jsons
       
       OUTPUT:
        
            :sepdates: (pandas dataframe) list of SEP dates for associated
                energy channels and thresholds
            
    """
    ek_id = keys.id_energy_channel
    thresh_id = keys.id_threshold_crossings
    id_pred_win = keys.id_prediction_window
    
    dict = {"Energy Channel": [],
            "Threshold": [],
            "Threshold Crossing Time": [],
            "Observation Window Start": [],
            "Observation Window End": []}
    
    obspath = outpath + "/opsep"
    onlyfiles = [obspath + "/" + f for f in os.listdir(obspath) if os.path.isfile(os.path.join(obspath, f))]
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
            if thresh_crossings == cfg.errval: continue
            for tc in thresh_crossings:
                threshold = {"threshold": tc["threshold"],
                            "threshold_units": tc["threshold_units"]}
                thresh_key = ccmc_json.threshold_to_key(threshold)
                
                str_cross_time = tc["crossing_time"]
                cross_time = ccmc_json.zulu_to_time(str_cross_time)
                if cross_time == None or cross_time == 0\
                    or cross_time == '':
                    continue
                
                dict["Energy Channel"].append(energy_key)
                dict["Threshold"].append(thresh_key)
                dict["Threshold Crossing Time"].append(cross_time)
                dict["Observation Window Start"].append(obs_st)
                dict["Observation Window End"].append(obs_end)
    
    df_sep = pd.DataFrame(dict)
    return df_sep




def read_approved_sep(target_dir):
    """ Read in the (manually) approved SEP events  in the
        list in the target directory.
        
        target_dir/approved_SEP.csv (csv of pandas DataFrame)
    """
    fname = target_dir + "/approved_SEP.csv"
    df_approved = pd.read_csv(fname)
    
    return df_approved



def check_target(target_dir):
    """ Check target_dir exists and other files.
    """
    if not os.path.isdir(target_dir):
        sys.exit("check_target: Specified target directory does not exist. "
            + target_dir + " Please create or check path and rerun. Exiting.")

    plot_dir = target_dir +"/plots"
    if not os.path.isdir(plot_dir):
        print("check_target: plots subdirectory not found. Creating. "
                + plot_dir)
        os.makedirs(plot_dir)
       

    fname = target_dir + "/observation_windows.csv"
    if not os.path.isfile(fname):
        print("check_target: Cannot find " + fname + ". Creating.")
        file = open(fname,"w+")
        file.write("#Observation Window Start, Observation Window End\n")
        file.close()

    fname = target_dir + "/approved_SEP.csv"
    if not os.path.isfile(fname):
        print("check_target: Cannot find " + fname + ". Creating.")
        
        dict = {"Energy Channel": [],
            "Threshold": [],
            "Threshold Crossing Time": [],
            "Observation Window Start": [],
            "Observation Window End": []}
        df = pd.DataFrame(dict)
        df.to_csv(fname)


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
    df_sep = check_for_sep()
    
    #Approved SEP in target directory
    df_approved = read_approved_sep(target_dir)
    print("Approved SEP in " + target_dir + "/approved_SEP.csv")
    print(df_approved)

    #Identify which files to move and move them
    obspath = outpath + "/opsep"
    allfiles = [f for f in os.listdir(obspath) if os.path.isfile(os.path.join(obspath, f))]
    
    pltpath = plotpath + "/opsep"
    allplots = [f for f in os.listdir(pltpath) if os.path.isfile(os.path.join(pltpath, f))]


    for i in range(len(new_obs_st)):
        obs_st = new_obs_st[i]
        obs_end = new_obs_end[i]

        print("move_output: Checking files with observations windows "
            + str(obs_st) + " to " + str(obs_end))

        ###SEP CHECKS
        #Check if there is an SEP event in this observations
        df = df_sep.loc[(df_sep["Observation Window Start"] == obs_st)]
        
        sep_time = None
        sep = None
        if not df.empty:
            sep_time = df["Threshold Crossing Time"]
            if df_approved.empty:
                if enforce_sep_stop:
                    sys.exit("move_obs: The observation with observation "
                        "window starting at " + str(obs_st) + " contains an "
                        "SEP event with threshold crossing time "
                        + str(sep_time) + ". "
                        "There are NO EVENTS in the approved SEP file "
                        + target_dir + "/approved_SEP.csv. Outputs cannot "
                        "be moved until the SEP event has been approved. "
                        "Exiting.")
            
            
            
            for j in range(len(df)):
                sep = df["Threshold Crossing Time"].iloc[j]
                df_check = df_approved.loc[(df_approved["Threshold Crossing Time"] == str(sep))]
                print(df_check)
                if df_check.empty:
                    if enforce_sep_stop:
                        sys.exit("move_obs: The observation with observation "
                            "window starting at " + str(obs_st) + " contains "
                            "an SEP event with threshold crossing time "
                            + str(sep_time) + ". "
                            "This SEP is not in the approved SEP file "
                            + target_dir + "/approved_SEP.csv. Outputs "
                            "cannot be moved until the SEP event has been "
                            "approved. Exiting.")

        ###MOVE FILES
        #Move the quiet and approved observations
        zulu_st = ccmc_json.make_ccmc_zulu_time(obs_st)
        zulu_st = zulu_st.replace(':','')
        zulu_sep = ccmc_json.make_ccmc_zulu_time(sep)
        
        selectfiles = [f for f in allfiles if zulu_st in f]
        selectplots = [f for f in allplots if zulu_st in f]
        
        if zulu_sep != None:
            zulu_sep = zulu_sep.replace(':','')
            selectplots2 = [f for f in allplots if zulu_sep in f]
            if selectplots2:
                selectplots.extend(selectplots2)
                
        success = False
        for file in selectfiles:
            fname = obspath + "/" + file
            targetname = target_dir + "/" + file
            try:
                shutil.move(fname, targetname)
                success = True
            except:
                print("move_output: Could not move " + fname
                    + ". Continuing.")

        #Append observation windows file in target directory with newly
        #added observations windows
        if success:
            target_win_file = open(target_dir + "/observation_windows.csv", "a")
            target_win_file.write(str(obs_st) + "," + str(obs_end) + "\n")
            target_win_file.close()

        
        for plot in selectplots:
            fname = pltpath + "/" + plot
            targetname = target_dir + "/plots/" + plot
            try:
                shutil.move(fname, targetname)
            except:
                print("move_output: Could not move " + fname
                    + ". Continuing.")
