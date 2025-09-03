from . import config as cfg
from . import tools
from ..json import ccmc_json_handler as ccmc_json
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.dates import DateFormatter
import math
from sklearn.utils.validation import check_consistent_length
from sklearn.utils.validation import check_array
import datetime
from scipy.stats import pearsonr
from math import log10
from pandas.plotting import register_matplotlib_converters
import warnings
import os

__author__ = "Phil Quinn, Kathryn Whitman"
__maintainer__ = "Phil Quinn"
__email__ = "philip.r.quinn@nasa.gov"

'''Contains functions for plotting data from forecasting
    models and observations.
    Written on 2020-07-17.

    Phil Quinn may be reached at philip.r.quinn@nasa.gov.
    Kathryn Whitman may be reached at kathryn.whitman@nasa.gov.
'''

#Changes in 0.5: correlation_plot subroutine added by K. Whitman
#2021-09-03, changes in 0.6: Set a limiting value or 1e-4 on the
#   axes in correlation_plot



def plot_time_profile(date, values, labels, dy=None, dyl=None,
                dyh=None, title=None, x_min=None, x_max=None,
                x_label="Date", y_min=None, y_max=None,
                y_label="Value", uselog_x = False, uselog_y = False,
                date_format="year", save="time_profile",
                showplot=False, closeplot=False, saveplot=False):
    """
    Plots multiple time profiles in same plot

    Parameters
    ----------
    date : array-like datetime objects, shape=(n profiles, n dates)
        Datetimes
        [[dates1,dates2,dates3,...],[dates1,dates2,dates3....],...]

    values : array-like float, shape=(n profiles, n dates)
        Metric values as function of datetime
        [[val1, val2, val3,...],[val1,val2,val3,...],...]

    labels : array-like string, shape=(n profiles, n dates)
        Labels for the different time profiles

    title : string
        Title for plot
        Optional

    x_min : datetime object
        Minimum datetime for x-axis
        Optional

    x_max : datetime object
        Maximum datetime for x-axis
        Optional

    x_label : string
        Label for x-axis
        Optional. Defaults to "Date"
        
    y_min : float
        Minimum for y-axis
        Optional

    y_max : float
        Maximum for y-axis
        Optional

    y_label : string
        Label for y-axis
        Optional. Defaults to "Metric"
        
    date_format : string
        May be "year" or "day" or "none"
        Default year
        Determines format of date on x-axis

    save : string
        Name to save PNG as (should not include ".png")
        Optional. Defaults to "metric_profile"

    showplot : boolean
        Indicator for displaying the plot on screen or not
        Optional. Defaults to False

    closeplot : boolean
        Indicator for clearing the figure from memory
        Optional. Defaults to False

    Returns
    -------
    fig, figname
    """

    register_matplotlib_converters()

    #check_consistent_length(date, metric)

    # checking if items in date are datetime objects
    #if not all(isinstance(x, datetime.datetime) for x in date):
    #    raise TypeError("Dates must be datetime objects.")

    plt.style.use('seaborn-whitegrid')

    fig = plt.figure(figsize=(13,8))
    ax = plt.subplot(111)
    
    #colors = plt.cm.tab10(np.linspace(0,1,len(values)+1))
    
    #color_metric = '#247afd'
    color_nans = '#05ffa6'
    markers = ["o","v","^","<",">","s","P","X","D","d","p","H",".","x","*","p"]

    for i in range(len(date)): #number of time profiles
        if 0 in values[i]:
            y_values = np.array(values[i])
            values[i] = np.ma.masked_where(y_values <= 0 , y_values)
#            values[i] = np.ma.masked_invalid(y_values)

        # getting indices where the metric is nan or +/-inf
        #indices = [j for j, arr in enumerate(metric[i]) if not np.isfinite(arr).all()]

        # replacing nan and +/-inf with None types so the results are still plottable
        #values[i] = [None if np.isnan(x) else x for x in values[i]]
        #values[i] = [None if x==np.inf else x for x in values[i]]
        if dy==None:
            if "REleASE" in labels[i]:
                ax.plot(date[i], values[i], label=labels[i], marker=".", linestyle=":")
            #elif "GOES" in labels[i]:
            #    ax.plot(date[i], values[i], label=labels[i], color="k")
            elif len(date[i]) == 1:
                ax.plot(date[i], values[i], label=labels[i], marker="D")#, fillstyle='none', #markeredgewidth=2)
            elif "IMP" in labels[i]:
                ax.plot(date[i], values[i], label=labels[i], linestyle="dashed")
            else:
                ax.plot(date[i], values[i], label=labels[i])
        else:
            if "REleASE" in labels[i]:
                ax.errorbar(date[i], values[i], label=labels[i], yerr=dy[i], marker=".", linestyle=":", elinewidth=2)
            #elif "GOES" in labels[i]:
            #    ax.plot(date[i], values[i], label=labels[i], color="k")
            elif len(date[i]) == 1:
                ax.errorbar(date[i], values[i], label=labels[i], yerr=dy[i], marker="D",  elinewidth=2, capsize=4)#, #fillstyle='none', markeredgewidth=2)
            elif "ASPECS" in labels[i]:
                ax.errorbar(date[i], values[i], label=labels[i], yerr=dy[i], linestyle="dashed", elinewidth=2)
            else:
                ax.errorbar(date[i], values[i], label=labels[i], yerr=dy[i], elinewidth=2)

        #ax.axhline(mean[i])

    # plotting vertical dashed lines where the metric is nan or +/-inf
    #for i in indices:
    #    plt.axvline(x=date[i], color=color_nans, linestyle='--')

    ax.grid(True, linestyle=':', alpha=0.5)

    if x_min != None and x_max != None:

        if not isinstance(x_min, datetime.datetime):
            raise TypeError("x_min must be datetime object.")
        if not isinstance(x_max, datetime.datetime):
            raise TypeError("x_max must be datetime object.")
        ax.set_xlim(x_min, x_max)
        
    if y_min != None and y_max != None:
        ax.set_ylim(y_min, y_max)

    ax.set(xlabel=x_label, ylabel=y_label)
    if date_format == "year" or date_format == "Year":
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    if date_format == "day" or date_format == "Day":
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d\n%H:%M'))
    
    plt.setp(ax.get_xticklabels(), rotation = 15)
    ax.set_title(title)
#    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
 #       item.set_fontsize(24)
    if uselog_x:
        ax.set_xscale('log')
    if uselog_y:
        ax.set_yscale('log')

    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.87,
                    chartBox.height])
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 0.95), fontsize='9', \
              framealpha=0.5)
        

    figname = save +'.png'
    if closeplot: plt.close(fig)

    return fig, figname




def correlation_plot(obs_values, model_values, plot_title, \
                        xlabel="Observations", ylabel="Model",  value="Value",
                        use_log = False, use_logx = False, use_logy = False):
    '''Make a correlation plot of two arrays.
        
        obs_values (1D array of floats) for x-axis
        
        model_values (1D array of floats) for y-axis)
        
        plot_title (string) is title for plot
        
        value (string) indicates which value you are comparing, e.g. Peak Flux
            (not used)
            
        returns plt
    '''
    corr = 0
    obs_np = []
    model_np = []
    slope = 0
    yint = 0
    
    obs_np = np.array(obs_values)
    model_np = np.array(model_values)
    
    if use_log or use_logx:
        log_obs = [log10(val) for val in obs_values]
        obs_np = np.array(log_obs)

    if use_log or use_logy:
        log_model = [log10(val) for val in model_values]
        model_np = np.array(log_model)
    

    #CORRELATION
    corr, _ = pearsonr(obs_np, model_np)

    #LINEAR REGRESSION
    slope, yint = np.polyfit(obs_np, model_np, 1)

    #1-to-1 Line
    mx = max(np.amax(obs_np),np.amax(model_np))
    mn = min(np.amin(obs_np), np.amin(model_np))
    if use_log or use_logx or use_logy:
        mn = max(-4,mn)
    step = (mx - mn)/10.
    x1to1 = np.arange(mn, mx, step).tolist()
    y1to1 = x1to1


    ######MAKE CORRELATION PLOT########
    plt.figure(figsize=(8,5))
    ax = plt.subplot(111)
    plt.style.use('default')
    plt.grid(which="both", axis="both")
    plt.title(plot_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    ax.set_ylim([mn-0.1*abs(mn),mx+0.1*mx])
    ax.set_xlim([mn-0.1*abs(mn),mx+0.1*mx])

    ax.plot(obs_np, model_np, 'bo', \
                label=(f'Pearsons Correlation \nCoefficient: ' \
                        + ' {0:.3f}'.format(corr)))
    ax.plot(np.sort(obs_np), slope*np.sort(obs_np) + yint,\
                color='red', label=(f'Linear Regression \nSlope: '+ \
                '{0:.3f} \ny-intercept: {1:.3f}'.format(slope, yint)))
    ax.plot(x1to1, y1to1, color='black', label="1:1 Line", linestyle="dashed")

    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.75,
                    chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.25, 0.95))

    return plt



def box_plot(values, labels, x_label="Model", y_label="Metric",
             title=None, save="boxes", uselog=False, showplot=False,
             closeplot=False):
    """
    Plots ratio or skill score for vary thresholds
    and for each model subtype

    Parameters
    ----------
    values : array-like float, shape=(n model subtypes, n values)
        Values to plot for each model subtype

    labels : array-like string, shape=(n model subtypes, n labels)
        Labels of the model subtype

    x_label : string
        Label for x-axis
        Optional. Defaults to "Model"

    y_label : string
        Label for y-axis
        Optional. Defaults to "Metric"

    title : string
        Title for plot
        Optional

    save : string
        Name to save PNG as (should not include ".png")
        Optional. Defaults to "boxes"

    showplot : boolean
        Indicator for displaying the plot on screen or not
        Optional. Defaults to False

    closeplot : boolean
        Indicator for clearing the figure from memory
        Optional. Defaults to False

    Returns
    -------
    None
    """
    if len(values) <= 4:
        fig = plt.figure(figsize=(9, 6))
    if len(values) > 4 and len(values) <= 7:
        fig = plt.figure(figsize=(12, 6))
    if len(values) > 7:
        fig = plt.figure(figsize=(16, 6))
    ax = fig.add_subplot(111)

    sns.boxplot(data=values, fliersize=0, meanline=True, showmeans=True, \
                medianprops = {'color': 'w', 'linewidth': 1},
                meanprops = {'color': 'k', 'linewidth': 1})

    means = [np.mean(list) for list in values]
    medians = [np.median(list) for list in values]

    for i in range(len(values)):
        if means[i] == max(means[i], medians[i]):
            vmean = 'bottom'
            vmed = 'top'
        else:
            vmean = 'top'
            vmed = 'bottom'
        ax.text(i, means[i], "\u03BC=" + str(np.round(means[i], 2)), size='large', \
                color='k', weight='semibold', horizontalalignment='center', \
                verticalalignment=vmean)
        ax.text(i, medians[i], "M=" + str(np.round(medians[i], 2)), size='large', \
                color='b', weight='semibold', horizontalalignment='center', \
                verticalalignment=vmed)

    sns.stripplot(data=values, linewidth=0.5)

    ax.set_title(title)
    #ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xticklabels(labels, rotation=45)

    if uselog:
        ax.set_yscale('log')

    return fig




############ OPSEP PLOTS ##############
def check_for_good_fluxes(flux):
    """ Check that there are some non-zero points. """
    is_good = False
    positive = [flx for flx in flux if not pd.isnull(flx)]
    positive = [flx for flx in positive if flx > 0]
    if len(positive) > 0:
        is_good = True
    
    return is_good
    
def opsep_plot_bgfluxes(unique_id, experiment, flux_type, options, user_name,
    fluxes, dates, energy_bins, means, sigmas, saveplot,
    spacecraft='', doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False):
    """Plot fluxes with time for all of the energy bins on the same plot. The
        estimated mean background levels are plotted as dashed lines.
        Zero values are masked, which is useful when make plots of the
        background and SEP flux separately.
        
        INPUTS:
        
        :experiment: (string) name of experiment
        :flux_type: (string) "integral" or "differential"
        :options: (string array) array of options that can be applied
        :user_name: (string) name of data source is user input file
        :fluxes: (float nxm array) flux time profiles for n energy channels and
            m time points
        :dates: (datetime 1xm array) m time points for flux time profile
        :energy_bins: (float nx2 array) energy bins for n energy channels
        :means: (float 1xn array) mean values of histogram for n energy channels
        :sigmas: (float 1xn array) sigma of histograms for n energy channels
        :saveplot: (bool) True to save plot automatically
        
        OUTPUTS:
        
        Plot to screen and plot saved to file
        
    """
    #All energy channels in specified date range with event start and stop
    #Plot all channels of user specified data
    #Additions to titles and filenames according to user-selected options
    suffix = f"{unique_id}_All_Bins"
    figname, subdir = tools.opsep_naming_scheme(dates[0], suffix, experiment, flux_type, user_name, options,
        spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
    exp_name = experiment
    if experiment == "user":
        exp_name = user_name

    fig = plt.figure(figname,figsize=(13.5,6))
    ax = plt.subplot(111)
    nbins = len(energy_bins)
    for i in range(nbins):
        #Check to make sure there are positive fluxes
        is_good = check_for_good_fluxes(fluxes[i])
        if not is_good:
            continue
        #Don't want to plot zero values, particularly in
        #background-subtracted plots
        maskfluxes = np.ma.masked_where(fluxes[i] <= 0, fluxes[i])
#        maskfluxes = np.ma.masked_invalid(maskfluxes)

        mean =  means[i]
        if isinstance(means[i], list):
            mean = np.nanmean(np.array(means[i]))
            
        legend_label = tools.setup_energy_bin_label(energy_bins[i])
        p = ax.plot_date(dates,maskfluxes, '-', label=legend_label)
        color = p[0].get_color()
        if i==0:
            plt.axhline(mean,color=color,linestyle=':', label="Mean Background")
        else:
            plt.axhline(mean,color=color,linestyle=':')


    if flux_type == "integral": flux_units = cfg.flux_units_integral
    if flux_type == "differential": flux_units = cfg.flux_units_differential

    fig.suptitle((f"{unique_id}\n {exp_name} {title_mod} {flux_type}"))
    
    #Formatting of axes
    flux_units = make_math_label(flux_units)
    ax.set_ylabel((f"Flux ({flux_units})"))
    ax.set_xlabel("Date")

    plt.yscale("log")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                     chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))
    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath, 'opsep', subdir, figname + '.png'))



def plot_weibull_fit(energy_bin, threshold, experiment, flux_type, user_name,
    options, sep_start_time, trim_times, trim_fluxes,
    best_pars, best_fit, max_time, max_val, max_meas_time, max_meas,
    max_curve_model_time, max_curve_model_peak, max_curve_meas_time, max_curve_meas_peak,
    saveplot, showplot, spacecraft='', doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False):
    """ Plot Weibull fit used to get onset peak """

    best_a = best_pars['alpha']
    best_b = best_pars['beta']
    best_Ip = best_pars['peak_intensity']

    suffix = f"Weibull_profile_fit_{energy_bin[0]}MeV"
    figname, subdir = tools.opsep_naming_scheme(sep_start_time, suffix, experiment, flux_type, user_name,
        options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
    exp_name = experiment
    if experiment == "user":
        exp_name = user_name

    fig = plt.figure(figname,figsize=(9,5))
    label = f"{energy_bin[0]} - {energy_bin[1]} MeV"
    if energy_bin[1] == -1:
        label = f">{energy_bin[0]} MeV"
    plt.plot(trim_times,trim_fluxes,label=label,marker='.', linestyle='none')
    label_fit = (f"Fit\n Ip: {best_Ip:.2f}"
                f"\n alpha: {best_a:.2f}"
                f"\n beta: {best_b:.2f}")
    plt.plot(trim_times,best_fit,label=label_fit)
    plt.plot(max_time, max_val,"o",label="max fit")
    plt.plot(max_meas_time, max_meas,">",label="measured peak near fit max")
    plt.plot(max_curve_model_time, max_curve_model_peak,"D",label="Min 2nd Derivative")
    plt.plot(max_curve_meas_time, max_curve_meas_peak,"^",label="measured peak near min 2nd Derivative")
    #plt.plot(onset_time, onset_peak[i],">",label="onset Weibull")
    plt.legend(loc='lower right')
    plt.title(f"Onset Peak Weibull Fit\n {exp_name} {title_mod} {flux_type}\n {sep_start_time}")
    plt.xlabel("Hours")
    plt.ylabel("Flux")
    plt.yscale("log")
    plt.ylim(1e-4,1e5)
    
    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath, 'opsep', subdir, figname + '.png'))
    if not showplot:
        plt.close(fig)


####PLOT STUFF FOR find_max_curvature in tools.py
#    if showplot or saveplot:
#        tzulu = ccmc_json.make_ccmc_zulu_time(crossing_time)
#        tzulu = tzulu.replace(":","")
#        figname = tzulu + "_" + "SecondDerivative_" + experiment + "_"\
#                + str(energy_threshold) + "MeV"
#        if spacecraft:
#            figname = tzulu + "_" + "SecondDerivative_" + experiment + "_"\
#                + spacecraft + "_" + str(energy_threshold) + "MeV"
#        fig = plt.figure(figname,figsize=(9,5))
#        plt.plot(xarr,yarr,label="orig")
#        plt.plot(xarr[max_k_idx+2], yarr[max_k_idx+2],"o",label="Min 2nd Derivative on Fit")
#        plt.plot(xarr[max_y_idx], yarr[max_y_idx],"o",label="max Fit")
#        plt.plot(xarr[2:],k_x,label="2nd Derivative")
#        plt.plot(xarr[max_k_idx+2],k_x[max_k_idx],"o",label="Min 2nd Derivative")
#        plt.legend(loc='lower left')
#        plt.xlabel("Hours")
#        plt.ylabel("Y")
#        #plt.yscale("log")
#        #plt.ylim(1e-4,1e6)
#        if saveplot:
#            fig.savefig(plotpath + '/' + figname + '.png')
#        if not showplot:
#            plt.close(fig)


def make_math_label(label):
    label = label.replace("^-1", "^{-1}")
    label = label.replace("^-2", "^{-2}")
    label = label.replace("*", "")
    return label

def opsep_plot_event_definitions(experiment, flux_type, user_name, options,
    evaluated_dates, evaluated_fluxes, evaluated_energy_bins, event_definitions,
    sep_start_times, sep_end_times, onset_peaks, onset_peak_times,
    max_fluxes, max_flux_times, showplot, saveplot, spacecraft='',
    doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False):
    """ Plot the fluxes used for event definitions with threshold,
        start and end times, onset peak and max flux.
        
        INPUT:
        
            :experiment: (str) e.g. "GOES-13" or "user" for user input
            :flux_type: (str) integral or differential
            :user_name: (str) if user input or if specified, will be used
                instead of experiments
            :options: (arr of str) any options applied to fluxes/energy bings
            :doBGSub: (bool) was background subtraction applied?
            :evaluated_dates: (arr) dates for the fluxes that were evaluated
                using event definitions, called evaluated_dates in Data obj
            :evaluated_fluxes: (2D arr, probably numpy) all fluxes that were
                evaluted with event definitions, called evaluated_fluxes in Data obj
            :event_definitions: (2D list of dict) Dictionaries of EnergBin and 
                Threshold objects that identify the applied event definitions
            :sep_start_dates: (list of datetime) start of SEP events
                in the same order of the event_definitions list 
            :sep_end_dates: (list of datetime) end of SEP events
                in the same order of the event_definitions list
            :spacecraft: (str) primary or secondary if specified
            
        OUTPUT:
        
            Multi-pane plot showing event definitions
    
    """
 
    suffix = "Event_Def"
    figname, subdir = tools.opsep_naming_scheme(evaluated_dates[0], suffix, experiment, flux_type, user_name,
        options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
 
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
    exp_name = experiment
    if experiment == "user":
        exp_name = user_name
 
    #Plot selected results
    #Event definition from integral fluxes
    if flux_type == "differential":
        print("Generating figure of estimated integral fluxes with threshold crossings.")
    if flux_type == "integral":
        print("Generating figure of integral fluxes with threshold crossings.")

    #plot integral fluxes (either input or estimated)
    nthresh = len(event_definitions)
    energy_units = event_definitions[0]['energy_channel'].units

    if nthresh > 4:
        #fig = plt.figure(figname,figsize=(12,12))
        fig, ax = plt.subplots(nthresh, 1, sharex=True, figsize=(12,12))
    else:
        #fig = plt.figure(figname,figsize=(12,9))
        fig, ax = plt.subplots(nthresh, 1, sharex=True, figsize=(12,9))
        if nthresh == 1: ax = [ax]

    plot_title = f"Event Definitions\n {exp_name} {title_mod} {flux_type} Fluxes"
    plt.suptitle(plot_title)

    for i in range(nthresh):
        #Get energy bin
        energy_bin = [event_definitions[i]['energy_channel'].min,
                    event_definitions[i]['energy_channel'].max]
        threshold = event_definitions[i]['threshold'].threshold
        flux_units = event_definitions[i]['threshold'].threshold_units
                    
        #Event definitions and fluxes are in the same order
        dates = evaluated_dates #for ease
        fluxes = np.array(evaluated_fluxes[i])
        
        is_good = check_for_good_fluxes(fluxes)
        if not is_good:
            continue
        
        #Create labels
        ylabel = f"Flux [${flux_units}$]"
        ylabel = make_math_label(ylabel)
        if energy_bin[1] == -1 and flux_type == "integral":
            data_label = f"{exp_name} >{energy_bin[0]} {energy_units}"
        elif energy_bin[1] == -1 and flux_type == "differential":
            data_label = f"{exp_name} Estimated >{energy_bin[0]} {energy_units}"
        else:
            data_label = f"{exp_name} {energy_bin[0]}-{energy_bin[1]} {energy_units}"

    
#        maskfluxes = np.ma.masked_invalid(fluxes)
        #Don't want to plot negative values, particularly in background-subtracted plots
#        if doBGSubOPSEP or doBGSubIDSEP:
        maskfluxes = np.ma.masked_where(fluxes <= 0, fluxes)
        ax[i].plot_date(dates,maskfluxes,'-',label=data_label,marker=".")
#        else:
#            ax[i].plot_date(dates,maskfluxes,'-',label=data_label)

        if not pd.isnull(sep_start_times[i]):
            ax[i].axvline(sep_start_times[i],color='black',linestyle=':')
            ax[i].axvline(sep_end_times[i],color='black',linestyle=':',
                        label="Start, End")
        ax[i].axhline(threshold,color='red',linestyle=':', label="Threshold")
        if not pd.isnull(onset_peaks[i]) and not pd.isnull(onset_peak_times[i]):
            ax[i].plot_date(onset_peak_times[i],onset_peaks[i],'o',color="black",
                    label="Onset Peak")
        if not pd.isnull(max_fluxes[i]) and not pd.isnull(max_flux_times[i]):
            ax[i].plot_date(max_flux_times[i],max_fluxes[i],'ro',mfc='none',
                    label="Max Flux")


        if i == nthresh-1: ax[i].set_xlabel('Date')
        ax[i].set_ylabel(ylabel)
#        if sum(maskfluxes) > 0: #If NaN present, returns False
        ax[i].set_yscale("log")
        #ymin = max(1e-6, min(integral_fluxes[i]))
        # plt.ylim(ymin, peak_flux[i]+peak_flux[i]*.2)
        ax[i].legend(loc='upper right')
#        for item in ([ax[i].title, ax[i].xaxis.label, ax[i].yaxis.label] + ax[i].get_xticklabels() + ax[i].get_yticklabels()):
#            item.set_fontsize(12)

    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath,'opsep', subdir, figname + '.png'))
    if not showplot:
        plt.close(fig)


def define_colors():
    colors = ['black','red','blue','green','cyan','magenta','violet',\
            'orange','brown','darkred','deepskyblue','mediumseagreen',
            'lightseagreen','purple','sandybrown','cadetblue','goldenrod',
            'navy','palevioletred','saddlebrown']
    return colors


def opsep_plot_all_bins(experiment, flux_type, user_name, options,
    all_dates, all_fluxes, all_energy_bins, event_definitions,
    sep_start_times, sep_end_times, showplot, saveplot, spacecraft='',
    doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False):
    """ Plot all energy bins with all event definitions """

    suffix = "All_Bins"
    figname, subdir = tools.opsep_naming_scheme(all_dates[0], suffix, experiment, flux_type, user_name,
        options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
 
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
    exp_name = experiment
    if experiment == "user":
        exp_name = user_name

    energy_units = event_definitions[0]['energy_channel'].units

    plot_title = f"All Energy Bins with Threshold Crossings\n {exp_name} {title_mod} {flux_type}"
    
    fig = plt.figure(figname,figsize=(12,6))
    ax = plt.subplot(111)
    colors = define_colors()
    #Plot the fluxes
    for j in range(len(all_energy_bins)):
        energy_bin = all_energy_bins[j]

        if energy_bin[1] == -1:
            legend_label = f">{energy_bin[0]} {energy_units}"
        else:
            legend_label = f"{energy_bin[0]}-{energy_bin[1]} {energy_units}"

        is_good = check_for_good_fluxes(all_fluxes[j])
        if not is_good:
            continue

#        maskfluxes = np.ma.masked_invalid(all_fluxes[j])
#        if doBGSubOPSEP or doBGSubIDSEP:
        maskfluxes = np.ma.masked_where(all_fluxes[j] <=0, all_fluxes[j])
        ax.plot_date(all_dates,maskfluxes,'-',label=legend_label, marker='.')
#        else:
#            ax.plot_date(all_dates,maskfluxes,'-',label=legend_label)

    #Plot the threshold crossing times
    for i in range(len(event_definitions)):
        energy_bin = [event_definitions[i]['energy_channel'].min,
                    event_definitions[i]['energy_channel'].max]
        threshold = event_definitions[i]['threshold'].threshold
        flux_units = event_definitions[i]['threshold'].threshold_units

        threshold_label = f"{threshold} ${flux_units}$"
        threshold_label = make_math_label(threshold_label)

        if energy_bin[1] == -1:
            line_label = f">{energy_bin[0]} {energy_units}, {threshold_label}"
        else:
            line_label = f"{energy_bin[0]}-{energy_bin[1]} {energy_units},\n{threshold_label}"

        if not pd.isnull(sep_start_times[i]):
            ax.axvline(sep_start_times[i],color=colors[i],linestyle=':',
                        label=line_label)
            ax.axvline(sep_end_times[i],color=colors[i],linestyle=':')

    plt.title(plot_title)
    ylabel = f"Flux [${flux_units}$]"
    ylabel = make_math_label(ylabel)
    plt.ylabel(ylabel)
    plt.xlabel('Date')
    #ax.xaxis_date()
    #ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d\n%H:%M'))
    plt.yscale("log")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                     chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))
#    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        #item.set_fontsize(12)
    if saveplot:
        fname = os.path.join(cfg.plotpath,'opsep',subdir, figname + '.png')
        fig.savefig(fname)
    if not showplot:
        plt.close(fig)
    
    


def opsep_plot_fluence_spectrum(experiment, flux_type, user_name, options,
    event_definitions, evaluated_dates, energy_bin_centers, fluence_spectra,
    fluence_spectra_units, showplot, saveplot, spacecraft='',
    doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False):
    """ Plot the fluence spectrum calculated by OpSEP. 
        
        fluence_spectra correspond to multiple event_definitions.
        energy_bin_centers are a single array defining the energy bin
            centers for all of the fluence spectra.
    
    """
    
    #Event-integrated fluence for energy channels
    print("Generating figure of event-integrated fluence spectrum.")
    
    suffix = "Fluence"
    figname, subdir = tools.opsep_naming_scheme(evaluated_dates[0], suffix, experiment, flux_type, user_name,
        options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
 
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
    exp_name = experiment
    if experiment == "user":
        exp_name = user_name


    #plot integral fluxes (either input or estimated)
    nthresh = len(event_definitions)
    energy_units = event_definitions[0]['energy_channel'].units
    fluence_units = fluence_spectra_units[0]
 
    ylabel = "Fluence"
 
    ncross = 0 #Check if any thresholds were crossed, if not no need plot

    fig = plt.figure(figname,figsize=(10,8))
    ax = plt.subplot(111)
    markers = ['o','P','D','v','^','<','>','*','d','+','8','p','h','1','X','x']
    colors = define_colors()
    for i in range(nthresh):
    
        if len(fluence_spectra[i]) == 0:
            continue
        ncross = ncross + 1

        energy_bin = [event_definitions[i]['energy_channel'].min, event_definitions[i]['energy_channel'].max]
        flspec_units = fluence_spectra_units[i]
        threshold_label = f"{event_definitions[i]['threshold'].threshold} ${event_definitions[i]['threshold'].threshold_units}$"
        threshold_label = make_math_label(threshold_label)

        #Create labels
        if energy_bin[1] == -1:
            legend_label = f">{energy_bin[0]} {energy_units}, {threshold_label}"
        else:
            legend_label = f"{energy_bin[0]}-{energy_bin[1]} {energy_units}, {threshold_label}"

        ax.plot(energy_bin_centers,fluence_spectra[i],markers[i],
                color=colors[i], mfc='none', label=legend_label)
    
    plt.grid(which="both", axis="both")
    plot_title = f"Event-Integrated Fluence Spectra\n {exp_name} {title_mod} {flux_type}"
    plt.title(plot_title)
    plt.xlabel(f"Energy [{energy_units}]")
    ylabel = f"Fluence [${fluence_units}$]"
    ylabel = make_math_label(ylabel)
    plt.ylabel(ylabel)

    plt.xscale("log")
    plt.yscale("log")
    ax.legend(loc='upper right')
#    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        #item.set_fontsize(14)

    if ncross == 0: plt.close(fig) #no thresholds crossed, empty plot

    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath,'opsep', subdir, figname + '.png'))
    if not showplot:
        plt.close(fig)



############ IDSEP PLOTS ##############
def setup_idsep_plot(figname, experiment, title_mod, unique_id, flux_units):
    """ Set up figure and axes for idsep 3 row plots.
    
    """
    nrow = 3 #number of rows of subplots
    
    plt.rcParams.update({'font.size': 16})
    fig, ax = plt.subplots(nrow, 1, sharex=True, figsize=(13,8), gridspec_kw={'height_ratios' : [1, 1, 1], 'hspace' : 0.4})
    fig.canvas.manager.set_window_title(figname)
    fig.suptitle((f"{experiment} {title_mod} {unique_id}"))
    
    #Formatting of axes
    ax[1].set_ylabel((f"Flux ({flux_units})"))
    ax[2].set_xlabel("Date")
    
    #Apply tight layout and suppress warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=UserWarning)
        plt.tight_layout()

    for iax in range(nrow):
        ax[iax].set_yscale('log')
        ax[iax].grid(axis='both')
    
    fig.autofmt_xdate(rotation=45)

    return fig, ax


def idsep_make_plots(unique_id, experiment, flux_type, exp_name, options, dates,
        fluxes, energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates,
        threshold, doBGSub, showplot, saveplot, disable_sigma=False, spacecraft="",
        close_plot=False):
    """ Make multiple plots with 3 vertical subplots representing individual energy
        channels.
        
        INPUTS:
            
            :unique_id: (string) used in filename and figname
            :experiment: (string)
            :flux_type: (string) integral or differential
            :exp_name: (string) if experiment is "user", then this holds
                the name of the satellite/model
            :options: (array of strings) options applied to data set, used
                to append to filenames and plot title
            :dates: (1xn array of datetime) n time steps
            :fluxes: (pxn array of float) p energy channels for n time steps
            :energy_bins: (px2 array of float) p energy channels
            :ave_dates: (1xm array of datetime) m time steps
            :ave_fluxes: (pxm array of float) average background fluxes for
                p energy channels and m time steps
            :ave_sigma: (px2xm array of float) +- sigma for p energy channels
                and m time steps
            :threshold_dates: (1xn array of datetime) each time step of fluxes
            :threshold: (pxn array of float) ave_fluxes + n*ave_sigma for
                n defined in the fetchsep.config for each time step of fluxes
            :doBGSub: (bool) do background subtraction
            :showplot: (bool)
            :saveplot: (bool)
            :disable_sigma: (bool) True = plot only ave_fluxes withou ave_sigma

        OUPTPUTS:
        
            None, except possibly to write out a plot to file (if showplot, saveplot)
    
    
    """
    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft)
    name = tools.idsep_naming_scheme(experiment, flux_type, exp_name, options, spacecraft=spacecraft)

    figname = (f"{name}_FluxWithThreshold_{unique_id}")
    
    #UNITS
    if flux_type == "integral": flux_units = cfg.flux_units_integral
    if flux_type == "differential": flux_units = cfg.flux_units_differential

    nbins = len(energy_bins)
    for i in range(nbins):
        if not i%3:
            exp = experiment
            if experiment == 'user' and exp_name != '':
                exp = exp_name
            fig, ax = setup_idsep_plot((f"{figname}_{i}"), exp, title_mod, unique_id, flux_units)
            iax = 0

        legend_label = tools.setup_energy_bin_label(energy_bins[i])
            
        #PLOT FLUXES
        maskfluxes = np.ma.masked_less_equal(fluxes[i], 0)
        ax[iax].plot_date(dates,maskfluxes,'.-',label=legend_label,color='tab:blue')
    
        #PLOT BACKGROUND
        if not disable_sigma:
            ax[iax].errorbar(ave_dates, ave_fluxes[i],fmt='.', yerr=ave_sigma[i], label=(f"ave bg {legend_label}"),zorder=100, color='tab:orange')
        if disable_sigma:
            ax[iax].errorbar(ave_dates, ave_fluxes[i],fmt='-',
                label=(f"ave bg {legend_label}"),zorder=100, color='tab:orange')
        
        #PLOT THRESHOLD = n*sigma (n in config file)
        ax[iax].plot_date(threshold_dates,threshold[i],'-',label=(f"threshold {legend_label}"), zorder=200, color='tab:green')

        ax[iax].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=3, mode="expand", borderaxespad=0., fontsize=11,
                      facecolor='lightgray', edgecolor='black')

        if saveplot and (iax ==2 or i == nbins-1):
            fig.savefig(os.path.join(cfg.plotpath,"idsep", name, (f"{figname}_{i}.png")))
            if not showplot:
                plt.close(fig)
            if close_plot:
                plt.close(fig)

        #increment to next axis
        iax += 1
 


def idsep_make_timeseries_plot(unique_id, experiment, flux_type, exp_name,
        options, dates, fluxes, energy_bins, doBGSub, showplot, saveplot,
        spacecraft="", close_plot=False):

    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft)
    name = tools.idsep_naming_scheme(experiment, flux_type, exp_name, options, spacecraft=spacecraft)


    figname = (f"{name}_{unique_id}")
 
    flux_units = ''
    if flux_type == "integral": flux_units = cfg.flux_units_integral
    if flux_type == "differential": flux_units = cfg.flux_units_differential

    nbins = len(energy_bins)
    for i in range(nbins):
        if not i%3:
            exp = experiment
            if experiment == 'user' and exp_name != '':
                exp = exp_name
            fig, ax = setup_idsep_plot((f"{figname}_{i}"), exp, title_mod, unique_id, flux_units)
            iax = 0

        legend_label = tools.setup_energy_bin_label(energy_bins[i])

        maskfluxes = np.ma.masked_less_equal(fluxes[i], 0)
        ax[iax].plot_date(dates,maskfluxes,'.-',label=legend_label)
    
        ax[iax].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=3, mode="expand", borderaxespad=0., fontsize=11,
                      facecolor='lightgray', edgecolor='black')
                      
        if saveplot and (iax ==2 or i == nbins-1):
            fig.savefig(os.path.join(cfg.plotpath,"idsep", name, (f"{figname}_{i}.png")))
            if not showplot:
                plt.close(fig)
            if close_plot:
                plt.close(fig)

        #increment to next axis
        iax += 1
    



def idsep_make_bg_sep_plot(unique_id, experiment, flux_type, exp_name, options,\
            dates, fluxes_bg, fluxes_sep, energy_bins, doBGSub,
            showplot, saveplot, spacecraft="", close_plot=False):
    
    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft)
    name = tools.idsep_naming_scheme(experiment, flux_type, exp_name, options, spacecraft=spacecraft)


    figname = (f"{name}_SEP_BG_{unique_id}")

    flux_units = ''
    if flux_type == "integral": flux_units = cfg.flux_units_integral
    if flux_type == "differential": flux_units = cfg.flux_units_differential

    nbins = len(energy_bins)
    for i in range(nbins):
        if not i%3:
            exp = experiment
            if experiment == 'user' and exp_name != '':
                exp = exp_name
            fig, ax = setup_idsep_plot((f"{figname}_{i}"), exp, title_mod, unique_id, flux_units)
            iax = 0

        legend_label = tools.setup_energy_bin_label(energy_bins[i])

        maskfluxes_bg = np.ma.masked_less_equal(fluxes_bg[i], 0)
        ax[iax].plot_date(dates,maskfluxes_bg,'.-',label=(f"Background {legend_label}"))
        maskfluxes_sep = np.ma.masked_less_equal(fluxes_sep[i], 0)
        ax[iax].plot_date(dates,maskfluxes_sep,'.-',label=(f"SEP {legend_label}"), zorder=100)

        ax[iax].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=3, mode="expand", borderaxespad=0., fontsize=11,
                      facecolor='lightgray', edgecolor='black')
        
        if saveplot and (iax ==2 or i == nbins-1):
            fig.savefig(os.path.join(cfg.plotpath,"idsep", name, (f"{figname}_{i}.png")))
            if not showplot:
                plt.close(fig)
            if close_plot:
                plt.close(fig)
        
        #increment to next axis
        iax += 1
    
    




#def idsep_make_diff_plot(unique_id, experiment, flux_type, exp_name, options, dates,\
#            diff_fluxes, ave_sigma, energy_bins, doBGSub, showplot, saveplot,
#            close_plot=False):
#    #NEEDS TO BE CLEANED UP
#    #Additions to titles and filenames according to user-selected options
#    modifier = ''
#    title_mod = ''
#    if "uncorrected" in options:
#        modifier = modifier + '_uncorrected'
#        title_mod = title_mod + 'uncorrected '
#    if doBGSub:
#        modifier = modifier + '_bgsub'
#        title_mod = title_mod + 'BG-subtracted '
#    if "S14" in options:
#        modifier = modifier + '_S14'
#        title_mod = title_mod + 'S14 '
#    if "Bruno2017" in options:
#        modifier = modifier + '_Bruno2017'
#        title_mod = title_mod + 'Bruno2017 '
#
#
#    figname = experiment + '_' + flux_type + modifier \
#            + '_' + 'Diff_' + unique_id
#    if experiment == 'user' and exp_name != '':
#        figname = exp_name + '_' + flux_type + modifier \
#                + '_' + 'Diff_' + unique_id
#    
#    fig = plt.figure(figname,figsize=(12,8))
#    plt.rcParams.update({'font.size': 16})
#    ax = plt.subplot(111)
#    nbins = len(energy_bins)
#    ifig = 0
#    for i in range(nbins):
#        thresh = np.multiply(ave_sigma[i][1],nsigma)
#        if i != 0 and not i%3:
#            if saveplot:
#                fig.savefig(cfg.plotpath + '/idsep/' +figname + '.png')
#            figname = figname + str(i)
#            fig = plt.figure(figname,figsize=(12,8))
#            ax = plt.subplot(111)
#            ifig = 0
#
#        ax = plt.subplot(min(3,nbins), 1, ifig+1)
#        ifig = ifig + 1
#        legend_label = ""
#        if energy_bins[i][1] != -1:
#            legend_label = str(energy_bins[i][0]) + '-' \
#                    + str(energy_bins[i][1]) + ' ' + cfg.energy_units
#        else:
#            legend_label = '>'+ str(energy_bins[i][0]) + ' ' + cfg.energy_units
#
#        ax.plot_date(dates,diff_fluxes[i],'.',label="diff " + legend_label)
#        ax.plot_date(dates,thresh,'-',label="threshold " + legend_label, zorder=100)
#        
#        flux_units = ''
#        if flux_type == "integral": flux_units = cfg.flux_units_integral
#        if flux_type == "differential": flux_units = cfg.flux_units_differential
#        
#        if i==0:
#            plt.title(experiment + ' '+ title_mod + ' ' + unique_id\
#                    + "\nDiff = Flux - Mean BG")
#            if experiment == 'user' and exp_name != '':
#                plt.title(exp_name + ' '+ title_mod + ' ' + unique_id\
#                    + "\nDiff = Flux - Mean BG")
#        plt.xlabel('Date')
#        #plt.ylabel('Flux [' + flux_units + ']')
#        plt.ylabel(r'Flux (MeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
#        fig.autofmt_xdate(rotation=45)
#        chartBox = ax.get_position()
##        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
##                         chartBox.height])
##        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05),fontsize=11)
# 
# 
#        if saveplot and i == nbins-1:
#            fig.savefig(cfg.plotpath + '/idsep/' +figname + '.png')
#            if not showplot:
#                plt.close(fig)
#            if close_plot:
#                plt.close(fig)
