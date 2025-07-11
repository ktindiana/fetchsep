from . import config as cfg
import numpy as np
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


def plot_marginals(y_true, y_pred, scale="linear",
        x_label="Observations", y_label="Forecast", thresh=None,
        save="marginal_plot", showplot=False, closeplot=False):
    """
    Plots model forecast against observations
    in a scatter plot with marginals. Includes an
    option for displaying thresholds use when discretizing
    into a categorical forecast

    Parameters
    ----------
    y_true : array-like
        Observed (true) values

    y_pred : array-like
        Forecasted (estimated) values

    scale : string
        Numeric scale to display results on
        Accepts "linear" or "log"
        Optional. Defaults to "linear"

    x_label : string
        Label for x-axis
        Optional. Defaults to "Observations"

    y_label : string
        Label for y-axis
        Optional. Defaults to "Forecast"

    thresh : float
        Value of threshold when displaying discretization
        Option. Defaults to None

    save : string
        Name to save PNG as (should not include ".png")
        Optional. Defaults to "marginal_plot"

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

    check_consistent_length(y_true, y_pred)

    y_true = check_array(y_true, force_all_finite=True, ensure_2d=False)
    y_pred = check_array(y_pred, force_all_finite=True, ensure_2d=False)

    if (scale == "log") and ((y_true < 0).any() or (y_pred < 0).any()):
        raise ValueError("Values cannot be negative on a logarithmic scale")

    if scale not in ("linear", "log"):
        raise ValueError("Scale must either be 'linear' or 'log'")

    #plt.style.use('dark_background')

    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(3, 3)
    ax_main = plt.subplot(gs[1:3, :2])
    ax_xDist = plt.subplot(gs[0, :2], sharex=ax_main)
    ax_yDist = plt.subplot(gs[1:3, 2], sharey=ax_main)

    color_data = '#0339f8'
    color_thresh = '#02c14d'
    color_unity_line = '#fc2647'
    color_mean = '#aa23ff'

    ax_main.scatter(y_true, y_pred, marker='.', color=color_data)
    ax_main.set(xlabel=x_label, ylabel=y_label)

    # getting max and min of x and y
    if scale == "log":
        xmax = 10**(math.ceil(np.log10(np.max(y_true))))
        ymax = 10**(math.ceil(np.log10(np.max(y_pred))))
        xymax = max(xmax, ymax)
        xmin = 10**(math.floor(np.log10(np.min(y_true))))
        ymin = 10**(math.floor(np.log10(np.min(y_pred))))
        xymin = min(xmin, ymin)
    elif scale == "linear":
        xmax = np.max(y_true)
        ymax = np.max(y_pred)
        xymax = max(xmax, ymax)
        xmin = np.min(y_true)
        ymin = np.min(y_pred)
        xymin = min(xmin, ymin)

    ax_main.set_xlim(xymin, xymax)
    ax_main.set_ylim(xymin, xymax)

    ax_main.grid(True, linestyle=':', alpha=0.5)

    ax_main.set_xscale(scale)
    ax_main.set_yscale(scale)

    # creating bins for histograms
    if scale == "log":
        hbins = np.logspace(np.log10(xymin), np.log10(xymax), 100)
    elif scale == "linear":
        hbins = np.linspace(xymin, xymax, 100)

    # histogram for x-axis
    ax_xDist.hist(y_true, bins=hbins, align='mid', color=color_data, zorder=0)
    ax_xDist.set(ylabel='Counts')

    # histogram for y-axis
    ax_yDist.hist(y_pred, bins=hbins, orientation='horizontal', align='mid', \
                  color=color_data, zorder=0)
    ax_yDist.set(xlabel='Counts')

    # drawing the unity line
    ax_main.plot([xymin, xymax], [xymin, xymax], linestyle='--', linewidth=1.5, \
                 color=color_unity_line, zorder=2)

    # plotting the mean value and line from the unity line to the mean value
    meanx = np.mean(y_true)
    meany = np.mean(y_pred)
    #ax_main.plot(meanx, meany, marker='o', color=color_mean, zorder=1)
    #ax_main.plot([meanx, meany], [meany, meany], linestyle='-', linewidth=2.5, \
    #             color=color_mean, zorder=1)

    if thresh != None:
        # drawing threshold lines for contingency table info
        ax_main.plot([thresh, thresh], [xymin, xymax], [xymin, xymax], \
                     [thresh, thresh], linestyle='-', linewidth=1.5, \
                     color=color_thresh, zorder=2)
        # printing contingency table labels
        plt.text(0.98, 0.98, 'Hits', fontsize=12, color=color_thresh, \
                 horizontalalignment='right', verticalalignment='top', \
                 transform=ax_main.transAxes)
        plt.text(0.02, 0.98, 'False Alarms', fontsize=12, color=color_thresh, \
                 horizontalalignment='left', verticalalignment='top', \
                 transform=ax_main.transAxes)
        plt.text(0.02, 0.02, 'Correct Negatives', fontsize=12, \
                 color=color_thresh, horizontalalignment='left', \
                 verticalalignment='bottom', transform=ax_main.transAxes)
        plt.text(0.98, 0.02, 'Misses', fontsize=12, color=color_thresh, horizontalalignment='right', verticalalignment='bottom', \
                 transform=ax_main.transAxes)

    if showplot: plt.show()

    fig.savefig('plots/MarginalPlots/'+save+'.png', dpi=300, bbox_inches='tight')

    if closeplot: plt.close(fig)

    return



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



def setup_modifiers(options, doBGSub, spacecraft=""):
    """ Add modifier strings according to options.
    
    """
    modifier = '' #for appending to filenames
    title_mod = '' #for appending to plot titles

    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '
    if spacecraft:
        modifier = modifier + '_' + spacecraft
        title_mod = title_mod + spacecraft + ' '

    return modifier, title_mod


def setup_energy_bin_label(energy_bin):
    """ Label for a single energy bin.
    
    """
    label = ""
    if energy_bin[1] != -1:
        label = (f"{energy_bin[0]}-{energy_bin[1]} {cfg.energy_units}")
    else:
        label = (f">{energy_bin[0]} {cfg.energy_units}")

    return label


def idsep_make_plots(unique_id, experiment, flux_type, exp_name, options, dates,
        fluxes, energy_bins, ave_dates, ave_fluxes, ave_sigma, threshold_dates,
        threshold, doBGSub, showplot, saveplot, disable_sigma=False, spacecraft=""):
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
    modifier, title_mod = setup_modifiers(options, doBGSub, spacecraft=spacecraft)


    figname = (f"{experiment}_{flux_type}{modifier}_FluxWithThreshold_{unique_id}")
    if experiment == 'user' and exp_name != '':
        figname = (f"{exp_name}_{flux_type}{modifier}_FluxWithThreshold_{unique_id}")
    
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

        legend_label = setup_energy_bin_label(energy_bins[i])
            
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
            fig.savefig(os.path.join(cfg.plotpath,"idsep",(f"{figname}_{i}.png")))
            if not showplot:
                plt.close(fig)

        #increment to next axis
        iax += 1
 



def idsep_make_timeseries_plot(unique_id, experiment, flux_type, exp_name,
        options, dates, fluxes, energy_bins, doBGSub, showplot, saveplot,
        spacecraft=""):

    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = setup_modifiers(options, doBGSub, spacecraft=spacecraft)


    figname = (f"{experiment}_{flux_type}{modifier}_FluxTimeseries_{unique_id}")
    if experiment == 'user' and exp_name != '':
        figname = (f"{exp_name}_{flux_type}{modifier}_FluxTimeseries_{unique_id}")
 
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

        legend_label = setup_energy_bin_label(energy_bins[i])

        maskfluxes = np.ma.masked_less_equal(fluxes[i], 0)
        ax[iax].plot_date(dates,maskfluxes,'.-',label=legend_label)
    
        ax[iax].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=3, mode="expand", borderaxespad=0., fontsize=11,
                      facecolor='lightgray', edgecolor='black')
                      
        if saveplot and (iax ==2 or i == nbins-1):
            fig.savefig(os.path.join(cfg.plotpath,"idsep",(f"{figname}_{i}.png")))
            if not showplot:
                plt.close(fig)

        #increment to next axis
        iax += 1
    



def idsep_make_bg_sep_plot(unique_id, experiment, flux_type, exp_name, options,\
            dates, fluxes_bg, fluxes_sep, energy_bins, doBGSub,
            showplot, saveplot, spacecraft=""):
    
    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = setup_modifiers(options, doBGSub, spacecraft=spacecraft)

    figname = (f"{experiment}_{flux_type}{modifier}_SEP_BG_{unique_id}")
    if experiment == 'user' and exp_name != '':
        figname = (f"{exp_name}_{flux_type}{modifier}_SEP_BG_{unique_id}")

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

        legend_label = setup_energy_bin_label(energy_bins[i])

        maskfluxes_bg = np.ma.masked_less_equal(fluxes_bg[i], 0)
        ax[iax].plot_date(dates,maskfluxes_bg,'.-',label=(f"Background {legend_label}"))
        maskfluxes_sep = np.ma.masked_less_equal(fluxes_sep[i], 0)
        ax[iax].plot_date(dates,maskfluxes_sep,'.-',label=(f"SEP {legend_label}"), zorder=100)

        ax[iax].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=3, mode="expand", borderaxespad=0., fontsize=11,
                      facecolor='lightgray', edgecolor='black')
        
        if saveplot and (iax ==2 or i == nbins-1):
            fig.savefig(os.path.join(cfg.plotpath,"idsep",(f"{figname}_{i}.png")))
            if not showplot:
                plt.close(fig)

        #increment to next axis
        iax += 1
    
    




def idsep_make_diff_plot(unique_id, experiment, flux_type, exp_name, options, dates,\
            diff_fluxes, ave_sigma, energy_bins, doBGSub, showplot, saveplot):
    #NEEDS TO BE CLEANED UP
    #Additions to titles and filenames according to user-selected options
    modifier = ''
    title_mod = ''
    if "uncorrected" in options:
        modifier = modifier + '_uncorrected'
        title_mod = title_mod + 'uncorrected '
    if doBGSub:
        modifier = modifier + '_bgsub'
        title_mod = title_mod + 'BG-subtracted '
    if "S14" in options:
        modifier = modifier + '_S14'
        title_mod = title_mod + 'S14 '
    if "Bruno2017" in options:
        modifier = modifier + '_Bruno2017'
        title_mod = title_mod + 'Bruno2017 '


    figname = experiment + '_' + flux_type + modifier \
            + '_' + 'Diff_' + unique_id
    if experiment == 'user' and exp_name != '':
        figname = exp_name + '_' + flux_type + modifier \
                + '_' + 'Diff_' + unique_id
    
    fig = plt.figure(figname,figsize=(12,8))
    plt.rcParams.update({'font.size': 16})
    ax = plt.subplot(111)
    nbins = len(energy_bins)
    ifig = 0
    for i in range(nbins):
        thresh = np.multiply(ave_sigma[i][1],nsigma)
        if i != 0 and not i%3:
            if saveplot:
                fig.savefig(cfg.plotpath + '/idsep/' +figname + '.png')
            figname = figname + str(i)
            fig = plt.figure(figname,figsize=(12,8))
            ax = plt.subplot(111)
            ifig = 0

        ax = plt.subplot(min(3,nbins), 1, ifig+1)
        ifig = ifig + 1
        legend_label = ""
        if energy_bins[i][1] != -1:
            legend_label = str(energy_bins[i][0]) + '-' \
                    + str(energy_bins[i][1]) + ' ' + cfg.energy_units
        else:
            legend_label = '>'+ str(energy_bins[i][0]) + ' ' + cfg.energy_units

        ax.plot_date(dates,diff_fluxes[i],'.',label="diff " + legend_label)
        ax.plot_date(dates,thresh,'-',label="threshold " + legend_label, zorder=100)
        
        flux_units = ''
        if flux_type == "integral": flux_units = cfg.flux_units_integral
        if flux_type == "differential": flux_units = cfg.flux_units_differential
        
        if i==0:
            plt.title(experiment + ' '+ title_mod + ' ' + unique_id\
                    + "\nDiff = Flux - Mean BG")
            if experiment == 'user' and exp_name != '':
                plt.title(exp_name + ' '+ title_mod + ' ' + unique_id\
                    + "\nDiff = Flux - Mean BG")
        plt.xlabel('Date')
        #plt.ylabel('Flux [' + flux_units + ']')
        plt.ylabel(r'Flux (MeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
        fig.autofmt_xdate(rotation=45)
        chartBox = ax.get_position()
#        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
#                         chartBox.height])
#        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05),fontsize=11)
 
 
        if saveplot and i == nbins-1:
            fig.savefig(cfg.plotpath + '/idsep/' +figname + '.png')
            if not showplot:
                plt.close(fig)

