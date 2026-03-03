from . import config as cfg
from . import tools
from . import date_handler as dh
from . import experiments as expts
from ..json import ccmc_json_handler as ccmc_json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.dates import DateFormatter
import math
from itertools import cycle
from sklearn.utils.validation import check_consistent_length
from sklearn.utils.validation import check_array
import datetime
from scipy.stats import pearsonr
from math import log10
from pandas.plotting import register_matplotlib_converters
import warnings
import os

__author__ = "Kathryn Whitman"
__maintainer__ = "Kathryn Whitman"
__email__ = "kathryn.whitman@nasa.gov"

def plot_time_profile(date, values, labels, dy=None, dyl=None,
                dyh=None, title=None, x_min=None, x_max=None,
                x_label="Date", y_min=None, y_max=None,
                y_label="Value", uselog_x = False, uselog_y = False,
                date_format="year", save="time_profile",
                showplot=False, closeplot=False, saveplot=False):
    """
    Plots multiple time profiles in same plot.
    Includes specialized formatting for specific models or datasets, 
    but can be called generically.

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
    OPSEPEnhancement=False, IDSEPEnhancement=False,
    color_scheme=1, no_goes_colors=False):
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
    flux_units = tools.get_flux_units(flux_type)
    flux_units = make_math_label(flux_units)
    energy_units = tools.get_energy_units()
    
    colors = define_colors(energy_bins, color_scheme=color_scheme, no_goes_colors=no_goes_colors)
    
    #All energy channels in specified date range with event start and stop
    #Plot all channels of user specified data
    #Additions to titles and filenames according to user-selected options
    suffix = f"All_Bins_{unique_id}"
    figname, subdir = tools.opsep_naming_scheme(dates[0], suffix, experiment, flux_type,
        user_name, options, spacecraft=spacecraft,
        doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft,
        doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    exp_name = experiment

    fig = plt.figure(figname,figsize=(13.5,8))
    ax = plt.subplot(111)
    nbins = len(energy_bins)
    for i in range(nbins):
        #Check to make sure there are positive fluxes
        is_good = check_for_good_fluxes(fluxes[i])
        if not is_good:
            continue
        #Don't want to plot zero values, particularly in background-subtracted plots
        maskfluxes = np.ma.masked_where(fluxes[i] <= 0, fluxes[i])

        mean =  means[i]
        if isinstance(means[i], list):
            mean = np.nanmean(np.array(means[i]))
            
        legend_label = tools.setup_energy_bin_label(energy_bins[i])
        p = ax.plot_date(dates,maskfluxes, '-', label=legend_label, color=colors[legend_label])
        color = p[0].get_color()
        if i==0:
            plt.axhline(mean,color=color,linestyle=':', label="Mean Background")
        else:
            plt.axhline(mean,color=color,linestyle=':')


    fig.suptitle((f"{unique_id}\n {exp_name} {title_mod} {flux_type}"))
    
    #Formatting of axes
    ax.set_ylabel((f"Intensity ({flux_units})"))
    ax.set_xlabel("Date")
    plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d\n%H:%M"))
    plt.xticks(rotation=45, ha="right")

    plt.yscale("log")
    plt.ylim(5e-5,1e5)
    plt.xlim(dates[0], dates[-1])
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                     chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(10)

    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath, 'opsep', subdir, figname + '.png'))



def plot_weibull_fit(energy_bin, threshold, experiment, flux_type, user_name,
    options, sep_start_time, trim_times, trim_fluxes,
    best_pars, best_fit, max_time, max_val, max_meas_time, max_meas,
    max_curve_model_time, max_curve_model_peak, max_curve_meas_time, max_curve_meas_peak,
    saveplot, showplot, spacecraft='', doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False):
    """ Plot Weibull fit used to get onset peak """

    flux_units = tools.get_flux_units_bin(energy_bin)
    flux_units = make_math_label(flux_units)
    energy_units = tools.get_energy_units()

    exp_name = experiment
    if experiment == "user":
        exp_name = user_name

    best_a = best_pars['alpha']
    best_b = best_pars['beta']
    best_Ip = best_pars['peak_intensity']

    suffix = f"Weibull_profile_fit_{energy_bin[0]}{energy_units}"
    figname, subdir = tools.opsep_naming_scheme(sep_start_time, suffix, experiment, flux_type, user_name,
        options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP,
        OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)

    fig = plt.figure(figname,figsize=(9,5))
    label = f"{energy_bin[0]} - {energy_bin[1]} {energy_units}"
    if energy_bin[1] == -1:
        label = f">{energy_bin[0]} {energy_units}"
    plt.plot(trim_times,trim_fluxes,label=label,marker='.', linestyle='none')
    label_fit = (f"Fit\n Ip: {best_Ip:.2f}"
                f"\n alpha: {best_a:.2f}"
                f"\n beta: {best_b:.2f}")
    plt.plot(trim_times,best_fit,label=label_fit)
    plt.plot(max_time, max_val,"o",label="max fit")
    plt.plot(max_meas_time, max_meas,">",label="measured peak near\nfit max")
    plt.plot(max_curve_model_time, max_curve_model_peak,"D",label="Min 2nd Derivative")
    plt.plot(max_curve_meas_time, max_curve_meas_peak,"^",label="measured peak near\nmin 2nd Derivative")
    #plt.plot(onset_time, onset_peak[i],">",label="onset Weibull")
    plt.legend(loc='lower right')
    plt.title(f"Onset Peak Weibull Fit\n {exp_name} {title_mod} {flux_type}\n {sep_start_time}")
    plt.xlabel("Hours")
    plt.ylabel(flux_units)
    plt.yscale("log")
    plt.ylim(1e-4,5*max_val)
    
    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath, 'opsep', subdir, figname + '.png'))
    if not showplot:
        plt.close(fig)


####PLOTTING FOR find_max_curvature in tools.py
#    if showplot or saveplot:
#        tzulu = dh.time_to_zulu(crossing_time)
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
#            fig.savefig(cfg.plotpath + '/opsep/' + figname + '.png')
#        if not showplot:
#            plt.close(fig)


def make_math_label(label):
    label = label.replace("^-1", "$^\\mathregular{-1}$")
    label = label.replace("^-2", "$^\\mathregular{-2}$")
    label = label.replace("*", " ")
    return label


def define_colors(energy_bins, event_definitions=None, color_scheme=1, no_goes_colors=False):
    """ Colors for flux time series. Colors are mapped to energy bins
        and event definitions (energy bin + threshold).
        
        Users can choose different color schemes according to their
        preference.
        
        The standard integral energy channels are plotted with 
        SWPC/SEP Scoreboard colors. Other energy channels are plotted
        with colors that are distinct from the integral channel colors.
        
        INPUTS:
        
            :energy_bins: (list) e.g. [[10,-1],[30,-1],[50,-1]]
            :event_definitions: (object) contains energy channel and threshold
            :color_scheme: (int) allows user to select different color schemes
            :no_goes_colors: (bool) True will turn off using the integral_colors
                for >5, >10, >30, etc MeV
        
    
    """
    #########
    #Integral channel GOES colors to match NOAA SWPC and SEP Scoreboards
    integral_colors = {
                '>1 MeV'     : '#bdbdbd', #'#b3b3b3',
                '>5 MeV'     : '#ffd480', # '#f0ad4e'
                '>10 MeV'    : '#ff0000',
                '>30 MeV'    : '#6b3d9a', #'#9467bd',
                '>50 MeV'    : '#0000ff',
                '>60 MeV'    : '#000000', #'#6e6e6e'
                '>100 MeV'   : '#00ff00',
                '>500 MeV'   : 'darkturquoise', #SEP Scoreboard choice '#00cfd1',
                '>1.0 MeV'     : '#bdbdbd',
                '>5.0 MeV'     : '#ffd480', #SEP Scoreboard choice,
                '>10.0 MeV'    : '#ff0000',
                '>30.0 MeV'    : '#6b3d9a',
                '>50.0 MeV'    : '#0000ff',
                '>60.0 MeV'    : '#000000',
                '>100.0 MeV'   : '#00ff00',
                '>500.0 MeV'   : 'darkturquoise', # SWPC
                }

    #########
    #Colors for other energy bins
 
    #### DEFAULT ####
    ######Tab 10-like bright, exciting, cosmic DEFAULT
    if color_scheme == 1:
        #Fluxes
        flux_colors = [
            "#1A3D8C",  # vivid cosmic blue
            "#FF6D2E",  # bright molten orange
            "#00A28A",  # strong teal / alien green
            "#E63946",  # vivid cosmic red
            "#8C4BDC",  # adjusted nebula purple (distinct from #6A4DC9)
            "#F0C27B",  # glowing sand / orange
            "#4DA6FF",  # bright sky blue
            "#FF6B91",  # adjusted coral: less orange, more pink/magenta
            "#3FA76D",  # forest green
            "#2EC4E6",  # electric plasma blue (replaces bright violet)
            "#B078F0",  # bright nebula violet
        ]

        threshold_colors = [
            #Go with Set 1-like bright, exciting, cosmic DEFAULT
            "#5B4F6D",  # muted cool purple (slate-like)
            "#B88628",  # warm mustard / starlight gold
            "#853841",  # muted wine red
            "#0E715B",  # deep teal / cosmic ocean
            "#BFB288",  # soft khaki / nebula dust
            "#6B8C2A",  # muted olive-green / planetary surface
        ]


    ######Tab 20-like exciting, dark, and soft space exploration palette
    if color_scheme == 2:
        flux_colors = [
            # Blue pair (space & deep sky)
            "#1D3B6A",  # deep vibrant cosmic blue
            "#5B8BC9",  # brightened light sky blue

            # Orange pair (Mars / rover dust)
            "#C26834",  # vivid burnt Mars orange
            "#E3A875",  # lively desert sand orange

            # Green / Teal pair (Earth, instruments)
            "#2E8E82",  # energetic dark teal
            "#6CC6B7",  # brighter pale teal

            # Red / Coral pair (Mars, spacecraft accents)
            "#A53433",  # bold brick red
            "#E07D79",  # warm dusty coral

            # Purple / Indigo pair (space nebula / shadows)
            "#4B49A8",  # more vibrant dark indigo
            "#8F78D4",  # brighter muted violet

            # Brown / Sand pair (lunar / asteroid surfaces)
            "#876048",  # rich dark brown
            "#D1B081",  # warm sand

            # Extra accent greens (Earth / habitat elements)
            "#538160",  # lively forest green
            "#9ACF9D",  # fresh soft green

            # Extra accent oranges (engine / solar panels)
            "#B5743B",  # brighter metallic burnt orange
            "#E8AC5A",  # muted gold-orange

            # Extra blues (space instruments)
            "#1B4463",  # deep navy
            "#5797D5",  # soft sky blue
                       ]

        threshold_colors = [
            "#4DA6FF",  # bright sky blue
            "#FF6D2E",  # vivid molten orange
            "#70E0C8",  # bright aqua / nebula gas
            "#FF7F91",  # bright coral / stellar highlights
            "#B44FD9",  # bright magenta-violet / nebula accent
            "#FFD460",  # glowing starlight gold

        ]


    ######Set 1-like color scheme inspired by Mars
    if color_scheme == 3:
        #Set 1-like Mars
        flux_colors = [
            "#1D3B6A",  # deep vibrant cosmic blue
            "#C26834",  # vivid burnt Mars orange
            "#2E8E82",  # energetic dark teal
            "#A53433",  # bold brick red
            "#4B49A8",  # more vibrant dark indigo
            "#876048",  # rich dark brown
            "#5B8BC9",  # brightened light sky blue
            "#E3A875",  # lively desert sand orange
            "#6CC6B7",  # brighter pale teal
            "#E07D79",  # warm dusty coral
                       ]

        threshold_colors = [
            "#4DA6FF",  # bright sky blue
            "#FF6D2E",  # vivid molten orange
            "#70E0C8",  # bright aqua / nebula gas
            "#FF7F91",  # bright coral / stellar highlights
            "#B44FD9",  # bright magenta-violet / nebula accent
            "#FFD460",  # glowing starlight gold
        ]


    ####Set 1-like inspired by nebulae and supernovae
    if color_scheme == 4:
        #Set 1-like nebulae
        flux_colors = [
            "#2C2D72",  # deep cosmic indigo (space background)
            "#4B6CB7",  # muted sky blue (stellar gas)
            "#D94F1A",  # fiery orange (supernova shock)
            "#FF9E3D",  # glowing amber (nebula dust)
            "#E03E5C",  # vivid cosmic red (hot star regions)
            "#B75CC2",  # nebula pink/purple (ionized gases)
            "#1F8F7A",  # teal (aurora-like gas clouds)
            "#70C9C4",  # soft aqua (diffuse gas)
            "#8B5C42",  # asteroid / rocky brown (planetary surfaces)
            "#E0B56B",  # soft gold (starlight glow)
        ]

        threshold_colors = [
            #With the Orion nebula scheme
            "#3B9CFF",  # bright sky blue
            "#FF5733",  # vivid orange-red
            "#1ABC9C",  # bright teal
            "#FF6F91",  # pink / coral
            "#9B59B6",  # strong magenta-violet
            "#FFD700",  # gold / starlight yellow
        ]


    ####Set 1-like, inspired by Ring Nebula and Orion Nebula
    if color_scheme == 5:
        flux_colors = [
            "#0B285D",  # deep cosmic blue
            "#355BA1",  # darker nebula blue
            "#1C7AA0",  # rich Orion blue (avoids darkturquoise)
            "#D9451A",  # fiery nebula orange
            "#E68C2F",  # muted glowing amber
            "#9B3FBF",  # vibrant magenta-violet
            "#2A8C80",  # dark teal (distinct from darkturquoise)
            "#D05F5B",  # dusty coral / stellar highlights
            "#3B7D9A",  # medium darker blue
            "#D1AA53",  # starlight gold / glowing edges
                       ]

        threshold_colors = [
            #With the Orion nebula scheme
            "#3B9CFF",  # bright sky blue
            "#FF5733",  # vivid orange-red
            "#1ABC9C",  # bright teal
            "#FF6F91",  # pink / coral
            "#9B59B6",  # strong magenta-violet
            "#FFD700",  # gold / starlight yellow

        ]


    ####Tab 20-like very bright color scheme
    if color_scheme == 6:
        #BRIGHT
        flux_colors = [
            # Blue pair
            "#3366CC",  # brighter dark blue
            "#66B2FF",  # brighter light blue

            # Orange pair
            "#FF7F33",  # brighter dark orange
            "#FFCC66",  # brighter light orange

            # Green / Teal pair
            "#33B2AA",  # brighter dark teal
            "#80E0D9",  # brighter light teal

            # Red / Coral pair
            "#E64C3C",  # brighter dark red-orange
            "#FF8A88",  # brighter light coral

            # Purple / Indigo pair
            "#6666FF",  # brighter dark indigo
            "#B366FF",  # brighter light violet

            # Brown / Sand pair
            "#C78C51",  # brighter dark brown
            "#E8C18D",  # brighter light sand

            # Extra accent greens
            "#73A46C",  # brighter dark forest green
            "#A8D18A",  # brighter soft green

            # Extra accent oranges
            "#D68B46",  # brighter brown-orange
            "#FFA333",  # brighter warm orange

            # Extra blues
            "#264C99",  # brighter deep navy
            "#66C2FF",  # brighter light azure
        ]

        threshold_colors = [
            "#6C757D",  # cool lunar gray
            "#C2185B",  # deep supernova magenta
            "#7B1FA2",  # dark cosmic plum
            "#A3B18A",  # sage / muted planetary green
            "#8D6E63",  # asteroid taupe
            "#B8860B",  # dark metallic gold
        ]



    #Standard Python Tab10 with red replaced with muted coral
    #Distinct scheme used for additional thresholds
    if color_scheme == 7:
        flux_colors = [
            "#1f77b4",  # tab:blue
            "#ff7f0e",  # tab:orange
            "#2ca02c",  # tab:green
            "#E65A60",  # reddish coral
            "#9467bd",  # tab:purple
            "#8c564b",  # tab:brown
            "#e377c2",  # tab:pink
            "#7f7f7f",  # tab:gray
            "#bcbd22",  # tab:olive
            "#17becf",  # tab:cyan
        ]

        
        threshold_colors = [
            #Set 2
            "#66c2a5",  # teal
            "#fc8d62",  # orange
            "#8da0cb",  # blue-violet
            "#e78ac3",  # pink
            "#a6d854",  # green
            "#ffd92f",  # yellow
            "#e5c494",  # tan
            "#b3b3b3",  # gray
        ]


    #Standard Python Set1 with red replaced with muted coral
    #Distinct scheme used for additional thresholds
    if color_scheme == 8:
        flux_colors = [
            "#E65A60",  # reddish coral
            "#377eb8",  # blue
            "#4daf4a",  # green
            "#984ea3",  # purple
            "#ff7f00",  # orange
            "#ffff33",  # yellow
            "#a65628",  # brown
            "#f781bf",  # pink
            "#999999",  # gray
        ]

        
        threshold_colors = [
            "#1F5F8B",  # cool – slate-teal blue
            "#D4A017",  # warm – mustard
            "#8E3B46",  # warm – wine
            "#0B6E4F",  # cool – deep sea green
            "#C2B280",  # warm – khaki
            "#4A4E69",  # cool – cool slate
            "#B8E620",  # warm – chartreuse
            "#2D6A4F",  # cool – forest green
            "#6B8E23",  # warm – olive-drab
        ]


    #Standard Python Tab20 with red replaced with coral
    #Set2 scheme used for additional thresholds
    if color_scheme == 9:
        flux_colors = [
            "#1f77b4",  # blue
            "#aec7e8",  # light blue
            "#ff7f0e",  # orange
            "#ffbb78",  # light orange
            "#2ca02c",  # green
            "#98df8a",  # light green
            "#E65A60",  # reddish coral
            "#ff9896",  # light red
            "#9467bd",  # purple
            "#c5b0d5",  # light purple
            "#8c564b",  # brown
            "#c49c94",  # light brown
            "#e377c2",  # pink
            "#f7b6d2",  # light pink
            "#7f7f7f",  # gray
            "#c7c7c7",  # light gray
            "#bcbd22",  # olive
            "#dbdb8d",  # light olive
            "#17becf",  # cyan
            "#9edae5",  # light cyan
        ]

        
        threshold_colors = [
            "#1F5F8B",  # cool – slate-teal blue
            "#D4A017",  # warm – mustard
            "#8E3B46",  # warm – wine
            "#0B6E4F",  # cool – deep sea green
            "#C2B280",  # warm – khaki
            "#4A4E69",  # cool – cool slate
            "#B8E620",  # warm – chartreuse
            "#2D6A4F",  # cool – forest green
            "#6B8E23",  # warm – olive-drab
        ]


    #Tab 20-like muted blues, olives, purples colors to browns, all distinct from GOES colors
    if color_scheme == 10:
        flux_colors = [
            # --- Distinct hues (cool + accent colors first) ---
            "#1B3A4B",  # deep slate blue
            "#3E7CB1",  # muted sky blue
            #"#5AA9E6",  # light azure
            "#2A9D8F",  # sea green
            "#588157",  # muted green
            "#A7C957",  # yellow-green
            "#28536B",  # steel teal
            #"#52796F",  # cool sage
            "#3A5A40",  # deep forest
            "#9F86C0",  # soft violet
            "#CDB4DB",  # pale lavender
            "#6D597A",  # smoky plum
            "#E9C46A",  # muted gold
            "#E56B6F",  # muted coral

            # --- Browns, tans, and grays at the end ---
            "#BC6C25",  # brown-orange
            "#9C6644",  # sienna
            "#7F5539",  # cocoa
            "#B08968",  # warm tan
            #"#D4A373",  # sand
            #"#8D99AE",  # blue-gray
        ]

        threshold_colors = [
            "#5AA9E6",  # light azure (cool)
            "#7FB069",  # soft leaf (cool)
            "#EAAC8B",  # soft peach (warm)
            "#4C5B5C",  # blue slate (cool-neutral)
            "#344E41",  # dark moss (cool)
            "#ADB5BD",  # light cool gray (neutral)
        ]


    #Tab 20-like autumn colors; Reds to blues and are not the same as GOES colors
    if color_scheme == 11:
        flux_colors = [
            "#C44536",  # strong red-orange
            "#E36414",  # burnt orange
            "#F77F00",  # warm orange
            "#FF9F1C",  # light orange
            
            "#BC6C25",  # brown-orange
            "#D4A373",  # sand
            "#E9C46A",  # muted gold
            "#F4D35E",  # soft yellow
            
            "#588157",  # forest
            "#6A994E",  # medium green
            "#7FB069",  # soft green
            "#A7C957",  # yellow-green
            
            "#2A9D8F",  # sea green
            "#1B9E77",  # teal
            "#20B2AA",  # light teal
            "#118AB2",  # blue-teal
            
            "#2C5282",  # dark blue
            "#3A86FF",  # bright blue (not #0000ff)
            "#5E60CE",  # indigo
            "#9D4EDD",  # light violet
        ]

        threshold_colors = [
            "#9E2A2B",  # deep brick
            "#A47148",  # warm brown
            "#386641",  # deep green
            "#1B4332",  # deep teal-green
            "#1D3557",  # deep navy (not pure blue)
            "#7209B7",  # violet
        ]


    #Tab 20-like bold, dark autumn colors with strong visual separation
    if color_scheme == 12:
        flux_colors = [
            "#8E2C2C",  # deep brick (not red)
            "#A23E48",  # wine
            "#B5651D",  # strong amber-brown
            "#C67C00",  # dark goldenrod (not orange)
            "#D4A017",  # bold gold
            "#6B8E23",  # olive drab
            "#4F772D",  # deep olive green
            "#3A5F0B",  # dark moss
            "#006400",  # dark green (not bright green)
            "#2F4F4F",  # dark slate teal (not cyan)
            "#0B3C49",  # deep blue-slate (not blue)
            "#1B4965",  # strong steel blue
            "#2E4057",  # dark desaturated blue
            "#5C415D",  # smoky plum-brown (not purple)
            "#7A306C",  # deep magenta-brown (not violet)
            "#9A031E",  # dark crimson-brown (not pure red)
            "#BB9457",  # bold tan
            "#99582A",  # saddle brown
            "#6C757D",  # strong cool gray
            "#495057",  # dark gray (not black)
        ]
        
        threshold_colors = [
            "#FF6B3A",  # vivid vermilion (brighter than prior rust tones)
            "#FFD23F",  # strong sunflower yellow
            "#B8E620",  # electric yellow-green (not neon green)
            "#00A86B",  # jade green (not bright green or cyan)
            "#4DA8DA",  # bright sky blue (not royal blue)
            "#7B2CBF",  # bright amethyst (not GOES purple)
        ]


    #Tab 20-like bright, soft rainbow primary cycle with high visibility
    if color_scheme == 13:
        flux_colors = [
            "#F94144",  # bright coral-red (not pure red)
            "#F3722C",  # vivid orange-red (distinct from GOES orange)
            "#F9C74F",  # bright golden yellow
            "#90BE6D",  # soft bright green (not neon)
            "#43AA8B",  # aqua-green (not cyan)
            "#4D908E",  # bright teal
            "#577590",  # strong steel blue
            "#4895EF",  # bright sky blue (not #0000ff)
            "#4361EE",  # clear blue-indigo
            "#B5179E",  # magenta (not GOES purple)
            "#F72585",  # bright pink
            "#FF8FAB",  # light rose
            "#FFB703",  # saturated yellow-orange
            "#A7C957",  # bright yellow-green
            "#06D6A0",  # vivid mint (not cyan)
            "#3A86FF",  # electric blue (not pure)
            "#8338EC",  # bright violet (distinct from GOES purple)
            "#FFBE0B",  # strong sunflower
            "#E76F51",  # bright terracotta
            "#2EC4B6",  # vibrant teal-green
        ]

        threshold_colors = [
            "#7B2CBF",  # deep amethyst
            "#1D3557",  # dark navy-slate (not pure blue)
            "#2D6A4F",  # deep forest green (not bright)
            "#9E2A2B",  # dark brick (not pure red)
            "#6B4226",  # dark brown
            "#3A0CA3",  # deep indigo (not GOES purple)
        ]


    #Flux energy bin colors
    color_map = {}
    #cycle through colors
    colorcycler = cycle(flux_colors)
    for i, bin in enumerate(energy_bins):
        label = tools.setup_energy_bin_label(bin)
        if label in integral_colors.keys() and not no_goes_colors:
            color_map.update({label: integral_colors[label]})
        else:
            color_map.update({label: next(colorcycler)})

    #Threshold colors
    threshcycler = cycle(threshold_colors)
    if event_definitions != None:
        for i, evdef in enumerate(event_definitions):
            bin = [evdef['energy_channel'].min, evdef['energy_channel'].max]
            label = tools.setup_energy_bin_label(bin)
            if label in color_map.keys():
                continue
            elif label in integral_colors.keys() and not no_goes_colors:
                color_map.update({label: integral_colors[label]})
            else:
                color_map.update({label: next(threshcycler)})


    return color_map


def vline_styles(event_definitions):
    """ Set line styles for vertical lines for different 
        event definitions.
        
    """
    lines = ["--","-.",":"]
    linecycler = cycle(lines)
    styles = []
    for i in range(len(event_definitions)):
        styles.append(next(linecycler))
        
    return styles
    


def plot_fluxes_basic(experiment, user_name, flux_type, dates, fluxes,
    energy_bins, showplot, saveplot, ylog=True, color_scheme=1, no_goes_colors=False,
    figname='basic.png'):
    """ Plot the fluxes for visualization only.
        
        INPUT:
        
            :experiment: (str) e.g. "GOES-13" or "user" for user input
            :user_name: (str) if user input or if specified, will be used
                instead of experiments
            :flux_type: (str) integral or differential
            :dates: (arr) dates for the fluxes that were evaluated
                using event definitions, called evaluated_dates in Data obj
            :fluxes: (2D arr, probably numpy) all fluxes that were
                evaluted with event definitions, called evaluated_fluxes in Data obj
            :color_scheme: (int) choices of 1 to 7
            :no_goes_colors: (bool) if True won't use SWPC colors for >5, >10, etc MeV

        OUTPUT:
        
            Multi-pane plot showing event definitions
    
    """
    flux_units = tools.get_flux_units(flux_type)
    flux_units = make_math_label(flux_units)
    energy_units = tools.get_energy_units()

    exp_name = experiment
    if experiment == "user":
        exp_name = user_name

    plot_title = f"{exp_name}"
    
    fig = plt.figure(figname,figsize=(16,8))
    ax = plt.subplot(111)
    colors = define_colors(energy_bins, color_scheme=color_scheme, no_goes_colors=no_goes_colors)
    
    #Plot the fluxes
    for j in range(len(energy_bins)):
        energy_bin = energy_bins[j]
        energy_label = tools.setup_energy_bin_label(energy_bin)
        is_good = check_for_good_fluxes(fluxes[j])
        if not is_good:
            continue
        
        #Don't plot integral channels on plots with differential units
        if flux_type == 'differential' and '>' in energy_label:
            continue

        maskfluxes = np.ma.masked_where(fluxes[j] <=0, fluxes[j])
        ax.plot(dates, maskfluxes,'.-', markersize=3, label=energy_label, color=colors[energy_label]) #, marker='.')

    plt.title(plot_title)
    ylabel = f"Intensity [{flux_units}]"
    plt.ylabel(ylabel)
    plt.xlabel('Date')
    plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d\n%H:%M"))
    plt.xticks(rotation=45, ha="right")

    plt.grid(axis="y")
    #ax.xaxis_date()
    #ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d\n%H:%M'))
    if 'counts' in flux_units: ylog=False
    if ylog: plt.yscale("log")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                     chartBox.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.01), fontsize=12)
    for item in ([ax.title, ax.yaxis.label] + ax.get_yticklabels()):
        item.set_fontsize(14)
    for item in ([ax.xaxis.label] + ax.get_xticklabels()):
        item.set_fontsize(10)

    if saveplot:
        print(f"plot_fluxes_basic: Saving plot {figname}")
        fig.savefig(figname)

    if not showplot:
        plt.close(fig)



def opsep_plot_event_definitions(experiment, flux_type, all_energy_bins, user_name, options,
    evaluated_dates, evaluated_fluxes, evaluated_energy_bins, event_definitions,
    sep_start_times, sep_end_times, onset_peaks, onset_peak_times,
    max_fluxes, max_flux_times, showplot, saveplot, spacecraft='',
    doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False,
    color_scheme=1, no_goes_colors=False):
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
 
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft, doBGSubOPSEP=doBGSubOPSEP, doBGSubIDSEP=doBGSubIDSEP, OPSEPEnhancement=OPSEPEnhancement, IDSEPEnhancement=IDSEPEnhancement)
    exp_name = experiment
    if experiment == "user":
        exp_name = user_name
 
    colors = define_colors(all_energy_bins, event_definitions=event_definitions,
        color_scheme=color_scheme, no_goes_colors=no_goes_colors)
 
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
        fig, ax = plt.subplots(nthresh, 1, sharex=True, figsize=(12,12), layout='constrained')
    else:
        #fig = plt.figure(figname,figsize=(12,9))
        fig, ax = plt.subplots(nthresh, 1, sharex=True, figsize=(12,9), layout='constrained')
        if nthresh == 1: ax = [ax]

    fig.canvas.manager.set_window_title(figname)
    plot_title = f"Event Definitions\n {exp_name} {title_mod} {flux_type} Fluxes"
    plt.suptitle(plot_title)

    for i in range(nthresh):
        #Get energy bin
        energy_bin = [event_definitions[i]['energy_channel'].min,
                    event_definitions[i]['energy_channel'].max]
        threshold = event_definitions[i]['threshold'].threshold
        flux_units = event_definitions[i]['threshold'].threshold_units

        threshold_label = f"{threshold} {flux_units}"
        threshold_label = make_math_label(threshold_label)

        #Event definitions and fluxes are in the same order
        dates = evaluated_dates #for ease
        fluxes = np.array(evaluated_fluxes[i])
        
        is_good = check_for_good_fluxes(fluxes)
        if not is_good:
            continue
 
        #Find energy bin to get index for colors
        energy_label = tools.setup_energy_bin_label(energy_bin)
 
        #Create labels
        ylabel = f"Intensity\n[{flux_units}]"
        ylabel = make_math_label(ylabel)
        data_label = f"{energy_label}"
        if energy_bin[1] == -1 and flux_type == "differential":
            data_label = f"Estimated {energy_label}"

    
        maskfluxes = np.ma.masked_where(fluxes <= 0, fluxes)
        ax[i].plot(dates, maskfluxes,'.-', markersize=3, label=data_label, color=colors[energy_label])

        if threshold != cfg.opsep_min_threshold:
            ax[i].axhline(threshold,color='black',linestyle='--', label=f"Threshold {threshold_label}")
        else:
            ax[i].axhline(max_fluxes[i],color='black',linestyle='--', lw=0, label=f"Above background")

        if not pd.isnull(sep_start_times[i]):
            ax[i].axvline(sep_start_times[i],color='black',linestyle=':', linewidth=2)
            ax[i].axvline(sep_end_times[i],color='black',linestyle=':',
                        label="Start, End", linewidth=2)


        if not pd.isnull(onset_peaks[i]) and not pd.isnull(onset_peak_times[i]):
            ax[i].plot_date(onset_peak_times[i],onset_peaks[i],'o',color='black',
                    label="Onset Peak")
        if not pd.isnull(max_fluxes[i]) and not pd.isnull(max_flux_times[i]):
            ax[i].plot_date(max_flux_times[i],max_fluxes[i],'s',mfc='none', color='black',
                    label="Max Flux", mew=2)


        if i == nthresh-1:
            ax[i].set_xlabel('Date')
            plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d\n%H:%M"))
            plt.xticks(rotation=45, ha="right")
        ax[i].set_ylabel(ylabel)
        ax[i].set_yscale("log")
        if 'counts' in flux_units:
            ax[i].set_yscale("linear")
        ax[i].legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
        for item in ([ax[i].title, ax[i].xaxis.label, ax[i].yaxis.label] + ax[i].get_xticklabels() + ax[i].get_yticklabels()):
            item.set_fontsize(10)

        if nthresh <= 2:
            for item in ([ax[i].title, ax[i].yaxis.label] + ax[i].get_yticklabels()):
                item.set_fontsize(14)

    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath,'opsep', subdir, figname + '.png'))
    if not showplot:
        plt.close(fig)



def opsep_plot_all_bins(experiment, flux_type, user_name, options,
    all_dates, all_fluxes, all_energy_bins, event_definitions,
    sep_start_times, sep_end_times, showplot, saveplot, spacecraft='',
    doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False, color_scheme=1,
    no_goes_colors=False):
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
    
    fig = plt.figure(figname,figsize=(13.5,8))
    ax = plt.subplot(111)
    colors = define_colors(all_energy_bins, color_scheme=color_scheme,
        event_definitions=event_definitions, no_goes_colors=no_goes_colors)
    vstyles = vline_styles(event_definitions)

    #Plot the fluxes
    for j in range(len(all_energy_bins)):
        energy_bin = all_energy_bins[j]
        energy_label = tools.setup_energy_bin_label(energy_bin)
    
        is_good = check_for_good_fluxes(all_fluxes[j])
        if not is_good:
            continue

        #Don't plot integral channels on plots with differential units
        if flux_type == 'differential' and '>' in energy_label:
            continue

        maskfluxes = np.ma.masked_where(all_fluxes[j] <=0, all_fluxes[j])
        ax.plot(all_dates, maskfluxes,'.-', markersize=3, label=energy_label, color=colors[energy_label]) #, marker='.')

    #Plot the threshold crossing times
    for i in range(len(event_definitions)):
        energy_bin = [event_definitions[i]['energy_channel'].min,
                    event_definitions[i]['energy_channel'].max]
        energy_label = tools.setup_energy_bin_label(energy_bin)
        
        threshold = event_definitions[i]['threshold'].threshold
        flux_units = event_definitions[i]['threshold'].threshold_units

        threshold_label = f"{threshold} {flux_units}"
        threshold_label = make_math_label(threshold_label)

        if threshold == cfg.opsep_min_threshold:
            threshold_label = "above background"

        line_label = f"{energy_label}, {threshold_label}"

        if not pd.isnull(sep_start_times[i]):
            ax.axvline(sep_start_times[i],color=colors[energy_label],linestyle=vstyles[i],
                        label=line_label, linewidth=2)
            ax.axvline(sep_end_times[i],color=colors[energy_label],linestyle=vstyles[i], linewidth=2)

    plt.title(plot_title)
    
    #Flux in original energy bins
    flux_units = tools.get_flux_units(flux_type)
    ylabel = f"Intensity [{flux_units}]"
    ylabel = make_math_label(ylabel)
    plt.ylabel(ylabel)
    plt.xlabel('Date')
    plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d\n%H:%M"))
    plt.xticks(rotation=45, ha="right")

    plt.grid(axis="y")
    plt.yscale("log")
    if 'counts' in flux_units:
        plt.yscale("linear")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                     chartBox.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.01), fontsize=10)
    for item in ([ax.title, ax.yaxis.label] + ax.get_yticklabels()):
        item.set_fontsize(14)
    for item in ([ax.xaxis.label] + ax.get_xticklabels()):
        item.set_fontsize(10)
    if saveplot:
        fname = os.path.join(cfg.plotpath,'opsep',subdir, figname + '.png')
        fig.savefig(fname)
    if not showplot:
        plt.close(fig)
    
    


def opsep_plot_fluence_spectrum(experiment, flux_type, user_name, options,
    event_definitions, evaluated_dates, energy_bin_centers, fluence_spectra,
    fluence_spectra_units, showplot, saveplot, spacecraft='',
    doBGSubOPSEP=False, doBGSubIDSEP=False,
    OPSEPEnhancement=False, IDSEPEnhancement=False,
    color_scheme=1, no_goes_colors=False):
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

    energy_bins = []
    if experiment != "user":
        expt = expts.experiment_info(experiment)
        energy_bins = expt[flux_type]['energy_bins']
    else:
        energy_bins = cfg.user_energy_bins

    #plot integral fluxes (either input or estimated)
    nthresh = len(event_definitions)
    energy_units = event_definitions[0]['energy_channel'].units
    fluence_units = fluence_spectra_units[0]
 
    ylabel = "Fluence"
 
    ncross = 0 #Check if any thresholds were crossed, if not no need plot

    fig = plt.figure(figname,figsize=(10,8))
    ax = plt.subplot(111)
    markers = ['o','P','D','v','^','<','>','*','d','+','8','p','h','1','X','x']
    colors = define_colors(energy_bins, event_definitions=event_definitions, color_scheme=color_scheme,
        no_goes_colors=no_goes_colors)
    
    for i in range(nthresh):
    
        if len(fluence_spectra[i]) == 0:
            continue
        ncross = ncross + 1

        energy_bin = [event_definitions[i]['energy_channel'].min, event_definitions[i]['energy_channel'].max]
        energy_label = tools.setup_energy_bin_label(energy_bin)
        
        flspec_units = fluence_spectra_units[i]
        threshold_label = f"{event_definitions[i]['threshold'].threshold} {event_definitions[i]['threshold'].threshold_units}"
        threshold_label = make_math_label(threshold_label)

        if event_definitions[i]['threshold'].threshold == cfg.opsep_min_threshold:
            threshold_label = "above background"

        #Create labels
        legend_label = f"{energy_label}, {threshold_label}"

        ax.plot(energy_bin_centers,fluence_spectra[i],markers[i],
                color=colors[energy_label], mfc='none', label=legend_label, linewidth=3)
    
    plt.grid(which="both", axis="both")
    plot_title = f"Event-Integrated Fluence Spectra\n {exp_name} {title_mod} {flux_type}"
    plt.title(plot_title)
    plt.xlabel(f"Energy [{energy_units}]")
    ylabel = f"Fluence [{fluence_units}]"
    ylabel = make_math_label(ylabel)
    plt.ylabel(ylabel)

    plt.xscale("log")
    plt.yscale("log")
    if 'counts' in fluence_units:
        plt.yscale("linear")
    ax.legend(loc='upper right')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(14)

    if ncross == 0: plt.close(fig) #no thresholds crossed, empty plot

    if saveplot:
        fig.savefig(os.path.join(cfg.plotpath,'opsep', subdir, figname + '.png'))
    if not showplot:
        plt.close(fig)



############ IDSEP PLOTS ##############
def setup_idsep_plot(figname, experiment, title_mod, unique_id, flux_units, nrow):
    """ Set up figure and axes for idsep nrow row plots.
    
    """
    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots(nrow, 1, sharex=True, figsize=(13,8), gridspec_kw={'height_ratios' : [1, 1, 1], 'hspace' : 0.1})
    fig.canvas.manager.set_window_title(figname)
    fig.suptitle((f"{experiment} {title_mod} {unique_id}"))
    
    #Formatting of axes
    flux_label = make_math_label(f"Intensity ({flux_units})")
    ax[1].set_ylabel(flux_label)
    ax[2].set_xlabel("Date")
    
    #Apply tight layout and suppress warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=UserWarning)
        plt.tight_layout()

    for iax in range(nrow):
        ax[iax].set_yscale('log')
        if 'counts' in flux_units:
            ax[iax].set_yscale('linear')
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
    flux_units = tools.get_flux_units(flux_type)
    flux_units = make_math_label(flux_units)
    energy_units = tools.get_energy_units()

    nbins = len(energy_bins)
    nrow = 3 #3 plots per page
    for i in range(nbins):
        if not i%3:
            exp = experiment
            if experiment == 'user' and exp_name != '':
                exp = exp_name
            fig, ax = setup_idsep_plot((f"{figname}_{i}"), exp, title_mod, unique_id, flux_units, nrow)
            ax_right = [0]*len(ax)
            iax = 0

        legend_label = tools.setup_energy_bin_label(energy_bins[i])
            
        #PLOT FLUXES
        maskfluxes = np.ma.masked_less_equal(fluxes[i], 0)
        if np.isnan(maskfluxes).all():
            print(f"idsep_make_plots: All values in flux array for {legend_label} are nan. Skipping.")
            iax += 1
            continue
        ax[iax].plot(dates,maskfluxes,'.-', markersize=3, label="Background Fluxes", color='tab:blue')
    
        #PLOT BACKGROUND
        if not disable_sigma:
            ax[iax].errorbar(ave_dates, ave_fluxes[i],fmt='.', yerr=ave_sigma[i], label=(f"Mean Background and Sigma"),zorder=100, color='tab:orange')
        if disable_sigma:
            ax[iax].errorbar(ave_dates, ave_fluxes[i],fmt='-',
                label=(f"Mean Background"),zorder=100, color='tab:orange')
        
        #PLOT THRESHOLD = n*sigma (n in config file)
        ax[iax].plot(threshold_dates,threshold[i],'-',label=(f"Threshold"), zorder=200, color='tab:green')

        if iax == 0:
            ax[iax].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=3, mode="expand", borderaxespad=0., fontsize=11,
                      facecolor='whitesmoke', edgecolor='black')

        #Put energy channel label on right y-axis
        ax_right[iax] = ax[iax].twinx()
        ax_right[iax].set_ylabel(f"{legend_label}\n ", rotation=270, labelpad=20)
        ax_right[iax].set_yticks([])
    
        #Set clean y-axis ranges
        ymin = 10 ** np.floor(np.log10(np.nanmin(maskfluxes)))
        ymax = 10 ** np.ceil(np.log10(np.nanmax(maskfluxes)))
        if 'counts' in flux_units:
            ymin = np.floor(np.nanmin(maskfluxes))
            ymax = np.ceil(np.nanmax(maskfluxes))
        ax[iax].set_ylim([ymin, ymax])


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
        spacecraft="", close_plot=False, subdir="idsep"):

    #Additions to titles and filenames according to user-selected options
    modifier, title_mod = tools.setup_modifiers(options, spacecraft=spacecraft)
    name = tools.idsep_naming_scheme(experiment, flux_type, exp_name, options, spacecraft=spacecraft)


    figname = (f"{name}_{unique_id}")
 
    #UNITS
    flux_units = tools.get_flux_units(flux_type)
    flux_units = make_math_label(flux_units)
    energy_units = tools.get_energy_units()
    
    nbins = len(energy_bins)
    nrow = 3 #3 plots per page
    for i in range(nbins):
        if not i%nrow:
            exp = experiment
            if experiment == 'user' and exp_name != '':
                exp = exp_name
            fig, ax = setup_idsep_plot((f"{figname}_{i}"), exp, title_mod, unique_id, flux_units, nrow)
            ax_right = [0]*len(ax)
            iax = 0

        legend_label = tools.setup_energy_bin_label(energy_bins[i])

        maskfluxes = np.ma.masked_less_equal(fluxes[i], 0)
        if np.isnan(maskfluxes).all():
            print(f"idsep_make_timeseries_plot: All values in flux array for {legend_label} are nan. Skipping.")
            iax += 1
            continue
        ax[iax].plot(dates,maskfluxes,'.-', markersize=3, label=legend_label)

        #Put energy channel label on right y-axis
        ax_right[iax] = ax[iax].twinx()
        ax_right[iax].set_ylabel(f"{legend_label}\n ", rotation=270, labelpad=20)
        ax_right[iax].set_yticks([])
    
        #Set clean y-axis ranges
        ymin = 10 ** np.floor(np.log10(np.nanmin(maskfluxes)))
        ymax = 10 ** np.ceil(np.log10(np.nanmax(maskfluxes)))
        if 'counts' in flux_units:
            ymin = np.floor(np.nanmin(maskfluxes))
            ymax = np.ceil(np.nanmax(maskfluxes))
        ax[iax].set_ylim([ymin, ymax])
        
        if saveplot and (iax ==2 or i == nbins-1):
            fig.savefig(os.path.join(cfg.plotpath,subdir, name, (f"{figname}_{i}.png")))
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

    #UNITS
    flux_units = tools.get_flux_units(flux_type)
    flux_units = make_math_label(flux_units)
    energy_units = tools.get_energy_units()

    nbins = len(energy_bins)
    nrow = 3
    for i in range(nbins):
        if not i%3:
            exp = experiment
            if experiment == 'user' and exp_name != '':
                exp = exp_name
            fig, ax = setup_idsep_plot((f"{figname}_{i}"), exp, title_mod, unique_id, flux_units, nrow)
            ax_right = [0]*len(ax)
            iax = 0

        legend_label = tools.setup_energy_bin_label(energy_bins[i])

        maskfluxes_bg = np.ma.masked_less_equal(fluxes_bg[i], 0)
        if np.isnan(maskfluxes_bg).all():
            print(f"idsep_make_bg_sep_plot: All values in flux array for {legend_label} are nan. Skipping.")
            iax += 1
            continue
        ax[iax].plot(dates,maskfluxes_bg,'.-', markersize=3, label=(f"Background"))
        maskfluxes_sep = np.ma.masked_less_equal(fluxes_sep[i], 0)
        ax[iax].plot(dates,maskfluxes_sep,'.-', markersize=3, label=(f"Enhanced"), zorder=100)

        if iax == 0:
            ax[iax].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=3, mode="expand", borderaxespad=0., fontsize=11,
                      facecolor='whitesmoke', edgecolor='black')


#        ax[iax].legend(loc='upper left', fontsize=11)


        #Put energy channel label on right y-axis
        ax_right[iax] = ax[iax].twinx()
        ax_right[iax].set_ylabel(f"{legend_label}\n ", rotation=270, labelpad=20)
        ax_right[iax].set_yticks([])
    
        #Set clean y-axis ranges
        ymin = 10 ** np.floor(np.log10(np.nanmin(maskfluxes_bg)))
        ymax = 10 ** np.ceil(np.log10(np.nanmax(maskfluxes_sep)))
        if 'counts' in flux_units:
            ymin = np.floor(np.nanmin(maskfluxes_bg))
            ymax = np.ceil(np.nanmax(maskfluxes_sep))
        ax[iax].set_ylim([ymin, ymax])


        if saveplot and (iax ==2 or i == nbins-1):
            fig.savefig(os.path.join(cfg.plotpath,"idsep", name, (f"{figname}_{i}.png")))
            if not showplot:
                plt.close(fig)
            if close_plot:
                plt.close(fig)
        
        #increment to next axis
        iax += 1



