import pandas as pd
import datetime
import os
import math
import numpy as np
import matplotlib.pylab as plt
import matplotlib.pyplot as mplt
import matplotlib.dates as mdates
from statistics import mode
from sklearn.utils.validation import check_consistent_length
from sklearn.utils.validation import check_array
from scipy.stats import pearsonr
from scipy.optimize import curve_fit 
from scipy.ndimage import uniform_filter1d
from lmfit import minimize, Parameters, fit_report
# from . import tools
# from ..opsep import calculate_fluence

def determine_time_resolution(dates):
    """ The time resolution is found by taking the difference between
        every consecutive data point. The most common difference is
        taken as the time resolution. Even if the data set has gaps,
        if there are enough consecutive time points in the observational
        or model output, the correct time resolution should be identified.
        
        INPUTS:
        
        :dates: (datetime 1xm array) - dates associated with flux time profile
        
        OUTPUTS:
        
        :time_resolution: (time delta object)
        
    """
    ndates = len(dates)
    time_diff = [str_to_datetime(a) - str_to_datetime(b) for a,b in zip(dates[1:ndates],dates[0:ndates-1])]
    if not time_diff:
        sys.exit("determine_time_resolution: Require more than 1 data point "
                "to determine time resolution. Please extend your "
                "requested time range and try again. Exiting.")
    time_resolution = mode(time_diff)
    return time_resolution

def str_to_datetime(date):
    """ String date to datetime
        
        INPUTS:
        
        :date: (string) - date as "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS"
        
        OUTPUTS:
        
        :dt: (datetime) - datetime conversion of date
    """
    if len(date) == 10: #only YYYY-MM-DD
        date = date  + ' 00:00:00'
    dt = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
    return dt

def calculate_fluence(dates, flux):
    """ This subroutine sums up all of the flux in the 1xn array "flux". The
        "dates" and "flux" arrays input here should reflect only the intensities
        between the SEP start and stop times, determined by the subroutine
        calculate_threshold_crossing. The subroutine does not differentiate between
        differential or integral flux.
        
        The extract_date_range subroutine is used prior to calling this one to
        make the dates and fluxes arrays covering only the SEP time period.
        The flux will be multiplied by time_resolution and summed for all of the
        data between the start and end times. Negative or bad flux values
        should have been set to None or interpolated with check_bad_data().
        None values will be skipped.
        
        The time resolution is found by taking the difference between
        every consecutive data point. The most common difference is
        taken as the time resolution. Even if the data set has gaps,
        if there are enough consecutive time points in the observational
        or model output, the correct time resolution should be identified.
        
        Fluence units will be 1/[MeV cm^2] for GOES or SEPEM differential
        fluxes or 1/[cm^2] for integral fluxes.
        
        INPUTS:
        
        :flux: (float 1xn array) - intensity time series for a single energy bin
            or single integral channel (1D array).
        :dates: (datetime 1xn array) - datetimes that correspond to the fluxes
        
        OUTPUTS:
        
        :fluence: (float) - sum of all the flux values in flux
        
    """
    ndates = len(dates)
    time_resolution = determine_time_resolution(dates)
    
#    print("calculate_fluence: Identified a time resolution of "
#            + str(time_resolution.total_seconds()) + " seconds.")
    
    fluence = 0
    for i in range(ndates):
        
        if flux[i] == None or flux[i] == 0.0 or math.isnan(flux[i]):
            continue

        if flux[i] >= 0:  #0 flux ok for models
                fluence = fluence + flux[i]*time_resolution.total_seconds()
        else:
            sys.exit('calculate_fluence: Bad flux data value of ' + str(flux[i]) +
                    ' found for bin ' + str(i) + ', '
                    + str(dates[i]) + '. This should not happen. '
                    + 'Did you call check_for_bad_data() first?')
                    
    fluence = fluence*4.0*math.pi #multiply 4pi steradians
    return fluence

def from_fetchsep(sep_dates, sep_fluxes, energy_bins, flux_type, plot_fit = False):
    """
    function that takes the input from fetchsep and passes it into 
    this new bit. This step should not be done for every event
    but only ones where a flag is hit for historical events
    and for all real-time FetchSEP runs
    
    
    needs:
    date array
    fluxes (all energies)
    energy bins
    
    output:
    start time of the secondary oneset

    """
    
    onset_times = function_two(sep_dates, sep_fluxes, energy_bins, flux_type, plot_fit)



    return onset_times

def function_two(dates, fluxes, energy_bins, flux_type, plot_fit = False):
    """
    the new bit

    takes the input fluxes and dates converts into a fluence
    spectra (~30 mins to an hour) to find when a new onset occurs when particle
    fluxes are already elevated


    passes new fluence spectrum to fitting_routine, assesses the fit parameters,
    then determines when to set the start time of the new event
    """


    n_dates = len(dates)
    m_energies = len(energy_bins)
    # create a new array to calculate the fluence from
    end_date = str_to_datetime(dates[0]) # start date to ensure we get to the while loop
    current_date = dates[6] # starting 30 mins into a flux file  


    fig = plt.figure(figsize=(20,10))
    ax = plt.subplot(111)
    if flux_type == 'differential':
        # lower_error = [4.2, 8.7, 15.0, 38.0, 84.0, 110.0, 330.0, 420.0, 510.0, 700.0]
        # upper_error =  [8.7, 14.5, 40.0, 38.0, 200.0, 900.0, 420.0, 510.0, 700.0, 700.0]
        lower_error = []
        upper_error = []
        for en in range(len(energy_bins)):
            if len(energy_bins[en].rsplit('-')) != 3:
                lower_error.append(float(energy_bins[en].rsplit('-')[0]))
                upper_error.append(float(energy_bins[en].rsplit('-')[1]))
            else:
                lower_error.append(float(energy_bins[en].rsplit('-')[0]))
                upper_error.append(1000.0)
        asymmetric_error = np.array(list(zip(lower_error, upper_error))).T
        bin_centers = []
        for m in range(m_energies):
            bin_centers.append((lower_error[m] + upper_error[m]) / 2)
    else:
        bin_centers = energy_bins
        
    


    # now creating the fluence spectrum

    current_index = dates.index[dates == current_date][0]
    current_datetime = str_to_datetime(current_date)
    foo = str_to_datetime(current_date)
    temp_date = foo - datetime.timedelta(minutes = 30) # 30 minutes backwards
    temp_index = dates.index[dates == str(temp_date)][0]
    end_date = temp_date
    low_energy_slope = []
    gammas = []
    breaks = []
    plot_times = []
    time_since_start = []
    chisqr = []
    time = 0
    while foo <= str_to_datetime(dates[n_dates-1]):
        
        
        current_index = dates.index[dates == current_date][0]
        # print(current_date, type(current_date), foo <= str_to_datetime(dates[n_dates-1]), end_date, current_index, temp_index)

        # maybe add something here to limit this loop - fitting can take a lot of time
        # only run when above background? 
        m = 0
        spectrum = []
        for m in range(m_energies):
    
            spectrum.append(calculate_fluence(dates[temp_index:current_index], fluxes[m][temp_index:current_index].tolist()))

        # now to do the fitting
        fit_parameters, fit_band, energy_array = fitting_routine(bin_centers, spectrum)

        if plot_fit:
            if flux_type == 'differential':
                ax.errorbar(bin_centers, spectrum, xerr = asymmetric_error, fmt = 'o', label = current_date)
            else:
                ax.plot(bin_centers, spectrum, fmt = 'o', label = current_date)
            ax.plot(energy_array, fit_band, linestyle = 'dashed', label = 'Band')


        low_energy_slope.append(fit_parameters['Gamma A'])
        gammas.append(fit_parameters['Gamma B'])
        breaks.append(fit_parameters['Break Energy'])
        plot_times.append(str_to_datetime(current_date))
        chisqr.append(fit_parameters['ChiSqr'])

        current_date = dates[current_index + 1]
        foo = str_to_datetime(current_date)
        temp_date = foo + datetime.timedelta(minutes = 5) # 30 minutes just for now
        time += 5
        time_since_start.append(time)
        if temp_date > str_to_datetime(dates[n_dates - 1]):
            break
        else:
            temp_index += 1
            # end_index = dates.index[dates == str(temp_date)][0]
            
            end_date = temp_date
    # ax.legend(loc = 'lower left', fontsize='10', framealpha=1.0) 
    plt.yscale('log')
    plt.xscale('log')
    savename = 'band_fitting_' + flux_type + '_running_avg.png'

    fig.savefig(savename)

    plt.close(fig)
        

   

    # running_avg_gamma = np.convolve(gammas, np.ones(6)/6, mode='same')
    # running_avg_energy = np.convolve(breaks, np.ones(6)/6, mode='same')
    # fig = plt.figure(figsize=(20,10))
    # ax = plt.subplot(111)
    # ax.plot(plot_times, running_avg_gamma, color = 'black')
    # ax.set_xlabel('Dates')
    # ax.set_ylabel('High Energy Spectral Slope')
    # axis_2 = ax.twinx()
    # axis_2.plot(plot_times, running_avg_energy, color = 'blue')
    # axis_2.set_ylabel('Break Energy (MeV)')
    # axis_2.set_yscale('log')
    # # Set major tick locator to place ticks at the beginning of each month
    # ax.xaxis.set_major_locator(mdates.DayLocator())

    # # Set major tick formatter to display month and year
    # ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %Y'))
    # savename = 'test_time_evo_' + flux_type + '_' + '20120101_running_avg.png'
    # fig.savefig(savename)
    # plt.close(fig)

    # input()
    # running_avg_energy = uniform_filter1d(breaks, size = 6)
    # # now to determine if there was a new onset
    # derivatives
    dgamma_b = np.gradient(gammas)
    dgamma_a = np.gradient(low_energy_slope)
    ddates = np.gradient(time_since_start)
    dbreak = np.gradient(breaks)

    deriv_gamma_b = dgamma_b / ddates
    deriv_gamma_a = dgamma_a / ddates
    deriv_break = dbreak / ddates


    ddbreak = np.gradient(deriv_break)
    
    sec_deriv_break = ddbreak / ddates

    high_energy = 100.0
    for e in range(m_energies):
        if flux_type == 'differential':
            if high_energy > lower_error[e] and high_energy < upper_error[e]:
                energy_index = e
        else:
            #tbd
            energy_index = 5
            continue

    onset_times = []
    r = 3
    for r in range(len(gammas)-3):
        current_gamma = gammas[r]
        next_gamma = gammas[r + 1]
        # what's a reasonable change in gammma to detect an onset?
        # 
        # print(plot_times[r], d_gamma[r], d_breaks[r], current_gamma, next_gamma, gammas[r+2], gammas[r+3])
        
        if next_gamma - current_gamma < -0.2 or current_gamma - gammas[r-1] < -0.2: # and np.abs(gammas[r+2] - next_gamma) < 0.1: # and np.abs(gammas[r+3] - gammas[r+2]) < 0.01:
            if np.abs(breaks[r] - breaks[r-1]) > 10.0 or np.abs(breaks[r] - breaks[r+1]) > 10.0:
                # if chisqr[r] / np.mean(chisqr) < 1.0:
                # if np.abs(low_energy_slope[r] - low_energy_slope[r-1]) > 0.2 or np.abs(low_energy_slope[r] - low_energy_slope[r+1]) > 0.2:
                    if np.abs(deriv_gamma_b[r+1] - deriv_gamma_b[r]) > 0.15 or np.abs(deriv_gamma_b[r-1] - deriv_gamma_b[r]) > 0.25 : #and deriv_gamma_b[r] < 0.1: # and np.abs(deriv_gamma_a[r]) > 0.02:  # Detecting large changes in gamma_b
                        if fluxes[energy_index][r] > fluxes[energy_index][r-1] and fluxes[energy_index][r] > fluxes[energy_index][r-2] and fluxes[energy_index][r] < fluxes[energy_index][r+1]: # and fluxes[energy_index][r] < fluxes[energy_index][r+2] and fluxes[energy_index][r] < fluxes[energy_index][r+3]:
                        # if breaks[r] < 700.0 or breaks[r+1] < 700.0 or breaks[r+2] < 700.0: # or np.abs(breaks[r] - breaks[r+1]) > 25.0 and deriv_break[r-1] > 10.0 :
                # further filter out points 
                            if onset_times == []:
                                # Gonna say this is an onset for now
                                # print('Found an Onset')
                                onset_times.append(plot_times[r])
                                # print(dates[r])
                                
                            elif plot_times[r] - onset_times[len(onset_times)-1] < datetime.timedelta(hours = 4):
                                # too close to another onset skip
                                continue
                            else:
                                # Gonna say this is an onset for now
                                # print('Found an Onset')
                                onset_times.append(plot_times[r])
                                # print(dates[r], chisqr[r], chisqr[r] / np.mean(chisqr) < 1.0)
                                # print(onset_times)
                    
            # if fluxes[energy_index][r] > fluxes[energy_index][r-1] and fluxes[energy_index][r] > fluxes[energy_index][r-2]:# and fluxes[energy_index][r] > fluxes[energy_index][r-3]:
                # low_energy_slope[r+1] - low_energy_slope[r] < -0.5 and
                # np.abs(low_energy_slope[r+1] - low_energy_slope[r]) > 0.75
                    
                
    if plot_fit:
        fig = plt.figure(figsize=(20,10))
        ax = plt.subplot(111)
        ax.plot(plot_times, low_energy_slope, color = '#BFDEAC', label = 'Gamma A')
        ax.plot(plot_times, gammas, color = 'red', label = 'Gamma B')
        # mplt.axvline(should_get_this_line, color = 'purple', linestyle = '--')
        for o in range(len(onset_times)):
            mplt.axvline(onset_times[o], color = 'black', linestyle = '--')
        
        ax.set_xlabel('Dates')
        ax.set_ylabel('Spectral Slopes')
        axis_2 = ax.twinx()
        axis_2.plot(plot_times, breaks, color = 'blue', label = 'Break Energy')
        axis_2.set_ylabel('Break Energy (MeV)')
        # axis_2.set_yscale('log')
        
        # Set major tick locator to place ticks at the beginning of each month
        ax.xaxis.set_major_locator(mdates.DayLocator())

        # Set major tick formatter to display month and year
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
        minor_locator = mdates.HourLocator(byhour=range(0, 24, 1)) # 
        ax.xaxis.set_minor_locator(minor_locator) 
        plt.legend()
        savename = 'time_evo_' + flux_type + '_' + str(dates[0]).rsplit(' ')[0] + '_running_avg.png'
        fig.savefig(savename)
        mplt.show()
        plt.close(fig)

            
        fig = plt.figure(figsize=(20,10))
        ax = plt.subplot(111)
        ax.plot(plot_times, deriv_gamma_a, color = '#BFDEAC', label = 'Gamma A')
        ax.plot(plot_times, deriv_gamma_b, color = 'red', label = 'Gamma B')
        # mplt.5e(should_get_this_line, color = 'purple', linestyle = '--')
        for o in range(len(onset_times)):
            mplt.axvline(onset_times[o], color = 'black', linestyle = '--')
        ax.set_xlabel('Dates')
        ax.set_ylabel('Spectral Slopes')
        axis_2 = ax.twinx()
        axis_2.plot(plot_times, deriv_break, color = 'blue', label = 'Break Energy')
        # axis_2.plot(plot_times, sec_deriv_break, color = '#BFDEAC', label = 'Second Derivative Break Energy')
        axis_2.set_ylabel('Break Energy (MeV)')
        
        
        # Set major tick locator to place ticks at the beginning of each month
        ax.xaxis.set_major_locator(mdates.DayLocator())

        # Set major tick formatter to display month and year
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
        minor_locator = mdates.HourLocator(byhour=range(0, 24, 1)) #  # Locates ticks at 00:00, 06:00, 12:00, 18:00
        ax.xaxis.set_minor_locator(minor_locator) 
        plt.legend()
        savename = 'derivs_' + flux_type + '_'+ str(dates[0]).rsplit(' ')[0] +'_running_avg.png'
        fig.savefig(savename)
        mplt.show()
        plt.close(fig)
    
    
    return onset_times


def fitting_routine(energy_bins, fluence_spectrum):
    """
    takes the fluence spectra passed from function_two
    and fits it to a standard (band/ellison-ramaty/powerlaw) form.
    """



    min_obs_band_params, fit_obs_band, transition_obs_band, energy = fit_band_function(energy_bins, fluence_spectrum)
    # print('Best Fit Band Parameters Obs: ', min_obs_band_params.params, min_obs_band_params.chisqr)#, min_obs_band.message)
    obs_band_params = min_obs_band_params.params.valuesdict() 
    obs_norm_band = obs_band_params['norm']
    obs_gamma_a_band = obs_band_params['gamma_a']
    obs_gamma_b_band = obs_band_params['gamma_b']
    obs_e_0_band = obs_band_params['break_energy']
    chisqr = min_obs_band_params.chisqr
    # print(fit_report(min_obs_band_params))
    

    fit_parameters ={'Normalization': obs_norm_band, 'Gamma A': obs_gamma_a_band, 'Gamma B': obs_gamma_b_band, 'Break Energy': obs_e_0_band, 'ChiSqr': chisqr}
    
    return fit_parameters, fit_obs_band, energy


###################### Band Function Subroutines ##############################################
def fit_band_function(energy_bins, fluence_array):
    """ Fitting the Band function to spectral data - see
    description of fit_ellisonramaty_function for full details of 
    inputs/outputs
    """
    from lmfit import minimize, Parameters
    import numpy as np
    energy = np.arange(1,1000,1) # Creating an array of energies across all bins (1 - 1000 MeV) to use for plotting the end fit
    params_band = Parameters()
    params_band.add('norm', value = 10000.0, min = 1.0, max = 10**13)
    params_band.add('gamma_b', value = 5.0, min = 0.1, max =20)
    params_band.add('break_energy', value = 7.0, min = 5.0, max = 1000)
    params_band.add('gamma_a', value = 0.70, min = 0.1, max = 4.0)
    minimize_band = minimize(residual_band, params_band, args= [energy_bins, fluence_array], nan_policy= 'propagate', max_nfev= 100000000)
    
    fit_band, transition_band = band_func(minimize_band.params, energy) #Band
    return minimize_band, fit_band, transition_band, energy

def residual_band(first_guess, *args):
    """
    Defining the Band function for minimization

    Inputs:
    First Guess. Array of length 5, first guess for the values
        the free parameters will take.
        Norm - normalization factor (A)
        Gamma_a - low energy spectral index
        Gamma_b - high energy spectral index
        E_0 -  break energy
       

    args. Array of length 2, contains any other arguments used in the function
        since least_squares requires functions of the form func(x, *args)
        args[0] - energy bins
        args[1] - fluence values       
    """
    parvals = first_guess.valuesdict()
    norm = parvals['norm']
    gamma_a = parvals['gamma_a']
    gamma_b = parvals['gamma_b']
    E_0 = parvals['break_energy']


    
    fluence = args[1]
    energy = args[0]

    band = []
    fit_error = []
    total_fit = []
    e = 0

    E_t = (gamma_b - gamma_a)*E_0

    for e in range(len(energy)-1):
        
        if energy[e] <= E_t:

            band.append(norm*energy[e]**(-gamma_a)*math.exp(-energy[e] / E_0))
            
        else:
            band.append(norm*energy[e]**(-gamma_b)*(((gamma_b-gamma_a)*E_0)**(gamma_b-gamma_a))*math.exp(gamma_a-gamma_b))   
        total_fit.append((band[e]))
        fit_error.append(np.abs(fluence[e] / total_fit[e])-1)
        e = e+1

    return fit_error

def band_func(first_guess, energy):
    """
    Defining the Band function for plotting the fit
    """
    

    parvals = first_guess.valuesdict()
    norm = parvals['norm']
    gamma_a = parvals['gamma_a']
    gamma_b = parvals['gamma_b']
    E_0 = parvals['break_energy']
        # print(norm, gamma_a, gamma_b, E_0, E_t, E_r)
    
    band = []
    e = 0
    E_transition = (gamma_b - gamma_a)*E_0
    # print('Transition Energy ', E_transition)

    for e in range(len(energy)):
        
        if energy[e] <= (gamma_b - gamma_a)*E_0:
            band.append(norm*energy[e]**(-gamma_a)*math.exp(-energy[e] / E_0))
        else:
            band.append(norm*energy[e]**(-gamma_b)*((gamma_b-gamma_a)*E_0)**(gamma_b-gamma_a)*math.exp(gamma_a-gamma_b))
        band[e] = band[e]
        e = e+1

    return band, E_transition



if __name__ == "__main__":
    # building the main function for testing with manual input files
    
    test_file = "C:\\Users\\cfalliso\\Documents\\Dev_SPHINX\\fetchsep\\dev_fetchsep\\data\\fluxes_GOES-07_differential_19891015_19891108.csv"
    #fluxes_GOES-13_differential_20120301_20120320.csv"
    #fluxes_GOES-13_differential_20120113_20120129.csv"
    #fluxes_GOES-18_differential_20240505_20240529.csv"
    #fluxes_GOES-07_differential_19891015_19891108.csv"
    # with open(test_file, 'r', newline='') as csv_file:
    #     csv_read = csv.reader(csv_file)
    #     for row in csv_read:
    #         if row[0] == '#'
    df = pd.read_csv(test_file)
    #dates,4.2-8.7,8.7-14.5,15.0-40.0,38.0-82.0,84.0-200.0,110.0-900.0,330.0-420.0,420.0-510.0,510.0-700.0,700.0--1
    dates = df['#dates']
    if 'differential' in test_file:
        fluxes = [df['4.2-8.7'], df['8.7-14.5'], df['15.0-44.0'], df['39.0-82.0'], df['84.0-200.0'], df['110.0-500.0']]#, df['330.0-420.0'], df['420.0-510.0'], df['510.0-700.0'], df['700.0--1']]
        energy_bins = ['4.2-8.7','8.7-14.5','15.0-44.0','39.0-82.0','84.0-200.0','110.0-500.0']#,'330.0-420.0','420.0-510.0','510.0-700.0','700.0--1']
        flux_type = 'differential'
        # fluxes = [df['4.2-8.7'], df['8.7-14.5'], df['15.0-40.0'], df['38.0-82.0'], df['84.0-200.0'], df['110.0-900.0'], df['330.0-420.0'], df['420.0-510.0'], df['510.0-700.0'], df['700.0--1']]
        # energy_bins = ['4.2-8.7','8.7-14.5','15.0-44.0','39.0-82.0','84.0-200.0','110.0-900.0','330.0-420.0','420.0-510.0','510.0-700.0','700.0--1']
        # fluxes = [df['1.02-1.86'], df['1.9-2.3'], df['2.31-3.34'], df['3.4-6.48'], df['5.84-11.0'],df['11.64-23.27'], df['24.9-38.1'], df['40.3-73.4'], df['83.7-98.5'], df['99.9-118.0'], df['115.0-143.0'], df['160.0-242.0'], df['276.0-404.0'], df['500.0--1']]
        # energy_bins = ['1.02-1.86','1.9-2.3','2.31-3.34','3.4-6.48','5.84-11.0','11.64-23.27','24.9-38.1','40.3-73.4','83.7-98.5','99.9-118.0','115.0-143.0','160.0-242.0','276.0-404.0','500.0--1']
    else:
        flux_type = 'integral'
        fluxes = [df['5.0--1'], df['10.0--1'], df['30.0--1'], df['50.0--1'], df['60.0--1'], \
            df['100.0--1'], df['700.0--1']]
        energy_bins = [5.0, 10.0, 30.0, 50.0, 60.0, 100.0, 700.0]
    onset_times = from_fetchsep(dates, fluxes, energy_bins, flux_type, True)
    
    
    datetime_dates = []
    for j in range(len(dates)):
        datetime_dates.append(str_to_datetime(dates[j]))


    fig = plt.figure(figsize=(20,20))
    ax = plt.subplot(111)
    for i in range(len(fluxes)):
        ax.plot(datetime_dates, fluxes[i], label = energy_bins[i])
    ax.set_xlabel('Date')
    ax.set_yscale('log')
    ax.set_ylabel('Flux (pfu)')
    ax.xaxis.set_major_locator(mdates.DayLocator())
    minor_locator = mdates.HourLocator(byhour=range(0, 24, 1)) # Locates ticks at 00:00, 06:00, 12:00, 18:00
    ax.xaxis.set_minor_locator(minor_locator) 

    # # Set major tick formatter to display month and year
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    for o in range(len(onset_times)):
        ax.vlines(onset_times[o], ymin = np.min(fluxes[0]), ymax = np.max(fluxes[0]), color = 'black', linestyle = '--')
    savename = str(dates[0]).rsplit(' ')[0] + '_event_' + flux_type + '.png'
    plt.legend()
    fig.savefig(savename)
    mplt.show()
    # plt.close(fig)