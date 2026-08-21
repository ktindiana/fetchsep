import os
import argparse
import warnings
import multiprocessing
from functools import partial
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import trapezoid
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def default_energy_bins():
    """Returns the default 100 target energy bins from GOES4DRM."""
    return np.array([
        1.000000E-02, 1.338088E-02, 1.790479E-02, 2.390681E-02, 3.183221E-02,
        4.210285E-02, 5.510735E-02, 7.111503E-02, 9.027458E-02, 1.124807E-01,
        1.375337E-01, 1.657185E-01, 1.968371E-01, 2.302813E-01, 2.675015E-01,
        3.081838E-01, 3.525475E-01, 4.010248E-01, 4.541620E-01, 5.125794E-01,
        5.769560E-01, 6.480289E-01, 7.263492E-01, 8.129234E-01, 9.085720E-01,
        1.014071E+00, 1.129672E+00, 1.258328E+00, 1.399567E+00, 1.556046E+00,
        1.728746E+00, 1.919146E+00, 2.129073E+00, 2.361070E+00, 2.617216E+00,
        2.899599E+00, 3.210916E+00, 3.554511E+00, 3.933012E+00, 4.350185E+00,
        4.810052E+00, 5.316827E+00, 5.875216E+00, 6.490437E+00, 7.167928E+00,
        7.914280E+00, 8.736185E+00, 9.641356E+00, 1.063789E+01, 1.173515E+01,
        1.294341E+01, 1.427364E+01, 1.573864E+01, 1.735183E+01, 1.912722E+01,
        2.108195E+01, 2.323460E+01, 2.560620E+01, 2.821697E+01, 3.109330E+01,
        3.426162E+01, 3.775280E+01, 4.160081E+01, 4.584282E+01, 5.051828E+01,
        5.567888E+01, 6.136987E+01, 6.765432E+01, 7.459531E+01, 8.226303E+01,
        9.074045E+01, 1.001212E+02, 1.105347E+02, 1.220314E+02, 1.348146E+02,
        1.489691E+02, 1.648166E+02, 1.823669E+02, 2.018336E+02, 2.239260E+02,
        2.484365E+02, 2.756300E+02, 3.060040E+02, 3.407294E+02, 3.793953E+02,
        4.231804E+02, 4.728099E+02, 5.290970E+02, 5.930393E+02, 6.664974E+02,
        7.505074E+02, 8.471435E+02, 9.587661E+02, 1.090824E+03, 1.243621E+03,
        1.417821E+03, 1.624766E+03, 1.869156E+03, 2.158058E+03, 2.500000E+03
    ])

def weibull(E, J0, k, a):
    return [J0 * k * a * En**(a - 1) * np.exp(-((k * En)**a)) for En in E]

def loginterp(x, y, xn, reduced=0):
    if reduced == 0:
        fl = interp1d(np.log(x), np.log(y), fill_value=(-50, -50), kind='slinear', bounds_error=False)
    else:
        fl = interp1d(np.log(x), np.log(y), fill_value='extrapolate', kind='linear', bounds_error=False)
    
    xn_log = np.log(xn)
    yn_log = fl(xn_log)
    return np.exp(yn_log)

def weibull_fit(x, y, xn):
    Z = (y[0], 1, 1)
    
    def err(params, E):
        J0, k, a = params
        e = 0
        y4 = weibull(E, J0, k, a)
        if any([np.isnan(i) for i in y4]):
            raise ValueError
        for i, _ in enumerate(y):
            e += (np.log(max(y[i], 1e-12)) - np.log(max(y4[i], 1e-12)))**2
            if i == 0 or i == len(y) - 1:
                e = e * 10
        return e

    try:
        m = minimize(err, Z, args=(x,), method='Nelder-Mead', tol=1E-12, options={"maxiter": 1000})
        return weibull(xn, *m.x)
    except Exception:
        return loginterp(x, y, xn, reduced=1)

def data_processor(y, x, xn, extrapolate=True):
    y1 = [i for i in y if i > 0]
    x1 = [i for n, i in enumerate(x) if y[n] > 0]
    
    if len(y1) == 0:
        return np.zeros_like(xn)

    if extrapolate:
        xext = xn[:-6]
        yext = weibull_fit(x1, y1, xext)
        y1 = [i for x_val, i in zip(xext, yext) if x_val < min(x1)] + y1 + [i for x_val, i in zip(xext, yext) if x_val > max(x1)]
        x1 = [i for i in xext if i < min(x1)] + x1 + [i for i in xext if i > max(x1)]
        
    return loginterp(x1, y1, xn)

def parse_header_energies(flux_cols):  
    """Current process for creating edata for interpolation"""
    edata = []
    for flux in flux_cols:
        avg = 0
        if "--" in flux:
            i = flux.index("-")
            f1 = float(flux[:i])
            avg = f1 + 100
        else:
            i = flux.index("-")
            f1 = float(flux[:i])
            f2 = float(flux[i+1:])
            avg = round((f1 + f2) / 2, 2)
        edata.append(avg)
    return np.array(edata)

def setup_eng_plot(ytype ="Differential Flux (cm$^{-2}$min$^{-1}$MeV$^{-1}$)"):
    plt.figure(dpi=150)
    plt.yscale("log");plt.xscale("log")
    plt.xlim(1,2500)
    plt.ylabel(ytype)
    plt.xlabel("Energy (MeV)")
    plt.grid()


def process_file(input_file, ebins, extrapolate, cores, plots=True):
    """Processes the fetchsep input file and returns the interpolated version of that data"""
    df = pd.read_csv(input_file)
    
    time_col = df.columns[0]
    flux_cols = df.columns[1:]

    # Parse dates for time-series plots
    dtime = pd.to_datetime(df[time_col]).dt.to_pydatetime()
    
    # Calculate time step interval in minutes from input data
    if len(dtime) > 1:
        tstep = (dtime[1] - dtime[0]).total_seconds() / 60.0
    else:
        tstep = 5.0
    
    dflux = df[flux_cols].values
    edata = parse_header_energies(flux_cols)
    
    
    with multiprocessing.Pool(cores) as P:
        interparray = P.map(partial(data_processor,x = edata,xn = ebins,extrapolate = extrapolate),dflux)
    interparray = np.array(interparray)*4*np.pi*60
    interparray[interparray < 1e-18] = 0

    # Grab the .csv file name
    base_name = os.path.splitext(input_file)[0]
    
    # Export full matrix matching standard 100-bin output format
    output_csv = f"{base_name}_interpolated.csv"
    out_df = pd.DataFrame(interparray, columns=[f"{e:.18e}" for e in ebins])
    out_df.insert(0, time_col, df[time_col])
    out_df.to_csv(output_csv, index=False)

    # Run each figure for the data
    if plots:
        integral_array = []
        ebins2 = [10, 30, 50, 60, 100, 500, 1000]
        
        for i in range(len(dtime)):
            integral = [trapezoid(interparray[i][n:], x=ebins[n:]) for n, _ in enumerate(ebins)]        
            integral2 = loginterp(ebins, integral, ebins2)
            integral_array.append(integral2)
        
        setup_eng_plot()
        spectrum = np.sum(interparray, axis=0) * tstep
        plt.step(ebins, spectrum)
        plt.savefig(f"{base_name}_sum-diff-fluence.png")
        plt.close()

        setup_eng_plot()
        for n, i in enumerate(interparray):
            if n == 0: 
                plt.step(ebins, i, "r-", lw=1, label=f"{input_file}, {int(tstep)}min")
            else:
                plt.step(ebins, i, lw=1)
        plt.legend()
        plt.savefig(f"{base_name}_differential-flux.png")
        plt.close()

        setup_eng_plot(ytype="Integral Flux (cm$^{-2}$min$^{-1}$)")
        for n, i in enumerate(integral_array):
            if n == 0: 
                plt.plot(ebins2, i, "r-", lw=1, label=f"{input_file}, {int(tstep)}min")
            else:
                plt.plot(ebins2, i, lw=1)
        plt.legend()
        plt.savefig(f"{base_name}_integral-flux.png")
        plt.close()

    print(f"Successfully processed {input_file} -> {output_csv}")

warnings.filterwarnings("ignore")
matplotlib.use('tkagg')
if __name__ == "__main__":
    # Example input: python fetchsep_post_processor.py GOES08July2000.csv 
    parser = argparse.ArgumentParser(description="Fetchsep post-processor for spectrum interpolation.")
    parser.add_argument("input_file", help="Path to fetchsep CSV output file.")
    parser.add_argument("-extrapolate", action="store_true", help="Enable Weibull extrapolation.")
    parser.add_argument("-ebin_file", type=str, default=None, help="File containing energy bins.")

    cores = multiprocessing.cpu_count()-1

    args = parser.parse_args()

    # Process ebin files
    if args.ebin_file and os.path.exists(args.ebin_file):
        ebins = np.loadtxt(args.ebin_file)
    else:
        ebins = default_energy_bins()

    # Pass arguments to processor for interpolated
    process_file(
        input_file=args.input_file,
        ebins=ebins,
        extrapolate=args.extrapolate,
        cores=cores
    )