#%pylab inline

import sys

sys.path.append("..")

from PythonLibrary.selectDomain import selectdomain
from PythonLibrary.residuals import *
from PythonLibrary.estError import est_error

import numpy as np
import matplotlib.pyplot as plt
from scipy import loadtxt
from scipy.optimize import curve_fit
from scipy.signal import argrelmax, argrelmin
from scipy.ndimage.filters import gaussian_filter

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def convert_time_freq():
    time, chan1 = np.genfromtxt("data/interferometer_cal.tsv", unpack=True, skip_header=1, usecols=(0,1))

    time = time

    t, ch1 = selectdomain(time, chan1, [.06, .16])

    ##### Fitting consine #####

    # [A, B, omega1, omega2, delta1, delta2, C]
    p = [1.8e-2, 6.68e-2, 7.264, 4.419e2, 1.773, -4.52, 0]

    popt, pcov = curve_fit(multi_cosine, t, ch1, p0 = p)

    yFit = multi_cosine(t, *popt)
    yuFit = multi_cosine(t, *p)

    # plt.plot(t, ch1)
    # plt.plot(t, yFit)

    ##### Getting exact peak locations #####
    def find_peaks():         # implement multiple runs?
        t_peak = []
        ch1_peak = []

        dt = []

        for number in t[argrelmax(yFit)]:
            tempx, tempy = selectdomain(t, ch1, [number-.007, number+.005])

            tempy = gaussian_filter(tempy, 7)

            p = [.14, .0027, tempx[argrelmax(tempy)[0]], -.05]

            # popt, pcov = curve_fit(gaussian, tempx, tempy)
            #
            # tempyFit = gaussian(tempx, *popt)
            # tempyuFit = gaussian(tempx, *p)
            #
            # plt.plot(tempx, tempy)
            # plt.plot(tempx, tempyFit, 'r')
            # plt.plot(tempx, tempyuFit, 'g')

            t_peak.append(tempx[argrelmax(tempy)[0]])
            ch1_peak.append(tempy[argrelmax(tempy)[0]])

            dt.append(p[1]/len(tempy))

        return t_peak, ch1_peak, dt

    t_peak, ch1_peak, dt = find_peaks()

    t_peak = np.ravel(t_peak)

    dt = np.absolute(dt)

    ##### Now we find the conversion from time to frequency #####

    def find_conv(nruns):
        parms_A = []
        parms_B = []
        deltat = []

        for j in range(0, nruns):
            t_pk = np.random.randn(len(t_peak)) * dt + t_peak

            l = len(t_pk)

            delt = []
            for i in range(0, l-1):
                delt.append(t_pk[i+1] - t_pk[i])
            delt = np.ravel(delt)

            ld = len(delt)

            deltat.append(delt)

        return deltat

    delt = find_conv(5000)
    dt = np.ravel(delt)

    df = 3e8/(2*(29.5e-2 - 10.5e-2))

    deltf = np.sqrt((.5e-2/29.5e-2)**2 + (.5e-2/10.5e-2)**2) * df

    f = np.ones_like(dt) * df

    # plt.figure(figsize = (10, 10))
    # plt.errorbar(deltat, f, fmt='o', label = "Data")
    # plt.xlabel("Change in Time (d(s))")
    # plt.ylabel("Change in Frequency (d(Hz))")
    # plt.title("Linear Fit to Find Conversion from Time to Frequency")
    # plt.savefig("plots/linear_conversion.png")
    # plt.show()

    # print np.average(deltat), np.std(deltat), df, deltf

    # interferometer_text = "Fit Form:\n \
    #                       $A\cos(\omega_1 t + \delta_1) + B\cos(\omega_2 t + \delta_2)$"
    #
    #
    # plt.figure(figsize = (10, 10))
    # plt.plot(t, ch1, 'o', label = "Raw Data")
    # plt.plot(t_peak, ch1_peak, 'ro', ms=15, label = "Peaks")
    # plt.plot(t, yFit, 'r', label = "Fitted Curve")
    #
    # plt.text(.08, .11, interferometer_text)
    # #plt.plot(t, yuFit, 'g', label = "Guesses")
    # plt.xlabel("Time (s)")
    # plt.ylabel("Absorbtion $W/(m^2)$")
    # plt.title("Finding Conversion Factors")
    # plt.legend()
    # plt.savefig("plots/conversion.pdf")
    # plt.show()

    return {"Time value":np.average(dt), "Time Uncertainty": np.std(dt), "Frequency value": df, "Frequency Uncertainty":deltf}

def plot_broadened_data(dataset):   #need to update this code to be frequency on x axis, also need to double-check converstion, might be off by about 3 orders of magnitude

    conversions = convert_time_freq()

    t, ch1, ch2 = loadtxt(dataset, unpack = True, skiprows=1, usecols=(0,1,3))

    k = conversions["Frequency Uncertainty"]/conversions["Time value"]

    ti, c1 = selectdomain(t, ch1, [.0825, .105])

    c1_filter = gaussian_filter(c1, 25)

    x0 = ti[argrelmax(c1_filter)]

    #print x0

    p = [4e-3, 1.5e-2, x0, -.07]

    popt, pcov = curve_fit(lorentzian, ti, c1, p0 = p)

    yFit = lorentzian(ti, *popt)
    yuFit = lorentzian(ti, *p)

    ch1err = est_error(ch1, 1e-3)
    c1err = est_error(c1, 1e-3)

    #print "Linewidth: %.2e Hz" % (k*popt[1])
    #print k

    f2 = "87 Rb F = 2"
    f3 = "85 Rb F = 3"
    f285 = "85 Rb F = 2"
    f1 = "87 Rb F = 1"

    plt.figure(figsize = (8, 8))
    plt.errorbar(t, ch1, yerr=ch1err, fmt='bo')
    plt.errorbar(ti, c1, c1err, fmt='mo')
    plt.plot(ti, yFit, 'r', label="Lorentzian Fit")
    #plt.plot(ti, yuFit, 'g')

    plt.annotate(f2, (.07, .14), (.06, .2), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate(f3, (.095, 0), (.085, -.2), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate(f285, (.12, .3), (.12, .4), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate(f1, (.19, .1), (.17, .2), arrowprops = dict(width=2, headwidth=4, facecolor="red"))

    plt.text(.05, -.25, "85 Rb F = 3 Linewidth: %.2e $\pm$ %.2e Hz" % (k*popt[1], np.sqrt((conversions["Frequency Uncertainty"]/conversions["Frequency value"])**2 + (conversions["Time Uncertainty"]/conversions["Time value"])**2) * k*np.sqrt(pcov[1,1])))

    plt.xlabel("Time (s)")
    plt.ylabel("Absorbtion $W/(m^2)$")
    plt.title("Determining the linewidth of the doppler-broadened peak")
    plt.savefig("plots/rb85_fwhm.pdf")
    plt.show()

    # need to point to each peak and point out which one it is

def plot_hyperfine(dataset):
    conversions = convert_time_freq()

    t, ch1, ch2 = loadtxt(dataset, unpack = True, skiprows=1, usecols=(0,1,3))

    k = conversions["Frequency value"]/conversions["Time value"]
    dk = np.sqrt((conversions["Time Uncertainty"]/conversions["Time value"])**2 + (conversions["Frequency Uncertainty"]/conversions["Frequency value"])**2)

    ti, c1 = selectdomain(t, ch1, [.0825, .105])

    x = [.08625, .0879, .0893, .0906, .0922, .0950]

    ts = []
    cs = []

    for i in range(0, len(x)):
        n = x[i]

        temp_ti, temp_c1 = selectdomain(ti, c1, [n-.0004, n+.0004])

        c1_filter = gaussian_filter(temp_c1, 4)

        c1_filter_minus_lin = c1_filter - ((c1_filter[-1] - c1_filter[0])/(temp_ti[-1] - temp_ti[0]) * temp_ti)

        # plt.plot(temp_ti, c1_filter)
        # plt.plot(temp_ti[argrelmin(c1_filter_minus_lin)], temp_c1[argrelmin(c1_filter_minus_lin)], 'o')

        ts.append(temp_ti[argrelmin(c1_filter_minus_lin)])
        cs.append(temp_c1[argrelmin(c1_filter_minus_lin)])

    diff = []
    ddiff = []
    for i in range(0, len(ts) - 1):
        temp = (ts[i+1] - ts[i]) * k * 1e-6
        diff.append(temp)
        ddiff.append((.0002 * k * 1e-6)/temp)

    dk = np.sqrt(dk**2 + np.average(ddiff)**2)

    print diff[0], dk * diff[0]
    print diff[1] + diff[2], dk * (diff[1] + diff[2])
    print diff[-1] + diff[-2], dk * (diff[-1] + diff[-2])

    ch1err = est_error(ch1, 1e-3)
    c1err = est_error(c1, 1e-3)



    plt.figure(figsize=(8,8))
    plt.errorbar(ti, c1, c1err, fmt='o')
    plt.plot(ts, cs, 'ro', label = "Peak Locations")
    #plt.plot(ti, yFit, 'r', label="Lorentzian Fit")
    #plt.plot(ti, yuFit, 'g')

    plt.annotate("Hyperfine, f = 0,1", (ts[0], cs[0]), (.085, -.05), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate("Crossover", (ts[1], cs[1]), (.083, .05), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate("Hyperfine, f = 1,2", (ts[2], cs[2]), (.091, -.05), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate("Crossover", (ts[3], cs[3]), (.083, .06), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate("Crossover", (ts[4], cs[4]), (.083, .07), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
    plt.annotate("Hyperfine, f = 2,3", (ts[5], cs[5]), (.097, -.05), arrowprops = dict(width=2, headwidth=4, facecolor="red"))

    plt.xlabel("Time (s)")
    plt.ylabel("Absorbtion $W/(m^2)$")
    plt.title("Determining the positions of the Hyperfine peaks")
    plt.legend()
    plt.savefig("plots/hyperfine.pdf")
    plt.show()

def FWHM_hyperfine(dataset):
    conversions = convert_time_freq()

    t, ch1, ch2 = loadtxt(dataset, unpack = True, skiprows=1, usecols=(0,1,3))

    k = conversions["Frequency value"]/conversions["Time value"]
    dk = np.sqrt((conversions["Time Uncertainty"]/conversions["Time value"])**2 + (conversions["Frequency Uncertainty"]/conversions["Frequency value"])**2)
    ti, c1 = selectdomain(t, ch1, [.09, .091])

    c1_filter = gaussian_filter(c1, 5)

    x0 = ti[argrelmin(c1_filter)]

    print x0

    p = [-9e-6, 2.35e-4, x0, 1e-2, 1.18e-1]

    popt, pcov = curve_fit(lorentzian_back, ti, c1, p0 = p)

    lw = 1.054e-34/(6.626e-34 * popt[1] * k)
    dlw = np.sqrt((np.sqrt(pcov[1,1])/popt[1])**2 + dk**2) * 1.054e-34/(6.626e-34 * popt[1] * k)

    print "Linewidth of Hyperfine State: %.2e +- %.2e" % (lw, dlw)

    yFit = lorentzian_back(ti, *popt)
    yuFit = lorentzian_back(ti, *p)

    ch1err = est_error(ch1, 1e-3)
    c1err = est_error(c1, 1e-3)

    hyperfine_text1 = r"$f(t) = \frac{A}{\pi}\frac{\frac{1}{2}\Gamma}{(t - t_0)^2 + (\frac{1}{2}\Gamma)^2} + Bt + C$"

    hyperfine_text2 = r"$f(t) = \frac{%.2e}{\pi}\frac{\frac{1}{2}%.2e}{(t - %.2e)^2  + (\frac{1}{2}%.2e)^2} + %.2et + %.2e$" % (popt[0], popt[1], popt[2], popt[1], popt[3], popt[4])

    hyperfine_text3 = "Linewidth of Hyperfine State: $%.2e\,\pm\,%.2e$" % (lw, dlw)

    plt.figure(figsize = (10, 10))
    plt.errorbar(ti, c1, c1err, fmt='o')
    plt.plot(ti, yFit, 'r', lw=2, alpha=.8, label="Lorentzian Fit")
    plt.text(.0899, .06, "%s\n%s\n%s" % (hyperfine_text1, hyperfine_text2, hyperfine_text3))
    plt.xlabel("Time (s)")
    plt.ylabel("Absorbtion $W/(m^2)$")
    plt.title("Determining the linewidth of a Hyperfine Feature peak")
    plt.legend()
    plt.savefig("plots/rb87_hyperfine.pdf")
    plt.show()



if __name__ == '__main__':
    #plot_broadened_data("data/doppler_free_all.tsv")
    #convert_time_freq()
    #plot_hyperfine("data/doppler_free_pk2.tsv")
    FWHM_hyperfine("data/doppler_free_pk2.tsv")
