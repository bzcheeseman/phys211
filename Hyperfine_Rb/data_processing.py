%pylab inline

import sys

sys.path.append("..")

from PythonLibrary.selectDomain import selectdomain
from PythonLibrary.residuals import lorentzian, multi_cosine, linear, gaussian
from PythonLibrary.estError import est_error

import numpy as np
import matplotlib.pyplot as plt
from scipy import loadtxt
from scipy.optimize import curve_fit
from scipy.signal import argrelmax
from scipy.ndimage.filters import gaussian_filter

def convert_time_freq():
    time, chan1 = np.genfromtxt("data/interferometer_cal.tsv", unpack=True, skip_header=1, usecols=(0,1))

    time = time*1e-3

    t, ch1 = selectdomain(time, chan1, [.06e-3, .16e-3])

    ##### Fitting consine #####

    # [A, B, omega1, omega2, delta1, delta2, C]
    p = [1.8e-2, 6.68e-2, 7.264e4, 4.419e5, 1.773, -4.52, 0]

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
            tempx, tempy = selectdomain(t, ch1, [number-.006e-3, number+.005e-3])

            tempy = gaussian_filter(tempy, 7)

            p = [.14, .0027e-3, tempx[argrelmax(tempy)[0]], -.05]

            # popt, pcov = curve_fit(gaussian, tempx, tempy)

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
        fs = []
        i = 0
        while i < nruns:
            t_pk = np.random.randn(len(t_peak)) * dt + t_peak
            delt = []
            for i in range(0, len(t_pk)-1):
                delt.append(t_pk[i+1] - t_pk[i])
            delt = np.ravel(delt)

            fp = []
            for i in range(0, len(delt)):
                fp.append((2.*3e8)/(delt[i]*(29.5e-2-10.5e-2)))
            fp = np.ravel(fp)

            #h = 6.62e-34*np.ones_like(f)

            #print len(f), len(deltat), len(dt)

            p_conv = [-1.569e15, 4.45e13]

            convert_popt, convert_pcov = curve_fit(linear, delt, fp, p0 = p_conv)

            #conv_yFit = linear(deltat, *convert_popt)
            #conv_yuFit = linear(deltat, *p_conv)

            deltat=delt
            fs=fp
            parms_A.append(convert_popt[0])
            parms_B.append(convert_popt[1])

            i += 1

        return parms_A, parms_B, deltat, fs

    parms_A, parms_B, delt, fs = find_conv(3)

    deltat = np.ravel(delt)
    f = np.ravel(fs)

    conv_popt = [np.average(parms_A), np.average(parms_B)]
    conv_yFit = linear(deltat, *conv_popt)

    redchi = np.sum(np.absolute(f - conv_yFit)/conv_yFit)/len(f)

    plt.figure(figsize = (10, 10))
    plt.errorbar(deltat, f, fmt='o', label = "Data")
    plt.plot(deltat, conv_yFit, 'r', label = "Fitted Curve")
    plt.text(.01365e-3,f[0],"f = kt + C\nk = %.2e $\pm$ %.2f \nC = %.2e $\pm$ %.2f" % (conv_popt[0], np.std(parms_A), conv_popt[1], np.std(parms_B)))
    plt.text(.01365e-3,f[0]-.0075e14, r"$\tilde{\chi}^2$ = %.2e" % redchi)
    plt.xlabel("Time (s)")
    plt.ylabel("Frequency (Hz)")
    plt.title("Linear Fit to Find Conversion from Time to Frequency")
    #plt.show()

    plt.figure(figsize = (10, 10))
    plt.plot(t, ch1, 'o', label = "Raw Data")
    plt.plot(t_peak, ch1_peak, 'ro', ms=15, label = "Peaks")
    plt.plot(t, yFit, 'r', label = "Fitted Curve")
    plt.plot(t, yuFit, 'g', label = "Guesses")
    plt.xlabel("Time (s)")
    plt.ylabel("Absorbtion")
    plt.title("Finding Conversion Factors")
    plt.legend()
    #plt.show()

def plot_data(dataset):   #need to update this code to be frequency on x axis, also need to double-check converstion, might be off by about 3 orders of magnitude
    t, ch1, ch2 = loadtxt(dataset, unpack = True, skiprows=1, usecols=(0,1,3))

    ti, c1 = selectdomain(t, ch1, [.122, .15])

    c1_filter = gaussian_filter(c1, 10)

    x0 = ti[argrelmax(c1_filter)]

    p = [1e-3, .02, x0, -.02]

    popt, pcov = curve_fit(lorentzian, ti, c1, p0 = p)

    print popt

    yFit = lorentzian(ti, *popt)
    yuFit = lorentzian(ti, *p)

    ch1err = est_error(ch1, 1e-3)
    c1err = est_error(c1, 1e-3)

    plt.figure(figsize = (10, 10))
    plt.errorbar(t, ch1, yerr=ch1err, fmt='bo')
    plt.errorbar(ti, c1, c1err, fmt='mo')
    plt.plot(ti, yFit, 'r')
    plt.plot(ti, yuFit, 'g')






if __name__ == '__main__':
    #plot_data("data/doppler_broadened.tsv")
    convert_time_freq()
