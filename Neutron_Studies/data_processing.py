#%pylab inline

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import sys

sys.path.append("..")

from PythonLibrary.residuals import *
from PythonLibrary.selectDomain import selectdomain
from PythonLibrary.estError import est_error

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def calibration(plot_cal):
    try:
        data = np.genfromtxt("data/run_3/na_22_cal_2.tsv", skip_header=26, usecols=(0,2))
    except Exception:
        data = np.genfromtxt("data/run_3/na_22_cal_2.tsv", skip_header=26)
    data_err = np.sqrt(data[:,1])

    peak_1_x, peak_1_y = selectdomain(data[:,0], data[:,1], [120, 180])

    p1 = [1e6, 50, 170, -10, 200]

    popt1, pcov1 = curve_fit(gaussian_back, peak_1_x, peak_1_y, p0 = p1)

    yFit1 = gaussian_back(peak_1_x, *popt1)
    yuFit1 = gaussian_back(peak_1_x, *p1)

    #print popt1[2]

    peak_2_x, peak_2_y = selectdomain(data[:,0], data[:,1], [330, 430])

    p2 = [1e6, 100, 350, 0]

    popt2, pcov2 = curve_fit(gaussian, peak_2_x, peak_2_y, p0 = p2)

    yFit2 = gaussian(peak_2_x, *popt2)
    yuFit2 = gaussian(peak_2_x, *p2)

    #print popt2[2]

    data_cs = np.genfromtxt("data/run_3/cs_137_cal_2.tsv", skip_header=26, usecols=(0,2))
    cs_err = np.sqrt(data_cs[:,1])

    peak_3_x, peak_3_y = selectdomain(data_cs[:,0], data_cs[:,1], [160, 230])

    p3 = [1e3, 50, 200, -10, 200]

    popt3, pcov3 = curve_fit(gaussian_back, peak_3_x, peak_3_y, p0 = p3)

    yFit3 = gaussian_back(peak_3_x, *popt3)
    yuFit3 = gaussian_back(peak_3_x, *p3)

    #print popt3[2]

    if plot_cal:
        plt.figure(figsize=(10,10))

        # plt.annotate("Na-22 511 keV", (popt1[2], 1.7e4), (popt1[2]+50, 1.7e4), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
        # plt.annotate("Cs-137 662 keV", (popt3[2], .98e4),  (popt3[2], 1.3e4), arrowprops = dict(width=2, headwidth=4, facecolor="red"))
        # plt.annotate("Na-22 1275 keV", (popt2[2], 3e3), (popt2[2], 4e3), arrowprops = dict(width=2, headwidth=4, facecolor="red"))

        plt.errorbar(data[:,0], data[:,1], data_err)
        # plt.plot(peak_1_x, yuFit1)
        plt.plot(peak_1_x, yFit1, alpha=.8, lw=3, label="Center: %.0f" % popt1[2])
        # plt.plot(peak_2_x, yuFit2)
        plt.plot(peak_2_x, yFit2, alpha=.8, lw=3, label="Center: %.0f" % popt2[2])
        #plt.figure(figsize=(10, 10))
        plt.errorbar(data_cs[:,0], data_cs[:,1], cs_err)
        # plt.plot(peak_3_x, yuFit3)
        plt.plot(peak_3_x, yFit3, alpha=.8, lw=3, label="Center: %.0f" % popt3[2])
        plt.xlabel("Channel")
        plt.ylabel("Counts")
        plt.title("Energy Calibration Using Known Sources")
        plt.legend()
        plt.savefig("plots/calibration_peaks.pdf")

    cal_line_x = np.array([popt1[2], popt3[2], popt2[2]])
    cal_line_y = np.array([511, 662, 1275])
    x_err = np.array([np.sqrt(pcov1[2,2]), np.sqrt(pcov3[2,2]), np.sqrt(pcov2[2,2])])

    p_lin = [2.0, 150.0]

    lin, lin_pcov = curve_fit(linear, cal_line_x, cal_line_y, p0 = p_lin)

    print lin

    # yFit = linear(cal_line_x, *lin)
    # yuFit = linear(cal_line_x, *p_lin)
    #
    # plt.figure(figsize=(10, 10))
    # plt.errorbar(cal_line_x, cal_line_y, x_err, fmt='o')
    # plt.plot(cal_line_x, yFit, alpha = .7, lw = np.sqrt(lin_pcov[0,0]), label="Slope Uncertainty: %.3f" % np.sqrt(lin_pcov[0,0]))
    # plt.title("%.3fx + %.3f" % (lin[0], lin[1]))
    # plt.legend()
    # plt.savefig("plots/tentative_cal.pdf")

    return lin
    # save figure, etc.

def spectrum(dataset, plot_full = False):
    conversion = calibration(False)

    try:
        data = np.genfromtxt(dataset, skip_header=26, usecols=(0,2))
    except ValueError:
        data = np.genfromtxt(dataset, skip_header=26, usecols=(0,1))

    data_x = linear(data[:,0], *conversion)

    domain = [2000, 2450]

    peak_x, peak_y = selectdomain(data_x, data[:,1], domain)

    back_x, back_y = selectdomain(data_x, data[:,1], [800, 3500], domain)
    back_x_full, back_y_full = selectdomain(data_x, data[:,1], [800, 3500])

    p_back = np.array([1e4, -1e-3, 6e2])

    back_popt, back_pcov = curve_fit(exponential, back_x, back_y, p0 = p_back, maxfev=int(1e4))

    back_yFit = exponential(back_x_full, *back_popt)
    back_yuFit = exponential(back_x_full, *p_back)

    to_subtract_x, to_subtract_y = selectdomain(back_x_full, back_yFit, domain)

    if plot_full:
        plt.figure(figsize=(10, 10))
        plt.errorbar(data_x, data[:,1], np.sqrt(data[:,1]), fmt='o', ms=1)
        plt.errorbar(peak_x, peak_y, np.sqrt(peak_y), fmt='o', ms=1, label = "Region of Interest")
        plt.plot(back_x_full, back_yFit, label = "Background Fit")
        plt.ylabel("Counts")
        plt.xlabel("Energy (keV)")
        plt.title("Isolating the Peak")
        plt.legend()
        plt.savefig("plots/peak_isolation_%s.pdf" % dataset.split("/")[2].split(".")[0])

    flat_peak = peak_y - to_subtract_y

    peak_p = [450, 18, 2200, 11]

    peak_popt, peak_pcov = curve_fit(gaussian, peak_x, flat_peak, p0 = peak_p)
    peak_yFit = gaussian(peak_x, *peak_popt)

    print peak_popt

    with open(dataset) as f:
        f.seek(298)
        try:
            livetime = float(f.read(6))
        except Exception:
            f.seek(404)
            livetime = float(f.read(6))
        f.close()
    print livetime

    plt.figure(figsize=(10, 10))
    plt.errorbar(peak_x, flat_peak, np.sqrt(flat_peak), fmt='o')
    plt.plot(peak_x, peak_yFit, label = "Gaussian Fit\nCenter: %.0f $\pm$ %.0f keV\nCountrate: %.2f counts/s" % (peak_popt[2], np.absolute(peak_popt[1]/np.sqrt(len(peak_x))), peak_popt[0]/livetime))
    plt.ylabel("Counts (after subtracting background)")
    plt.xlabel("Energy (keV)")
    plt.title("Finding the Peak Center")
    plt.legend()
    plt.savefig("plots/peak_center_%s.pdf" % dataset.split("/")[2].split(".")[0])

def cross_section(dataset):

    def calibrate():
        data = np.genfromtxt("data/cross_section/na_compton.tsv", skip_header=26)
        xdata, ydata = selectdomain(data[:,0], data[:,1], [100, 2048])
        plt.figure(figsize=(10, 10))
        plt.plot(xdata, ydata)
        plt.savefig("plots/calibration.pdf")

        return 500

    calibrate()
    data = np.genfromtxt("data/cross_section/cu/%s" % dataset, skip_header=26)

    xdata, ydata = selectdomain(data[:,0], data[:,1], [1500, 2048])

    with open("data/cross_section/cu/%s" % dataset) as f:
        f.seek(300)
        try:
            livetime = float(f.read(6))
        except Exception:
            f.seek(404)
            livetime = float(f.read(6))
        f.close()
    print livetime

    countrate = np.trapz(ydata, xdata)/livetime

    with open("data/cross_section/countrates.tsv", 'a+') as f:
        f.write(str(countrate))
        f.write("\t")
        f.write(dataset)
        f.write("\n")
        f.close()

    print np.trapz(ydata, xdata)/livetime

    plt.plot(xdata, ydata)


if __name__ == '__main__':
    # calibration(True)
    # datasets = ["direct_spec_port_closed.tsv", "direct_spec_port_open.tsv", "shielded_carbon_port_open.tsv", "shielded_paraffin_port_open.tsv", "shielded_spec_port_closed.tsv", "shielded_spec_port_open.tsv"]
    # for dataset in datasets:
    #     spectrum("data/%s" % dataset, True) #convert channel axis to energy
    spectrum("data/run_3/shielded_paraffin.tsv", True)
    # cross_section("cu_1_25cm.tsv")
