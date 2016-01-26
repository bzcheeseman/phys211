%pylab inline

from matplotlib import rc
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def selectdomain(xdata,ydata,domain):
    ind=np.searchsorted(xdata,domain)
    return xdata[ind[0]:ind[1]],ydata[ind[0]:ind[1]]

def larmor_fit(dataset):
    data = np.genfromtxt(dataset, skip_header = 1, usecols = (0,1))

    dset = dataset.split(".")

    xdata = data[:,0]
    ydata = data[:,1]

    sel_xdata, sel_ydata = selectdomain(xdata, ydata, [0.00702, 0.00735])

    def sin_function(x, A, w, d, c):
        return A*np.sin(w*x + d) + c

    p = [.01, 150e3, np.pi/2, .035]

    popt, pcov = curve_fit(sin_function, sel_xdata, sel_ydata, p0 = p)

    yFit = sin_function(sel_xdata, *popt)
    yuFit = sin_function(sel_xdata, *p)

    parms = {"A":popt[0], "freq":popt[1], "shift":popt[2], "vertical":popt[3]}

    print parms["freq"]/(2*np.pi), "Hz"

    plt.figure(figsize = (10,8))
    plt.plot(sel_xdata, sel_ydata, 'o')
    plt.plot(sel_xdata, yuFit, 'g')
    plt.plot(sel_xdata, yFit, 'r', lw = 4, alpha = .7)
    plt.xlabel("time (s)")
    plt.ylabel("signal amplitude (-mV)")
    plt.savefig("plots/({}.{}).png".format(dset[0].split("/")[1], dset[1]))

    with open('data/fitparms_{}.{}.tsv'.format(dset[0].split("/")[1], dset[1]), 'w')\
      as f:
        for key in sorted(parms.keys()):
            f.write("{}\t{}\n".format(key, parms[key]))

def gyromagnetic_ratio_fit():
    data20 = np.genfromtxt("data/fitparms_larmor_y0.020A.tsv", usecols=(1))
    data33 = np.genfromtxt("data/fitparms_larmor_y0.033A.tsv", usecols=(1))
    data40 = np.genfromtxt("data/fitparms_larmor_y0.040A.tsv", usecols=(1))
    data45 = np.genfromtxt("data/fitparms_larmor_y0.045A.tsv", usecols=(1))
    data55 = np.genfromtxt("data/fitparms_larmor_y0.055A.tsv", usecols=(1))
    data65 = np.genfromtxt("data/fitparms_larmor_y0.065A.tsv", usecols=(1))

    Is = [.02, .033, .04, .045, .055, .065]
    dI = .001 * np.ones_like(Is)

    freqs = [data20[1], data33[1], data40[1], data45[1], data55[1], data65[1]]\
    # propagate df through previous section - larmor fit

    def lin(x, A, B):
        return A*x + B

    popt, pcov = curve_fit(lin, Is, freqs, p0 = np.polyfit(Is, freqs, 1))
    yFit = np.poly1d(popt)

    print yFit

    plt.errorbar(Is, freqs, xerr = dI, fmt = 'o')
    plt.plot(Is, yFit(Is), 'r')
    plt.xlabel("Current (A)")
    plt.ylabel("Frequency (Hz)")
    plt.savefig("plots/gyromagnetic_ratio_fit.png")

def earth_field_fit():
    data = np.genfromtxt("data/earth_B_field_data.csv")

    xdata = data[:,0]
    ydata = data[:,1]

    def fitform(x, Bhe, Bhapp, Bve, C):
        return 2.895e-9*np.sqrt((Bhe - Bhapp)**2 + (Bve - x)**2) + C

    p = [2e6, 2e6, 2e6, 300]

    popt, pcov = curve_fit(fitform, xdata, ydata, p0 = p)

    yFit = fitform(xdata, *popt)
    yuFit = fitform(xdata, *p)

    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, yuFit, 'g')
    plt.plot(xdata, yFit, 'r')
    plt.ylabel("Frequency (kHz)")
    plt.xlabel("")
    plt.show()


if __name__ == '__main__':

    earth_field_fit()
