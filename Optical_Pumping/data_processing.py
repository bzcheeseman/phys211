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

    with open('fitparms_{}.{}.tsv'.format(dset[0].split("/")[1], dset[1]), 'w')\
      as f:
        for key in sorted(parms.keys()):
            f.write("{}\t{}\n".format(key, parms[key]))

def gyromagnetic_ratio_fit():
    data20 = np.genfromtxt("fitparms_larmor_y0.020A.tsv", usecols=(1))
    data33 = np.genfromtxt("fitparms_larmor_y0.033A.tsv", usecols=(1))
    data40 = np.genfromtxt("fitparms_larmor_y0.040A.tsv", usecols=(1))
    data45 = np.genfromtxt("fitparms_larmor_y0.045A.tsv", usecols=(1))
    data55 = np.genfromtxt("fitparms_larmor_y0.055A.tsv", usecols=(1))
    data65 = np.genfromtxt("fitparms_larmor_y0.065A.tsv", usecols=(1))

    Is = [.02, .033, .04, .045, .055, .065]

    freqs = [data20[1], data33[1], data40[1], data45[1], data55[1], data65[1]]

    def lin(x, A, B):
        return A*x + B

    popt, pcov = curve_fit(lin, Is, freqs, p0 = np.polyfit(Is, freqs, 1))

    yFit = np.poly1d(popt)

    print yFit

    plt.plot(Is, freqs, 'o')
    plt.plot(Is, yFit(Is), 'r')


if __name__ == '__main__':
    gyromagnetic_ratio_fit()
