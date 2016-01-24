#%pylab inline

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

    xdata = data[:,0]
    ydata = data[:,1]

    sel_xdata, sel_ydata = selectdomain(xdata, ydata, [0.00702, 0.00735])

    def sin_function(x, A, w, d, c):
        return A*np.sin(w*x + d) + c

    p = [.01, 950e2, np.pi/4, .04]

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
    plt.savefig("plots/%s" % dataset)
    plt.show()



if __name__ == '__main__':
    larmor_fit("data/larmor_y0.033A.tsv")
