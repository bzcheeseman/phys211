import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
import pandas as pd

import sys
sys.path.append("..")

from PythonLibrary import *

# try:
#     %pylab inline
# except SyntaxError:
#     pass

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

class data_manage(object):
    """docstring for data_manage"""
    def __init__(self):
        self.data_folder = "data"

        self.primary_data = ["no_xray_filter.tsv", "with_xray_filter.tsv"]

        self.ss_301 = pd.read_table(self.data_folder + "/SS_302_001.tsv", delimiter = "\t", skiprows=24)
        self.velocity = pd.read_table(self.data_folder + "/velocity_iron_116passes.tsv", delimiter = "\t", skiprows=24)

    def calibrate(self):
        xdata, ydata = selectDomain.selectdomain(self.velocity["Channel"].values, self.velocity["Counts"].values, domain = [150, 450])

        yerr = np.sqrt(ydata)

        def abs(x, A, C):
            return A*np.absolute(x-C)

        p = [1000, 260]

        self.popt_counts, self.pcov_counts = curve_fit(abs, xdata, ydata, p0 = p)

    def cal(self, x):
        counts_chan = self.popt_counts[0] * (x - self.popt_counts[1])

        lamda = 632.8e-9
        t = 300e-6
        number_passes = 116

        v = (counts_chan * lamda)/(2 * t * number_passes)

        self.Ep = (1 + v/3e8)*14.4e3 # in eV

        return self.Ep

    '''
    The fit has a really hard time if we convert to energy first, so keep stuff in channels and then run it through the calibration.
    '''
    def SS_301(self):
        xdata = self.ss_301["Channel"].values
        ydata = self.ss_301["Counts"].values
        yerr = np.sqrt(ydata)

        def lorentzian(x, I0, Gamma, x0, C):
            return (-I0*(((.5*Gamma)**2)/((x-x0)**2+(.5*Gamma)**2))) + C

        #p = [5000, 421, 5e-7 + 1.44e4, 6000]

        popt, pcov = curve_fit(lorentzian, xdata, ydata)
        print popt

        yFit = lorentzian(xdata, *popt)
        #yuFit = lorentzian(xdata, *p)

        plt.figure(figsize = (10, 10))
        plt.errorbar(xdata, ydata, yerr = yerr, fmt = 'o', ms = 1)
        #plt.plot(xdata, yuFit)
        plt.plot(xdata, yFit, 'r')



if __name__ == '__main__':
    obj = data_manage()
    obj.calibrate()
    obj.SS_301()
