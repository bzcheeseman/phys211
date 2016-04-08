import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
import scipy.signal as signal

import sys
sys.path.append("/users/aman/desktop/PHYS211")

from PythonLibrary import *

try:
    %pylab inline
except SyntaxError:
    pass

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


def compton(E, theta, m0):
    return E/(1+(E/(m0 * 3e8**2))(1-np.cos(theta)))

class data_manage(object):
    """docstring for data_manage"""
    def __init__(self, dataset):
        try:
            data = np.genfromtxt(dataset, skip_header=25)

            self.xdata = data[:,0]
            self.ydata = data[:,1]
        except IndexError:
            data = np.genfromtxt(dataset, skip_header=25, delimiter=",", usecols=(0,2))

            self.xdata = data[:,0]
            self.ydata = data[:,1]

    def plot(self, rng = [0,0], style = 'o'):
        if rng == [0,0]:
            plt.figure(figsize=(10, 10))
            plt.plot(self.xdata, self.ydata, style)
        else:
            plt.figure(figsize=(10, 10))
            xd, yd = selectDomain.selectdomain(self.xdata, self.ydata, rng)
            plt.plot(xd, yd, style)

    def find_peaks(self, step = 100):
        length = len(self.ydata)

        peaks = np.zeros(length)

        for i in range(0, length, step):
            tempx, tempy = selectDomain.selectdomain(self.xdata, self.ydata, [i, i+step])
            yd = gaussian_filter(self.ydata, 12)
            peaks[signal.argrelmax(yd)] = 1

        print np.where(peaks == 1)





    def peak_fit(self):
        pass

if __name__ == '__main__':
    # mine = co_57
    obj = data_manage("data/co_57.tsv")
    obj.plot([0, 200])
    obj.find_peaks(100)
