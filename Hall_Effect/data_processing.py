import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
import pandas as pd

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

class data_manage(object):
    """
       docstring for data_manage
       Deals with plotting, peak finding, etc. as internal methods
    """
    def __init__(self, v_calib_dataset):
        try:
            data = pd.read_table(v_calib_dataset)
            print data
            print "\n"
            print data["Vert leads 1-2 +"].values

        except IndexError:
            pass


    def plot(self):
        plt.figure(figsize=(10, 10))
        plt.plot(self.xdata, self.ydata, 'o')


if __name__ == '__main__':
    obj = data_manage("data/room_temp_calibration.tsv")
