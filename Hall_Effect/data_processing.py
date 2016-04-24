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
    def __init__(self):
        v_calib_dataset = "data/room_temp_calibration.tsv"

        v_calib_data = pd.read_table(v_calib_dataset)

        self.headers = ["Voltage Across 100ohm", "Vert leads 1-2", "Vert leads 3-4", "Horz leads 1-3", "Horz leads 2-4", "Thermo Voltage", "Temperature"]

        val = .5 * (v_calib_data[self.headers[-1] + " +"] + v_calib_data[self.headers[-1] + " -"])

        print val


    def plot(self):
        plt.figure(figsize=(10, 10))
        plt.plot(self.xdata, self.ydata, 'o')


if __name__ == '__main__':
    obj = data_manage()
