import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
import pandas as pd

import sys
sys.path.append("..")

from PythonLibrary import *

try:
    %pylab inline
except SyntaxError:
    pass

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

class data_manage(object):
    """docstring for data_manage"""
    def __init__(self):
        self.data_folder = "data"

        self.primary_data = ["no_xray_filter.tsv", "with_xray_filter.tsv"]

        self.ss_301 = "SS_302_001.tsv"
        self.ss_velocity = "velocity_ss.tsv"

        
