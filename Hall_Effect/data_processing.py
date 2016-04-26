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
    """
       docstring for data_manage
       Deals with plotting, peak finding, etc. as internal methods

       DON'T USE V_24 - bad data
    """
    def __init__(self):
        v_calib_dataset = "data/room_temp_calibration.tsv"

        v_calib_data = pd.read_table(v_calib_dataset)

        self.headers = ["Voltage Across 100ohm", "Vert leads 1-2", "Vert leads 3-4", "Horz leads 1-3", "Horz leads 2-4", "Thermo Voltage", "Temperature"]

        val_0 = .5 * (v_calib_data[self.headers[-1] + " +"].values + v_calib_data[self.headers[-1] + " -"].values)

        self.dT = np.std(val_0)


    def V_Hall(self):
        data = pd.read_table("data/77K_B.tsv")
        v_h_12_0 = .5 * (data[self.headers[1] + " +"].values - data[self.headers[1] + " -"].values)

        v_h_34_0 = .5 * (data[self.headers[2] + " +"].values - data[self.headers[2] + " -"].values)

        v_h_12 = np.average(v_h_12_0)
        dv_h_12 = np.std(v_h_12_0)
        v_h_34 = np.average(v_h_34_0)
        dv_h_34 = np.std(v_h_34_0)

        self.v_hall = np.average([v_h_12, v_h_34])
        self.dv_hall = np.average([dv_h_12, dv_h_34])

        print "v_hall", self.v_hall
        print "dv_hall", self.dv_hall

    def v_fluctuations(self):
        data = pd.read_table("data/77K_B.tsv")

        v_13 = .5 * (data[self.headers[3] + " +"].values + data[self.headers[3] + " -"].values)


        self.v_parallel = np.average(v_13)

        dv_13 = np.std(v_13)

        self.dv = dv_13

    def hall_coeff_and_resistivity(self):
        self.V_Hall()
        self.v_fluctuations()

        self.t = 1.22e-1 #cm
        self.dt = .01e-1 #cm
        self.w = 4.03e-1 #cm
        self.dw = .01e-1 #cm
        self.I = 4e-3 #A
        self.dI = .1e-3 #A
        self.B = 900 * 1e-8 #V.s/cm^2
        self.dB = (20/1e3) * 1e-8 * self.B

        self.R_H = (np.absolute(self.v_hall) * self.t)/(self.I * self.B)
        self.dR_H = np.sqrt((self.dv_hall/self.v_hall)**2 + (self.dt/self.t)**2 + (self.dI/self.I)**2 + (self.dB/self.B)**2) * self.R_H

        print "R_H", self.R_H
        print "dR_H", self.dR_H

        print "N_A", 1/(self.R_H * 1.602e-19)

        print "v_parallel", self.v_parallel
        print "dv_parallel", self.dv

        self.rho_B = (self.v_parallel*self.w)/self.I

        print "rho_B", self.rho_B


    def resistivity_B_0(self):
        data = pd.read_table("data/77K_NoB.tsv")

        val_13 = .5 * (data[self.headers[3]+ " +"].values + data[self.headers[3] + " -"].values)

        self.rho_0 = (np.average(val_13) * self.w)/self.I

        print "rho_0", self.rho_0
        print "mu_h, hall/resist.", 4*self.R_H/self.rho_0

        print "mu_h, magnetoresistance", np.sqrt((self.rho_B - self.rho_0)/(self.rho_0 * self.B**2))




if __name__ == '__main__':
    obj = data_manage()

    obj.hall_coeff_and_resistivity()
    obj.resistivity_B_0()
