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
        self.l = 6.0e-1 #cm
        self.dl = 0.2e-1 #cm

        self.R_H = (np.absolute(self.v_hall) * self.t)/(self.I * self.B)
        self.dR_H = np.sqrt((self.dv_hall/self.v_hall)**2 + (self.dt/self.t)**2 + (self.dI/self.I)**2 + (self.dB/self.B)**2) * self.R_H

        print "R_H", self.R_H
        print "dR_H", self.dR_H

        print "N_A", 1/(self.R_H * 1.602e-19)

        print "v_parallel", self.v_parallel
        print "dv_parallel", self.dv

        self.rho_B = (self.v_parallel * (self.w*self.t))/(self.I*self.l)

        print "rho_B", self.rho_B


    def resistivity_B_0(self):
        data = pd.read_table("data/77K_NoB.tsv")

        val_13 = .5 * (data[self.headers[3]+ " +"].values + data[self.headers[3] + " -"].values)

        self.rho_0 = (np.average(val_13) * (self.w*self.t))/(self.I*self.l)

        self.drho_0 = np.absolute(np.sqrt( (self.dv/self.v_parallel)**2 + (self.dw/self.w)**2 + (self.dt/self.t)**2 + (self.dI/self.I)**2 + (self.dl/self.l)**2 ) * self.rho_0)

        print "rho_0", self.rho_0
        print "drho_0", self.drho_0

        mu_h = np.absolute(self.R_H/self.rho_B)
        dmu_h = np.sqrt(np.absolute( (self.dR_H/self.R_H)**2 + (self.drho_0/self.rho_0) )) * mu_h

        print "mu_h, hall/resist.", mu_h
        print "dmu_h, hall/resist.", dmu_h

        m_h = np.sqrt(np.absolute((self.rho_B - self.rho_0)/(self.rho_0 * self.B**2)))
        dm_h = np.sqrt(np.absolute( (self.dB/self.B)**2 + (self.drho_0/self.rho_B)**2 + (self.drho_0/self.rho_0)**2 )) * m_h

        print "mu_h, magnetoresistance", m_h
        print "dmu_h, magnetoresistance", dm_h

    def temp_scan(self):
        data = pd.read_table("data/temp_scan3.tsv")

        data = data.drop([0,22])

        v_h_12 = .5 * (data[self.headers[1] + " +"].values - data[self.headers[1] + " -"].values)
        v_h_34 = .5 * (data[self.headers[2] + " +"].values - data[self.headers[2] + " -"].values)

        R_H_12 = -(v_h_12 * self.t)/(self.I * self.B)
        R_H_34 = -(v_h_34 * self.t)/(self.I * self.B)

        T0 = np.where(R_H_12 < 2000)[0]

        val_13 = .5 * (data[self.headers[3]+ " +"].values + data[self.headers[3] + " -"].values)

        rho = -(val_13 * self.w * self.t)/(self.I * self.l)

        T = data[self.headers[-1] + " +"].values

        crit_temp = (.5 * (T[T0[0]] + T[T0[1]]))

        print "T_0", crit_temp

        dT = self.dT * np.ones_like(T)
        dR_12 = self.dR_H * np.ones_like(R_H_12)
        dR_34 = self.dR_H * np.ones_like(R_H_34)
        drho = self.drho_0 * np.ones_like(rho)

        # plt.figure(figsize=(10, 10))
        # plt.errorbar(T, R_H_12, xerr = dT, yerr = dR_12, fmt = '.', ms = 1, label = "$R_H$ 1-2")
        # plt.errorbar(T, R_H_34, xerr = dT, yerr = dR_34, fmt = '.', ms = 1, label = "$R_H$ 3-4")
        # plt.xlabel("Temperature (K)")
        # plt.ylabel("$R_H$ (1/C)")
        # plt.legend()
        # plt.savefig("plots/T_vs_R_H.pdf")
        #
        # plt.figure(figsize=(10, 10))
        # plt.errorbar(T, rho, xerr = dT, yerr = drho, fmt='.', ms=1, label = r"$\rho$")
        # plt.xlabel("Temperature (K)")
        # plt.ylabel(r"$\rho$ ($\Omega\cdot$cm)")
        # plt.legend()
        # plt.savefig("plots/rho_vs_T.pdf")

        self.values = {"R_H_12": R_H_12, "R_H_34": R_H_34, "dR_12": dR_12, "dR_34": dR_34, "rho": rho, "drho": drho, "T": T, "dT": dT, "crit_temp": crit_temp, "T0 range": T0}

    def fit_R_H(self):
        xdata = self.values["T"]

        begin = self.values["T0 range"][0]
        end = self.values["T0 range"][-1]

        xdata, ydata_12 = selectDomain.selectdomain(xdata, self.values["R_H_12"], indices = [begin, end])
        _, ydata_34 = selectDomain.selectdomain(xdata, self.values["R_H_34"], indices = [begin, end])

        def exponential(x, A, E):
            return -A * np.exp(E/(2*8.617e-5 * x))

        p = [2e-5, .76]

        popt, pcov = curve_fit(exponential, xdata, ydata_12, p0 = p, maxfev=int(2e6))

        print popt

        yFit = exponential(xdata, *popt)
        yuFit = exponential(xdata, *p)

        plt.plot(xdata, ydata_12, 'o')
        plt.plot(xdata, yuFit, 'r')
        plt.plot(xdata, yFit, 'g')

    def fit_resist(self):
        xdata = self.values["T"]

        begin = self.values["T0 range"][0]
        end = self.values["T0 range"][-1]

        xdata, ydata = selectDomain.selectdomain(xdata, self.values["rho"], indices = [begin, end])
        xerr, yerr = selectDomain.selectdomain(self.values["dT"], self.values["drho"], indices = [begin, end])

        def exponential(x, A, E):
            return A * np.exp(E/(2*8.617e-5 * x))

        p = [2e-5, .76]

        popt, pcov = curve_fit(exponential, xdata, ydata, sigma = yerr, p0 = p, maxfev=int(2e6))

        print "E_g", popt[1], "$\pm$", pcov[1,1]

        yFit = exponential(xdata, *popt)
        yuFit = exponential(xdata, *p)

        redchi = np.sum( (ydata - yFit)**2/(yFit) )
        print "redchi",  redchi

        text = r"$-A e^{\frac{E_g}{2 k_B T}}$" + "\n \
                $A = %.2e \pm %.1e$ \n \
                $E_g = %.6f \pm %.6f$ \n \
                $\chi^2 = %.3f$" % (popt[0], pcov[0,0], popt[1], pcov[1,1], redchi)

        plt.figure(figsize=(10, 10))
        plt.errorbar(xdata, ydata, xerr = xerr, yerr = yerr, fmt = 'o', ms = 1)
        plt.plot(xdata, yFit, 'r')
        plt.text(360, 12, text)
        plt.xlabel("Temperature (K)")
        plt.ylabel(r"$\rho$ ($\Omega\cdot$cm)")
        plt.title("Fitting Resistivity to find the Gap Energy, $E_g$")
        plt.savefig("plots/resistivity_eg.pdf")



if __name__ == '__main__':
    obj = data_manage()

    obj.hall_coeff_and_resistivity()
    obj.resistivity_B_0()
    obj.temp_scan()
    #obj.fit_R_H()
    #obj.fit_resist()
