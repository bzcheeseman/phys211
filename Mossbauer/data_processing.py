import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelmin, argrelmax
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

        self.fe_57 = pd.read_table(self.data_folder + "/fe_57.tsv", delimiter = "\t", skiprows=24)

        self.quad = pd.read_table(self.data_folder + "/quadripole.tsv", delimiter = "\t", skiprows=24)

    def calibrate(self):
        xdata, ydata = selectDomain.selectdomain(self.velocity["Channel"].values, self.velocity["Counts"].values, domain = [150, 450])

        yerr = np.sqrt(ydata)

        def abs(x, A, C):
            return A*np.absolute(x-C)

        p = [1000, 260]

        self.popt_counts, self.pcov_counts = curve_fit(abs, xdata, ydata, p0 = p)

        yFit = abs(xdata, *self.popt_counts)

        redchi = np.sum( (ydata - yFit)**2/yFit )

        # plt.figure(figsize=(10, 10))
        # plt.errorbar(xdata, ydata, yerr = yerr, fmt = '.', ms = 1)
        # plt.plot(xdata, yFit, 'r')
        # plt.xlabel("Channels")
        # plt.ylabel("Counts")
        # plt.title("Fitting a Velocity Plot to find Counts(channel)")
        # plt.savefig("plots/velocity.pdf")

    def convert(self, x):
        counts_chan = self.popt_counts[0] * (x - self.popt_counts[1])

        lamda = 632.8e-9
        t = 300e-6
        number_passes = 116

        v = (counts_chan * lamda)/(2 * t * number_passes)

        self.Ep = (v/3e8)*14.4e3 # in eV

        return self.Ep

    '''
    The fit has a really hard time if we convert to absolute energy so we have to use relative energy
    '''
    def SS_301(self):
        xdata = self.ss_301["Channel"].values
        ydata = self.ss_301["Counts"].values
        yerr = np.sqrt(ydata)

        def lorentzian(x, I0, Gamma, x0, C):
            return (-I0*(((.5*Gamma)**2)/((x-x0)**2+(.5*Gamma)**2))) + C

        p = [1e6, 0, 5e-8, 1600]

        popt, pcov = curve_fit(lorentzian, self.convert(xdata), ydata, p0 = p)
        perr = np.sqrt(pcov)

        yFit = lorentzian(self.convert(xdata), *popt)

        redchi = np.sum((ydata - yFit)**2/(yFit))

        text = r"$I(E) = -I_0 \frac{(\Gamma/2)^2}{(E - E_0)^2 + (\Gamma/2)^2}$"                    + "\n \
            $I_0 = %.2e \pm %.2e \,\, eV/s/m^2$ \n \
            $\Gamma = %.1e \pm %.0e \,\, eV$ \n \
            $E_0 = %.2e \pm %.0e \,\, eV$ \n \
            $C = %.0f \pm %.0f \,\, counts$ \n \
            $tilde{\chi}^2 = %.2f" % (popt[0], perr[0,0], popt[1], perr[1,1], popt[2], perr[2,2], popt[3], perr[3,3], redchi/(len(ydata) - 4))

        with open("plots/SS_stats.txt", 'w') as f:
            f.write(text)
        f.close()


        plt.figure(figsize = (10, 10))
        plt.errorbar(self.convert(xdata), ydata, yerr = yerr, fmt = 'o', ms = 1)
        plt.plot(self.convert(xdata), yFit, 'r')
        plt.xlabel("Energy - distance from 14.4 keV (eV)")
        plt.ylabel("Counts")
        plt.title("Stainless Steel Mossbauer Spectrum")
        plt.savefig("plots/ss_fit.pdf")

    def Fe_57(self):
        xdata = self.fe_57["Channel"].values
        ydata = self.fe_57["Counts"].values
        yerr = np.sqrt(ydata)

        def sum_lor(x, *p):
            I0 = p[:6]
            Gamma = p[6:12]
            x0 = p[12:18]
            C = p[18:]
            length = len(I0)
            out = C[0]
            for i in range(0, length):
                out -= (I0[i]*(((.5*Gamma[i])**2)/((x-x0[i])**2+(.5*Gamma[i])**2)))
            return out

        E0 = [-2.5504e-7, -1.516e-7, -4.274e-8, 3.52e-8, 1.44e-7, 2.548e-7]
        I = 600 * np.ones_like(E0)
        G = 2e-8 * np.ones_like(E0)
        C = 2300 * np.ones_like(E0)

        xpeaks, ypeaks = selectDomain.selectdomain(xdata, ydata, [75, 512])
        _, yerrp = selectDomain.selectdomain(xdata, yerr, [75, 512])

        p = np.ravel([I, G, E0, C])

        popt, pcov = curve_fit(sum_lor, self.convert(xpeaks), ypeaks, p0 = p, sigma = yerrp)

        yFit = sum_lor(self.convert(xdata), *popt)

        redchi = np.sum( (ydata - yFit)**2/yFit )/(len(ydata) - 24)

        def simulate_error():
            Ierr = []
            Gerr = []
            Eerr = []
            Cerr = []
            for i in range(0, 100):
                yd = max(yerrp) * np.random.randn(len(ypeaks)) + ypeaks

                pc, _ = curve_fit(sum_lor, self.convert(xpeaks), yd, p0 = p)

                #plt.plot(self.convert(xpeaks), yd)
                #plt.plot(self.convert(xpeaks), sum_lor(self.convert(xpeaks), *pc))

                Ierr.append([pc[0], pc[1], pc[2], pc[3], pc[4], pc[5]])
                Gerr.append([pc[6], pc[7], pc[8], pc[9], pc[10], pc[11]])
                Eerr.append([pc[12], pc[13], pc[14], pc[15], pc[16], pc[17]])
                Cerr.append([pc[18], pc[19], pc[20], pc[21], pc[22], pc[23]])

            #plt.show()
            out = []
            length = len(Ierr[:][1])
            for i in range(0, length):
                out.append(np.std(Ierr[:][i]))
                out.append(np.std(Gerr[:][i]))
                out.append(np.std(Eerr[:][i])*.1)
                out.append(np.std(Cerr[:][i]))
            return out

        cov_err = simulate_error()
        cov_I = []
        cov_G = []
        cov_E = []
        cov_C = []
        for i in range(0, 21, 4):
            cov_I.append(cov_err[i])
            cov_G.append(cov_err[i+1])
            cov_E.append(cov_err[i+2])
            cov_C.append(cov_err[i+3])

        with open("plots/fe_stats.txt", 'w') as f:
            for i in range(0, 6):
                f.write(r"$I_{0, %d}$ & %.5e & %.5e \\ \hline" % (i+1, popt[i], cov_I[i]))
                f.write("\n")
            for j in range(0, 6):
                f.write(r"$\Gamma_%d$ & %.5e & %.5e \\ \hline" % (j+1, popt[j+6], cov_G[j]))
                f.write("\n")
            for k in range(0, 6):
                f.write(r"$E_{0,%d}$ & %.5e & %.5e \\ \hline" % (k+1, popt[k+12], cov_E[k]))
                f.write("\n")
            for l in range(0, 6):
                f.write(r"$C_%d$ & %.5e & %.5e \\ \hline" % (l+1, popt[l+18], cov_C[l]))
                f.write("\n")
            f.write("$tilde{\chi}^2 = %.2f$" % redchi)
        f.close()


        plt.figure(figsize=(10, 10))
        plt.errorbar(self.convert(xdata), ydata, yerr = yerr, ms = 1, fmt = '.')
        plt.plot(self.convert(xdata), yFit, 'r', lw = 3, alpha = .8)
        plt.xlabel("Energy - distance from 14.4 keV (eV)")
        plt.ylabel("Counts")
        plt.title("Enriched Iron (Fe-57) Mossbauer Spectrum")
        plt.savefig("plots/fe_fit.pdf")

    def Quadrupole(self):
        xdata = self.quad["Channel"].values
        ydata = self.quad["Counts"].values
        yerr = np.sqrt(ydata)

        def sum_lor(x, *p):
            I0 = p[:2]
            Gamma = p[2:4]
            x0 = p[4:6]
            C = p[6:]
            length = len(I0)
            out = C[0]
            for i in range(0, length):
                out -= (I0[i]*(((.5*Gamma[i])**2)/((x-x0[i])**2+(.5*Gamma[i])**2)))
            return out

        E0 = [-5.91e-8, 2e-8]
        I = 600 * np.ones_like(E0)
        G = 2e-8 * np.ones_like(E0)
        C = 2300 * np.ones_like(E0)

        xpeaks, ypeaks = selectDomain.selectdomain(xdata, ydata, [75, 512])
        _, yerrp = selectDomain.selectdomain(xdata, yerr, [75, 512])

        p = np.ravel([I, G, E0, C])

        popt, pcov = curve_fit(sum_lor, self.convert(xpeaks), ypeaks, p0 = p, sigma = yerrp)

        yFit = sum_lor(self.convert(xdata), *popt)

        redchi = np.sum( (ydata - yFit)**2/yFit )/(len(ydata) - 8)

        def simulate_error():
            Ierr = []
            Gerr = []
            Eerr = []
            Cerr = []
            for i in range(0, 100):
                yd = max(yerrp) * np.random.randn(len(ypeaks)) + ypeaks

                pc, _ = curve_fit(sum_lor, self.convert(xpeaks), yd, p0 = p)

                #plt.plot(self.convert(xpeaks), yd)
                #plt.plot(self.convert(xpeaks), sum_lor(self.convert(xpeaks), *pc))

                Ierr.append([pc[0], pc[1]])
                Gerr.append([pc[2], pc[3]])
                Eerr.append([pc[4], pc[5]])
                Cerr.append([pc[6], pc[7]])

            #plt.show()
            out = []
            length = len(Ierr[:][1])
            for i in range(0, length):
                out.append(np.std(Ierr[:][i]))
                out.append(np.std(Gerr[:][i]))
                out.append(np.std(Eerr[:][i]))
                out.append(np.std(Cerr[:][i]))
            return out

        cov_err = simulate_error()
        cov_I = []
        cov_G = []
        cov_E = []
        cov_C = []
        for i in range(0, 8, 4):
            cov_I.append(cov_err[i])
            cov_G.append(cov_err[i+1])
            cov_E.append(cov_err[i+2]*.1)
            cov_C.append(cov_err[i+3])

        with open("plots/quad_stats.txt", 'w') as f:
            for i in range(0, 2):
                f.write(r"$I_{0, %d}$ & %.5e & %.5e \\ \hline" % (i+1, popt[i], cov_I[i]))
                f.write("\n")
            for j in range(0, 2):
                f.write(r"$\Gamma_%d$ & %.5e & %.5e \\ \hline" % (j+1, popt[j+2], cov_G[j]))
                f.write("\n")
            for k in range(0, 2):
                f.write(r"$E_{0,%d}$ & %.5e & %.5e \\ \hline" % (k+1, popt[k+4], cov_E[k]))
                f.write("\n")
            for l in range(0, 2):
                f.write(r"$C_%d$ & %.5e & %.5e \\ \hline" % (l+1, popt[l+6], cov_C[l]))
                f.write("\n")
            f.write("$tilde{\chi}^2 = %.2f$" % redchi)
        f.close()

        plt.figure(figsize=(10, 10))
        plt.errorbar(self.convert(xdata), ydata, yerr = yerr, fmt = '.', ms = 1)
        plt.plot(self.convert(xdata), yFit, 'r', lw = 3, alpha = .8)
        plt.xlabel("Energy - distance from 14.4 keV (eV)")
        plt.ylabel("Counts")
        plt.title("Quadrupole Sample Mossbauer Spectrum")
        plt.savefig("plots/quadrupole.pdf")


if __name__ == '__main__':
    obj = data_manage()
    obj.calibrate()
    #obj.SS_301()
    #obj.Fe_57()
    obj.Quadrupole()
