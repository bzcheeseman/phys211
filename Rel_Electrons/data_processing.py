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
    return E/(1+(E/(m0 * 3e8**2))*(1-np.cos(theta)))

class data_manage(object):
    """
       docstring for data_manage
       Deals with plotting, peak finding, etc. as internal methods
    """
    def __init__(self, dataset):
        name = dataset.split("/")[-1].split(".")[0]
        self.name = "%s%s-%s" % (name.split("_")[0][0].upper(), name.split("_")[0][1], name.split("_")[1])
        try:
            data = np.genfromtxt(dataset, skip_header=25)

            self.xdata = data[:,0]
            self.ydata = data[:,1]
            self.yerr = np.sqrt(self.ydata)
        except IndexError:
            data = np.genfromtxt(dataset, skip_header=25, delimiter=",", usecols=(0,2))

            self.xdata = data[:,0]
            self.ydata = data[:,1]

    def calibrate(self, data_folder = "data/", plot = False):
        data_files = ["co_57.tsv", "na_22.tsv", "bi_207.tsv", "cs_137.csv", "na_22.tsv", "in_116.tsv", "bi_207.tsv"]
        E = [.122, .511, .5696, .662, 1.2745, 1.2933, 1.7697] # MeV


        ch = [45, 394, 450, 522, 997, 1019,  1393]

        ch_real = []
        dch_real = []
        for i in range(0, len(ch)-1):
            try:
                data = np.genfromtxt(data_folder+data_files[i], skip_header=25)

                xdata = data[:,0]
                ydata = data[:,1]
            except IndexError:
                data = np.genfromtxt(data_folder+data_files[i], skip_header=25, delimiter=",", usecols=(0,2))

                xdata = data[:,0]
                ydata = data[:,1]

            p = [ydata[ch[i]], 5, ch[i], 0]
            popt, pcov = curve_fit(residuals.gaussian, xdata, ydata, p0 = p)

            # yuFit = residuals.gaussian(xdata, *p)
            # yFit = residuals.gaussian(xdata, *popt)
            # plt.figure(figsize=(10, 10))
            # plt.plot(xdata, ydata)
            # plt.plot(xdata, yFit)
            # plt.plot(xdata, yuFit)
            # plt.show()

            ch_real.append(popt[2])
            dch_real.append(pcov[2,2])

        ch_real.append(1393)
        dch_real.append(dch_real[0])
        del popt


        popt = np.polyfit(E, ch_real, deg = 1)

        def find_pcov(nruns):
            pcov = [0,0]
            p0 = []
            p1 = []

            for i in range(0, nruns):
                chfake = np.random.randn(len(ch)) * 1 + ch
                poptfake = np.polyfit(E, chfake, deg = 1)

                p0.append(poptfake[0])
                p1.append(poptfake[1])

                i+=1

            pcov[0] = np.std(p0)
            pcov[1] = np.std(p1)

            return pcov


        pcov = find_pcov(100)

        self.A = 1/popt[0]
        self.B = -popt[1]/popt[0]
        self.dA = np.absolute(pcov[0]/popt[0]) * self.A
        self.dB = np.sqrt((pcov[1]/popt[1])**2 + (pcov[0]/popt[0])**2) * self.B

        if plot:
            ch_E = np.poly1d(popt)

            redchi = np.sum(np.power((ch_real - ch_E(E)), 2)/ch_E(E))

            plt.figure(figsize = (10, 10))
            plt.title("Calibration")
            plt.errorbar(E, ch_real, yerr = dch_real, fmt = 'o', label = "Data")
            plt.plot(E, ch_E(E), 'r-', label = "A E + B")
            plt.text(1, 400, "$(%.2e \pm %.2e)E + (%.2e \pm %.2e)$ \n\
            $\chi^2 = %.2f$" % (self.A, self.dA, self.B, self.dB, redchi))
            plt.legend()
            plt.xlabel("Energy (MeV)")
            plt.ylabel("Channel")
            plt.savefig("plots/calibration.pdf")

    def convert(self, x):
        return residuals.linear(x, self.A, self.B)

    def plot(self, rng = [0,0], style = 'o', xlabel = None, ylabel = None, title = None, save = False):

        def pl(x, y, xerr, yerr, save = False):
            plt.figure(figsize = (10, 10))
            plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=style)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
            if save:
                plt.savefig("plots/%s.pdf" % self.name)

        if rng == [0,0]:
            pl(self.convert(self.xdata), self.ydata, None, np.sqrt(self.ydata))
        else:
            xd, yd = selectDomain.selectdomain(self.convert(self.xdata), self.ydata, rng)
            ye = np.sqrt(yd)
            pl(xd, yd, None, ye)

    def find_peaks(self):
        x = self.convert(self.xdata)
        y = gaussian_filter(self.ydata, 6)
        yd = self.ydata[signal.argrelmax(y, order=3)]

        yd = yd[yd >= 2]

        peaks = []
        sigmas = []
        As = []
        for yitem in yd:
            cent = np.where(self.ydata == yitem)

            p = [yitem, 1e-3, x[cent[0][0]], 0]

            popt, pcov = curve_fit(residuals.gaussian, x, self.ydata, p0 = p)

            if np.absolute(popt[1]) < 1e-2 and popt[0] > 0:
                peaks.append(popt[2])
                sigmas.append(popt[1])
                As.append(popt[0])

        self.peaks = {"locations":peaks, "sigmas":sigmas, "As":As}
        print self.peaks["locations"]

    def peak_fit(self, save = False):
        l = len(self.peaks["locations"])

        for i in range(0, l):
            p = [self.peaks["As"][i], self.peaks["sigmas"][i], self.peaks["locations"][i], 0]

            x, y = selectDomain.selectdomain(self.convert(self.xdata), self.ydata, [p[2] - 6*p[1], p[2] + 6*p[1]])

            popt, pcov = curve_fit(residuals.gaussian, x, y, p0 = p)

            xFit = np.linspace(x[0], x[-1], 100)

            yFit = residuals.gaussian(xFit, *popt)
            yF = residuals.gaussian(x, *popt)

            redchi = np.sum((y - yF)**2/yF)/(len(x)-4)

            plt.figure(figsize=(10,10))
            plt.title(self.name)
            plt.errorbar(x, y, yerr=np.sqrt(y), fmt='.', ms= 1, label = "Data")
            plt.plot(xFit, yFit, label = "Fitted Curve")
            plt.xlabel("Energy (MeV)")
            plt.ylabel("Counts")
            plt.text(x[0], popt[0]/2, "Fit Form: \n\
            $C(E) = Ae^{-(E-\mu)^2/2\sigma^2} + D$ \n\
            $A = %.2e \pm %.2e$ \n\
            $\sigma = %.2e \pm %.2e$ \n\
            $\mu = %.2e \pm %.2e$ \n\
            $D = %.2e \pm %.2e$ \n\
            $\overline{\chi}^2 = %.2f$" % (popt[0], pcov[0,0], popt[1], pcov[1,1], popt[2], pcov[2,2], popt[3], pcov[3,3], redchi))
            plt.legend()
            plt.savefig("plots/%s_%.4f.pdf" % (self.name, popt[2]))

    def notebook_process(self):
        T = [310, 675, 375, 160, 175, 830, 259, 200, 475, 700, 850]
        E = [450, 840, 522, 278, 301, 997, 394, 328, 649, 862, 1019]
        dT = 4*np.ones_like(T)
        dE = np.ones_like(E)

        T = self.convert(np.sort(T))
        dT = self.convert(dT)
        E = self.convert(np.sort(E))
        dE = self.convert(dE)

        # print E

        pc = 2*E - T

        mnr = pc**2/(2*T)
        dmnr = np.sqrt((dT/T)**2 + (dE/E)**2) * mnr

        p0 = [.5, .511]

        popt, pcov = curve_fit(residuals.linear, T, mnr, p0 = p0)
        b = popt[0]*2
        m = popt[1]/.511

        yFit = residuals.linear(T, *popt)
        yuFit = residuals.linear(T, popt[0], .511)

        plt.figure(figsize=(10, 10))
        plt.title("Finding $m_e^*$ and the Relativistic Dispersion Relation")
        plt.errorbar(T, mnr, xerr=dT, yerr=dmnr, fmt='.')
        plt.plot(T, yuFit, label="$m = 511 \,\, keV/c^2$")  #fix this
        plt.plot(T, yFit, label = "$m = (%.4f \pm %.4f) \, MeV/c^2$" % (popt[1], pcov[1,1]))  #fix this
        plt.ylabel("Effective Mass (eV)")  #not effective mass anymore
        plt.xlabel("Kinetic Energy (kg)")  #sure
        plt.text(T[0], 1, "$p^2c^2/2T = BT + C$ \n\
        $B = 1/2$ \n\
        $C = mc^2$ \n\
        $2B \, (unitless) = %.3f \pm %.3f$ \n\
        $C/.511 \, (unitless) = %.4f \pm %.4f$" % (b, pcov[0,0]/popt[0] * b, m, pcov[1,1]/popt[1] * m))
        plt.legend()
        plt.savefig("plots/dispersion.pdf") # gotta see if this actually works or not - effective mass is kinda iffy

    def plot_spectra(self):
        self.plot([0, 2.3], xlabel = "Energy (MeV)", ylabel = "Counts", title = self.name)
        plt.annotate("Full Energy Peak", xy=(.17, 4500), xytext=(.2, 4900), arrowprops=dict(facecolor='black'))
        plt.annotate("Compton Edge", xy=(.3, 1000), xytext=(.33, 4200), arrowprops=dict(facecolor='black'))
        plt.annotate("Full Energy Peak", xy=(.44, 3900), xytext=(.6, 4000), arrowprops=dict(facecolor='black'))
        plt.annotate("Full Energy Peak", xy=(.85, 1000), xytext=(.55, 2500), arrowprops=dict(facecolor='black'))
        plt.annotate("Compton Edge", xy=(.92, 600), xytext=(.65, 3000), arrowprops=dict(facecolor='black'))
        plt.annotate("Compton Edge", xy=(1.05, 500), xytext=(.65, 3000), arrowprops=dict(facecolor='black'))
        plt.annotate("Full Energy Peak", xy=(1.12, 1900), xytext=(1.3, 3000), arrowprops=dict(facecolor='black'))
        plt.annotate("Full Energy Peak", xy=(1.3, 2200), xytext=(1.3, 3000), arrowprops=dict(facecolor='black'))
        plt.annotate("Full Energy Peak", xy=(1.5, 300), xytext=(1.5, 1500), arrowprops=dict(facecolor='black'))
        plt.annotate("Compton Edge", xy=(1.88, 100), xytext=(1.6, 1000), arrowprops=dict(facecolor='black'))
        plt.annotate("Full Energy Peak", xy=(2.08, 300), xytext=(1.7, 1200), arrowprops=dict(facecolor='black'))
        plt.savefig("plots/%s.pdf" % self.name)



if __name__ == '__main__':

    # for dset in ["ba_133.tsv", "bi_207.tsv", "cs_137.csv", "in_116.tsv", "na_22.tsv"]:
    #     print dset
    #     obj = data_manage("data/%s" % dset)
    #     obj.calibrate(plot=False)
    #     obj.find_peaks()
    #     pks.append(obj.peak_fit())

    obj = data_manage("data/in_116.tsv")
    obj.calibrate(plot=False)
    # obj.find_peaks()
    # obj.peak_fit()
    # obj.notebook_process()
