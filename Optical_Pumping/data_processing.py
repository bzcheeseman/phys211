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

def larmor_fit(dataset, freq):
    data = np.genfromtxt(dataset, skip_header = 1, usecols = (0,1))

    dset = dataset.split(".")

    xdata = data[:,0]
    ydata = data[:,1]

    sel_xdata, sel_ydata = selectdomain(xdata, ydata, [0.00702, 0.00735])

    def sin_function(x, A, w, d, c):
        return A*np.sin(w*x + d) + c

    p = [.01, freq * 1e3, np.pi/2, .035]

    popt, pcov = curve_fit(sin_function, sel_xdata, sel_ydata, p0 = p)

    yFit = sin_function(sel_xdata, *popt)
    yuFit = sin_function(sel_xdata, *p)

    parms = {"A":popt[0], "freq":popt[1], "df":np.sqrt(pcov[1,1]), "shift":popt[2], "vertical":popt[3]}

    print parms["freq"]/(2*np.pi), "Hz"

    plt.figure(figsize = (10,8))
    plt.plot(sel_xdata, sel_ydata, 'o')
    plt.plot(sel_xdata, yuFit, 'g')
    plt.plot(sel_xdata, yFit, 'r', lw = 4, alpha = .7)
    plt.xlabel("time (s)")
    plt.ylabel("signal amplitude (-mV)")
    plt.savefig("plots/({}.{}).png".format(dset[0].split("/")[1], dset[1]))

    with open('data/fitparms_{}.{}.tsv'.format(dset[0].split("/")[1], dset[1]), 'w')\
      as f:
        for key in sorted(parms.keys()):
            f.write("{}\t{}\n".format(key, parms[key]))

def gyromagnetic_ratio_fit():
    data20 = np.genfromtxt("data/fitparms_larmor_y0.020A.tsv", usecols=(1))
    data33 = np.genfromtxt("data/fitparms_larmor_y0.033A.tsv", usecols=(1))
    data40 = np.genfromtxt("data/fitparms_larmor_y0.040A.tsv", usecols=(1))
    data45 = np.genfromtxt("data/fitparms_larmor_y0.045A.tsv", usecols=(1))
    data55 = np.genfromtxt("data/fitparms_larmor_y0.055A.tsv", usecols=(1))
    data65 = np.genfromtxt("data/fitparms_larmor_y0.065A.tsv", usecols=(1))

    Is = np.array([.02, .033, .04, .045, .055, .065])
    dI = .001 * np.ones_like(Is)

    Bs = (8*1.257e-6*9)/(np.sqrt(125)*17e-2)*Is  # y field - 9 turns

    dB = np.sqrt(np.average(dI/Is)**2 + (0.5e-2/17e-2)**2) * np.ones_like(Bs) #propagate errors through, da = .5cm

    freqs = [data20[2], data33[2], data40[2], data45[2], data55[2], data65[2]]
    df = np.array([data20[1], data33[1], data40[1], data45[1], data55[1], data65[1]])
    # propagate df through previous section - larmor fit

    def lin(x, A, B):
        return A*x + B

    popt, pcov = curve_fit(lin, Bs, freqs, p0 = np.polyfit(Bs, freqs, 1))
    yFit = np.poly1d(popt)

    parms_gamma = []
    parms_omega0 = []

    def mock_data():
        mock_freqs = []
        mock_Bs = []

        for j in range(0, len(freqs)):
            mock_freqs.append(np.random.normal(freqs[j], df[j]/2*2.5348))
            mock_Bs.append(np.random.normal(Bs[j], Bs[j] * dB[j]/2*2.5348))

        #print mock_Bs, dB

        mock_popt, mock_pcov = curve_fit(lin, mock_Bs, mock_freqs, p0 = [3.742e10, 4.313e4])

        ymockFit = lin(mock_freqs, *mock_popt)

        #plt.plot(mock_Bs, mock_freqs, 'o')

        return mock_popt[0], mock_popt[1]

    for i in range(0, 2):  #write c++ function to do this, python takes FOREVER
        tempg, tempo = mock_data()
        parms_gamma.append(tempg)
        parms_omega0.append(tempo)

    print yFit
    print np.average(parms_gamma), np.std(parms_gamma)

    plt.errorbar(Bs, freqs, xerr = Bs*dB, yerr = df, fmt = 'o')
    plt.plot(Bs, yFit(Bs), lw = np.std(parms_gamma), alpha = .6)
    plt.text(Bs[0], 1.4e5, "$\omega = \gamma B + \omega_0$\n \
                              %f B + %f" % (popt[0], popt[1]))
    plt.xlabel("Vertical B Field (T)")
    plt.ylabel("Frequency (Hz)")
    plt.savefig("plots/gyromagnetic_ratio_fit.png")
    #plt.show()

def earth_field_fit():
    data = np.genfromtxt("data/earth_B_field_data.csv")

    xdata = data[:,0]
    ydataprime = data[:,1]
    ydata = (8*1.257e-6*36)/(np.sqrt(125)*23.5e-2) * data[:,1]

    Bhapp = (8*1.257e-6*36)/(np.sqrt(125)*25e-2) * 0.179

    dydata = np.sqrt(np.average(.001/ydataprime)**2 + (0.5e-2/17e-2)**2)

    def fitform(x, Bhe, Bve):
        return 2.895e-9*np.sqrt((Bhe - Bhapp)**2 + (Bve - x)**2)

    p = [1e-3, .325]

    popt, pcov = curve_fit(fitform, xdata, ydata, p0 = p)

    print popt

    yFit = fitform(xdata, *popt)
    yuFit = fitform(xdata, *p)

    plt.plot(xdata, ydata, 'o')
    #plt.plot(xdata, yuFit, 'g')
    plt.plot(xdata, yFit, 'ro')
    plt.ylabel("Frequency (kHz)")
    plt.xlabel("Vertical Applied B Field (T)")
    plt.show()


if __name__ == '__main__':

    #larmor_fit("data/larmor_y0.065A.tsv", 150)
    gyromagnetic_ratio_fit()
    #earth_field_fit()
