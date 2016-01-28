%pylab inline

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
    #plt.plot(sel_xdata, yuFit, 'g')
    plt.plot(sel_xdata, yFit, 'r', lw = 4, alpha = .7)
    plt.xlabel("Time (s)")
    plt.ylabel("Signal Amplitude (-mV)")
    plt.title("Fitting to find the Larmor Frequency")
    #plt.text(0.00705, 0.0325, "f = %.3f $\pm$ %.3f kHz" % ((popt[1]*1e-3)/(2*np.pi), (np.sqrt(pcov[1,1])*1e-3)/(2*np.pi)))
    plt.savefig("plots/({}.{}).png".format(dset[0].split("/")[1], dset[1]))

    with open('data/fitparms_{}.{}.tsv'.format(dset[0].split("/")[1], dset[1]), 'w+')\
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

    freqs = (1/(2*np.pi))*np.array([data20[2], data33[2], data40[2], data45[2], data55[2], data65[2]])
    df = np.array([data20[1], data33[1], data40[1], data45[1], data55[1], data65[1]])
    # propagate df through previous section - larmor fit

    def lin(x, A, B):
        return A*x + B

    popt, pcov = curve_fit(lin, Bs, freqs, p0 = [87e6, .04])
    yFit = lin(Bs, *popt)

    parms_gamma = []
    parms_omega0 = []

    def mock_data():
        mock_freqs = []
        mock_Bs = []

        fsigma = (df/2*2.5348)
        Bsigma = (Bs * dB/2*2.5348)

        mock_freqs = np.random.randn(len(freqs))*fsigma + freqs
        mock_Bs = np.random.randn(len(freqs))*Bsigma + Bs

        #print mock_Bs, dB

        mock_popt, mock_pcov = curve_fit(lin, mock_Bs, mock_freqs, p0 = [3.742e10, 4.313e4])

        ymockFit = lin(mock_freqs, *mock_popt)

        #plt.plot(mock_Bs, mock_freqs)

        return mock_popt[0], mock_popt[1]

    for i in range(0, 250):  #write c++ function to do this, python takes FOREVER
        tempg, tempo = mock_data()
        parms_gamma.append(tempg)
        parms_omega0.append(tempo)

    #print yFit
    #print np.average(parms_gamma), np.std(parms_gamma)

    devg = np.std(parms_gamma)
    dev0 = np.std(parms_omega0)

    chi_squared = np.sum((freqs - yFit)**2/(df)**2)/(len(freqs))

    plt.figure(figsize = (10, 10))
    plt.errorbar(Bs, freqs, xerr = Bs*dB, yerr = df, fmt = 'o')
    plt.plot(Bs, yFit, alpha = .6)
    plt.text(Bs[0], 2.2e4, "$\omega = \gamma B + \omega_0$\n \
                              $\gamma$ = %.2f $\pm$ %.2f \, GHz/T\n \
                              $\omega_0$ = %.4f $\pm$ %.4f \, GHz \n \
                              $\chi^2/n = %.2f" % ((popt[0]*1e-9), (devg*1e-9), (popt[1]*1e-6), (dev0*1e-6), chi_squared))
    plt.xlabel("Vertical B Field (T)")
    plt.ylabel("Frequency (Hz)")
    plt.title("Finding the Gyromagnetic Ratio")
    plt.savefig("plots/gyromagnetic_ratio_fit.png")
    #plt.show()

def earth_field_fit():
    data = np.genfromtxt("data/earth_B_field_data.csv")

    ydata = data[:,1]
    xdataprime = data[:,0]
    xdata = (8*1.257e-6*36)/(np.sqrt(125)*23.5e-2) * data[:,0]

    #Bhapp = (8*1.257e-6*36)/(np.sqrt(125)*25e-2) * 0.179

    dxdata = np.sqrt(np.average(.001/xdataprime)**2 + (0.5e-2/17e-2)**2)

    def fitform(x, Bhe, Bve, C):
        return 2.895e-9*np.sqrt(pow((Bhe - ((8*1.257e-6*36)/(np.sqrt(125)*25e-2) * 0.179)),2) + pow((Bve - x),2)) + C

    p = [19000e-9, 500000e-9, 350]

    popt, pcov = curve_fit(fitform, xdata, ydata, p0 = p)

    print popt

    yFit = fitform(xdata, *popt)
    yuFit = fitform(xdata, *p)

    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, yuFit, 'g')
    plt.plot(xdata, yFit, 'r')
    plt.ylabel("Frequency (kHz)")
    plt.xlabel("Vertical Applied B Field (T)")
    plt.show()

def plot_trace(dataset):
    data = np.genfromtxt(dataset, skip_header = 1, usecols = (0,1,3))

    time = data[:,0]
    signal = data[:,1]
    sq_wave = .005*data[:,2]

    t, s = selectdomain(time, signal, [-.01, .1])
    t, sq_wave = selectdomain(time, sq_wave, [-.01, .1])

    tm, si = selectdomain(time, signal, [0.00025, .095])

    def exp_func(x, A, B, C):
        return A*np.exp(-B*x) + C

    p = [.03, 31.3, -.02]

    popt, pcov = curve_fit(exp_func, tm, si, p0 = p)

    print popt

    yFit = exp_func(tm, *popt)
    yuFit = exp_func(tm, *p)

    depump = "Depump cycle"

    pump = "Pump cycle"

    plt.figure(figsize = (8,8))
    plt.plot(t, s, label = "Data")
    #plt.plot(tm, yuFit, 'g')
    plt.plot(tm, yFit, 'r', lw = 4, alpha = .6)
    plt.xlabel("Time, (s)")
    plt.ylabel("Intensity, (-mV)")
    plt.title("The Pump-Depump Cycle")
    plt.annotate(depump, (0, 0), (-.018, 0), arrowprops = dict(width=2, headwidth=4, frac = .125, facecolor="red"))
    plt.annotate(pump, (0.04, -0.005), (0.03, 0.005), arrowprops = dict(width=2, headwidth=4, frac = .125, facecolor="red"))
    plt.text(0.05, .005, "$I = Ae^{-Bx} + C$ \n" \
                            r"$\frac{1}{e}$ time of Pump is %.3f $\pm$ %.4f sec" % (1/popt[1], 1/popt[1] * np.sqrt(pcov[1,1])/popt[1]))
    plt.plot(t, sq_wave, label = r"Square wave ($\times$ 0.001 for scale)")
    plt.legend()
    plt.savefig("plots/trace_example.png")

    #plt.show()

def plot_reduced_trace(dataset):
    data = np.genfromtxt(dataset, skip_header=1, usecols=(0,1,3))

    time = data[:,0]
    signal = data[:,1] + 0.07
    sq_wave = .001*data[:,2]

    plt.plot(time, signal, label = "Trace")
    plt.plot(time, sq_wave, label = r"Square Wave ($\times$ 0.001 for scale)")
    plt.xlabel("Time, (s)")
    plt.ylabel("Intensity, (-mV)")
    plt.title("The Pump-Depump Cycle, Reduced Wave")
    plt.legend()
    plt.savefig("plots/reduced_wave.png")

if __name__ == '__main__':

    #larmor_fit("data/larmor_y0.020A.tsv", 80)
    #gyromagnetic_ratio_fit() #DON"T RUN ANYMORE, GOT A GOOD VALUE
    #earth_field_fit() # need to figure out wtf - try to get some shape at least
    #plot_trace("data/s4_3_scope_trace.tsv")
    #plot_reduced_trace("data/s4_4_scope_trace.tsv")
