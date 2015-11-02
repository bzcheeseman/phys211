import numpy as np
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def main(output_plot):

	########rmg plots############

	if output_plot == "rmg":
		ydata = np.multiply([1e-6, .5, 1, 1.5, 2, 2.5, 3, 3.5],(1.410e-3*9.8)/100)
		xdata = np.multiply(1.36e-3, [1.45, 1.52, 1.69, 1.8, 1.92, 2, 2.12, 2.25])
		xdata_err = np.multiply(np.sqrt( np.power(np.multiply(np.ones_like(xdata), .2),2)
	                         +np.power(np.divide(.05, [1.45, 1.52, 1.69, 1.8, 1.92, 2, 2.12, 2.25]), 2)), xdata)
		ydata_err = np.multiply( np.sqrt(np.power(np.divide(.05e-2, np.array([1e-6, .5, 1, 1.5, 2, 2.5, 3, 3.5])/100), 2) 
	                                 + (1e-6/1.410e-3)**2), ydata )

		p, V = np.polyfit(xdata, ydata, 1, cov=True)

		lin = np.poly1d(p)

		data = np.transpose(np.vstack((xdata, xdata_err, ydata, ydata_err)))

		head = "muB \t delta_muB \t rmg \t delta_rmg"

		os.chdir("/users/aman/desktop/phys211/esr/data")

		np.savetxt("magnetic_moment_run1.tsv", data, delimiter="\t", header = head)

		os.chdir("/users/aman/desktop")

		chisq = np.sum( (lin(xdata) - ydata)**2/ydata_err**2 )/len(ydata)

		plt.errorbar(xdata, ydata, yerr=ydata_err, fmt='o',ms = .1, label = "Raw Data")
		plt.plot(xdata, lin(xdata), label = "Linear Fit")
		plt.xlabel("$B \, (T)$")
		plt.ylabel("$rmg \, (N \cdot m)$")
		plt.title("Finding $\mu$")
		plt.text(.002, .0003, "$-rmg = \mu B$"\
	         "\n"
	        "$\mu = %.4f\pm %.4f \, A \cdot m^2$" % (p[0], V[0,0]))
		plt.text(.0024, 0, "$\chi^2/N = %.2f$" % chisq)
		plt.legend(loc = 0)
		plt.savefig("/users/aman/desktop/phys211/esr/plots/rmgvsmub.pdf")
		plt.show()

	#########L plots############
	if output_plot == "L":
		omega = np.divide(2*np.pi, [np.mean([12.8, 14.15, 12.13, 14.08, 13.60]), 11.41, 7.93, 5.48, 4.5, 3.75])
		B = np.multiply(1.36e-3, [.85, 1, 1.5, 2, 2.5, 3])

		sigma = np.std([12.8, 14.15, 12.13, 14.08, 13.60])
		mu = np.mean([12.8, 14.15, 12.13, 14.08, 13.60])
		T = np.array([mu, 11.41, 7.93, 5.48, 4.5, 3.75])
		T_err = np.array([sigma/np.sqrt(5), T[1]*(sigma/mu), T[2]*(sigma/mu), T[3]*(sigma/mu), T[4]*(sigma/mu), T[5]*(sigma/mu)])

		omega_err = (T_err/T)*omega

		B_err = np.multiply(np.sqrt(np.power(np.multiply(np.ones_like(B), .2),2)
	                         +np.power(np.divide(.05, [.85, 1, 1.5, 2, 2.5, 3]), 2)), B)

		L = .00115
		dL = .00005

		def lin(x, *p):
			return p[0]*x

		p, V = curve_fit(lin, B, omega, p0 = [.2])

		yFit = lin(B, *p)

		data = np.transpose(np.vstack((omega, omega_err, B, B_err)))

		head = "omega \t delta_omega \t B \t delta_B"

		os.chdir("/users/aman/desktop/phys211/esr/data")

		np.savetxt("magnetic_moment_run2.tsv", data, delimiter="\t", header = head)

		os.chdir("/users/aman/desktop")

		chisq = np.sum( (yFit - omega)**2/omega_err**2 )/len(omega)

		plt.errorbar(B, omega, yerr = omega_err, fmt = 'o', ms = .1, label = "Raw Data")
		plt.plot(B, yFit, label = "Linear Fit")
		plt.ylabel("$\omega \, (Hz)$")
		plt.xlabel("$B \, (T)$")
		plt.title("Finding $\mu$")
		plt.text(0.0033,0.6, "$\omega = \gamma B$"\
	         "\n"\
	        "$\gamma = %.0f\pm %.0f \, Hz/T$"\
	         "\n"\
	        "$\mu = %.3f \pm %.3f \, A \cdot m^2$"% (p[0], V[0,0], p[0]*L, (V[0,0]/p[0] + dL/L)*p[0]*L))
		plt.text(.0015, 1.1, "$\chi^2/N = %.3f$" % chisq)
		plt.legend(loc = 0)
		plt.savefig("/users/aman/desktop/phys211/esr/plots/omegavsB.pdf")
		plt.show()

	#########Quantum Electron plots############
	if output_plot == "electron":
		frequency = np.array([20.91, 22.93, 24.01, 25.65, 26.5, 27.55, 29.8, 31.07, 31.83]) ## MHz
		vplus = np.array([35.2, 39.2, 40.0, 44.0, 44.8, 45.6, 48.8, 51.2, 52.0]) ## mV
		vminus = np.array([-26.4, -28.8, -29.6, -32.0, -33.6, -34.4, -38.2, -40.0, -41.6]) ## mV

		Bplus = vplus * 0.48e-6 ## T
		Bminus = vminus * 0.48e-6 ## T

		vplus_err = .8/vplus
		vminus_err = .8/vminus

		Bplus_err = vplus_err*Bplus
		Bminus_err = vminus_err*Bminus

		freq_err = np.ones_like(frequency) * .02

		def lin(x, *p):
		    return p[0]*x + p[1]

		pplus, Vplus = np.polyfit(Bplus, frequency, 1, cov=True)

		pminus, Vminus = np.polyfit(Bminus, frequency, 1, cov=True)

		print Vplus[0,0], Vminus[0,0]

		###chi squared

		chisqplus = np.sum( (lin(Bplus, *pplus) - frequency)**2/freq_err**2 )/len(frequency)
		chisqminus = np.sum( (lin(Bminus, *pminus) - frequency)**2/freq_err**2 )/len(frequency)

		print chisqplus, chisqminus

		##plots
		plt.figure(figsize=(10,10))
		plt.errorbar(Bplus, frequency, xerr = Bplus_err, yerr = freq_err, fmt = 'o', ms = .5, label = "Raw Data")
		plt.errorbar(Bminus, frequency, xerr = Bminus_err, yerr = freq_err, fmt = 'o', ms = .5, label = "Raw Data")
		plt.plot(Bplus, lin(Bplus, *pplus), label = "Fit (+)")
		plt.plot(Bminus, lin(Bminus, *pminus), label = "Fit (-)")
		plt.text(-.00001,3.2e1, "$\omega = \gamma B$"\
		         "\n"\
		        "$\gamma_+ = 1.37e\!+\!05 \,MHz/T$" % (pplus[0]))
		plt.text(-.00001,3.16e1, "$\gamma_- = -1.48e\!+\!05 \,MHz/T$")
		plt.text(-.00001, 3e1, "$\chi_+^2/N = %.0f$" % chisqplus)
		plt.text(-.00001, 2.96e1, "$\chi_-^2/N = %.0f$" % chisqminus)
		plt.xlabel("B (T)")
		plt.ylabel("$\omega$ (MHz)")
		plt.title("Finding $\gamma$")
		plt.legend(loc = 0)
		plt.savefig("/users/aman/desktop/phys211/esr/plots/findinggamma.pdf")
		#plt.show()

		plt.figure()
		plt.errorbar(np.linspace(frequency[0], frequency[-1], len(frequency)),
			frequency - lin(Bplus, *pplus), yerr = freq_err, fmt="o", ms = .5, 
			label = "Data Points (+)")
		plt.errorbar(np.linspace(frequency[0], frequency[-1], len(frequency)),
			frequency - lin(Bminus, *pminus), yerr = freq_err, fmt="o", ms = .5, 
			label = "Data Points (-)")
		plt.plot(np.linspace(frequency[0], frequency[-1], len(frequency)), 
			np.zeros_like(frequency), label = "Fit")
		plt.fill_between(np.linspace(frequency[0], frequency[-1], len(frequency)), 
			np.sqrt(np.sqrt(Vplus[1,1]/np.sqrt(2*np.log(2)))), -np.sqrt(np.sqrt(Vplus[1,1]/np.sqrt(2*np.log(2)))),
			label = r"$(.05 \times \sigma)^{1/4}$", alpha = .35)
		plt.ylabel("Difference (MHz)")
		plt.title("Cross Sectional Differences")
		plt.legend(loc = 0)
		plt.savefig("/users/aman/desktop/phys211/esr/plots/gamma_diff.pdf")
		plt.show()

	return 0

if __name__ == '__main__':
	import sys
	main(sys.argv[1])
