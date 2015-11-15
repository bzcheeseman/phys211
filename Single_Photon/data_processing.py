import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
from scipy import loadtxt

def main(file):
	
	voltages, incident, coinc12, coinc13 = loadtxt("data/"+file, unpack=True, skiprows=1)
	print file

	##Sinusoidal Fit, coinc12##

	def sin(x, *p):
		return p[0]*np.sin(p[1]*x+p[2]) + p[0]*np.cos(p[1]*x+p[3]) + p[4]

	p0 = [353, .6, 6, -.7, 903]

	popt, pcov = curve_fit(sin, voltages, coinc12, p0 = p0, maxfev=int(2e6))

	print popt

	yFit = sin(voltages, *popt)
	yuFit = sin(voltages, *p0)

	chisq=np.sum((coinc12-yFit)**2/np.sqrt(coinc12)**2)/(len(coinc12)-len(p0))
	print chisq

	##Linear Fit, coinc12##

	p, V = np.polyfit(voltages, coinc13, 1, cov=True)

	lin=np.poly1d(p)

	lin_chisq = np.sum((coinc13-lin(voltages))**2/np.sqrt(coinc13)**2)/(len(coinc13)-2)
	print lin_chisq

	##save chi squared, fit parameters and error, etc. to data file##
	parameter_names=["A", "$\theta$","$\delta_1$", "$\delta_2$", "C"]
	data = np.transpose(np.vstack([popt, np.sqrt(np.absolute(pcov))]))
	linear = np.zeros_like(data)
	linear[0,0] = p[0]
	linear[1,0] = p[1]
	linear[0,1] = V[0,0]
	linear[0,2] = V[0,1]
	linear[1,1] = V[1,0]
	linear[1,2] = V[1,1]
	data = np.vstack([data, linear])

	for i in range(0, 10):
		try:
			if data[i,0] == 0.0:
				data = np.delete(data, i, 0)
		except IndexError:
			pass

	data = np.delete(data, -1, 0)

	os.remove("data/fitparms.csv")

	np.savetxt("data/fitparms.csv", data)

	##Plotting##

	plt.figure(figsize=(10,10))
	plt.errorbar(voltages, coinc12, xerr=0.1*np.ones_like(voltages), 
		yerr=np.sqrt(coinc12), fmt='o', ms=0.1, label="Raw Data, %s" % os.path.splitext(file)[0])
	plt.errorbar(voltages, coinc13, xerr=0.1*np.ones_like(voltages), 
		yerr=np.sqrt(coinc13), fmt='o', ms=0.1, label="Raw Data, %s" % os.path.splitext(file)[0])
	plt.plot(voltages, yFit, linewidth=2, alpha=0.8, color='r', 
		label="Sinusoidal Fit to Interference Pattern")
	#plt.plot(voltages, yuFit, label="guesses")
	plt.plot(voltages, lin(voltages), linewidth=np.sqrt(V[1,1]), alpha=.1, 
		color='r', label="Variance in Linear Fit")
	plt.plot(voltages, lin(voltages), 'r', label="Linear Fit")
	#sinwave
	plt.text(31,2000,r"$A\sin(\theta V + \delta_1) + A\cos(\theta V + \delta_2) + C$")
	plt.text(31,1750,r"$%.1f\sin(%.1f V + (%.1f))$"\
					"\n"\
					r"$ + %.1f\cos(%.1f V + (%.1f)) + %.1f$" % 
		(popt[0], popt[1], popt[2], popt[0], popt[1], popt[3], popt[4]))
	plt.text(31, 1600, r"$\chi^2/\nu = %.2f, \, \nu = %.0f$" % (chisq, len(coinc12)-len(p0)))
	#linear
	plt.text(31, 3300, r"$AV + B$")
	plt.text(31, 3200, r"$%.2fV + %.2f$" % (p[0], p[1]))
	plt.text(31, 3090, r"$\chi^2/\nu = %.2f, \, \nu = %.0f$" % (lin_chisq, len(coinc13)-2))

	plt.legend(loc=0)
	plt.title("Piezo Modulation of Interference Patterns")
	plt.xlabel("Piezo Voltage (V)")
	plt.ylabel("Number of Counts")
	plt.savefig("plots/"+os.path.splitext(file)[0]+".pdf")
	plt.show()


if __name__ == '__main__':
	main("eraser.tsv")

