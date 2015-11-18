import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
from scipy import loadtxt

def main(file):
	
	voltages, incident, coinc12, coinc13 = loadtxt("data/"+file, unpack=True, skiprows=1)
	print file

	##Sinusoidal Fit, coinc12 - all##

	def sin(x, *p):
		return p[0]*np.sin(p[1]*x+p[2]) + p[0]*np.cos(p[1]*x+p[3]) + p[4]

	p0 = [353, .6, 6, -.7, 3100]

	popt, pcov = curve_fit(sin, voltages, coinc12, p0 = p0, maxfev=int(2e6))

	print popt

	yFit = sin(voltages, *popt)
	yuFit = sin(voltages, *p0)

	chisq=np.sum((coinc12-yFit)**2/np.sqrt(coinc12)**2)/(len(coinc12)-len(p0))
	print chisq

	##Sinusoidal Fit, coinc13 - distinguishable + indistinguishable##

	p01 = [353, .6, 6, -.7, 3100]

	popt1, pcov1 = curve_fit(sin, voltages, coinc13, p0 = p0, maxfev=int(2e6))

	print popt1

	yFit1 = sin(voltages, *popt1)
	yuFit1 = sin(voltages, *p01)

	chisq1=np.sum((coinc13-yFit1)**2/np.sqrt(coinc13)**2)/(len(coinc13)-len(p01))
	print chisq1

	##Linear Fit, coinc12##

	p, V = np.polyfit(voltages, coinc12, 1, cov=True)

	lin=np.poly1d(p)

	lin_chisq = np.sum((coinc12-lin(voltages))**2/np.sqrt(coinc12)**2)/(len(coinc12)-2)
	print lin_chisq

	##Linear Fit, coinc13##

	p1, V1 = np.polyfit(voltages, coinc13, 1, cov=True)

	lin1=np.poly1d(p1)

	lin_chisq1 = np.sum((coinc13-lin1(voltages))**2/np.sqrt(coinc13)**2)/(len(coinc13)-2)
	print lin_chisq1

	##save chi squared, fit parameters and error, etc. to data file##
	data = np.transpose(np.vstack([popt, np.sqrt(np.absolute(pcov))]))
	data1 = np.transpose(np.vstack([popt1, np.sqrt(np.absolute(pcov1))]))

	splitter = np.zeros(len(data[0]))
	#splitter[:] = '%'

	linear = np.zeros_like(data)
	linear1 = np.zeros_like(data1)
	linear1[0,0] = p[0]
	linear1[1,0] = p[1]
	linear1[0,1] = V[0,0]
	linear1[0,2] = V[0,1]
	linear1[1,1] = V[1,0]
	linear1[1,2] = V[1,1]

	linear[0,0] = p1[0]
	linear[1,0] = p1[1]
	linear[0,1] = V1[0,0]
	linear[0,2] = V1[0,1]
	linear[1,1] = V1[1,0]
	linear[1,2] = V1[1,1]

	data = np.vstack([data, splitter, data1, splitter, linear, splitter, linear1])

	#np.savetxt("data/fitparms_"+os.path.splitext(file)[0]+".csv", 
	#	data)

	##Plotting##

	plt.figure(figsize=(15,10))
	plt.errorbar(voltages, coinc12, xerr=0.1*np.ones_like(voltages), 
		yerr=np.sqrt(coinc12), fmt='o', ms=0.1, label="Raw Data, %s" % "Coincidences 1-2")
	plt.errorbar(voltages, coinc13, xerr=0.1*np.ones_like(voltages), 
		yerr=np.sqrt(coinc13), fmt='o', ms=0.1, label="Raw Data, %s" % "Coincidences 1-3")
	plt.plot(voltages, yFit, linewidth=2, alpha=0.8, color='b', 
		label="Sinusoidal Fit to Interference Pattern")
	#plt.plot(voltages, lin(voltages), linewidth=np.sqrt(V[1,1]), alpha=.1, 
	#	color='r', label="Variance in Linear Fit")
	plt.plot(voltages, lin(voltages), 'r', label="Linear Fit")

	plt.plot(voltages, yFit1, linewidth=2, alpha=0.8, color='b')
	#plt.plot(voltages, lin1(voltages), linewidth=np.sqrt(V1[1,1]), alpha=.1, 
	#	color='r', label="Variance in Linear Fit")
	plt.plot(voltages, lin1(voltages), 'r')
	#sinwave
	plt.text(43,1000,r"$A\sin(\theta V + \delta_1) + A\cos(\theta V + \delta_2) + C$")
	plt.text(43,500,r"$%.1f\sin(%.1f V + (%.1f))$"\
					"\n"\
					r"$ + %.1f\cos(%.1f V + (%.1f)) + %.1f$" % 
		(popt[0], popt[1], popt[2], popt[0], popt[1], popt[3], popt[4]))
	plt.text(44, 200, r"$\chi^2/\nu = %.2f, \, \nu = %.0f$" % (chisq, len(coinc12)-len(p0)))
	##1-3
	plt.text(43,6800,r"$A\sin(\theta V + \delta_1) + A\cos(\theta V + \delta_2) + C$")
	plt.text(43,6300,r"$%.1f\sin(%.1f V + (%.1f))$"\
					"\n"\
					r"$ + %.1f\cos(%.1f V + (%.1f)) + %.1f$" % 
		(popt1[0], popt1[1], popt1[2], popt1[0], popt1[1], popt1[3], popt1[4]))
	plt.text(44, 6000, r"$\chi^2/\nu = %.2f, \, \nu = %.0f$" % (chisq1, len(coinc13)-len(p0)))
	#linear
	plt.text(31, 1000, r"$AV + B$")
	plt.text(31, 800, r"$%.2fV + %.2f$" % (p[0], p[1]))
	plt.text(31, 500, r"$\chi^2/\nu = %.2f, \, \nu = %.0f$" % (lin_chisq, len(coinc12)-2))
	##1-3
	plt.text(31, 6800, r"$AV + B$")
	plt.text(31, 6600, r"$%.2fV + %.2f$" % (p1[0], p1[1]))
	plt.text(31, 6300, r"$\chi^2/\nu = %.2f, \, \nu = %.0f$" % (lin_chisq1, len(coinc13)-2))

	plt.legend(loc=0)
	plt.title("Piezo Modulation of Interference Patterns")
	plt.xlabel("Piezo Voltage (V)")
	plt.ylabel("Number of Counts")
	plt.savefig("plots/"+os.path.splitext(file)[0]+".pdf")
	plt.show()


if __name__ == '__main__':
	main("indistinguishable.tsv")

