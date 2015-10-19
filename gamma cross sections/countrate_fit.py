import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import loadtxt
import os

# will have to change the input dataset anyway - watch out!!

def choose_domain(xdata, ydata, domain, upper = None):
	import numpy as np

	if upper != None:
		locs = np.where(ydata >= upper)[0]
		for index in locs:
			ydata[index] = upper

	lower_bound = np.where(xdata == domain[0])[0]
	upper_bound = np.where(xdata == domain[1])[0]

	return xdata[lower_bound:upper_bound], ydata[lower_bound:upper_bound]

def main(dataset, show = False):
	xdata, n1data, n2data, n3data, n1err, n2err, n3err = loadtxt(dataset, unpack=True, usecols=[0,1,2,3,4,5,6], skiprows = 1)
	fitxdata = np.linspace(xdata[0], xdata[-1], 500)

	n1err = np.multiply(n1data, n1err)
	n2err = np.multiply(n2data, n2err)
	n3err = np.multiply(n3data, n3err)

	def expand_exponential(x, *p):
		return (p[0] * np.exp(p[1] * x)) + (p[2] * np.exp(p[3] * x)) + p[4]
	def exponential(x, *p):
		return p[0] * np.exp(p[1] * x) + p[2]

	p = [129.716, -.00236, -4.8764]
	exp_p1 = [129, -.00236, 300, -.09, -4]
	exp_p2 = [85, -.0236, -4]

	lamda = -.00236
	scale = 450

	#lambda x, a, c, d, e: expand_exponential(x, a, lamda, c, d, e)

	popt1, pcov1 = curve_fit(expand_exponential, xdata, n1data, p0 = exp_p1, sigma = n1err, maxfev = int(3e8)) #expanded exponential
	popt2, pcov2 = curve_fit(exponential, xdata, n2data, p0 = exp_p2, sigma = n2err, maxfev = int(3e8))
	popt3, pcov3 = curve_fit(exponential, xdata, n3data, p0 = p, sigma = n3err)

	pcov1 = np.absolute(pcov1)
	pcov2 = np.absolute(pcov2)
	pcov3 = np.absolute(pcov3)

	yFit1 = expand_exponential(xdata, *popt1)
	yuFit1 = expand_exponential(xdata, *exp_p1)

	yFit2 = exponential(xdata, *popt2)
	yuFit2 = exponential(xdata, *exp_p2)

	yFit3 = exponential(xdata, *popt3)
	yuFit3 = exponential(xdata, *p)

	chisq1 = np.sum(((yFit1-n1data)/n1err)**2)
	chisq2 = np.sum(((yFit2-n2data)/n2err)**2)
	chisq3 = np.sum(((yFit3-n3data)/n3err)**2)

	##linear best fit for finding lambda
	'''def lin(x, *p):
		return p[0]*x + p[1]

	def expanded_lin(x, *p):
		return p[0]*x + p[1]*x + p[2]

	p1, _ = curve_fit(lin, xdata, np.log(yFit1/popt1[0]), p0 = [.5, -.04])
	p2, _ = curve_fit(lin, xdata, np.log(yFit2/popt2[0]), p0 = [.5, -.04])
	p3, _ = curve_fit(expanded_lin, xdata, np.log(yFit3/popt3[0]), p0 = [.022, ])

	plt.figure(figsize = (14, 9))
	plt.plot(xdata, np.log(yFit1/popt1[0]), 'o', label = "Al Absorber/31 keV Energy")
	plt.plot(xdata, np.log(yFit2/popt2[0]), 'o', label = "Al Absorber/81 MeV Energy")
	plt.plot(xdata, np.log(yFit3/popt3[0]), 'o', label = "Al Absorber/356 keV Energy")
	plt.plot(fitxdata, lin(fitxdata, *p1), label = "511 keV Energy Fit")
	plt.plot(fitxdata, lin(fitxdata, *p2), label = "1.27 MeV Energy Fit")
	plt.plot(fitxdata, lin(fitxdata, *p3), label = "356 keV Energy Fit")
	plt.annotate(r"$\lambda_1 = %f \, mm^{-1}$" % np.absolute(p1[0]), xy = (xdata[3],lin(xdata[3], *p1)), 
		xytext = (xdata[3]+3,lin(xdata[3], *p1)), arrowprops = {"width":2, "frac":.3, "headwidth":7})
	plt.annotate(r"$\lambda_2 = %f \, mm^{-1}$" % np.absolute(p2[0]), xy = (xdata[3],lin(xdata[3], *p2)), 
		xytext = (xdata[3]+3,lin(xdata[3], *p2)), arrowprops = {"width":2, "frac":.3, "headwidth":7})
	plt.annotate(r"$\lambda_3 = %f \, mm^{-1}$" % np.absolute(p3[0]), xy = (xdata[3],lin(xdata[3], *p3)), 
		xytext = (xdata[3]+3,lin(xdata[3], *p3)), arrowprops = {"width":2, "frac":.3, "headwidth":7})
	plt.title("Fitting to find the Linear Attenuation Coefficient")
	plt.xlabel(r"$Thickness\,(mm)$")
	plt.ylabel(r"$ln(\frac{R}{R_0})$")
	plt.legend(loc = 0)
	plt.savefig("/users/aman/desktop/phys211/gamma cross sections/plots/linearcoeff_Na.pdf")
	plt.show()'''
	##back to the exponential fit

	yFit1 = expand_exponential(fitxdata, *popt1)
	yuFit1 = expand_exponential(fitxdata, *exp_p1)

	yFit2 = exponential(fitxdata, *popt2)
	yuFit2 = exponential(fitxdata, *exp_p2)

	yFit3 = exponential(fitxdata, *popt3)
	yuFit3 = exponential(fitxdata, *p)

	fig = plt.figure(figsize = (15, 10))
	fig.add_subplot(131)
	plt.errorbar(xdata, n1data, n1err, fmt = 'o', label = "Raw Data")
	plt.plot(fitxdata, yFit1, 
		linewidth = 2, alpha = .9, label = "Falling Exponential Fit")
	#plt.plot(fitxdata, yuFit1, 
	#	linewidth = 2, alpha = .9, label = "guesses")
	plt.text(20, 125,  r"$R(x) = R_0e^{-\lambda x} + R_0'e^{-\tau x} + C$"\
						"\n"\
					r"$R_0 = %.0f \pm %.1g \, counts/s$"\
						"\n"\
					r"$\lambda = %.4f \pm %.1g \, mm^{-1}$"\
						"\n"\
					r"$R_0' = %.0f \pm %.1g \, counts/s$"\
						"\n"\
					r"$\tau = %.3f \pm %.1g \, mm^{-1}$"\
						"\n"\
					r"$C = %.0f \pm %.1g \, counts/s$"\
						"\n"\
					r"$\chi^2 = %.2f $"\
						"\n"\
					r"$\frac{\chi^2}{\nu} = %.2f $"\
						"\n"\
					r"$\nu = %.0f$"\
					% (popt1[0], np.sqrt(pcov1[0,0]), np.absolute(popt1[1]), np.sqrt(pcov1[1,1]), 
						popt1[2], np.sqrt(pcov1[2,2]), np.absolute(popt1[3]), np.sqrt(pcov1[3,3]), popt1[4], np.sqrt(pcov1[4,4]),
						chisq1, chisq1/(len(xdata) - len(exp_p1)), len(xdata)-len(exp_p1)))
	plt.xlabel("Al Thickness")
	plt.ylabel("Countrate (counts/s)")
	plt.title("31 keV Transmission Intensity")
	plt.legend()

	fig.add_subplot(132)
	plt.errorbar(xdata, n2data, n2err, fmt = 'o', label = "Raw Data")
	plt.plot(fitxdata, yFit2, 
		linewidth = 2, alpha = .9, label = "Falling Exponential Fit")
	#plt.plot(fitxdata, yuFit2, 
	#	linewidth = 2, alpha = .9, label = "guesses")
	plt.text(20, 70,  r"$R(x) = R_0e^{-\lambda x} + C$"\
						"\n"\
					r"$R_0 = %.0f \pm %.1g \, counts/s$"\
						"\n"\
					r"$\lambda = %.4f \pm %.1g \, mm^{-1}$"\
						"\n"\
					r"$C = %.0f \pm %.1g \, counts/s$"\
						"\n"\
					r"$\chi^2 = %.2f $"\
						"\n"\
					r"$\frac{\chi^2}{\nu} = %.2f $"\
						"\n"\
					r"$\nu = %.0f$"\
					% (popt2[0], np.sqrt(pcov2[0,0]), np.absolute(popt2[1]), np.sqrt(pcov2[1,1]), 
						popt2[2], np.sqrt(pcov2[2,2]), chisq2, chisq2/(len(xdata) - len(exp_p2)), len(xdata)-len(exp_p2)))

	plt.xlabel("Al Thickness (mm)")
	plt.ylabel("Countrate (counts/s)")
	plt.title("81 keV Transmission Intensity")
	plt.legend()

	fig.add_subplot(133)
	plt.errorbar(xdata, n3data, n3err, fmt = 'o', label = "Raw Data")
	plt.plot(fitxdata, yFit3, 
		linewidth = 2, alpha = .9, label = "Falling Exponential Fit")
	plt.text(25, 100,  r"$R(x) = R_0e^{-\lambda x} + C$"\
						"\n"\
					r"$R_0 = %.0f \pm %.1g \, counts/s$"\
						"\n"\
					r"$\lambda = %.4f \pm %.1g \, mm^{-1}$"\
						"\n"\
					r"$C = %.0f \pm %.1g \, counts/s$"\
						"\n"\
					r"$\chi^2 = %.2f $"\
						"\n"\
					r"$\frac{\chi^2}{\nu} = %.2f $"\
						"\n"\
					r"$\nu = %.0f$"\
					% (popt3[0], np.sqrt(pcov3[0,0]), np.absolute(popt3[1]), np.sqrt(pcov3[1,1]), 
						popt3[2], np.sqrt(pcov3[2,2]), chisq3, chisq3/(len(xdata) - len(p)), len(xdata)-len(p)))

	plt.xlabel("Al Thickness (mm)")
	plt.ylabel("Countrate (counts/s)")
	plt.title("356 keV Transmission Intensity")
	plt.legend()

	plt.savefig("/users/aman/desktop/phys211/gamma cross sections/plots/%s_Ba.pdf" % dataset[0:dataset.find('.')])
	if show:
		plt.show()

if __name__ == '__main__':
	os.chdir("/users/aman/desktop/phys211/gamma cross sections/data/Ba_133")
	main("countrates.tsv", show = True)


