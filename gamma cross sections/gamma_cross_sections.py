import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import loadtxt
import os

#write for loop to iterate over each peak?  Apply gaussian filter to find peaks and then just fit each one?
#adjust count rate for duration of sweep also

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
	with open(dataset, 'r') as f:
		f.seek(270)
		time = float(f.read(6))

		f.close()

	xdata, ydata = loadtxt(dataset, unpack=True, usecols=[0,1], skiprows=25)

	sel_xdata_1, sel_ydata_1 = choose_domain(xdata, ydata, [215, 307])
	sel_xdata_2, sel_ydata_2 = choose_domain(xdata, ydata, [610, 700])
	#sel_xdata_3, sel_ydata_3 = choose_domain(xdata, ydata, [610, 760])

	error = np.sqrt(ydata)

	err1 = np.ones_like(sel_ydata_1) #ask about these tomorrow
	err2 = np.ones_like(sel_ydata_2)
	#err3 = np.ones_like(sel_ydata_3)

	for i in range(0, len(sel_ydata_1)):
		if sel_ydata_1[i] == 0:
			err1[i] = 1
		else:
			err1[i] = np.sqrt(sel_ydata_1[i])

	for i in range(0, len(sel_ydata_2)):
		if sel_ydata_2[i] == 0:
			err2[i] = 1
		else:
			err2[i] = np.sqrt(sel_ydata_2[i])

	#for i in range(0, len(sel_ydata_3)):
	#	if sel_ydata_3[i] == 0:
	#		err3[i] = 1
	#	else:
	#		err3[i] = np.sqrt(sel_ydata_3[i])


	def gaussian(x, *p):
		return p[0]/(p[1]*np.sqrt(2*np.pi)) * np.exp(-(x-p[2])**2/(2*p[1]**2)) + p[3] + p[4]*x

	def red_gaussian(x, *p):
		return p[0]/(p[1]*np.sqrt(2*np.pi)) * np.exp(-(x-p[2])**2/(2*p[1]**2)) + p[3]

	p1 = [48810, 9, 266, 433, -1]
	p2 = [9303, 14, 650, 2]
	#p3 = [400**2, 19, 675, 0, 0]

	popt1, pcov1 = curve_fit(gaussian, sel_xdata_1, sel_ydata_1, p0 = p1, sigma = err1)
	popt2, pcov2 = curve_fit(red_gaussian, sel_xdata_2, sel_ydata_2, p0 = p2, sigma = err2)
	#popt3, pcov3 = curve_fit(gaussian, sel_xdata_3, sel_ydata_3, p0 = p3, sigma = err3)

	yFit_1 = gaussian(sel_xdata_1, *popt1)
	yuFit_1 = gaussian(sel_xdata_1, *p1)
	
	yFit_2 = red_gaussian(sel_xdata_2, *popt2)
	yuFit_2 = red_gaussian(sel_xdata_2, *p2)

	#yFit_3 = gaussian(sel_xdata_3, *popt3)
	#yuFit_3 = gaussian(sel_xdata_3, *p3)

	chisq1 = np.sum(((yFit_1-sel_ydata_1)/err1)**2)
	chisq2 = np.sum(((yFit_2-sel_ydata_2)/err2)**2)
	#chisq3 = np.sum(((yFit_3-sel_ydata_3)/err3)**2)

	plt.figure(figsize = (13, 10))
	plt.errorbar(xdata, ydata, error, fmt = 'o', label = "Raw Data")
	plt.plot(sel_xdata_1, yFit_1, 
		linewidth = 2, alpha = .9, label = "511 keV Fit")
	plt.plot(sel_xdata_2, yFit_2, linewidth = 2, alpha = .9, label = "1270 keV Fit")
	#plt.plot(sel_xdata_3, yFit_3, linewidth = 2, alpha = .9, label = "356 keV Fit")
	plt.text(500, 700,  r"$f_1(x) = \frac{N_1}{\sigma_1\sqrt{2\pi}}e^{\frac{-(x-\mu_1)^2}{2\sigma_1^2}} + C + Bx$"\
						"\n"\
					r"$N_1 = %.0f \pm %.1g \, counts$"\
						"\n"\
					r"$\sigma_1 = %.2f \pm %.1g$"\
						"\n"\
					r"$\mu_1 = %.2f \pm %.1g$"\
						"\n"\
					r"$C = %.0f \pm %.1g \, counts$"\
						"\n"\
					r"$B = %.0f \pm %.1g \, counts/channel$"\
						"\n"\
					r"$\chi^2 = %f $"\
						"\n"\
					r"$\frac{\chi^2}{\nu} = %f $"\
						"\n"\
					r"$\nu = %d$"\
					% (popt1[0], np.sqrt(pcov1[0,0]), popt1[1], np.sqrt(pcov1[1,1]), popt1[2], np.sqrt(pcov1[2,2]), 
						popt1[3], np.sqrt(pcov1[3,3]), popt1[4], np.sqrt(pcov1[4,4]), chisq1, chisq1/(len(sel_xdata_1) - len(p1)), len(sel_xdata_1) - len(p1)))
	plt.text(800, 700,  r"$f_2(x) = \frac{N_2}{\sigma_2\sqrt{2\pi}}e^{\frac{-(x-\mu_2)^2}{2\sigma_2^2}} + C$"\
						"\n"\
					r"$N_2 = %.0f \pm %.1g \, counts$"\
						"\n"\
					r"$\sigma_2 = %.1f \pm %.1g$"\
						"\n"\
					r"$\mu_2 = %.1f \pm %.1g$"\
						"\n"\
					r"$C = %.0f \pm %.1g \, counts$"\
						"\n"\
					r"$\chi^2 = %f $"\
						"\n"\
					r"$\frac{\chi^2}{\nu} = %f $"\
						"\n"\
					r"$\nu = %d$"\
					% (popt2[0], np.sqrt(pcov2[0,0]), popt2[1], np.sqrt(pcov2[1,1]), popt2[2], np.sqrt(pcov2[2,2]), 
						popt2[3], np.sqrt(pcov2[3,3]), chisq2, chisq2/(len(sel_xdata_2) - len(p2)), len(sel_xdata_2) - len(p2)))
	#plt.text(850, 600,  r"$f_3(x) = \frac{N_3}{\sigma_3\sqrt{2\pi}}e^{\frac{-(x-\mu_3)^2}{2\sigma_3^2}} + C + Bx$"\
	#					"\n"\
	#				r"$N_3 = %.0f \pm %.1g \, counts$"\
	#					"\n"\
	#				r"$\sigma_3 = %.1f \pm %.1g$"\
	#					"\n"\
	#				r"$\mu_3 = %.1f \pm %.1g$"\
	#					"\n"\
	#				r"$C = %.0f \pm %.1g \, counts$"\
	#					"\n"\
	#				r"$B = %.0f \pm %.1g \, counts/channel$"\
	#					"\n"\
	#				r"$\chi^2 = %f $"\
	#					"\n"\
	#				r"$\frac{\chi^2}{\nu} = %f $"\
	#					"\n"\
	#				r"$\nu = %d$"\
	#				% (popt3[0], np.sqrt(pcov3[0,0]), popt3[1], np.sqrt(pcov3[1,1]), popt3[2], np.sqrt(pcov3[2,2]), 
	#					popt3[3], np.sqrt(pcov3[3,3]), popt3[4], np.sqrt(pcov3[4,4]), chisq3, chisq3/(len(sel_xdata_3) - len(p3)), len(sel_xdata_3) - len(p3)))
	plt.text(823, 350, r"$\frac{N_1}{t} = %f \, counts/sec$"\
							"\n"\
						r"$\frac{N_2}{t} = %f \, counts/sec$"\
						#	"\n"\
						#r"$\frac{N_3}{t} = %f \, counts/sec$"\
						% (popt1[0]/time, popt2[0]/time))

	#dist = dataset[2:dataset.find("m")]

	#with open("countrates.tsv", "a+") as f:
	#	f.write(str(dist)+"\t"+str(popt1[0]/time)+"\t"+str(popt2[0]/time)+"\t"+str(popt3[0]/time)+"\t"+
	#		str(np.sqrt(pcov1[0,0])/popt1[0])+"\t"+str(np.sqrt(pcov2[0,0])/popt2[0])+"\t"\
	#		+str(np.sqrt(pcov3[0,0])/popt3[0])+"\n")
	#f.close()

	plt.xlabel("Channel")
	plt.ylabel("Counts")
	plt.title("PHA Spectrum of Na-22 Energy Decay ~ %s" % dataset[2:dataset.find('.')])
	plt.legend(loc = 4)
	plt.savefig("/users/aman/desktop/phys211/gamma cross sections/plots/%s.pdf" % dataset[0:dataset.find('.')])
	if show:
		plt.show()

if __name__ == '__main__':
	os.chdir("/users/aman/desktop/phys211/gamma cross sections/data/Na_22/run_2")
	#with open("countrates.tsv", "w+") as f:
	#	f.write("Thickness in mm"+"\t"+"Countrate N1"+"\t"+"Countrate N2"+"\t"+"Countrate N3"\
	#		+"\t"+"% Error in N1"+"\t"+"% Error in N2"+"\t"+"% Error in N3"+"\n")
	#f.close()

	#for i in range(0, 60):
	try:
		main("Na%dmmAl.tsv" % 8, show = True)
	except IOError:
		pass





