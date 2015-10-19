import matplotlib.pyplot as plt
import numpy as np
from scipy import loadtxt
from scipy.special import erf
import sys

def main(dataset):
	#xdata, ydata = loadtxt(dataset, unpack = True, usecols=[0,1], skiprows = 25)

	xdata = [31, 81, 356, 511, 1270]
	ydata1 = [2.94, .513, .24, .19, .087]
	ydata2 = np.multiply(2.643,[1.128E+00, 2.018E-01, ((1.042E-01+9.276E-02)/2.0), 8.445E-02, 5.496E-02])
	err1 = [.1, .06, .04, .03, .05]
	err2 = [.03, .03, .1, .05, .07]

	delta = np.absolute(np.subtract(ydata1, ydata2))
	sigma = np.sqrt(np.multiply(err1, err1) + np.multiply(err2, err2))

	t = np.divide(delta, sigma)

	found_erf = erf(t)

	confidence = np.subtract(1, found_erf) * 100

	plt.figure(figsize = (7,7))
	plt.errorbar(xdata, ydata1, err1, fmt = 'o', label = "$\lambda_{exp}$")
	plt.errorbar(xdata, ydata2, err2, fmt = 'o', label = "$\lambda_{calc}$")
	plt.text(31, 3.2, " %.1f%%" % confidence[0])
	plt.text(81, .7, " %.1f%%" % confidence[1])
	plt.text(356, .5, " %.1f%%" % confidence[2])
	plt.text(511, .4, " %.1f%%" % confidence[3])
	plt.text(1270, .3, " %.1f%%" % confidence[4])
	plt.xlabel("Energy (keV)")
	plt.ylabel(r"$\lambda \, (cm^{-1}) $")
	plt.title("Comparison to Literature")
	plt.legend()
	plt.savefig("/users/aman/desktop/phys211/gamma cross sections/plots/confidence.pdf")
	plt.show()

if __name__ == '__main__':
	import sys

	main(sys.argv[1])