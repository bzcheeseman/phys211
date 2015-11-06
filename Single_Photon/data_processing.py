import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def main():
	xdata = np.linspace(0, 1, 35)

	plt.plot(xdata, 2*np.ones_like(xdata), label = "$R_{acc}$")
	plt.fill_between(xdata, 0.36+2*np.ones_like(xdata), -0.36+2*np.ones_like(xdata), alpha = .3, 
		label = "$\delta R_{acc}$")
	plt.errorbar([.3], [0.9], [0.1], label = "$R_1$")
	plt.errorbar([.7], [1.6], [0.1], label = "$R_2$")
	plt.legend(loc = 3)
	plt.savefig("plots/overlaps.pdf")
	plt.show()














if __name__ == '__main__':
	main()

