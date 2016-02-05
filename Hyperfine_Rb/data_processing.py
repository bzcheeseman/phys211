%pylab inline

import numpy as np
import matplotlib.pyplot as plt
from scipy import loadtxt

def selectdomain(xdata,ydata,domain):
    ind=np.searchsorted(xdata,domain)
    return xdata[ind[0]:ind[1]],ydata[ind[0]:ind[1]]

def plot_data(dataset):
    t, ch1, ch2 = loadtxt(dataset, unpack = True, skiprows=1, usecols=(0,1,3))



    plt.plot(t, ch1)





if __name__ == '__main__':
    plot_data("data/doppler_broadened.tsv")
