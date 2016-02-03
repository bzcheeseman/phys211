%pylab inline

import numpy as np
import matplotlib.pyplot as plt
from scipy import loadtxt

def plot_data(dataset):
    t, ch1, ch2 = loadtxt(dataset, skiprows=1, usecols=(0,1,3))

    plt.plot(time, ch1)





if __name__ == '__main__':
    plot_data("data/doppler_broadened.tsv")
