%pylab inline

import numpy as np
import matplotlib.pyplot as plt

def main():
    data = np.genfromtxt("data/doppler_broadened.tsv")

    xdata = data[:,0]
    ydata = data[:,1]

    plt.plot(xdata, ydata)





if __name__ == '__main__':
    main()
