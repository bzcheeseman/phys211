%pylab inline

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import sys

sys.path.append("..")

from PythonLibrary.residuals import *
from PythonLibrary.selectDomain import selectdomain

def calibration(plot_cal):
    data = np.genfromtxt("data/na_22_cal.tsv", skip_header=26, usecols=(0,2))

    peak_1_x, peak_1_y = selectdomain(data[:,0], data[:,1], [230, 400])

    p1 = [1e6, 50, 315, -10, 200]

    popt1, pcov1 = curve_fit(gaussian_back, peak_1_x, peak_1_y, p0 = p1)

    yFit1 = gaussian_back(peak_1_x, *popt1)
    yuFit1 = gaussian_back(peak_1_x, *p1)

    print popt1[2]

    peak_2_x, peak_2_y = selectdomain(data[:,0], data[:,1], [675, 875])

    p2 = [1e6, 100, 777, 0]

    popt2, pcov2 = curve_fit(gaussian, peak_2_x, peak_2_y, p0 = p2)

    yFit2 = gaussian(peak_2_x, *popt2)
    yuFit2 = gaussian(peak_2_x, *p2)

    print popt2[2]

    data_cs = np.genfromtxt("data/cs_137_cal.tsv", skip_header=26, usecols=(0,2))

    peak_3_x, peak_3_y = selectdomain(data_cs[:,0], data_cs[:,1], [325, 500])

    p3 = [1e3, 50, 500, -10, 200]

    popt3, pcov3 = curve_fit(gaussian_back, peak_3_x, peak_3_y, p0 = p3)

    yFit3 = gaussian_back(peak_3_x, *popt3)
    yuFit3 = gaussian_back(peak_3_x, *p3)

    print popt3[2]

    if plot_cal:
        plt.figure(figsize=(10,10))
        plt.plot(data[:,0], data[:,1], 'o')
        # plt.plot(peak_1_x, yuFit1)
        plt.plot(peak_1_x, yFit1)
        # plt.plot(peak_2_x, yuFit2)
        plt.plot(peak_2_x, yFit2)
        plt.figure(figsize=(10, 10))
        plt.plot(data_cs[:,0], data_cs[:,1])
        # plt.plot(peak_3_x, yuFit3)
        plt.plot(peak_3_x, yFit3)

    cal_line_x = sort([popt1[2], popt2[2], popt3[2]])
    cal_line_y = sort([511, 1275, 662])



    plt.plot(cal_line_x, cal_line_y, 'o')












if __name__ == '__main__':
    calibration(False)
