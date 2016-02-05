__author__ = "Aman LaChapelle"

import numpy as np

# Searches through the x data to find a region to keep and returns the same size
# slice in both x and y
def selectdomain(xdata, ydata, domain):
    indices = np.searchsorted(xdata,domain)
    return xdata[indices[0]:indices[1]],ydata[indices[0]:indices[1]]
