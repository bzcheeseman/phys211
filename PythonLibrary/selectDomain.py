__author__ = "Aman LaChapelle"

import numpy as np

# Searches through the x data to find a region to keep and returns the same size
# slice in both x and y
def selectdomain(xdata, ydata, domain = None, indices = None, remove = None):
    if indices == None and domain != None:
        indices = np.searchsorted(xdata,domain)

    if remove is not None:
        indices_first = np.searchsorted(xdata, [domain[0], remove[0]])
        indices_second = np.searchsorted(xdata, [remove[1], domain[1]])

        return np.append(xdata[indices_first[0]:indices_first[1]], xdata[indices_second[0]:indices_second[1]]),np.append(ydata[indices_first[0]:indices_first[1]], ydata[indices_second[0]:indices_second[1]])

    return xdata[indices[0]:indices[1]], ydata[indices[0]:indices[1]]
