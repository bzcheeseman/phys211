import numpy as np
import matplotlib.pyplot as plt
import os

import sys
sys.path.append("/users/aman/desktop/PHYS211")

from PythonLibrary import *

try:
    %pylab inline
except Exception:
    pass


def compton(E, theta, m0):
    return E/(1+(E/(m0 * 3e8**2))(1-np.cos(theta)))

def plot(dataset, rng = [0,0]):
    try:
        data = np.genfromtxt(dataset, skip_header=25)

        xdata = data[:,0]
        ydata = data[:,1]
    except IndexError:
        data = np.genfromtxt(dataset, skip_header=25, delimiter=",", usecols=(0,2))

        xdata = data[:,0]
        ydata = data[:,1]

    if rng != [0,0]:
        
