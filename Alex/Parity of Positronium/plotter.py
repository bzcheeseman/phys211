from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

def read_data():
    datafile = r'data/calibration1.csv'
    data = np.genfromtxt(datafile, delimiter=',', skiprows=1)
    return data
    
data = read_data()