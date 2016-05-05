from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

def read_data():
    datafile = r'data/calibration1.csv'
    data = np.genfromtxt(datafile, delimiter=',', skiprows=1)
    return data
    
def read_data2():
    datafile = r'data/calibration2.csv'
    data = np.genfromtxt(datafile, delimiter=',', skiprows=1)
    return data
    
data = read_data()
data2 = read_data2()

d1data0 = data[:6]
d1data100 = data[-6:]
d2data0 = data2[-6:]
d2data100 = data2[:6]

v0 = (np.mean(d1data0[:,2]) + np.mean(d2data0[:,2]))/2
v0_err = (np.sum(d1data0[:,3]) + np.sum(d2data0[:,3]))/12

v100 = (np.mean(d1data100[:,2]) + np.mean(d2data100[:,2]))/2
v100_err = (np.sum(d1data100[:,3]) + np.sum(d2data100[:,3]))/12

print v0
print v0_err

print v100
print v100_err
