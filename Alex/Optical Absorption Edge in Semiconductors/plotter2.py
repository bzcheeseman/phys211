from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt


V0 = 0.66
V0_err = 0.02
V100 = 24.4
V100_err = 0.06
def read_data(datafile):
    data = np.genfromtxt(datafile, delimiter=',', skiprows=1)
    return data
    
data = read_data(r'data/Ge.csv')
x = (6.626*pow(10,-34)*3*pow(10,8))/(data[:,0]*pow(10,-9))
T = (data[:,1] - V0)/(V100 - V0)
T_err = (data[:,2] + V0_err)/(data[:,1] - V0) + (V100_err + V0_err)/(V100 - V0)
T_err *= T
print T
print T_err

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1 = plt.axes()
ax1.errorbar(x, T*100, yerr=T_err*100, fmt='k.', label = 'Data')
ax1.set_title('Transmittance vs Energy for Ge')
ax1.set_xlabel('Energy (J)')
ax1.set_ylabel('Transmittance %')
plt.savefig('Ge.png')
plt.show()