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
    
data = read_data(r'data/Si_1.14mm.csv')
data2 = read_data(r'data/Si_2.09mm.csv')
data3 = read_data(r'data/Si_3.15mm.csv')
x = (6.626*pow(10,-34)*3*pow(10,8))/(data[:,0]*pow(10,-9))
x2 = (6.626*pow(10,-34)*3*pow(10,8))/(data2[:,0]*pow(10,-9))
x3 = (6.626*pow(10,-34)*3*pow(10,8))/(data3[:,0]*pow(10,-9))

T = (data[:,1] - V0)/(V100 - V0)
T_err = (data[:,2] + V0_err)/(data[:,1] - V0) + (V100_err + V0_err)/(V100 - V0)
T_err *= T

T2 = (data2[:,1] - V0)/(V100 - V0)
T_err2 = (data2[:,2] + V0_err)/(data2[:,1] - V0) + (V100_err + V0_err)/(V100 - V0)
T_err2 *= T2

T3 = (data3[:,1] - V0)/(V100 - V0)
T_err3 = (data3[:,2] + V0_err)/(data3[:,1] - V0) + (V100_err + V0_err)/(V100 - V0)
T_err3 *= T3

fig1 = plt.figure()
ax1 = fig1.add_subplot(131)
ax1 = plt.axes()
ax1.errorbar(x, T*100, yerr=T_err*100, fmt='k.', label = '1.14mm')
ax1.set_title('Transmittance vs Wavelength for Si - 1.14mm')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Transmittance Factor')

ax2 = fig1.add_subplot(132)
ax2 = plt.axes()
ax2.errorbar(x2, T2*100, yerr=T_err2*100, fmt='b.', label = '2.09mm')
ax2.set_title('Transmittance vs Wavelength for Si')
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Transmittance Factor')

ax3 = fig1.add_subplot(133)
ax3 = plt.axes()
ax3.errorbar(x3, T3*100, yerr=T_err3*100, fmt='r.', label = '3.15mm')
ax3.set_title('Transmittance vs Wavelength for Si')
ax3.set_xlabel('Energy (J)')
ax3.set_ylabel('Transmittance %')

plt.legend(loc=0)


plt.savefig('Si.png')
plt.show()