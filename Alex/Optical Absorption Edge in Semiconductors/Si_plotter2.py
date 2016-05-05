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

d = 0.00114
d_err = 0.000005
d2 = 0.00209
d3 = 0.00315
R = 0.33
R_err = 0.01
R2 = 0.34
R3 = 0.35

def alpha(R,T,d):
    return (1/d) * np.log(((1-R)*(1-R) + np.sqrt(pow((1-R),4) + 4*T*T*R*R))/(2*T))
    
def alpha_err(a,R,T,d,R_err,T_err,d_err):
    return np.sqrt( pow((((1-R)*(1-R)*T_err)/(T*d*np.sqrt(pow((1-R),4) + 4*R*R*T*T))),2) + pow(((2*(T*np.exp(a*d)+1-R)*R_err)/(R*d*np.sqrt(pow((1-R),4)+4*R*R*T*T))),2) + pow((a*d_err/d),2))

a = alpha(R, T, d)
a2 = alpha(R2, T2, d2)
a3 = alpha(R3, T3, d3)

a_err = alpha_err(a, R, T, d, R_err, T_err, d_err)
a_err2 = alpha_err(a2, R2, T2, d2, R_err, T_err2, d_err)
a_err3 = alpha_err(a3, R3, T3, d3, R_err, T_err3, d_err)

fig1 = plt.figure()
ax1 = fig1.add_subplot(131)
ax1 = plt.axes()
ax1.errorbar(x, a/1000, yerr=a_err/1000, fmt='k.', label = '1.14mm')
ax1.set_title('Transmittance vs Wavelength for Si - 1.14mm')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Transmittance Factor')

ax2 = fig1.add_subplot(132)
ax2 = plt.axes()
ax2.errorbar(x2, a2/1000, yerr=a_err2/1000, fmt='b.', label = '2.09mm')
ax2.set_title('Transmittance vs Wavelength for Si')
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Transmittance Factor')

ax3 = fig1.add_subplot(133)
ax3 = plt.axes()
ax3.errorbar(x3, a3/1000, yerr=a_err3/1000, fmt='r.', label = '3.15mm')
ax3.set_title('Absorption Coefficient against Energy for Si')
ax3.set_xlabel('Energy (J)')
ax3.set_ylabel('alpha /mm')

plt.legend(loc=0)


plt.savefig('Si_alpha.png')
plt.show()