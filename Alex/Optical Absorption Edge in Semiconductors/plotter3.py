from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt


V0 = 0.66
V0_err = 0.02
V100 = 24.4
V100_err = 0.06

R = 0.37
R_err = 0.01
d = 0.00111
d_err = 0.000005

def read_data(datafile):
    data = np.genfromtxt(datafile, delimiter=',', skiprows=1)
    return data
    
data = read_data(r'data/Ge.csv')
x = (6.626*pow(10,-34)*3*pow(10,8))/(data[:,0]*pow(10,-9))
T = (data[:,1] - V0)/(V100 - V0)
T_err = (data[:,2] + V0_err)/(data[:,1] - V0) + (V100_err + V0_err)/(V100 - V0)
T_err *= T

def alpha(R,T,d):
    return (1/d) * np.log(((1-R)*(1-R) + np.sqrt(pow((1-R),4) + 4*T*T*R*R))/(2*T))
    
def alpha_err(a,R,T,d,R_err,T_err,d_err):
    return np.sqrt( pow((((1-R)*(1-R)*T_err)/(T*d*np.sqrt(pow((1-R),4) + 4*R*R*T*T))),2) + pow(((2*(T*np.exp(a*d)+1-R)*R_err)/(R*d*np.sqrt(pow((1-R),4)+4*R*R*T*T))),2) + pow((a*d_err/d),2))

a = alpha(R, T, d)
a_err = alpha_err(a, R, T, d, R_err, T_err, d_err)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1 = plt.axes()
ax1.errorbar(x, a/1000, yerr=a_err/1000, fmt='k.', label = 'Data')
ax1.set_title('Absorption Coefficient against Energy for Ge')
ax1.set_xlabel('Energy (J)')
ax1.set_ylabel('alpha /mm')
plt.savefig('Ge_alpha.png')
plt.show()