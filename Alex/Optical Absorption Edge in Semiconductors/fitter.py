from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt


V0 = 0.66
V0_err = 0.02
V100 = 24.4
V100_err = 0.06

R = 0.31
R_err = 0.01
d = 0.000002
d_err = 0.0

def read_data(datafile):
    data = np.genfromtxt(datafile, delimiter=',', skiprows=1)
    return data
    
data = read_data(r'data/GaAs.csv')
x = (6.626*pow(10,-34)*3*pow(10,8))/(data[:,0]*pow(10,-9))
T = (data[:,1] - V0)/(V100 - V0)
T_err = (data[:,2] + V0_err)/(data[:,1] - V0) + (V100_err + V0_err)/(V100 - V0)
T_err *= T
print(zip(x,T,T_err))[11]

def alpha(R,T,d):
    return (1/d) * np.log(((1-R)*(1-R) + np.sqrt(pow((1-R),4) + 4*T*T*R*R))/(2*T))
    
def alpha_err(a,R,T,d,R_err,T_err,d_err):
    return np.sqrt( pow((((1-R)*(1-R)*T_err)/(T*d*np.sqrt(pow((1-R),4) + 4*R*R*T*T))),2) + pow(((2*(T*np.exp(a*d)+1-R)*R_err)/(R*d*np.sqrt(pow((1-R),4)+4*R*R*T*T))),2) + pow((a*d_err/d),2))

a = alpha(R, T, d)
a_err = alpha_err(a, R, T, d, R_err, T_err, d_err)
a_errf = alpha_err(a, R, T, d, R_err, T_err, d_err)/a
a = a*a
a_err = a_errf * a

a = np.sort(a)[-15:]
a_err = np.sort(a_err)[-15:]
x = np.sort(x)[-15:]

a = np.sort(a)[:12]
a_err = np.sort(a_err)[:12]
x = np.sort(x)[:12]

def linear(p, x):
    return p[0] * (x - p[1])
def residual(p, x, y, err):
    return (linear(p, x) - y) / err

a /= pow(1000000,2)
a_err /= pow(1000000,2)

p1 = [2.*pow(10,19), 0.]
pf1, cov1, info1, mesg1, success1 = optimize.leastsq(residual, p1, args=(x, a, a_err), full_output=1)
chisq1 = sum(info1["fvec"]*info1["fvec"])
dof1 = len(x)-len(pf1)
pferr1 = [np.sqrt(cov1[i,i]) for i in range(len(pf1))]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1 = plt.axes()
ax1.errorbar(x, a, yerr=a_err, fmt='k.', label = 'Data')
T = np.linspace(x.min(), x.max(), 5000)
ax1.plot(T, linear(pf1, T), 'r-', label = 'Fit')
ax1.set_title('Absorption Coefficient against Energy for GaAs')
ax1.set_xlabel('Energy (J)')
ax1.set_ylabel('alpha /um')
textfit = '$alpha^2 = A * (E - E_g)$ \n' \
          '$A = (%.2f \pm %.2f) * 10^{19} J$ \n' \
          '$E_g = (%.3f \pm %.3f) * 10^{-19} J$ \n' \
          '$\chi^2= %.2f$ \n' \
          '$N = %i$ (dof) \n' \
          '$\chi^2/N = % .2f$' \
           % (pf1[0]/pow(10,19), pferr1[0]/pow(10,19), pf1[1]*pow(10,19), pferr1[1]*pow(10,19), chisq1, dof1,
              chisq1/dof1)
ax1.text(0.05, .8, textfit, transform=ax1.transAxes, fontsize=12,
         verticalalignment='top')

plt.legend(loc=0)
plt.show()
plt.savefig('GaAs_alpha_fit.png')