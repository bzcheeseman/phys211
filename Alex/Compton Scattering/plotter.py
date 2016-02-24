from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt


def read_data():
    datafile = r'data/data.csv'
    data = np.genfromtxt(datafile, delimiter=',', skiprows=1)
    return data
    
data = read_data()

def linear(p, x):
    return p[0] * x + p[1]
def residual(p, x, y, err):
    return (linear(p, x) - y) / err

    
    
#constants
mc2 = 511. #keV
E = 662. #keV, energy of initial gamma photons

#Day 1 
angles1 = data[:,2] #x
cent1 = data[:,1] #y
#Calibration
e1 = np.array([81., 356., 662.])
c1 = np.array([118., 505., 921.]) # +/- 3
c1err = np.array([3., 3., 3.])
#Fitting
# Channel = A1 * Energy + A2
p1 = [0., 0.]
pf1, cov1, info1, mesg1, success1 = optimize.leastsq(residual, p1, args=(e1, c1, c1err), full_output=1)
chisq1 = sum(info1["fvec"]*info1["fvec"])
dof1 = len(e1)-len(pf1)
pferr1 = [np.sqrt(cov1[i,i]) for i in range(len(pf1))]
print 'pf1', pf1, '\n', cov1
A1 = 1 / pf1[0]
A1err = (pferr1[0] / pf1[0]) * (1 / pf1[0])
A2 = - pf1[1] / pf1[0]
A2err = np.sqrt(pferr1[0] * (2 * pf1[1] / pf1[0]*pf1[0]) + pferr1[1] * (1 / pf1[0]))
A = np.array([A1, A2])
Aerr = np.array([A1err, A2err])
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1 = plt.axes()
ax1.errorbar(c1, e1, xerr=c1err, yerr=0.0, fmt='k.', label = 'Data')
T = np.linspace(c1.min(), c1.max(), 5000)
ax1.plot(T, linear(A, T), 'r-', label = 'calibration')

ax1.set_title('Channel to Energy calibration - Day 1')
ax1.set_xlabel('Channel Number')
ax1.set_ylabel('Energy (keV)')
ax1.legend(loc=(0.7,0.5))
    
textfit = '$E = A_1 * channel + A_2$ \n' \
          '$A_1 = %.2f \pm %.2f$ keV \n' \
          '$A_2 = %.0f \pm %.0f$ keV \n' \
          '$\chi^2= %.2f$ \n' \
          '$N = %i$ (dof) \n' \
          '$\chi^2/N = % .2f$' \
           % (pf1[0], pferr1[0], pf1[1], pferr1[1], chisq1, dof1,
              chisq1/dof1)
ax1.text(0.15, .8, textfit, transform=ax1.transAxes, fontsize=12,
         verticalalignment='top')
             

plt.show()

energies = linear(A, cent1)
enerr = (A1err/A1) * cent1 + A2err

print angles1
print energies
print enerr

#m2, b2 = np.polyfit(c2, e1, 1)
#
#def fit2(x):
#    return m2 * x + b2
#en2 = fit2(c2)
#
#
#
### DATA is angles1, angles2, en1, en2
#y = np.append(angles1, angles2)
#x = np.append(en1, en2)
#yerr = [1 * np.pi / 180] * len(y)
#
#
#def func(p,x):
#    return np.arccos(p[0]/x)
#
#p = [1 - (E/(1 + E/mc2))] #expected
#
#popt, pcov = curve_fit(func, x, y, p, yerr, maxfev=int(2e6))