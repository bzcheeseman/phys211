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
e1 = np.array([81., 356., 662.])
c1err = np.array([3., 3., 3.])

##Day 2
angles2 = data[:,6]
angles2 = np.array(np.delete(angles2, len(angles2)-1))
cent2 = data[:,5]
cent2 = np.array(np.delete(cent2, len(cent2)-1))
##Calibration
c21 = [117., 499., 909.]
c22 = [116., 498., 903.]
c2 = [0]*3
for i in range(len(c21)):
    c2[i] = (c21[i] + c22[i]) / 2
c2 = np.asarray(c2)

#Fitting
# Channel = A1 * Energy + A2
p2 = np.array([0., 0.])
pf2, cov2, info2, mesg2, success2 = optimize.leastsq(residual, p2, args=(e1, c2, c1err), full_output=1)
chisq2 = sum(info2["fvec"]*info2["fvec"])
dof2 = len(e1)-len(pf2)
pferr2 = [np.sqrt(cov2[i,i]) for i in range(len(pf2))]
print 'pf2', pf2, '\n', cov2
B1 = 1 / pf2[0]
B1err = (pferr2[0] / pf2[0]) * (1 / pf2[0])
B2 = - pf2[1] / pf2[0]
B2err = np.sqrt(pferr2[0] * (2 * pf2[1] / pf2[0]*pf2[0]) + pferr2[1] * (1 / pf2[0]))
B = np.array([B1, B2])
Berr = np.array([B1err, B2err])
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2 = plt.axes()
ax2.errorbar(c2, e1, xerr=c1err, yerr=0.0, fmt='k.', label = 'Data')
T = np.linspace(c2.min(), c2.max(), 5000)
ax2.plot(T, linear(B, T), 'r-', label = 'calibration')

ax2.set_title('Channel to Energy calibration - Day 2')
ax2.set_xlabel('Channel Number')
ax2.set_ylabel('Energy (keV)')
ax2.legend(loc=(0.7,0.5))
    
textfit = '$E = B_1 * channel + B_2$ \n' \
          '$B_1 = %.2f \pm %.2f$ keV \n' \
          '$B_2 = %.0f \pm %.0f$ keV \n' \
          '$\chi^2= %.2f$ \n' \
          '$N = %i$ (dof) \n' \
          '$\chi^2/N = % .2f$' \
           % (pf2[0], pferr2[0], pf2[1], pferr2[1], chisq2, dof2,
              chisq2/dof2)
ax2.text(0.15, .8, textfit, transform=ax2.transAxes, fontsize=12,
         verticalalignment='top')
         
plt.show()

energies = linear(B, cent2)
enerr = (B1err/B1) * cent2 + B2err
print angles2
print energies
print enerr