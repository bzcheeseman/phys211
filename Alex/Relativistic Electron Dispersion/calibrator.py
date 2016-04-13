from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt


cal = np.array([1.7697, 1.2933, 0.031, 0.437])
for i in range(len(cal)):
    cal[i] = cal[i] * 1000 #switch to KeV
ch = np.array([1259, 922, 46, 305])
err = np.array([2,2,2,2])

def linear(p, x):
    return p[0] * x + p[1]
def residual(p, x, y, err):
    return (linear(p, x) - y) / err
    
p0 = np.array([3., 10.])

pf, cov, info, mesg, success = optimize.leastsq(residual, p0, args=(cal, ch, err), full_output=1, maxfev=100)
chisq = sum(info["fvec"]*info["fvec"])
dof = len(cal)-len(pf)
pferr = [np.sqrt(cov[i,i]) for i in range(len(pf))]

fig = plt.figure()
ax = plt.axes()
ax.errorbar(cal, ch, xerr=0., yerr=err, fmt='k.', label = 'Data')
T = np.linspace(cal.min(), cal.max(), 5000)
ax.plot(T, linear(pf, T), 'r-', label = 'calibration')

ax.set_title('Channel to Energy calibration')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Channel Number')
ax.legend(loc=(0.7,0.5))
    
textfit = '$Ch = A * energy + B$ \n' \
          '$A = %.3f \pm %.3f$ keV \n' \
          '$B = %.0f \pm %.0f$ keV \n' \
          '$\chi^2= %.2f$ \n' \
          '$N = %i$ (dof) \n' \
          '$\chi^2/N = % .2f$' \
           % (pf[0], pferr[0], pf[1], pferr[1], chisq, dof,
              chisq/dof)
ax.text(0.15, .8, textfit, transform=ax.transAxes, fontsize=12,
         verticalalignment='top')
         
plt.show()

A = 1/pf[0]
Aerr = pferr[0]/pf[0] * (1 / pf[0])
B = 1/pf[1]
Berr = pferr[1]/pf[1] * (1 / pf[1])
print A, Aerr
print B, Berr
