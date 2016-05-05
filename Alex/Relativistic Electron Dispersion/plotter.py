from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

a = 1.42372210086
aerr = 0.00295712984228
b = 0.0770992785753
berr = 0.00969212354148

cedges = np.array([338,234,749,279,614,186,634,773])
cedgerr = np.array([3,3,3,3,3,4,3,4])
cpeaks = np.array([469,353,900,405,758,296,782,922])
cperr = np.array([2] * len(cpeaks))

def c2e(cs, errs): #cs is an array of channels
    es = a * cs + b
    errs = (aerr/a + errs/cs) * cs + berr
    return es, errs
    
T, Terrs = c2e(cedges, cedgerr)
peaks, peakerr = c2e(cpeaks, cperr) #peak = initial photon Energy


pc = np.array(2*peaks - T)
pcerr = np.array(2*peakerr + Terrs)

ys = np.array(pc*pc / (2 * T))
yserr = np.array(np.sqrt((2*pc/(2*T))*(2*pc/(2*T))*pcerr*pcerr + (pc*pc*2/(T*T))*(pc*pc*2/(T*T))*Terrs*Terrs))
print ys, T, yserr


### FIT AND PLOT (PC)^2/2T (ys) against T

def poly(p, x):
    return p[0]*(x) + p[1]
    
def residual(p, x, y, err):
    return (poly(p, x) - y) / err

p0 = np.array([1.,1.])

pf, cov, info, mesg, success = optimize.leastsq(residual, p0, args=(T, ys, yserr), full_output=1, maxfev=1000)
print pf
chisq = sum(info["fvec"]*info["fvec"])
print chisq
dof = len(ys)-len(pf)
pferr = [np.sqrt(cov[i,i]) for i in range(len(pf))]

fig = plt.figure()
ax = plt.axes()
ax.errorbar(T, ys, xerr=0., yerr=yserr, fmt='k.', label = 'Data')
xs = np.linspace(T.min(), T.max(), 5000)
ax.plot(xs, poly(pf, xs), 'r-', label = 'fit')

ax.set_title('Energy - Momentum Relation')
ax.set_xlabel('T')
ax.set_ylabel('$(pc)^2/2T$')
ax.legend(loc=(0.77,0.65))

textfit = '$f(T) = A T + B$ \n' \
          '$A = %.2f \pm %.2f$ \n' \
          '$B = %.1f \pm %.1f$ keV \n' \
          '$\chi^2= %.2f$ \n' \
          '$N = %i$ (dof) \n' \
          '$\chi^2/N = % .2f$' \
           % (pf[0], pferr[0], pf[1], pferr[1], chisq, dof,
              chisq/dof)
ax.text(0.1, .9, textfit, transform=ax.transAxes, fontsize=12,
         verticalalignment='top')
plt.savefig('plots/energy_momentum.png')
plt.show()



restmass = 2*peaks*(peaks - T)/T
rmerr = np.sqrt((4*peaks/T)*(4*peaks/T)*(peakerr)*(peakerr) + (2*peaks*peaks/(T*T))*(2*peaks*peaks/(T*T))*Terrs*Terrs)
rmerr += 2*peakerr
meanmass = np.mean(restmass)
print meanmass
massx = [meanmass] * len(restmass)
massx1 = [meanmass + 7] * 5000
massx2 = [meanmass - 7] * 5000
masserr = np.std(restmass)
print masserr

fig2 = plt.figure()
ax2 = plt.axes()
ax2.errorbar(T, restmass, xerr=0., yerr=rmerr, fmt='k.', label = 'Data')
ex = np.linspace(np.min(T), np.max(T), 5000)
ax2.plot(T, massx, 'r-', label = 'mean = 512.1keV')
ax2.plot(ex, massx1, 'b-.')
ax2.plot(ex, massx2, 'b-.')

ax2.set_title('Rest Mass Calculation')
ax2.set_xlabel('Compton Edge (keV)')
ax2.set_ylabel('Rest Mass')
ax2.legend(loc=(0.6,0.85))
plt.savefig('plots/restmass_T.png')
plt.show()

### BETA STUFF ###

beta = T * (2*peaks - T)/(T*T - 2*peaks*T + 2*peaks*peaks)
betaerr = (4*peaks*(peaks-T)*np.sqrt((peaks*Terrs)*(peaks*Terrs)+(T*peakerr)*(T*peakerr)))
betaerr /= (T*T - 2*T*peaks + 2*peaks*peaks)*(T*T - 2*T*peaks + 2*peaks*peaks)

fig3 = plt.figure()
ax3 = plt.axes()
ax3.errorbar(beta, pc, xerr=betaerr, yerr=pcerr, fmt='k.', label = 'Data')
ax3.set_title('Momentum vs Beta')
ax3.set_xlabel('Beta (v/c)')
ax3.set_ylabel('Momentum (pc)')
plt.savefig('plots/momentum_beta.png')
plt.show()

fig4 = plt.figure()
ax4 = plt.axes()
ax4.errorbar(beta, T, xerr=betaerr, yerr=Terrs, fmt='k.', label = 'Data')
ax4.set_title('Kinetic Energy vs Beta')
ax4.set_xlabel('Beta (v/c)')
ax4.set_ylabel('Kinetic Energy (T)')
plt.savefig('plots/T_beta.png')
plt.show()




