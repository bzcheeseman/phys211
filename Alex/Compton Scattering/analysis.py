from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt


#constants
mc2 = 511. #keV
E = 662. #keV, energy of initial gamma photons


angles1 = np.array([0.17453293,0.26179939,0.38397244,0.52359878,0.61086524,1.04719755,1.22173048,1.30899694])
angles2 = np.array([1.30899694,1.3962634,1.57079633,1.74532925,2.0943951,2.35619449,0.78539816])


en1 = np.array([645.94535489,638.70777021,609.75743148,572.12199113,558.37058023,404.93378497,361.50827687,342.6905567])
en2 = np.array([339.09451519,322.89846155,291.97872279,274.31030064,227.9306925,206.58134907,482.65044515])
  
en1err = np.array([6.3184713,6.26564529,6.05434124,5.77964599,5.67927657,4.55936515,4.24240908,4.10506146])
en2err = np.array([4.10100687,3.98279443,3.75711613,3.6281571, 3.28963966,3.13381417,5.14879897])

angles = np.append(angles1, angles2)
en = np.append(en1, en2)
enerr = np.append(en1err, en2err)

def fitfunc(p, x):
    y=(p[0]/p[1])*(1-np.cos(x))
    return p[0]/(1+y)
def residual(p, x, y, dy):
    return (fitfunc(p, x)-y)/dy
    
p = [E, mc2] #expected

pf, cov, info, mesg, success = optimize.leastsq(residual, p, args=(angles, en, enerr), full_output=1)
chisq = sum(info["fvec"]*info["fvec"])
dof = len(en)-len(pf)
pferr = [np.sqrt(cov[i,i]) for i in range(len(pf))]
print pf
print cov

fig1 = plt.figure()
ax1 = plt.axes()
ax1.errorbar(angles, en, xerr=0., yerr=enerr, fmt='k.', label = 'Data')
T = np.linspace(angles.min(), angles.max(), 5000)
ax1.plot(T, fitfunc(pf, T), 'r-', label = 'Fit')
ax1.plot(T, fitfunc(p, T), 'b-', label = 'Theory')

ax1.set_title('Compton Scattering Energy Fit')
ax1.set_xlabel('Angle (radians), x')
ax1.set_ylabel('Energy (keV), E\'')
ax1.legend()
    
textfit = '$E\' = E / ((1 + (E/m_e c^2)(1 - cos(x)))$ \n' \
          '$E = %.2f \pm %.2f$ keV \n' \
          '$m_e c^2 = %.2f \pm %.2f$ keV \n' \
          '$\chi^2= %.2f$ \n' \
          '$N = %i$ (dof) \n' \
          '$\chi^2/N = % .2f$' \
           % (pf[0], pferr[0], pf[1], pferr[1], chisq, dof,
              chisq/dof)
ax1.text(0.05, .33, textfit, transform=ax1.transAxes, fontsize=12,
         verticalalignment='top')
             

plt.show()