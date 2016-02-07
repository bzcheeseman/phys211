__author__ = "Aman LaChapelle"

import numpy as np

def lorentzian(x, A, Gamma, x0, C):
    return A/np.pi*((.5*Gamma)/((x-x0)**2+(.5*Gamma)**2)) + C

def multi_lorentzian(x, *p):
    lor = P[-1]
    for i in range(0, len(p)-1, 3):
        p[i] = A
        p[i+1] = Gamma
        p[i+2] = x0
        lor +=  A/np.pi*((.5*Gamma)/((x-x0)**2+(.5*Gamma)**2))
    return lor

def gaussian(x, A, sigma, center, C):
    return A*np.exp(-(x-center)**2/(2*sigma**2)) + C

def linear(x, A, B):
    return A*x + B

def quadratic(x, A, B, C):
    return A*x**2 + B*x + C

def exponential(x, A, B, C):
    return A*np.exp(B*x) + C

def cosine(x, A, omega, delta, C):
    return A*np.cos(omega*x + delta) + C

def multi_cosine(x, A, B, omega1, omega2, delta1, delta2, C):
    return A*np.cos(omega1*x + delta1) + B*np.cos(omega2*x + delta2) + C
