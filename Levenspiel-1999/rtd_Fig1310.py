"""
Dispersion model curves as E(theta) and E(t) for open vessels having large 
deviation from plug flow where D/uL > 0.01. Refer to Ch. 13, pg 300-302 for 
Eq. 14 and Figure 13.10 in Levenspiel 1999 book.

Reference:
Levenspiel 1999. Chemical Reaction Engineering, 3rd Edition, Wiley & Sons, Inc.
"""

import numpy as np
import matplotlib.pyplot as py

# Functions for residence time distribution
# -----------------------------------------------------------------------------

def eth(DLu, theta):
    """
    E curve dispersion model for open vessels w/ large deviation from plug flow
    where D/uL > 0.01, see Eq 14 E(theta)00 on pg 301 in Levenspiel 1999 book.
    EXAMPLE
    theta = np.linspace(0, 2, 400)
    Eth = eth(0.1, theta)
    PARAMETERS
    DLu = D/Lu = dimensionless vessel dispersion number, (-)
    theta = dimensionless time, (-)
    RETURN
    Eth = E(theta)00 = dimensionless exit-age distribution, (-)
    """
    tm1 = 1/np.sqrt(4*np.pi*DLu*theta)
    tm2 = ((1-theta)**2)/(4*theta*DLu)
    Eth = tm1*np.exp(-tm2)
    return Eth
    
    
def eth2(D, L, u, theta):
    """
    E curve dispersion model for open vessels w/ large deviation from plug flow
    where D/uL > 0.01, see Eq 14 E(theta)00 on pg 301 in Levenspiel 1999 book.
    EXAMPLE
    theta = np.linspace(0, 2, 400)
    Eth = eth(0.1, theta)
    PARAMETERS
    D = dispersion coefficient, m^2/s
    L = length of vessel, m
    u = inlet velocity, m/s
    theta = dimensionless time range, (-)
    RETURN
    Eth = E(theta)00 = dimensionless exit-age distribution, (-)
    """
    tm1 = 1/(np.sqrt(4*np.pi*(D/(u*L))*theta))
    tm2 = ((1-theta)**2)/(4*theta*(D/(u*L)))
    Eth = tm1*np.exp(-tm2)
    return Eth


def et(D, L, u, t):
    """
    E curve dispersion model for open vessels w/ large deviation from plug flow
    where D/uL > 0.01, see Eq 14 E(t) on pg 301 in Levenspiel 1999 book.
    EXAMPLE
    t = np.linspace(0, 4, 100)
    Et = et(0.5, 1, 1, t)
    PARAMETERS
    D = dispersion coefficient, m^2/s
    L = length of vessel, m
    u = inlet velocity, m/s
    t = time range, s
    RETURN
    Et = E(t)00 = exit-age distribution, 1/s
    """
    tm1 = u/(np.sqrt(4*np.pi*D*t))
    tm2 = ((L-u*t)**2)/(4*D*t)
    Et = tm1*np.exp(-tm2)
    return Et
    
    
# Calculate E-curve for Dispersion Model
# -----------------------------------------------------------------------------

# Replicate E(theta),00 curves in Fig 13.10, pg 302 in Levenspiel 1999 book
theta = np.linspace(1e-6, 2, 400)  # dimensionless time range, (-)
Eth10 = eth(10, theta)          # E(theta),00 for D/Lu = 10
Eth1 = eth(1, theta)            # E(theta),00 for D/Lu = 1
Eth01 = eth(0.1, theta)         # E(theta),00 for D/Lu = 0.1
Eth001 = eth(0.01, theta)       # E(theta),00 for D/Lu = 0.01

# Data points from Fig 13.10 in Levenspiel 1999 book
x10, y10 = np.loadtxt('DuL10.csv', delimiter=",", unpack=True)
x1, y1 = np.loadtxt('DuL1.csv', delimiter=",", unpack=True)
x01, y01 = np.loadtxt('DuL01.csv', delimiter=",", unpack=True)
x001, y001 = np.loadtxt('DuL001.csv', delimiter=",", unpack=True)

# E(theta),00 curve as function of D, L, u, and theta
Eth = eth2(0.5, 1, 1, theta) 

# E(t) curve as function of D, L, u, and t
t = np.linspace(1e-6, 4, 400)  # time range, s
Et = et(0.5, 1, 1, t)

# Plot
# -----------------------------------------------------------------------------

py.close('all')

py.figure(1)
py.plot(theta, Eth10, lw=2, label='D/uL = 10')
py.plot(theta, Eth1, lw=2, label='D/uL = 1')
py.plot(theta, Eth01, lw=2, label='D/uL = 0.1')
py.plot(theta, Eth001, lw=2, label='D/uL = 0.01')
py.plot(x10, y10, 'o', mec='b', mew=2, mfc='none')
py.plot(x1, y1, 'o', mec='g', mew=2, mfc='none')
py.plot(x01, y01, 'o', mec='r', mew=2, mfc='none')
py.plot(x001, y001, 'o', mec='c', mew=2, mfc='none')
py.title('E-curves for Fig 13.10 in Levenspiel 1999')
py.xlabel('$\Theta$ (-)')
py.ylabel('E($\Theta$)')
py.ylim([0, 2])
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(theta, Eth, lw=2)
py.xlabel('$\Theta$ (-)')
py.ylabel('E($\Theta$)')
py.grid()

py.figure(3)
py.plot(t, Et, lw=2)
py.xlabel('Time (s)')
py.ylabel('E(t)')
py.grid()