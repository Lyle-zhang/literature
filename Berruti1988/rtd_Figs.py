"""
Compare Vusse 1962 RTD model to Figures 5, 6, and 7 in Berruti 1988 paper.
Bubbling bed experiment with bed material as silica sand at dp = 710 um and 
density of rho = 2470 kg/m^3.
"""

import numpy as np
import matplotlib.pyplot as py

# Function for residence time distribution
#------------------------------------------------------------------------------

def rtd(n, tau, t, q=1, r=1):
    """
    Residence time distribution function for stirred tank reactor from paper by 
    Vusse 1962. Assumes circulating flow, pumping action of stirrer, isotropic 
    homogeneous turbulence in the circulating fluid.
    INPUTS: 
    n = number of mixing stages, (-)
    tau = circulation time, s
    t = time vector, s
    q = feed rate, m^3/s
    r = circulation rate, m^3/s
    OUTPUT:
    rt = residence time distribution, 1/s
    EXAMPLE:
    r = rtd(n, tau, t) or rtd(n, tau, t, q=2, r=3)
    """
    a = (n/tau)*(r/(r+q))**(1/n)    # term from Eq. 30
    s = 0                           # initialize summation for g(at) Eq. 31
    
    for k in range(0, n):
        # summation loop for number of mixing stages, see g(at) from Eq. 31
        tm = np.exp(a*t*(np.cos(2*np.pi*k/n) + 1j*np.sin(2*np.pi*k/n)) + 1j*(2*np.pi*k)/n)
        s += tm
    
    # g as g(at) from Eq. 31 and rt as R(t) from Eq. 30
    g = (1/n)*s
    # vector of complex values for residence time distribution
    rt = (q/(r+q))*((r+q)/r)**((n-1)/n)*(n/tau)*np.exp(-(n*t)/tau)*g
    return rt
    
   
def weibull(x, lam, k):
    """
    Weibull distribution function.
    x = time parameter
    k = shape parameter
    lam = lambda as scale parameter
    """
    w = (k/lam)*((x/lam)**(k-1))*np.exp(-(x/lam)**k)
    return w
    
    
# Vusse 1962 RTD model
# -----------------------------------------------------------------------------

t = np.linspace(0, 20, 200)     # time range, s

r5 = rtd(10, 1.8, t)            # Fig. 5, Exp. 4 in Berruti 1988
r6 = rtd(13, 2.8, t)            # Fig. 6, Exp. 6 in Berruti 1988
r7 = rtd(10, 3.2, t)            # Fig. 7, Exp. 9 in Berruti 1988

# Compare Vusse 1962, CSTR series, and Weibull distribution
# -----------------------------------------------------------------------------

r5v = rtd(9, 2.2, t)            # Vusse model for Figure 5
r5c = rtd(3, 2.9, t, r=1e-4)    # CSTR series for Figure 5
r5w = weibull(t, 3, 2)          # Weibull distribution for Figure 5

# Data from Berruti 1988 paper
# -----------------------------------------------------------------------------

x5, y5 = np.loadtxt('fig5.csv', delimiter=",", unpack=True) # Figure 5
x6, y6 = np.loadtxt('fig6.csv', delimiter=",", unpack=True) # Figure 6
x7, y7 = np.loadtxt('fig7.csv', delimiter=",", unpack=True) # Figure 7

# Plot
# -----------------------------------------------------------------------------

py.close('all')

py.figure(1)
py.plot(t, abs(r5.real), 'b-', lw=2, label='Vusse model')
py.plot(x5, y5, 'g--', lw=2, label='Berruti exp')
py.xlabel('Time (s)')
py.ylabel('RTD function (1/s)')
py.title('Berruti 1988, Figure 5')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, abs(r6.real), 'b-', lw=2, label='Vusse model')
py.plot(x6, y6, 'g--', lw=2, label='Berruti exp')
py.xlabel('Time (s)')
py.ylabel('RTD function (1/s)')
py.title('Berruti 1988, Figure 6')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(3)
py.plot(t, abs(r7.real), 'b-', lw=2, label='Vusse model')
py.plot(x7, y7, 'g--', lw=2, label='Berruti exp')
py.xlabel('Time (s)')
py.ylabel('RTD function (1/s)')
py.title('Berruti 1988, Figure 7')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(4)
py.plot(t, abs(r5v.real), 'b-', lw=2, label='Vusse model')
py.plot(t, abs(r5c.real), 'k-', lw=2, label='CSTR series')
py.plot(t, abs(r5w.real), 'm-', lw=2, label='Weibull distribution')
py.plot(x5, y5, 'o', mew=2, mec='g', mfc='none', label='Berruti experiment')
py.xlabel('Time (s)')
py.ylabel('RTD (1/s)')
py.title('Berruti 1988, Figure 5')
py.legend(loc='best', numpoints=1)
py.grid()

py.show()