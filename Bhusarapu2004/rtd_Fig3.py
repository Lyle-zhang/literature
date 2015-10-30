"""
Compare RTD data from Figure 3a in Bhusarapu 2004 paper to the RTD results from 
Vusse 1962.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
#------------------------------------------------------------------------------

#tau = 0.5                   # circulation time, s
#t = np.linspace(0, 2, 100)  # time range, s

# Function for residence time distribution
#------------------------------------------------------------------------------

def rtd(n, tau, t, q=1, r=1):
    """
    Residence time distribution function for stirred tank reactor.
    Assumes circulating flow, pumping action of stirrer, isotropic homogeneous 
    turbulence in the circulating fluid.
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
    
    
# Calculate residence time distribution for n and tau
#------------------------------------------------------------------------------

n3a = 3
tau3a = 42.8
theta3a = np.linspace(0, 80, 100)/42.8  # dimensionless time range
rt3a = rtd(n3a, tau3a, theta3a)

n3b = 3
tau3b = 15.5
theta3b = np.linspace(0, 80, 100)/15.5  # dimensionless time range
rt3b = rtd(n3b, tau3b, theta3b)

# Data from paper to compare to rtd() results
#------------------------------------------------------------------------------

x3a, y3a = np.loadtxt('fig3a.csv', delimiter=",", unpack=True)
x3b, y3b = np.loadtxt('fig3b.csv', delimiter=",", unpack=True)

# Plot
#------------------------------------------------------------------------------

py.close('all')

py.figure(1)
py.plot(theta3a, abs(rt3a.real)/2, 'b-', lw=2, label='vusse')
py.plot(x3a, y3a, 'o', mec='g', mew=2, mfc='none', label='exp')
py.xlabel('$\Theta$')
py.ylabel('E($\Theta$)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(theta3b, abs(rt3b.real), 'b-', lw=2, label='vusse')
py.plot(x3b, y3b, 'o', mec='g', mew=2, mfc='none', label='exp')
py.xlabel('$\Theta$')
py.ylabel('E($\Theta$)')
py.legend(loc='best', numpoints=1)
py.grid()
