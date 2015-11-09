"""
Residence time distribution function for Figure 6 in Vusse 1962 paper for a 
stirred tank reactor. Assumes circulating flow by pumping action of stirrer and 
an isotropic homogeneous turbulence in the circulating fluid.
"""

import numpy as np
import matplotlib.pyplot as py

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

tau = 0.5                   # circulation time, s
t = np.linspace(0, 2, 100)  # time range, s

r1 = rtd(1, tau, t)         # one mixing stage, (-)
r2 = rtd(2, tau, t)         # two mixing stages, (-)
r4 = rtd(4, tau, t)         # four mixing stages, (-)
r10 = rtd(10, tau, t)       # ten mixing stages, (-)

# Data from paper to compare to rtd() results
#------------------------------------------------------------------------------

x1, y1 = np.loadtxt('n1.csv', delimiter=",", unpack=True)
x2, y2 = np.loadtxt('n2.csv', delimiter=",", unpack=True)
x4, y4 = np.loadtxt('n4.csv', delimiter=",", unpack=True)
x10, y10 = np.loadtxt('n10.csv', delimiter=",", unpack=True)

# Plot
#------------------------------------------------------------------------------

py.close('all')

py.figure(1)
py.plot(t, r1.real, 'r-', lw=2, label='n=1')
py.plot(x1, y1, 'o', mec='r', mew=2, mfc='none', label='paper')
py.plot(t, r2.real, 'g-', lw=2, label='n=2')
py.plot(x2, y2, 'o', mec='g', mew=2, mfc='none', label='paper')
py.plot(t, abs(r4.real), 'b-', lw=2, label='n=4')
py.plot(x4, y4, 'o', mec='b', mew=2, mfc='none', label='paper')
py.plot(t, abs(r10.real), 'm-', lw=2, label='n=10')
py.plot(x10, y10, 'o', mec='m', mew=2, mfc='none', label='paper')
py.xlabel('Time (s)')
py.ylabel('RTD function (1/s)')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()
