"""
Residence time distribution function for Figure 6 in Vusse 1962 paper for a 
stirred tank reactor. Assumes circulating flow by pumping action of stirrer and 
an isotropic homogeneous turbulence in the circulating fluid.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
#------------------------------------------------------------------------------

tau = 0.5                   # circulation time, s
t = np.linspace(0, 2, 100)  # time range, s

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

r1 = rtd(1, tau, t)
r2 = rtd(2, tau, t)
r4 = rtd(4, tau, t)
r10 = rtd(10, tau, t)

# Data from paper to compare to rtd() results
#------------------------------------------------------------------------------

file_n1 = 'n1.csv'
x1,y1 = np.loadtxt(file_n1, delimiter=",", unpack=True)

file_n2 = 'n2.csv'
x2, y2 = np.loadtxt(file_n2, delimiter=",", unpack=True)

file_n4 = 'n4.csv'
x4, y4 = np.loadtxt(file_n4, delimiter=",", unpack=True)

file_n10 = 'n10.csv'
x10, y10 = np.loadtxt(file_n10, delimiter=",", unpack=True)

# Plot
#------------------------------------------------------------------------------

py.close('all')

py.figure(1)
py.plot(t, r1.real, 'r-', lw=2, label='n=1')
py.plot(x1, y1, 'ro', mec='r')
py.plot(t, r2.real, 'g-', lw=2, label='n=2')
py.plot(x2, y2, 'go', mec='g')
py.plot(t, abs(r4.real), 'b-', lw=2, label='n=4')
py.plot(x4, y4, 'bo', mec='b')
py.plot(t, abs(r10.real), 'm-', lw=2, label='n=10')
py.plot(x10, y10, 'mo', mec='m')
py.xlabel('Time (s)')
py.ylabel('Distribution function R(t) (1/s)')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()
