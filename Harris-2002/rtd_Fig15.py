"""
Compare RTD data from Figure 15 in Harris 2002 paper to the RTD results from 
Vusse 1962.
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

n = 8                               # number of mixing stages, (-)
tau = 1.3                           # circulation time, s
t15a = np.linspace(0, 6, 100)       # time range, s
r15a = rtd(n, tau, t15a)            # RTD for Figure 15a

n = 4                               # number of mixing stages, (-)
tau = 3.6                           # circulation time, s
t15b = np.linspace(0, 40, 100)      # time range, s
r15b = rtd(n, tau, t15b)            # RTD for Figure 15b

n = 3                               # number of mixing stages, (-)
tau = 2.6                           # circulation time, s
t15c = np.linspace(0, 35, 100)      # time range, s
r15c = rtd(n, tau, t15c)            # RTD for Figure 15c

n = 8                               # number of mixing stages, (-)
tau = 1.8                           # circulation time, s
t15d = np.linspace(0, 14, 100)      # time range, s
r15d = rtd(n, tau, t15d)            # RTD for Figure 15d

# Data from paper to compare to rtd() results
#------------------------------------------------------------------------------

x15a, y15a = np.loadtxt('fig15a.csv', delimiter=",", unpack=True)
x15b, y15b = np.loadtxt('fig15b.csv', delimiter=",", unpack=True)
x15c, y15c = np.loadtxt('fig15c.csv', delimiter=",", unpack=True)
x15d, y15d = np.loadtxt('fig15d.csv', delimiter=",", unpack=True)

# Plot
#------------------------------------------------------------------------------

py.close('all')

py.figure(1)
py.plot(t15a, abs(r15a.real), 'b-', lw=2, label='vusse')
py.plot(x15a, y15a, 'o', mec='g', mew=2, mfc='none', label='exp')
py.xlabel('Time (s)')
py.ylabel('Distribution function (1/s)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t15b, abs(r15b.real), 'b-', lw=2, label='vusse')
py.plot(x15b, y15b, 'o', mec='g', mew=2, mfc='none', label='exp')
py.xlabel('Time (s)')
py.ylabel('Distribution function (1/s)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(3)
py.plot(t15c, abs(r15c.real), 'b-', lw=2, label='vusse')
py.plot(x15c, y15c, 'o', mec='g', mew=2, mfc='none', label='exp')
py.xlabel('Time (s)')
py.ylabel('Distribution function (1/s)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(4)
py.plot(t15d, abs(r15d.real), 'b-', lw=2, label='vusse')
py.plot(x15d, y15d, 'o', mec='g', mew=2, mfc='none', label='exp')
py.xlabel('Time (s)')
py.ylabel('Distribution function (1/s)')
py.legend(loc='best', numpoints=1)
py.grid()
