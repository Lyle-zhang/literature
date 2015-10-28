"""
Compare Vusse 1962 RTD model to Figure 6 in Smolders 2000 paper
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
    
   
# RTD model from Vusse 1962 
# -----------------------------------------------------------------------------

t = np.linspace(0, 25, 200)     # time range, s

r6a = rtd(4, 3.2, t)            # Fig. 6a in Smolders 2000
r6b = rtd(6, 2.8, t)            # Fig. 6b in Smolders 2000
r6c = rtd(5, 2.8, t)            # Fig. 6c in Smolders 2000

# Data from Smolders 2000 paper
# -----------------------------------------------------------------------------

x6a, y6a = np.loadtxt('fig6a.csv', delimiter=",", unpack=True) # Figure 6a
x6b, y6b = np.loadtxt('fig6b.csv', delimiter=",", unpack=True) # Figure 6b
x6c, y6c = np.loadtxt('fig6c.csv', delimiter=",", unpack=True) # Figure 6c

# Plot
# -----------------------------------------------------------------------------

py.close('all')

py.figure(1)
py.plot(t, abs(r6a.real), 'b-', lw=2, label='model')
py.plot(x6a, y6a, 'g--', lw=2, label='exp')
py.xlabel('Time (s)')
py.ylabel('Distribution function R(t) (1/s)')
py.title('Figure 6a')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, abs(r6b.real), 'b-', lw=2, label='model')
py.plot(x6b, y6b, 'g--', lw=2, label='exp')
py.xlabel('Time (s)')
py.ylabel('Distribution function R(t) (1/s)')
py.title('Figure 6b')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(3)
py.plot(t, abs(r6c.real), 'b-', lw=2, label='model')
py.plot(x6c, y6c, 'g--', lw=2, label='exp')
py.xlabel('Time (s)')
py.ylabel('Distribution function R(t) (1/s)')
py.title('Figure 6c')
py.legend(loc='best', numpoints=1)
py.grid()

py.show()