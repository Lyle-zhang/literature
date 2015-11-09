"""
Residence time distribution (RTD) of Figure 3 from pg. 339 in Kunii 1991 book
"""

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# Parameters
# -----------------------------------------------------------------------------

tau = 1.8                   # residence time of solids in bed, s
t = np.linspace(0, 10, 200) # time vector, s

# Residence Time Distribution, (RTD)
# -----------------------------------------------------------------------------

def etd(n, tau, t):
    """
    Exit age distribution (RTD) for solids in multistaged fluidized beds from 
    Eq 5 in Kunii 1991 book from pg 339.
    INPUTS
    n = number of stages, number of equal-sized beds in series
    tau = total solids residence time, s
    t = time vector, s
    OUTPUT
    et = exit age distribution or RTD of solids as a whole
    EXAMPLE
    
    """
    # solids residence time for each stage from Eq. 4
    ti = tau/n
    # RTD of the solids for the beds as a whole from Eq. 5
    et = 1/(np.math.factorial(n-1)*ti)*(t/ti)**(n-1)*np.exp(-t/ti)
    return et

# RTD for solids in a single bed, see Eq. 3
et = (1/tau)*np.exp(-t/tau)

# RTD for multistage bed operations, see Eq. 5
# n = number of stages, equal-sized beds in series
n = 1
et1 = etd(n, tau, t)    # should be equal to et above from Eq. 3

n = 2
et2 = etd(n, tau, t)

n = 3
et3 = etd(n, tau, t)

n = 5
et5 = etd(n, tau, t)

n = 10
et10 = etd(n, tau, t)

n = 20
et20 = etd(n, tau, t)

# check area under curves, should all be similar
print('et1 ', np.trapz(et1))
print('et2 ', np.trapz(et2))
print('et3 ', np.trapz(et3))
print('et5 ', np.trapz(et5))
print('et10 ', np.trapz(et10))
print('et20 ', np.trapz(et20))

# Plot
# -----------------------------------------------------------------------------

py.figure(1)
py.plot(t, et, 'b-', lw=2, label='Eq. 3')
py.plot(t, et1, 'r--', lw=2, label='Eq. 5')
py.ylabel('E(t)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, et1, lw=2, label='n=1')
py.plot(t, et2, lw=2, label='n=2')
py.plot(t, et3, lw=2, label='n=3')
py.plot(t, et5, lw=2, label='n=5')
py.plot(t, et10, lw=2, label='n=10')
py.plot(t, et20, lw=2, label='n=20')
py.ylabel('E(t)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.grid()
