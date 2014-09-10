# Correlations from Churchill1974 paper
# see section "Functions Which Cross One Limiting Solution", pg. 41
# plot should be same as plot in Fig. 10, pg. 41

# use Python 3 print function
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# parameters
# -----------------------------------------------------------------------------

# range of x = 0-6
x = np.linspace(np.spacing(1), 6)

# vector to store zinf values
zinf = np.zeros(len(x))

# vector to store z at n=4 values
z4 = np.zeros(len(x))

# vector to store z at n=2 values
z2 = np.zeros(len(x))

# vector to store z at n=1 values
z1 = np.zeros(len(x))

# vector to store z at n=1/2 values
z12 = np.zeros(len(x))


# functions
# -----------------------------------------------------------------------------

# where zi{x} = x,  z{inf} = 1,  xA = 5,  alpha = 2

# zinf{x} function as seen in Fig. 9 and Eq. 8
# zinf{x} = 1 + (5/x)^2
def z_inf(x):
    return 1 + (5.0 / x)**2

    
# z{x} function as seen in Fig. 9 and Eq. 9
# z{x} = x / ( [1 + t^n]^(1/n) ) where t = x / (1 + (5/x)^a)
def z(x, n, a):
    t = x / (1 + (5.0/x)**a)
    return x / (( 1 + t**n )**(1/n))    
    
    
# calculations
# -----------------------------------------------------------------------------

# evaluate zi{x} = x
zi = x

# evaluate zinf{x}
k = 0

for i in x:
    zinf[k] = z_inf(i)
    k+=1

# evaluate z{x} for n = 2
# and alpha = 4, 2, 1, 1/2
k = 0

for i in x:
    z4[k] = z(i, 2, 4)
    z2[k] = z(i, 2, 2)
    z1[k] = z(i, 2, 1)
    z12[k] = z(i,2, 0.5)
    k+=1


# plot results
# -----------------------------------------------------------------------------

py.figure(1)
py.plot(x, zinf, 'g--', label=r'z$_\infty${x}')
py.plot(x, zi, 'b--', label=r'z$_0${x}')
py.plot(x, z4, 'y', label='z{x} a=4')
py.plot(x, z2, 'r', label='z{x} a=2')
py.plot(x, z1, 'm', label='z{x} a=1')
py.plot(x, z12, 'c', label='z{x} a=1/2')
py.axhline(y=1, color='k', label=r'z{$\infty$}=1')
py.ylim([0, 4])
py.legend(loc='best', numpoints=1)
py.xlabel('X')
py.ylabel('Z')
py.title('Churchill1974 - Figure 10\n'+r'z$_0${x}=x, z{$\infty$}=1, x$_A$=5, n=2')
py.grid()
py.show()
