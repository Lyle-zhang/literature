"""
Function for 1D transient heat conduction within a solid sphere, cylinder, or 
slab shape. Convection at the surface and symmetry at center.

References:
1) Ozisik, M. Necati, 1994. Finite Difference Methods in Heat Transfer.
2) Bergman, Lavine, Incropera, Dewitt, 2011. Fundamentals of Heat and Mass 
   Transfer, 7th Edition.
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np
import scipy.linalg as sp

# Function
# -----------------------------------------------------------------------------

def hc(m, dr, b, dt, h, Tinf, g, T, i, r, pbar, cpbar, kbar):
    """
    1D transient heat conduction within a solid sphere, cylinder, or slab shape 
    with convection at the surface and symmetry at center. Returns an array of 
    temperatures [T] at each intraparticle node point.
    
    Solves system of equations [A]*[x] = [b] where:
    A = known coefficent matrix, tridiagonal for 1-D problem
    x = unknown vector, next temperature at each node
    b = known vector, current temperature at each node
    
    Example:
    T = hc(m, dr, b, dt, h, Tinf, g, T, i, r, pbar, cpbar, kbar)
    
    where:
    m = number of nodes from center (m=0) to surface (m)
    dr = radius step, m
    b = shape factor where 2 is sphere, 1 is cylinder, 0 is slab
    dt = time step, s
    h = heat transfer coefficient, W/m^2*K
    Tinf = ambient temperature, K
    g = heat generation
    T = temperature at node
    i = time index
    r = radius of particle, m
    pbar = effective density or concentration
    cpbar = effective heat capacity, J/kg*K
    kbar = effective thermal conductivity, W/m*K
    """
    
    ab = np.zeros((3, m))   # banded array from the tridiagonal matrix
    bb = np.zeros(m)        # column vector
    
    k = np.arange(1, m-1)
    ri = (k * dr)**b
    rminus12 = ((k-0.5) * dr)**b
    rplus12 = ((k+0.5) * dr)**b
    
    v = dt / (pbar[0] * cpbar[0])
    
    # create internal terms
    kminus12 = (kbar[k] + kbar[k-1])/2
    kplus12 = (kbar[k] + kbar[k+1])/2
    w = dt / (pbar[k] * cpbar[k] * ri * (dr**2))
    z = dt / (pbar[k] * cpbar[k])
    
    # create surface terms
    ww = dt / (pbar[m-1] * cpbar[m-1])
    krminus12 = (kbar[m-1] + kbar[m-2])/2
    
    # upper diagonal
    ab[0, 1] = -(2 * v * kbar[0] * (1+b)) / (dr**2)     # center node T1
    ab[0, 2:] = -w * rplus12 * kplus12                  # internal nodes Tm+1
    
    # center diagonal
    ab[1, 0] = 1 + (2* v * kbar[0] * (1+b)) / (dr**2)                         # center node T0
    ab[1, 1:m-1] = 1 + w * rminus12 * kminus12 + w * rplus12 * kplus12        # internal nodes Tm
    ab[1, m-1] = 1 + (2*ww/(dr**2)) * krminus12 + ww * ((2/dr) + (b/r))*h     # surface node Tr
    
    # lower diagonal
    ab[2, 0:m-2] = -w * rminus12 * kminus12     # internal nodes Tm-1
    ab[2, m-2] = -(2*ww/(dr**2)) * krminus12    # surface node Tr-1
    
    # column vector
    bb[0] = T[i-1, 0] + v*g[0]                                       # center node T0
    bb[1:m-1] = T[i-1, k] + z*g[k]                                   # internal nodes Tm
    bb[m-1] = T[i-1, m-1] + ww*((2/dr)+(b/r))*h*Tinf + ww*g[m-1]     # surface node Tr    
    
    # temperatures
    T = sp.solve_banded((1, 1), ab, bb)
    
    return T