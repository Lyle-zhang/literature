# Sadhukhan2009 paper
# 1-D transient heat conduction

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries and packages
import numpy as np
import matplotlib.pyplot as py

#---- parameters

rhow = 682      # density of wood, kg/m^3
d = 0.02        # wood particle diameter, m
Ti = 285        # initial particle temp, K
Tinf = 683      # ambient temp, K
h = 48          # heat transfer coefficient, W/m^2*K
H = -220000     # heat of reaction, J/kg where (-) = exothermic, (+) = endothermic

#---- kinetic parameters

A1 = 168.4      # pre-factor, 1/s
E1 = 51.965     # activation energy, kJ/mol
A2 = 13.2
E2 = 45.960
A3 = 5.7e6
E3 = 92.4
R = 0.008314    # universal gas constant, kJ/mol*K

#---- initial calculations

b = 1           # run model as a cylinder (b = 1) or as a sphere (b = 2)

nt = 2000       # number of time steps
tmax = 800     # max time, s
dt = tmax/nt    # time step, s
t = np.arange(0,tmax+dt,dt)

nr = 19     # number or radius steps
r = d/2     # radius of particle, m
dr = r/nr   # radius step, delta r
m = nr+1    # nodes from center m=0 to surface m=steps+1

#---- initial temperatures, densities, mass fraction

T = np.zeros((m,1))     # column vector of temps
TT = np.zeros((1,m))    # row vector to store temps
pw = np.zeros((1,m))    # wood density at each node
pc = np.zeros((1,m))    # char density at each node
pg = np.zeros((1,m))    # gas density at each node

T[:,0] = Ti     # initial temp
TT[0,:] = Ti
pw[0,:] = rhow  # initial wood density

Ys = np.zeros(1)
Ys[0] = 1       # initial mass fraction, Ys=1 for all wood

#---- initial thermal properties 

cpw = 1112.0 + 4.85*(T.T-273.15)    # wood heat capacity, J/(kg*K) 
kw = 0.13 + (3e-4)*(T.T-273.15)     # wood thermal conductivity, W/(m*K)
cpc = 1003.2 + 2.09*(T.T-273.15)    # char heat capacity, J/(kg*K)
kc = 0.08 - (1e-4)*(T.T-273.15)     # char thermal conductivity, W/(m*K)

Yw = pw[0,:]/rhow               # wood fraction, Yw = 1 all wood, Yw = 0 all char
cpbar = Yw*cpw + (1-Yw)*cpc     # effective heat capacity
kbar = Yw*kw + (1-Yw)*kc        # effective thermal conductivity
pbar = pw[0,:] + pc[0,:]        # effective density

#---- initial kinetic reactions (used for initial heat generation term, g)
K1 = A1 * np.exp(-E1 / (R*T.T));

rp = -K1 * pw[0,:]
g = H * rp

#---- initial A and C arrays

A = np.zeros((m,m))
C = np.zeros((m,1))

#---- solve system of equations [A]{T}={C} where T = A\C for each time step

for i in range(1,nt+1):
    
    #-- heat transfer
    v = dt / (pbar[0]*cpbar[0,0])
    
    A[0,0] = 1 + (2*v*kbar[0,0]*(1+b)) / (dr**2)
    A[0,1] = -(2*v*kbar[0,0]*(1+b)) / (dr**2)
    C[0,0] = T[0] + v*g[0,0]
    
    for k in range(1,m-1):
        
        ri = (k*dr)**b
        rminus12 = ((k-0.5)*dr)**b
        rplus12 = ((k+0.5)*dr)**b
        kminus12 = (kbar[0,k]+kbar[0,k-1])/2
        kplus12 = (kbar[0,k]+kbar[0,k+1])/2
        
        w = dt/(pbar[k]*cpbar[0,k]*ri*(dr**2))
        z = dt/(pbar[k]*cpbar[0,k])
        
        A[k,k-1] = -w*rminus12*kminus12
        A[k,k] = 1 + w*rminus12*kminus12 + w*rplus12*kplus12
        A[k,k+1] = -w*rplus12*kplus12
        C[k,0] = T[k] + z*g[0,k]
    
    ww = dt/(pbar[m-1]*cpbar[0,m-1])
    krminus12 = (kbar[0,m-1]+kbar[0,m-2])/2
    
    A[m-1,m-2] = -(2*ww/(dr**2))*krminus12
    A[m-1,m-1] = 1 + (2*ww/(dr**2))*krminus12 + ww*((2/dr)+(b/r))*h
    C[m-1,0] = T[m-1] + ww*((2/dr)+(b/r))*h*Tinf + ww*g[0,m-1]
    
    T = np.linalg.solve(A,C)
    TT = np.vstack((TT, T.T))
    
    #-- kinetic reactions
    K1 = A1*np.exp(-E1/(R*T.T))
    K2 = A2*np.exp(-E2/(R*T.T))
    K3 = A3*np.exp(-E3/(R*T.T))
    
    rw = -(K1+K2)*pw[i-1,:]
    rg1 = K1*pw[i-1,:] - K3*pg[i-1,:]*pc[i-1,:]
    rc1 = K2*pw[i-1,:] - K3*pg[i-1,:]*pc[i-1,:]
    rg2 = K3*pg[i-1,:]*pc[i-1,:]
    rc2 = K3*pg[i-1,:]*pc[i-1,:]
    
    pw = np.vstack((pw,  pw[i-1,:]+rw*dt))
    pc = np.vstack((pc, pc[i-1,:] + (rc1+rc2)*dt))
    pg = np.vstack((pg, pg[i-1,:] + (rg1+rg2)*dt))
    
    #-- update thermal properties
    cpw = 1112.0 + 4.85*(T.T-273.15)
    kw = 0.13 + (3e-4)*(T.T-273.15)
    cpc = 1003.2 + 2.09*(T.T-273.15)
    kc = 0.08 - (1e-4)*(T.T-273.15)
    
    Yw = pw[i,:]/(pw[i,:]+pc[i,:])
    cpbar = Yw*cpw + (1-Yw)*cpc
    kbar = Yw*kw + (1-Yw)*kc
    pbar = pw[i,:] + 1.3*pc[i,:]
    Ys = np.vstack((Ys, np.mean(pbar)/rhow))
    
    rp = -K1*pw[i,:]
    g = H*rp

#---- plot results

# grab data from csv file, data taken from Fig 1 (cylinder) in Sadhukhan2009
Tcyl = np.loadtxt('Fig2_Tcylinder.csv', delimiter=',')
mcyl = np.loadtxt('Fig2_Mcylinder.csv', delimiter=',')

# time vs temperature
py.figure(1)
py.plot(t, TT[:, 0], '-b', label='center')
py.plot(t, TT[:, m-1], '-r', label='surface')
py.plot(Tcyl[:,0], Tcyl[:,1]+273.15, 'ob', label='expt')
py.axhline(Tinf, color='k', linestyle='-.')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title('Fig 2 - Cylinder')
py.legend(loc='best', numpoints=1)
py.ylim([Ti-20,Tinf+50])
py.grid()
py.show()

# time vs mass fraction 
py.figure(2)
py.plot(t, Ys, '-g', label='model')
py.plot(mcyl[:,0], mcyl[:,1], 'ob', label='expt')
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (-)')
py.title('Fig 2 - Cylinder')
py.legend(loc='best', numpoints=1)
py.ylim([0,1.05])
py.grid()
py.show()