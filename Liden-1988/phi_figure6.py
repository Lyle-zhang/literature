"""
Calculate tar yield as phi for Table 3 and Figure 6 in Liden 1988 paper
"""

import numpy as np
import matplotlib.pyplot as py

# Functions
# -----------------------------------------------------------------------------

def phi(ts, tg, T):
    """
    Observed tar yield
    ts = solids residence time, s
    tg = gas residence time, s
    T = temperature, K
    """
    R = 8.314       # univeral gas constant, J/mol*K
    phis = 0.703    # ultimate tar yeild or max theoretical tar yield Eq. 18, (-)
    
    A = 1e13                    # pre-factor, 1/s
    E = 183.3e3                 # activation energy, J/mol
    k = A*np.exp(-E/(R*T))      # overall wood decomposition rate constant, 1/s
    
    A2 = 4.28e6                 # pre-factor, 1/s
    E2 = 107.5e3                # activation energy, J/mol
    k2 = A2*np.exp(-E2/(R*T))   # tar decomposition rate constant, 1/s
    
    # observed tar yield Eq. 18, (-)
    tar = phis*(1-np.exp(-k*ts))*((1-np.exp(-k2*tg))/(k2*tg))
    return tar
    
    
def phik(ts, tg, k, k2, T):
    phis = 0.703    # ultimate tar yeild or max theoretical tar yield
    tar = phis*(1-np.exp(-k*ts))*((1-np.exp(-k2*tg))/(k2*tg))
    
    return tar
    
    
# Run No. 30 Table 3
# -----------------------------------------------------------------------------

T = 698.15      # temperature, K
tg = 0.5        # gas residence time, s
ts = 8.4        # solids residence time, s

tar = phi(ts, tg, T)    # observed tar yield with calculated k, k2 from Eq. 18
tark = phik(ts, tg, 0.194, 0.0387, T)   # observed tar yield with Table 3 k, k2

print('Run No. 30')
print('phi is', tar)    # shoud equal 0.559
print('phi is', tark)
print('---')


# Run No. 33 Table 3
# -----------------------------------------------------------------------------

T = 723.15      # temperature, K
tg = 0.69       # gas residence time, s
ts = 2.9        # solids residence time, s

tar = phi(ts, tg, T)    # observed tar yield with calculated k, k2 from Eq. 18
tark = phik(ts, tg, 0.578, 0.0735, T)   # observed tar yield with Table 3 k, k2

print('Run No. 33')
print('phi is', tar)    # shoud equal 0.558
print('phi is', tark)
print('---')

# Run No. 52 Table 3
# -----------------------------------------------------------------------------

T = 728.15      # temperature, K
tg = 0.51       # gas residence time, s
ts = 3.2        # solids residence time, s

tar = phi(ts, tg, T)    # observed tar yield with calculated k, k2 from Eq. 18
tark = phik(ts, tg, 0.660, 0.0831, T)   # observed tar yield with Table 3 k, k2

print('Run No. 52')
print('phi is', tar)    # shoud equal 0.6017
print('phi is', tark)
print('---')

# Plot Tar Yield for Gas Residence Times in Figure 6
#------------------------------------------------------------------------------
py.close('all')

tg = np.linspace(0.2, 0.8, 100)
tar773 = phi(3, tg, 773)
tar873 = phi(3, tg, 873)

py.figure(1)
py.plot(tg, tar773, lw=2, label='773 K')
py.plot(tg, tar873, lw=2, label='873 K')
py.xlabel('Gas Residence Time (s)')
py.ylabel('Tar Yield (-)')
py.title('Figure 6')
py.yticks(np.arange(0, 1.1, 0.1))
py.legend(loc='best', numpoints=1)
py.grid()
py.show()

# Plot Tar Yield for Solids Residence Times
#------------------------------------------------------------------------------

ts = np.linspace(1, 4, 100)
tar = phi(ts, 0.51, 728)
tark = phik(ts, 0.51, 0.660, 0.0831, 728)

py.figure(2)
py.plot(ts, tar, lw=2, label='k-calc')
py.plot(ts, tark, lw=2, label='k-given')
py.xlabel('Solids Residence Time (s)')
py.ylabel('Tar Yield (-)')
py.title('Calculated vs Table 3 k-values, Run No. 52')
py.legend(loc='lower right', numpoints=1)
py.grid()
py.show()

# Plot Tar Yield at Various Temperatures
# -----------------------------------------------------------------------------

tar698 = phi(ts, 0.5, 698)
tar728 = phi(ts, 0.5, 728)
tar773 = phi(ts, 0.5, 773)
tar800 = phi(ts, 0.5, 800)

py.figure(3)
py.plot(ts, tar698, lw=2, label='698 K')
py.plot(ts, tar728, lw=2, label='728 K')
py.plot(ts, tar773, lw=2, label='773 K')
py.plot(ts, tar800, lw=2, label='800 K')
py.xlabel('Solids Residence Time (s)')
py.ylabel('Tar Yield (-)')
py.title('Effect of Temp. on Tar Yield, Gas Res. Time = 0.5 s')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()
