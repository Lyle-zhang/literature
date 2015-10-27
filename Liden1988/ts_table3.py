"""
Calculate values k, k2, ts for Table 3 in Liden 1988 paper
"""

import numpy as np

# Functions
# -----------------------------------------------------------------------------

def k(T):
    """
    Overall wood decomposition rate constant
    T = temperature, K
    """
    A = 10**13                  # pre-factor, 1/s
    E = 183.3e3                 # activation energy, J/mol
    R = 8.314                   # univeral gas constant, J/mol*K
    k = A*np.exp(-E/(R*T))      # overall wood decomposition rate constant, 1/s
    return k
    
    
def k2(T):
    """
    Tar decomposition rate constant
    T = temperature, K
    """
    A2 = 4.28e6                 # pre-factor, 1/s
    E2 = 107.5e3                # activation energy, J/mol
    R = 8.314                   # univeral gas constant, J/mol*K
    k2 = A2*np.exp(-E2/(R*T))   # tar decomposition rate constant, 1/s
    return k2


# Calculate k, k2 for ts-value in Table 3
# -----------------------------------------------------------------------------

# Run No. 30 Table 3

T = 698.15      # temperature, K
tg = 0.5        # gas residence time, s
phi = 0.559     # observed tar yield Table 2 Run no. 30
phis = 0.703    # ultimate tar yeild or max theoretical tar yield

k_30 = k(T)     # overall wood decomposition rate constant, 1/s
k2_30 = k2(T)   # tar decomposition rate constant, 1/s

tm1 = phi*k2_30*tg
tm2 = phis*(1-np.exp(-k2_30*tg))
ts = np.log(1-tm1/tm2)/-k_30   # solids residence time from Eq. 21, s

print('Run No. 30')
print('k is', k_30)     # should equal 0.194
print('k2 is', k2_30)   # should equal 0.0387
print('ts is', ts)      # should equal 8.4
print('---')

# Run No. 33 Table 3

T = 723.15     # temperature, K
tg = 0.69      # gas residence time Table 2 Run no. 33, s
phi = 0.558    # observed tar yield Table 2 Run no. 33
phis = 0.703   # ultimate tar yeild or max theoretical tar yield

k_33 = k(T)    # overall wood decomposition rate constant, 1/s
k2_33 = k2(T)  # tar decomposition rate constant, 1/s

tm1 = phi*k2_33*tg
tm2 = phis*(1-np.exp(-k2_33*tg))
ts = np.log(1-tm1/tm2)/-k_33   # solids residence time from Eq. 21, s

print('Run No. 33')
print('k is', k_33)      # should equal 0.578
print('k2 is', k2_33)    # should equal 0.0735
print('ts is', ts)       # should equal 2.9
print('---')

# Run No. 52 Table 3

T = 728.15      # temperature, K
tg = 0.51       # gas residence time Table 2 Run no. 52, s
phi = 0.6017    # observed tar yield Table 2 Run no. 52
phis = 0.703    # ultimate tar yeild or max theoretical tar yield

k_52 = k(T)     # overall wood decomposition rate constant, 1/s
k2_52 = k2(T)   # tar decomposition rate constant, 1/s

tm1 = phi*k2_52*tg
tm2 = phis*(1-np.exp(-k2_52*tg))
ts = np.log(1-tm1/tm2)/-k_52   # solids residence time from Eq. 21, s

print('Run No. 52')
print('k is', k_52)     # should equal 0.660
print('k2 is', k2_52)   # should equal 0.0831
print('ts is', ts)      # should equal 3.2
print('***')

# Use k, k2, phi directly for ts-value in Table 3
#------------------------------------------------------------------------------

phis = 0.703    # ultimate tar yeild or max theoretical tar yield
tg = 0.5        # gas residence time, s

# Run No. 30 Table 3

k = 0.194       # overall wood decomposition rate constant, 1/s
k2 = 0.0387     # tar decomposition rate constant, 1/s
phi = 0.559     # observed tar yield

tm1 = phi*k2*tg
tm2 = phis*(1-np.exp(-k2*tg))
ts = np.log(1-tm1/tm2)/-k   # solids residence time from Eq. 21, s

print('Run No. 30')
print('ts is', ts)
print('---')

# Run No. 52 Table 3

k = 0.660       # overall wood decomposition rate constant, 1/s
k2 = 0.0831     # tar decomposition rate constant, 1/s
phi = 0.6017    # observed tar yield

tm1 = phi*k2*tg
tm2 = phis*(1-np.exp(-k2*tg))
ts = np.log(1-tm1/tm2)/-k   # solids residence time from Eq. 21, s

print('Run No. 52')
print('ts is', ts)
print('---')

