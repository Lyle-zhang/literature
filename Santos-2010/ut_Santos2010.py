# Calculate terminal velocity, Ut, for a near spherical particle
# from Santos2010 book, pg 86, Eqs 4.12 - 4.15

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

#---- Parameters

# air properties at T = 300K, P = 1 atm
rhog = 1.17     # density (kg/m^3)
ug = 1.85e-5    # dynamic viscosity (kg/ms)
g = 9.81        # gravity (m/s^2)

# particle properties (near spherical particle)
dp = 0.000207   # diameter of particle (m)
rhop = 2500      # density of particle (kg/m^3)

#---- Calculations

# particle terminal velocity, Ut (m/s)
# particle Reynolds number, Re (-)

# for Re <= 2
Ut1 = (g*(dp**2)*(rhop - rhog)) / (18*ug)
Re1 = (dp*rhog*Ut1) / ug

print('Re <= 2 ....... Ut =', Ut1, 'Re =', Re1)

# for 2 < Re <= 500
Ut2 = ( (g*(dp**1.6)*(rhop-rhog)) / (13.9*(rhog**0.4)*(ug**0.6)) )**0.71
Re2 = (dp*rhog*Ut2) / ug

print('2 < Re <= 500 . Ut =', Ut2, 'Re =', Re2)

# for Re > 500
Ut3 = ( (3.03*g*dp*(rhop - rhog)) / rhog )**0.5
Re3 = (dp*rhog*Ut3) / ug

print('Re > 500 ...... Ut =', Ut3, 'Re =', Re3)

# if statement that checks Re for the Ut
Ut = (g*(dp**2)*(rhop - rhog)) / (18*ug)
Re = (dp*rhog*Ut) / ug

if Re <= 2:
    Utt = Ut
    Ree = Re
elif Re > 2 and Re <= 500:
    Utt = ( (g*(dp**1.6)*(rhop-rhog)) / (13.9*(rhog**0.4)*(ug**0.6)) )**0.71
    Ree = (dp*rhog*Utt) / ug
elif Re > 500:
    Utt = ( (3.03*g*dp*(rhop - rhog)) / rhog )**0.5
    Ree = (dp*rhog*Utt) / ug

print('--- Santos2010 ---')
print('Re =', Ree)
print('Ut =', Utt)