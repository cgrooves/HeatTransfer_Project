# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 12:06:43 2018

@author: pizzaslayer
"""

# Solar Panel Project

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# Known variables
Tinf = 300.
Tsur = Tinf
L = 10e-2
k = 177.
t = 1e-3
T0 = 320.

V = 3.

Beta = 0.95
qdp_solar = 800.
emissivity = 0.2

boltzmann = 5.67e-8


#%% Convection analysis
Tavg = 329.4

Tref = (Tinf + Tavg)/2.

# Properties - calculate fluid properties at Tref and Patm using CoolProp
rho = CP.PropsSI('D','T',Tref,'P',101325,'Air')
mu = CP.PropsSI('V','T',Tref,'P',101325,'Air')
kair = CP.PropsSI('conductivity','T',Tref,'P',101325,'Air')
Pr = CP.PropsSI('PRANDTL','T',Tref,'P',101325,'Air')
nu = mu/rho


Lc = 10.*L # Assume Lc is entire length of plate
ReL = V*Lc/nu

print "Reynold's Number: %e" % ReL
if ReL < 5e5:
    print "Laminar flow"
else:
    print "Turbulent flow"
    
if Pr < 0.6:
    print "Prandtl number is < 0.6 - Nusselt correlation may not be reliable"
    
Nu = 0.664*ReL**.5*Pr**(1/3.)

hconv = Nu*kair/Lc

# Hrad - using Tavg as approximation
hrad = emissivity*boltzmann*(Tavg + Tsur)*(Tavg**2 + Tsur**2)

#%% Analysis - Model Building and Verification
h = hconv + hrad

m = np.sqrt(h/(k*t))
H = Beta*qdp_solar/h + Tinf

# Simplified temperature profile using hyperbolic functions
def Temp_Profile(x):
    return (T0-H)/np.cosh(m*L)*np.cosh(m*(x-L)) + H
    
# Plot from x = 0:2L
n_pts = 100
x_range = np.linspace(0,2*L,n_pts)
T = [0]*n_pts

for x in range(n_pts):
    T[x] = Temp_Profile(x_range[x])

plt.plot(x_range,T)

Tavg2 = (T0-H)/(m*L)*np.tanh(m*L) + H

print "Initial Tavg guessed: %f K" % Tavg
print "Analytical Tavg: %f K" % Tavg2

#%% Calculations using model
qdp_w = k*m*(H-T0)*np.tanh(m*L)

efficiency = qdp_w/qdp_solar

print "Estimated efficiency: %f " % efficiency