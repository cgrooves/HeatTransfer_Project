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
g = 9.81


#%% Convection analysis
def convection_analysis(Tavg):
        
    Tref = (Tinf + Tavg)/2.
    
    # Properties - calculate fluid properties at Tref and Patm using CoolProp
    rho = CP.PropsSI('D','T',Tref,'P',101325,'Air')
    mu = CP.PropsSI('V','T',Tref,'P',101325,'Air')
    cp_air = CP.PropsSI('Cp0mass','T',Tref,'P',101325,'Air')
    kair = CP.PropsSI('conductivity','T',Tref,'P',101325,'Air')
    Pr = CP.PropsSI('PRANDTL','T',Tref,'P',101325,'Air')
    alpha_air = kair/rho/cp_air
    nu = mu/rho
    
    Ltot = 10.*L # Assume Lc is entire length of plate
    Lc = Ltot/2.
    Ra = g*(Tavg-Tinf)*Lc**3/(Tref*nu*alpha_air)
    
    print "Rayleigh Number: %e" % Ra
    if Ra < 9e7 and Ra > 1e4 and Pr > 0.7:
        print "Using correlation (9.30)"
        Nu = 0.54*Ra**.25
    elif Ra < 9e11 and Ra > 1e7:
        print "Using correlation (9.31)"
        Nu = 0.15*Ra**(1/3.)
    else:
        print "No correlation selected for Rayleigh and Pr #"
        return
        
    hconv = Nu*kair/Lc
    
    # Hrad - using Tavg as approximation
    hrad = emissivity*boltzmann*(Tavg + Tsur)*(Tavg**2 + Tsur**2)
    
    return (hconv, hrad)

#%% Analysis - Model Building and Verification
Tavg1 = 329.7
(hconv, hrad) = convection_analysis(Tavg1)

h = hconv + hrad

m = np.sqrt(h/(k*t))
H = Beta*qdp_solar/h + Tinf

# Simplified temperature profile using hyperbolic functions
def Temp_Profile(x):
    return (T0-H)/np.cosh(m*L)*np.cosh(m*(x-L)) + H
    
# Plot from x = 0:2L
n_pts = 100
x_range = np.linspace(0,2*L,n_pts)
T1 = [0]*n_pts

for x in range(n_pts):
    T1[x] = Temp_Profile(x_range[x])

Tmax1 = np.max(T1)

Tavg2 = (T0-H)/(m*L)*np.tanh(m*L) + H

print "Initial Tavg guessed: %f K" % Tavg1
print "Analytical Tavg: %f K" % Tavg2

# Calculations using model
qdp_w = k*m*(H-T0)*np.tanh(m*L)

efficiency = qdp_w/qdp_solar

print "Estimated efficiency without convection shield: %f " % efficiency

#%%
print "\n\nAnalysis with convection shield:"
Tavg1 = 333.3
(hconv, hrad) = convection_analysis(Tavg1)
# Add cover-plate to absorber plate design
h = hrad

# Re-do calculations
m = np.sqrt(h/(k*t))
H = Beta*qdp_solar/h + Tinf

# Simplified temperature profile using hyperbolic functions
def Temp_Profile(x):
    return (T0-H)/np.cosh(m*L)*np.cosh(m*(x-L)) + H
    
# Plot from x = 0:2L
n_pts = 100
x_range = np.linspace(0,2*L,n_pts)
T2 = [0]*n_pts

for x in range(n_pts):
    T2[x] = Temp_Profile(x_range[x])

Tmax2 = np.max(T2)

Tavg2 = (T0-H)/(m*L)*np.tanh(m*L) + H

print "Initial Tavg guessed: %f K" % Tavg1
print "Analytical Tavg: %f K" % Tavg2

qdp_w = k*m*(H-T0)*np.tanh(m*L)

efficiency2 = qdp_w/qdp_solar

print "Estimated efficiency: %f " % efficiency2

#%% SUMMARY AND ANALYSIS
eta_gain = (efficiency2/efficiency - 1)*100
print "\n\nSUMMARY:\n"
print "Increase in efficiency due to convective shield is %f %%" % eta_gain
print "Max temp. without shield: %f K\nMax temp. with shield: %f K" % (Tmax1, Tmax2)
print "Increase in max temperature is %f %%" % (100*(Tmax2/Tmax1 - 1))

# Plots of temp. profile for shielded, unshielded (separate and together)
x1 = np.concatenate((x_range, x_range+2*L, x_range+4*L, x_range+6*L, x_range+8*L))
T1 = np.concatenate((T1,T1,T1,T1,T1))

plt.plot(x1,T1)