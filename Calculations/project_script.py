# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 12:06:43 2018

@author: pizzaslayer
"""

# Solar Panel Project

import numpy as np
import matplotlib.pyplot as plt

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


#%% Convection analysis
hconv = 10.
hrad = 10.

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