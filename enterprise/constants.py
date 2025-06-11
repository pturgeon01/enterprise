# constants.py
"""Declares physical constants for use in enterprise.
Depends on numpy for base mathematical constants, and
scipy.constants for physical constants.
"""

import numpy as np
import scipy.constants as sc
# Cosmological parameters used
hc = 0.6736 #1807.06209 TT,TE,EE+lowE+lensing
Om_Mat = 0.3153 #1807.06209
Om_Rad = 2.47 * 10**(-5) / hc**2 #1801.04268 Look further 
H_0 = 67.36 #In kms-1Mpc-1 1807.06209
T_0 = 2.35* 10**(-4) #In eV 
A_s = 2.1/(10**9) #1801.04268
M_PL = 1.2*10**(28) #in eV
h_bar = 6.582*10**(-16) #in eVs
# mathematical constants from numpy
# the log constants are useful for converting between log base 10 and e
pi = np.pi
e = np.e
log10e = np.log10(np.e)
ln10 = np.log(10.0)

# physical constancts in mks
c = sc.speed_of_light
G = sc.gravitational_constant
h = sc.Planck

# astronomical times in sec (and frequencies in Hz)
yr = sc.Julian_year
day = sc.day
fyr = 1.0 / yr

# astronomical distances in meters
AU = sc.astronomical_unit
ly = sc.light_year
pc = sc.parsec
kpc = pc * 1.0e3
Mpc = pc * 1.0e6
Gpc = pc * 1.0e9

# solar mass in kg and m,s natural units
GMsun = 1.327124400e20  # measured more precisely than Msun alone!
Msun = GMsun / G
Rsun = GMsun / (c**2)
Tsun = GMsun / (c**3)

# other useful conversions in mks
erg = sc.erg

# other things
DM_K = 2.41e-16  # for DM variation design matrix

# relative angle between the Earth's ecliptic and the galactic equator
e_ecl = 23.43704 * np.pi / 180.0

# unit vector pointing in direction of angle between Earth's ecliptic and the galactic equator
M_ecl = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(e_ecl), -np.sin(e_ecl)], [0.0, np.sin(e_ecl), np.cos(e_ecl)]])

# Other used parameters
f_pl = 1.854 * 10**(43) #In Hz
T_BBN = 10**(5) #In eV 1801.04268
f_eq = np.sqrt(2)*H_0*10**3/Mpc * Om_Mat/np.sqrt(Om_Rad) / (2*np.pi)  #Derived 
f_ref = c / (2*np.pi) * 0.05/Mpc #1511.05994
z_eq = 3402 #1807.06209
T_eq = T_0*(1 + z_eq)
f_LVK = 25 # In Hz
Om_LVK = 1.7 * 10**(-8) #2101.12130 flat GWB (prior)
del_N = 0.2 #1801.04268 Compare to #1807.06209
f_0 = H_0*1000/(2*np.pi*Mpc) #In Hz derived
eta_0 = 2*(H_0*1000/Mpc)**(-1)*1/(np.sqrt(Om_Rad + Om_Mat) + np.sqrt(Om_Rad)) #In sec derived
