# constants.py
"""Declares physical constants for use in enterprise.
Depends on numpy for base mathematical constants, and
scipy.constants for physical constants.
"""

import numpy as np
import scipy.constants as sc
# Cosmological parameters used
hc = 0.67 #1801.04268
Om_Mat = 0.344 #2303.10095
Om_Rad = 2.47 * 10**(-5) / hc**2 #1801.04268
H_0 = 68.4 #In kms-1Mpc-1 2412.13045
T_0 = 2.35* 10**(-4) #In eV 
A_s = (np.e**(3.053))/(10**10) #Planck data PL21+BK18+LV21 2208
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
f_BBN = 1.5 * 10 **(-12) #1801.04268 in Hz
f_eq = 1.3*10**(-2)*c/Mpc #1801.04268 in Hz
z_eq = 3402 #Baumann
T_eq = T_0*(1 + z_eq)
f_LVK = 25 # In Hz
Om_LVK = 3.4 * 10**(-9)
del_N = 0.4 #2210.14159
eta_0 = 2 / (H_0*1000*(np.sqrt(Om_Rad) + np.sqrt(Om_Rad + Om_Mat))) * Mpc #In sec
