# file pars.py
# constants in cgs units
import math

gN = 6.67408e-08 # Newton's gravitational constant
mSun = 1.9884754153381438e+33 # mass of the sun in grams
mEarth = 5.97219e+27 # mass of the earth
au = 1.495978707e13 # astronomical unit in cm
yr = 2*math.pi / (gN*mSun/au**3)**0.5 # 1 year in seconds
# ...

# Problem Parameters
Np = 2 # no. of particles
