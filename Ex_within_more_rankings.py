# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 14:21:12 2020
@author: Jordan
Using the Ex_within to try out new rankings and incorperate a grb prob
"""


import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d; from astropy_healpix import HEALPix; 
from astropy.coordinates import ICRS, SkyCoord; from astropy import units as u;
from timeit import default_timer as timer

plt.close('all')



N = 1000  ##Change to alter the number of loops the code runs for

'''
I would like to know the relative ratio of of the number of GRBs which occur within
250 MPc to those outside, up to z = 1 (we assume that everything past this point
is just considered to be at z = 1). To first order this value is just given by 
the ratio of volume.
'''

R_min = 250
R_z1 = 3550.7

Ratio = (R_min/R_z1)**3

"""
Factor to increase the volume ratio is (z+1) as the arrivial rates of GRBs should 
be taken into account when considering the volume out to Z = 1
"""
z_max = 1
factor = (1 + z_max)
Ratio = Ratio * factor

placement = np.zeros(N)
placement2 = np.zeros(N)
placement3 = np.zeros(N)
placement4 = np.zeros(N)
placement5 = np.zeros(N)
placement6 = np.zeros(N)
placement7 = np.zeros(N)
placement8 = np.zeros(N)
placement9 = np.zeros(N)
placement10 = np.zeros(N)
placement11 = np.zeros(N)
placement12 = np.zeros(N)
placement13 = np.zeros(N)
placement14 = np.zeros(N)
placement15 = np.zeros(N)
placement16 = np.zeros(N)
placement17 = np.zeros(N)
placement18 = np.zeros(N)
placement19 = np.zeros(N)
placement20 = np.zeros(N)
placement21 = np.zeros(N)
placement22 = np.zeros(N)
placement23 = np.zeros(N)
placement24 = np.zeros(N)
placement25 = np.zeros(N)
placement26 = np.zeros(N)
placement27 = np.zeros(N)
placement28 = np.zeros(N)
placement29 = np.zeros(N)

percentages = np.zeros(N)
percentages2 = np.zeros(N)
percentages3 = np.zeros(N)
percentages4 = np.zeros(N)
percentages5 = np.zeros(N)
percentages6 = np.zeros(N)
percentages7 = np.zeros(N)
percentages8 = np.zeros(N)
percentages9 = np.zeros(N)
percentages10 = np.zeros(N)
percentages11 = np.zeros(N)
percentages12 = np.zeros(N)
percentages13 = np.zeros(N)
percentages14 = np.zeros(N)
percentages15 = np.zeros(N)
percentages16 = np.zeros(N)
percentages17 = np.zeros(N)
percentages18 = np.zeros(N)
percentages19 = np.zeros(N)
percentages20 = np.zeros(N)
percentages21 = np.zeros(N)
percentages22 = np.zeros(N)
percentages23 = np.zeros(N)
percentages24 = np.zeros(N)
percentages25 = np.zeros(N)
percentages26 = np.zeros(N)
percentages27 = np.zeros(N)
percentages28 = np.zeros(N)
percentages29 = np.zeros(N)

no_se_func = []
ras_dex = np.zeros(shape = (N, 2))
test_case = np.zeros(shape = (N, 2))


def Ang_Dist(ra1, ra2, dec1, dec2):## Calculates the angular distance between apparent position and galaxy
    
    ra1 *= (np.pi/180); ra2 *= (np.pi/180)
    dec1 *= (np.pi/180); dec2 *= (np.pi/180)
    
    return (180/np.pi) * np.arccos(np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
################################################################# David's ranks


def rank(theta, sigma, d_lum, luminosity, luminosity_probability): ## Normal
    ## Implements a ranking statistic defined in report
    return np.exp(-(theta**2/(2 * (sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability #* Colour_factor

def rank2(theta, sigma, d_lum, luminosity, luminosity_probability): ## Luminosity
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity**2)[:, 0] * luminosity_probability**2 

def rank3(theta, sigma, d_lum, luminosity, luminosity_probability): ## Luminosity Distance
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity**2)[:, 0] * luminosity_probability**2 

def rank4(theta, sigma, d_lum, luminosity, luminosity_probability): ## Lum_prob
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity**3)[:, 0] * luminosity_probability**2

def rank5(theta, sigma, d_lum, luminosity, luminosity_probability): ## Lum_prob, Lum
    return np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2 

def rank6(theta, sigma, d_lum, luminosity, luminosity_probability): ## D_Lum, Lum_prob
    return np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2 

def rank7(theta, sigma, d_lum, luminosity, luminosity_probability): ## D_lum, Lum
    return np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability 

def rank8(theta, sigma, d_lum, luminosity, luminosity_probability): ## All
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2

def rank9(theta, sigma, d_lum, luminosity, luminosity_probability): ## Angular Distance
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability

def rank10(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, D_Lum
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity)[:, 0] * luminosity_probability  

def rank11(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, Lum
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity**2)[:, 0] * luminosity_probability 

def rank12(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, Lum_Prob
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability**2

def rank13(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Ang_Dist
    return np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2 

def rank14(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Lum
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2

def rank15(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except d_lum
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2

def rank16(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Lum_prob
    return np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability

def rank17(theta, sigma, d_lum, luminosity, luminosity_probability): ## No angular Distance
    return np.exp(0 * -(theta**2/(2 *(sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability 

def rank18(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity Distance
    return np.exp(-(theta**2/(2 * (sigma)**2))) *(1/d_lum**0 * luminosity)[:, 0] * luminosity_probability

def rank19(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity
    return np.exp(-(theta**2/(2 * (sigma)**2))) *(1/d_lum * luminosity**0)[:, 0] * luminosity_probability**2

def rank20(theta, sigma, d_lum, luminosity, luminosity_probability): ## 23 of daves old functions
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2

def rank21(theta, sigma, d_lum, luminosity, luminosity_probability): ## Optimise 1
    return np.exp(-(2 * theta**2*((sigma)**2))/10) *(1/d_lum**0 * luminosity**5)[:, 0] * luminosity_probability

def rank22(theta, sigma, d_lum, luminosity, luminosity_probability): ## angular /4, dlum**2, lum**2 and lum prob **0.5
    return np.exp(-((theta**2) * (sigma**2))/(4)) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**(0.5)

def rank23(theta, sigma, d_lum, luminosity, luminosity_probability): ## dL**1, theta sigma/100
    return np.exp(-((theta**2)*(2 * (sigma)**2))/100) *(1/d_lum**1 * luminosity)[:, 0] * luminosity_probability
'''

###########################################################################################

def rank(theta, sigma, d_lum, luminosity, luminosity_probability): ## Normal
    ## Implements a ranking statistic defined in report
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability #* Colour_factor

def rank2(theta, sigma, d_lum, luminosity, luminosity_probability): ## Luminosity
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability 

def rank3(theta, sigma, d_lum, luminosity, luminosity_probability): ## Luminosity Distance
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability 

def rank4(theta, sigma, d_lum, luminosity, luminosity_probability): ## Lum_prob
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability**2

def rank5(theta, sigma, d_lum, luminosity, luminosity_probability): ## Lum_prob, Lum
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2 

def rank6(theta, sigma, d_lum, luminosity, luminosity_probability): ## D_Lum, Lum_prob
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2 

def rank7(theta, sigma, d_lum, luminosity, luminosity_probability): ## D_lum, Lum
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability 

def rank8(theta, sigma, d_lum, luminosity, luminosity_probability): ## All
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2

def rank9(theta, sigma, d_lum, luminosity, luminosity_probability): ## Angular Distance
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability

def rank10(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, D_Lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability  

def rank11(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, Lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability 

def rank12(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, Lum_Prob
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability**2

def rank13(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Ang_Dist
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2 

def rank14(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2

def rank15(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except d_lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2

def rank16(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Lum_prob
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability

def rank17(theta, sigma, d_lum, luminosity, luminosity_probability): ## No angular Distance
    return np.exp(0 * -(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability 

def rank18(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity Distance
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability

def rank19(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum * luminosity**0)[:, 0] * luminosity_probability**2

def rank20(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity Probability
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability**0

def rank21(theta, sigma, d_lum, luminosity, luminosity_probability): ## Optimise 1
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(4) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2

def rank22(theta, sigma, d_lum, luminosity, luminosity_probability): ## Optimise 2
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2

def rank23(theta, sigma, d_lum, luminosity, luminosity_probability): ## Optimise 2
    return np.exp(-((theta**2)**1000/(2 * (sigma)**2))) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2
'''
################################################################################################
# My Ranks 

def rank24(theta, sigma, d_lum, luminosity, luminosity_probability):## square on angle term, no luminosity and 4 on prob
    return np.exp(-(4*(theta**2) * (sigma**2))/(2)) *(1/d_lum**8 * luminosity**2)[:, 0] * luminosity_probability

def rank25(theta, sigma, d_lum, luminosity, luminosity_probability):##no luminosity or distance dependence 
    return np.exp(-(theta**2)/(2*(sigma**2))) * ((1/d_lum**8 * luminosity**(1/2))[:, 0])**0 * luminosity_probability 

def rank26(theta, sigma, d_lum, luminosity, luminosity_probability):
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity)[:, 0] * luminosity_probability

def rank27(theta, sigma, d_lum, luminosity, luminosity_probability):
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity)[:, 0] * luminosity_probability**2

def rank28(theta, sigma, d_lum, luminosity, luminosity_probability):
    return (np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity**3)[:, 0] * luminosity_probability**2)**2 + np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity)[:, 0] * luminosity_probability**2

def rank29(theta, sigma, d_lum, luminosity, luminosity_probability):
    return np.exp(-(100*theta**2)/(2*sigma**2)) * (1/d_lum * luminosity)[:, 0] * luminosity_probability * (abs(theta - sigma)) + np.exp(-(theta**2/(2*sigma**2)))*(1/d_lum**3 * luminosity**2)[:, 0] * luminosity_probability
"""
#################################################################
#Daves old functions before I fixed them
def rank(theta, sigma, d_lum, luminosity, luminosity_probability): ## Normal
    ## Implements a ranking statistic defined in report
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability #* Colour_factor

def rank2(theta, sigma, d_lum, luminosity, luminosity_probability): ## Luminosity
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability 

def rank3(theta, sigma, d_lum, luminosity, luminosity_probability): ## Luminosity Distance
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability 

def rank4(theta, sigma, d_lum, luminosity, luminosity_probability): ## Lum_prob
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability**2

def rank5(theta, sigma, d_lum, luminosity, luminosity_probability): ## Lum_prob, Lum
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2 

def rank6(theta, sigma, d_lum, luminosity, luminosity_probability): ## D_Lum, Lum_prob
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2 

def rank7(theta, sigma, d_lum, luminosity, luminosity_probability): ## D_lum, Lum
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability 

def rank8(theta, sigma, d_lum, luminosity, luminosity_probability): ## All
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2

def rank9(theta, sigma, d_lum, luminosity, luminosity_probability): ## Angular Distance
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability

def rank10(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, D_Lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability  

def rank11(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, Lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability 

def rank12(theta, sigma, d_lum, luminosity, luminosity_probability): ## Ang_Dist, Lum_Prob
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability**2

def rank13(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Ang_Dist
    return np.exp(-(theta**2/(2 *(sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2 

def rank14(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2

def rank15(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except d_lum
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2

def rank16(theta, sigma, d_lum, luminosity, luminosity_probability): ## All except Lum_prob
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability

def rank17(theta, sigma, d_lum, luminosity, luminosity_probability): ## No angular Distance
    return np.exp(0 * -(theta**2/(2 *(sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability 

def rank18(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity Distance
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability

def rank19(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum * luminosity**0)[:, 0] * luminosity_probability**2

def rank20(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity Probability
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability**0

def rank21(theta, sigma, d_lum, luminosity, luminosity_probability): ## Optimise 1
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(4) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2

def rank22(theta, sigma, d_lum, luminosity, luminosity_probability): ## Optimise 2
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2

def rank23(theta, sigma, d_lum, luminosity, luminosity_probability): ## Optimise 2
    return np.exp(-((theta**2)**100/(2 * (sigma)**2))) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2
"""


#making a blank txt file to keep track of what the statistics were that produced each ranking
file = open("Statistics.txt", "w")
file.write("1 np.exp(-(theta**2/(2 * (sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability\n")
file.write("\n2 np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity**2)[:, 0] * luminosity_probability**2\n")
file.write("\n3 np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity**2)[:, 0] * luminosity_probability**2\n")
file.write("\n4 np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity**3)[:, 0] * luminosity_probability**2\n")
file.write("\n5 np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2\n")
file.write("\n6 np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2\n")
file.write("\n7 np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability\n")
file.write("\n8 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2\n")
file.write("\n9 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability\n")
file.write("\n10 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity)[:, 0] * luminosity_probability\n")
file.write("\n11 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity**2)[:, 0] * luminosity_probability\n")
file.write("\n12 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability**2\n")
file.write("\n13 np.exp(-(theta**2/(2 *(sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**2\n")
file.write("\n14 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity)[:, 0] * luminosity_probability**2\n")
file.write("\n15 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum * luminosity**2)[:, 0] * luminosity_probability**2\n")
file.write("\n16 np.exp(-(theta**2/((sigma)**2))) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability\n")
file.write("\n17 np.exp(0 * -(theta**2/(2 *(sigma)**2))) *(1/d_lum * luminosity)[:, 0] * luminosity_probability\n")
file.write("\n18 np.exp(-(theta**2/(2 * (sigma)**2))) *(1/d_lum**0 * luminosity)[:, 0] * luminosity_probability\n")
file.write("\n19 np.exp(-(theta**2/(2 * (sigma)**2))) *(1/d_lum * luminosity**0)[:, 0] * luminosity_probability**2\n")
file.write("\n20 np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability**2\n")
file.write("\n21 np.exp(-(2 * theta**2*((sigma)**2))/10) *(1/d_lum**0 * luminosity**5)[:, 0] * luminosity_probability\n")
file.write("\n22 np.exp(-((theta**2) * (sigma**2))/(4)) *(1/d_lum**2 * luminosity**2)[:, 0] * luminosity_probability**(0.5)\n")
file.write("\n23 np.exp(-((theta**2)*(2 * (sigma)**2))/100) *(1/d_lum**1 * luminosity)[:, 0] * luminosity_probability\n")
file.write("\n24 np.exp(-(4*(theta**2) * (sigma**2))/(2)) *(1/d_lum**8 * luminosity**2)[:, 0] * luminosity_probability\n")
file.write("\n25 np.exp(-(theta**2)/(2*(sigma**2))) * ((1/d_lum**8 * luminosity**(1/2))[:, 0])**0 * luminosity_probability\n")
file.write("\n26 np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity)[:, 0] * luminosity_probability\n")
file.write("\n27 np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity)[:, 0] * luminosity_probability**2\n")
file.write("\n28 (np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**0 * luminosity**3)[:, 0] * luminosity_probability**2)**2 + np.exp(-(theta**2/(2 * (sigma)**2)))**(sigma**4) * (1/d_lum**1 * luminosity)[:, 0] * luminosity_probability**2\n")
file.write("\n29 np.exp(-(100*theta**2)/(2*sigma**2)) * (1/d_lum * luminosity)[:, 0] * luminosity_probability * (abs(theta - sigma)) + np.exp(-(theta**2/(2*sigma**2)))*(1/d_lum**3 * luminosity**2)[:, 0] * luminosity_probability\n")
file.close()

#################################################################
def convert(h, m, s): #Hours minutes seconds to degrees (More for applied code than here)
    return h + (m/60) + (s/3600)
#################################################################
def Luminosity_Handling(magnitude): ##Converts Absolute B Magnitude to Luminosity
    
    solar_b = 4.74
    solar_l = 1 #3.846e26 W
    
    return solar_l * 10**(0.4 * (solar_b - magnitude)) ## Gives an array in terms of solar luminosity
###########################################################
def spherical_convert(ra, dec): ##Test  ##Converts ra and dec to an xyz array
    r = 1
    #ra = phi
    #dec = theta
    ##Convert to radians

    ra = ra * np.pi/180
    dec = dec * np.pi/180
    
    x = np.cos(ra) * np.cos(dec)
    y = np.sin(ra) * np.cos(dec)
    z = np.sin(dec)
    
    return np.array([x, y, z])
############################################################
def rotation(x, angle):##Test  #Rotation about the z axis
    #need angle in radians
    
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                          [np.sin(angle), np.cos(angle), 0],
                          [0, 0, 1]])
    
    return x * rotation
############################################################
def back_convert(axyz): ##Test  ## Converts xyz coordinates to ra and dec
    
    x = axyz[0]
    y = axyz[1]
    z = axyz[2]
    r = modulus(axyz)
       
    arg1 = float(y/x)
    arg2 = float(z/r)
    
    phi = np.arctan(arg1)
    theta = np.arccos(arg2)
        
    return (180/np.pi) * phi, (90 - theta * (180/np.pi))## Returns ra, dec in that order in degrees
#################################################################
def modulus(array):  ##Test  ##Finds the modulus of a matrix/array
    
    return np.sqrt(array[0]**2 + array[1]**2 + array[2]**2)
#################################################################
def find_nearest(array, value): #Kind of a hash and not exactly interpolation, but for this point, should be okay
    
    array = np.asarray(array) - value
    truey = [i for i, val in enumerate(array) if val >= 0]
    idx = truey[0]#(np.abs(array - value)).argmin()

    return idx    
#################################################################
def reduction(RA_dec, Dec_dec, df_master):  ##Reduces the df_master by considering angular distance
    
    #host = df_master.iloc[current_i]
    #RA_dec = ra_prime[0]#host[["RA"]].values.tolist()[0]
    #Dec_dec = dec_prime[0]#host[["dec"]].values.tolist()[0]    
    ## Testing purposes only (hashed out lines)
    
    RA = df_master[["RA"]].values.tolist()    
    ra_arry = np.isclose(RA, RA_dec, atol = error_radius)
    res_ra = [i for i, val in enumerate(ra_arry) if val == False] ##Something up here - removing too many items
    
    DEC = df_master[["dec"]].values.tolist()
    dec_arry = np.isclose(DEC, Dec_dec, atol = error_radius)
    res_dec = [i for i, val in enumerate(dec_arry) if val == False]
    
    indices_to_keep = set(range(df_master.shape[0])) - set(res_ra) - set(res_dec)
    df_sliced = pd.DataFrame.take(df_master, list(indices_to_keep), axis = 0)
        
    ra = df_sliced[["RA"]].values
    dec = df_sliced[["dec"]].values

    return np.array(ra[:, 0]), np.array(dec[:, 0]), df_sliced
#################################################################
def Luminosity_back_convert(L_given, d_L): #  ##Converts luminosity to luminosity at source
    #L = L0/4 *np.pi * d_l**2
    return (L_given) * (4 * np.pi * (3.086e22 * d_L)**2)

def Luminosity_for_convert(L_given, d_L): #  ##Converts luminosity at source to apparent luminosity
    return(L_given)/(4 * np.pi * (3.086e22 * d_L)**2)
#################################################################
def L_func(L_test, c, d_L): ##   ##Takes an input and returns a probability based on the broken power law
    
    L_star = np.log10(4.61e51 * 1e7)        ##All from Guetta/Piran 2005
    del_1 = 30
    del_2 = 10
    alpha = 0.5
    beta = 1.5
    
    L = np.zeros(len(d_L))
    SGR_test = np.zeros(len(d_L))
    
    for j in range(len(d_L)):  ## Slightly inefficient, but on the scales of reduced catalog, not too drastic
        L[j]  = np.log10(Luminosity_back_convert(L_test, d_L[j]))
        
    L_prob = np.zeros(len(L))
    
    for i in range(len(L)):
        if L[i] < L_star and (L_star/del_1) < L[i]:
            L_prob[i] = c * (L[i]/L_star)**-alpha
    
        elif L[i] > L_star and L[i] < del_2 * L_star:
             L_prob[i] = c * (L[i]/L_star)**-beta
    
        elif L[i] < (L_star/del_1):   
            L_prob[i] = 0 ## What to do when the values fall outside the range that defines the power law?
            SGR_test[i] = 1  ##Creates a flag for if the luminosity at source would be low enough to be considered an SGR
        
        else:
            L_prob[i] = 0
            
    return L_prob, SGR_test
#################################################################
def L_func1(L): ##  ##Builds the broken power law based on a log scale from 52 to 59
    
    L_star = np.log10(4.61e51 * 1e7)
    del_1 = 30
    del_2 = 10
    alpha = 0.5
    beta = 1.5
    
    N = len(L)
    L2 = np.zeros(N)
    summ = 0
    sum1 = np.zeros(N)
    
    for i in range(N):
        
        if L[i] < L_star and (L_star/del_1) < L[i]:
            L2[i] = (L[i]/L_star)**-alpha
        
        elif L[i] > L_star and L[i] < del_2 * L_star:
                L2[i] = (L[i]/L_star)**-beta
        
        else:
            L2[i] = L_star
        summ += L2[i]
    c = 1/(summ)
    sum1[i] = summ
    
    L2 *= c
    
    return L2, c
#################################################################
def cumulative(array): ###   #Builds cumulative distributions
    
    N = array.shape[0]
    summing = np.zeros(N + 1)
    #array = L2
    
    for i in range(1, N + 1):
        df = pd.DataFrame(array[:i])
        summing[i] = df.sum().values.tolist()[0]
        
    return summing# /= summing[-1]

##If you have N galaxies
##########################################################################################
def axis_rotation(axis, point, angle):  ## Rotation about an axis function
    
    init_matrix = np.array([[0, -1 *  axis[2], axis[1]],
                           [axis[2], 0, -1 * axis[0]],
                           [-1 * axis[1], axis[0], 0]])
    
    matrix_2 = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]])
    
    term_2 = np.sin(angle) * init_matrix
    
    rot_matrix = (1 - np.cos(angle)) * np.dot(init_matrix, init_matrix) + term_2 + matrix_2
    
    rotated_point = np.dot(rot_matrix, point)
    
    return rotated_point

def Sector_find_reduced(RA_grb, Dec_grb, err_radius):
    '''
    Give coordinates of the grb location and an error in the position, this function
    will use cone_search to find all sky sectors that the cone intersects and 
    will read the corresponding csv files and compile them into one dataframe
    '''
    
    #corrects for if the rotations of the galaxy coords puts the GRB in an invalid position
    if abs(Dec_grb) > 90:
        x = RA_grb
        parity = Dec_grb/abs(Dec_grb)
        
        Dec_grb = (180 - abs(Dec_grb))*parity
        
        RA_grb = RA_grb + 180
        
        if RA_grb > 360:
            RA_grb = x - 180
    
    elif RA_grb < 0:
        RA_grb = 360 + RA_grb
        
    #making the sky coordinates
    coords = SkyCoord(RA_grb, Dec_grb, unit = "deg")
    
    #finding intersecting sectors 
    sectors = hp.cone_search_skycoord(coords, radius = err_radius*u.degree)
    
    #making the empty dataframe
    df_container = pd.DataFrame()
    
    for i in sectors:
        '''
        loop over the intersecting sectors to read the files and append to 
        the df_container
        '''
        
        name = name = str("Sector_{}".format(i))
        holder = pd.read_csv("Data Files/GLADE_Sectioned_Reduced/{}.csv".format(name),\
                             delimiter = ",", index_col = 0)
        
        df_container = df_container.append(holder)
    
    return df_container

def Sector_find(RA_grb, Dec_grb, err_radius):
    '''
    Give coordinates of the grb location and an error in the position, this function
    will use cone_search to find all sky sectors that the cone intersects and 
    will read the corresponding csv files and compile them into one dataframe
    '''
    
    #corrects for if the rotations of the galaxy coords puts the GRB in an invalid position
    if abs(Dec_grb) > 90:
        x = RA_grb
        parity = Dec_grb/abs(Dec_grb)
        
        Dec_grb = (180 - abs(Dec_grb))*parity
        
        RA_grb = RA_grb + 180
        
        if RA_grb > 360:
            RA_grb = x - 180
    
    elif RA_grb < 0:
        RA_grb = 360 + RA_grb
        
    #making the sky coordinates
    coords = SkyCoord(RA_grb, Dec_grb, unit = "deg")
    
    #finding intersecting sectors 
    sectors = hp.cone_search_skycoord(coords, radius = err_radius*u.degree)
    
    #making the empty dataframe
    df_container = pd.DataFrame()
    
    for i in sectors:
        '''
        loop over the intersecting sectors to read the files and append to 
        the df_container
        '''
        
        name = name = str("Sector_{}".format(i))
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(name),\
                             delimiter = ",", index_col = 0)
        
        df_container = df_container.append(holder)
    
    return df_container
#########################################################################################
#########################################################################################
#%%
start = timer()

df_master = pd.read_csv("Data Files/GLADE_Master_comma_100Mpc.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

L1 = np.linspace(56, 59, 101) #In J now
L2, c = L_func1(L1) #  ##Builds broken power law

cumuL = cumulative(L2) ##Luminosity Distribution
df_cumLum =   pd.read_csv("Data Files/Cumulative_Luminosity_100Mpc.csv")
df_cumLum.columns = ["NaN", "Cumulative Luminosity"]
normal_c = df_cumLum[["Cumulative Luminosity"]].values[-1][0]
L_rank = df_cumLum[["Cumulative Luminosity"]].values * 1/normal_c
df_cumLum = df_cumLum[["Cumulative Luminosity"]].values#                ## This is all to do with building a usable and callable power law
                                                                         

lum_N = np.linspace(0, df_cumLum.shape[0], df_cumLum.shape[0])
df_dL = df_master[["Luminosity Distance"]]

#using HEALPix to split the sky into equal area sectors
hp = HEALPix(nside=16, order='ring', frame=ICRS())

tests = randint(0, 2, size = N) ## If tests[i] = 0, use test galaxy, or if = 1, choose random point beyond the catalog
dummies = random(N)
RandL = random(N)

gals = np.zeros(N)          ## Picks out a luminosity
gal_index = np.zeros(N)
"""
aa = np.zeros(shape = (N, 5))  # Storing Angular distance
ab = np.zeros(shape = (N, 5))  # Storing Luminosity Distance
ac = np.zeros(shape = (N, 5))  # Storing B Luminosity
ad = np.zeros(shape = (N, 5))  # Storing Luminosity Probability
"""
lum_atsource = np.zeros(N)
received_luminosity = np.zeros(N)
cumul_N = np.zeros(N)

lum_list = list(L_rank)

df_dL = df_dL.values.tolist()   ## Luminosity distance values for use
a = np.zeros(N)                 ## For storing temporary and unimportant values
b = np.zeros(N)                 ## For storing temporary and unimportant values

test_ra = df_master[["RA"]]
test_dec = df_master[["dec"]]

indices = list(np.arange(df_master.shape[0]))

error_radius = 2 * (2.62)  ## Change as necessary - this is an example value from HEARSCH 

percentages = np.zeros(N)
distances = np.zeros(N)
luminosity_i = np.zeros(N)    
rank_host = np.zeros(N)
faulty = np.zeros(shape = (N, 5))       ## All of this used to store values

phi = 2 * np.pi * random(N) * (180/np.pi)               ## Random positions for rotations
theta = np.arccos(2 * random(N) - 1) * (180/np.pi)      

thph = spherical_convert(theta, phi)
mod = np.zeros(N)

for i in range(N):
    mod[i] = modulus(thph[:, i])
    thph[:, i] /= mod[i]

xyz = np.zeros(shape = (N, 3))
m = np.zeros(shape = (N, 3))
ra_prime = np.zeros(N); dec_prime = np.zeros(N)
rotation_angle = error_radius * normal(size = N) * (np.pi/180) 


#I want to try and keep the summation of the ranks and the largest rank so that I can work with them later

names = np.array(["Max Rank", "Sum Ranks", "Top 5 Avg"], dtype = str)
number_ranks = np.arange(2, 30, 1)

for i in number_ranks:
    name_hold = np.array([names[0] + str(i), names[1] + str(i), names[2] + str(i)])
    names = np.append(names, name_hold)
    
Ranking_holder = pd.DataFrame(columns = names, index = range(N))

for i in range(N):
    
    gals[i] = find_nearest(L_rank, dummies[i])  ## Picks out galaxies from the cumulative luminosity distribution
       
    a[i] = (find_nearest(cumuL, (RandL[i])))
    if a[i] == len(L1):
        a[i] = len(L1) - 1
    
    b[i] = 10**(L1[int(a[i])])
    received_luminosity[i] = Luminosity_for_convert((b[i]), df_dL[int(gals[i])][0])
    ## Takes dummy luminosity and converts it to luminosity at source by using the luminosity distance of
    ## the host galaxy
    
    current_i = indices.index(gals[i])
    
    testr = np.array(test_ra.iloc[[current_i]].values.tolist())
    testd = np.array(test_dec.iloc[[current_i]].values.tolist())
    ## Extracting data about the host
    
    ##Rotation of test ra and dec
####################################################################################################    
    xyz[i, :] = spherical_convert((50), (10))
    
    m[i, :] = np.cross(xyz[i, :], thph[:, i])#Defines an orthogonal axis
    m_mod = modulus(m[i, :])
    m[i, :] /= m_mod #Normalises orthoganal axis
        
    x_prime = axis_rotation(m[i, :], xyz[i, :], rotation_angle[i]) ##Rotates about an axis
    
    xmod = modulus(x_prime)
    x_prime /= xmod
    ra_prime[i], dec_prime[i] = back_convert(x_prime)
    
    ra_prime[i] = testr[0][0] + (ra_prime[i] - 50)
    dec_prime[i] = testd[0][0] + (dec_prime[i] - 10)

###################################################################################################

    #ident = np.zeros(df_master.shape[0])
      
    print(str(i + 1), "out of " + str(N))
    print("Test galaxy: ", str(gals[i]))
    
    #ident[current_i] = 1
    #df_master["Identifier"] = ident  ## Creates a mask for identifying the host galaxy
    
    
    #q, t, df_sliced = reduction(abs(ra_prime[i]), dec_prime[i], df_master) ## Reduces the catalog by RA and dec
    
    ''''My new function'''
    #selects the corresponding sectors to look through
    df_sliced = Sector_find_reduced(ra_prime[i], dec_prime[i], error_radius)
    df_sliced = df_sliced.rename(columns = {"Unnamed: 0.1": "Unnamed: 0"})
    
    #creates a mask to identify the host galaxy, the host having an identifier of 1
    ident = np.zeros(df_sliced.shape[0])
    df_sliced["Identifier"] = ident
    df_sliced.at[current_i, "Identifier"] = 1
    
    #if statement resolves an issue where sometimes the host galaxy has its info corrupted
    if math.isnan(df_sliced.loc[current_i][ "RA"]) == True:
        '''
        checks if the position data is corrupted, if so then it retrives the information
        from the master file. The only thing that isn't recovered is the sector but 
        that won't really matter, plus I can grab that if it is needed
        '''
        common = df_sliced.columns & df_master.columns
        x = df_master.loc[current_i]
        df_sliced.at[current_i, common] = list(x)
    
    
    
    
    ra = np.array(df_sliced[["RA"]].values.tolist())[:, 0]
    dec = np.array(df_sliced[["dec"]].values.tolist())[:, 0]
    
    Luminosity = np.array(df_sliced[["B Luminosity"]].values.tolist()) #Luminosity_Handling(np.array(df_sliced[["Absolute B Magnitude"]].values.tolist())) ## Converts A    
    dl = np.array(df_sliced[["Luminosity Distance"]].values.tolist())    
    
    lum_prob, SGR_test = L_func(received_luminosity[i], c, dl) ##Uses the luminosity function to calculate probabilities
    df_sliced["Luminosity Probability"] = lum_prob
    df_sliced["SGR flag"] = SGR_test
    
    angular_distaance = np.zeros(df_sliced.shape[0])
    
    for k in range(df_sliced.shape[0]):
        angular_distaance[k] = Ang_Dist(ra[k], ra_prime[i], dec[k], dec_prime[i])
  
    id_check = [i for i, val in enumerate(angular_distaance) if math.isnan(val) == True]
    
    for k in range(len(id_check)):
        angular_distaance[int(id_check[k])] = Ang_Dist(ra_prime[i], testr, dec_prime[i], testd)
  
    angular_distance = Ang_Dist(ra, testr[0][0], dec, testd[0][0])
    
    # Spit out comparison ra and dec
    # Sky position and true luminosity 
    # We might find that knowing the details might help better interpret the results
    # Test revisions
    
    df_sliced["Angular Distance"] = angular_distaance

    ranking = rank(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank"] = ranking

    Ranking_holder.loc[i, "Max Rank"] = max(ranking)
    Ranking_holder.loc[i, "Sum Ranks"] = np.sum(ranking)
    
    x = -np.sort(-ranking)
    Ranking_holder.loc[i, "Top 5 Avg"] = np.mean(x[:5])
    
    ranking2 = rank2(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank2"] = ranking2
    
    Ranking_holder.loc[i, "Max Rank2"] = max(ranking2)
    Ranking_holder.loc[i, "Sum Ranks2"] = np.sum(ranking2)
        
    x = -np.sort(-ranking2)
    Ranking_holder.loc[i, "Top 5 Avg2"] = np.mean(x[:5])
    
    ranking3 = rank3(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank3"] = ranking3
        
    Ranking_holder.loc[i, "Max Rank3"] = max(ranking3)
    Ranking_holder.loc[i, "Sum Ranks3"] = np.sum(ranking3)
            
    x = -np.sort(-ranking3)
    Ranking_holder.loc[i, "Top 5 Avg3"] = np.mean(x[:5])
    
    ranking4 = rank4(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank4"] = ranking4
        
    Ranking_holder.loc[i, "Max Rank4"] = max(ranking4)
    Ranking_holder.loc[i, "Sum Ranks4"] = np.sum(ranking4)
            
    x = -np.sort(-ranking4)
    Ranking_holder.loc[i, "Top 5 Avg4"] = np.mean(x[:5])
    
    ranking5 = rank5(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank5"] = ranking5
        
    Ranking_holder.loc[i, "Max Rank5"] = max(ranking5)
    Ranking_holder.loc[i, "Sum Ranks5"] = np.sum(ranking5)
            
    x = -np.sort(-ranking5)
    Ranking_holder.loc[i, "Top 5 Avg5"] = np.mean(x[:5])
    
    ranking6 = rank6(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank6"] = ranking6
        
    Ranking_holder.loc[i, "Max Rank6"] = max(ranking6)
    Ranking_holder.loc[i, "Sum Ranks6"] = np.sum(ranking6)
            
    x = -np.sort(-ranking6)
    Ranking_holder.loc[i, "Top 5 Avg6"] = np.mean(x[:5])
    
    ranking7 = rank7(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank7"] = ranking7
        
    Ranking_holder.loc[i, "Max Rank7"] = max(ranking7)
    Ranking_holder.loc[i, "Sum Ranks7"] = np.sum(ranking7)
            
    x = -np.sort(-ranking7)
    Ranking_holder.loc[i, "Top 5 Avg7"] = np.mean(x[:5])
    
    ranking8 = rank8(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank8"] = ranking8
    
    Ranking_holder.loc[i, "Max Rank8"] = max(ranking8)
    Ranking_holder.loc[i, "Sum Ranks8"] = np.sum(ranking8)
            
    x = -np.sort(-ranking8)
    Ranking_holder.loc[i, "Top 5 Avg8"] = np.mean(x[:5])
    
    ranking9 = rank9(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank9"] = ranking9
        
    Ranking_holder.loc[i, "Max Rank9"] = max(ranking9)
    Ranking_holder.loc[i, "Sum Ranks9"] = np.sum(ranking9)
            
    x = -np.sort(-ranking9)
    Ranking_holder.loc[i, "Top 5 Avg9"] = np.mean(x[:5])
    
    ranking10 = rank10(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank10"] = ranking10
        
    Ranking_holder.loc[i, "Max Rank10"] = max(ranking10)
    Ranking_holder.loc[i, "Sum Ranks10"] = np.sum(ranking10)
            
    x = -np.sort(-ranking10)
    Ranking_holder.loc[i, "Top 5 Avg10"] = np.mean(x[:5])
    
    ranking11 = rank11(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank11"] = ranking11
        
    Ranking_holder.loc[i, "Max Rank11"] = max(ranking11)
    Ranking_holder.loc[i, "Sum Ranks11"] = np.sum(ranking11)
            
    x = -np.sort(-ranking11)
    Ranking_holder.loc[i, "Top 5 Avg11"] = np.mean(x[:5])
    
    ranking12 = rank12(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank12"] = ranking12
        
    Ranking_holder.loc[i, "Max Rank12"] = max(ranking12)
    Ranking_holder.loc[i, "Sum Ranks12"] = np.sum(ranking12)
            
    x = -np.sort(-ranking12)
    Ranking_holder.loc[i, "Top 5 Avg12"] = np.mean(x[:5])
    
    ranking13 = rank13(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank13"] = ranking13
        
    Ranking_holder.loc[i, "Max Rank13"] = max(ranking13)
    Ranking_holder.loc[i, "Sum Ranks13"] = np.sum(ranking13)
            
    x = -np.sort(-ranking13)
    Ranking_holder.loc[i, "Top 5 Avg13"] = np.mean(x[:5])
    
    ranking14 = rank14(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank14"] = ranking14
    
    Ranking_holder.loc[i, "Max Rank14"] = max(ranking14)
    Ranking_holder.loc[i, "Sum Ranks14"] = np.sum(ranking14)
            
    x = -np.sort(-ranking14)
    Ranking_holder.loc[i, "Top 5 Avg14"] = np.mean(x[:5])
    
    ranking15 = rank15(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank15"] = ranking15
        
    Ranking_holder.loc[i, "Max Rank15"] = max(ranking15)
    Ranking_holder.loc[i, "Sum Ranks15"] = np.sum(ranking15)
            
    x = -np.sort(-ranking15)
    Ranking_holder.loc[i, "Top 5 Avg15"] = np.mean(x[:5])
    
    ranking16 = rank16(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank16"] = ranking16
        
    Ranking_holder.loc[i, "Max Rank16"] = max(ranking16)
    Ranking_holder.loc[i, "Sum Ranks16"] = np.sum(ranking16)
            
    x = -np.sort(-ranking16)
    Ranking_holder.loc[i, "Top 5 Avg16"] = np.mean(x[:5])
    
    ranking17 = rank17(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank17"] = ranking17
        
    Ranking_holder.loc[i, "Max Rank17"] = max(ranking17)
    Ranking_holder.loc[i, "Sum Ranks17"] = np.sum(ranking17)
            
    x = -np.sort(-ranking17)
    Ranking_holder.loc[i, "Top 5 Avg17"] = np.mean(x[:5])
    
    ranking18 = rank18(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank18"] = ranking18
        
    Ranking_holder.loc[i, "Max Rank18"] = max(ranking18)
    Ranking_holder.loc[i, "Sum Ranks18"] = np.sum(ranking18)
            
    x = -np.sort(-ranking18)
    Ranking_holder.loc[i, "Top 5 Avg18"] = np.mean(x[:5])
    
    ranking19 = rank19(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank19"] = ranking19
        
    Ranking_holder.loc[i, "Max Rank19"] = max(ranking19)
    Ranking_holder.loc[i, "Sum Ranks19"] = np.sum(ranking19)
            
    x = -np.sort(-ranking19)
    Ranking_holder.loc[i, "Top 5 Avg19"] = np.mean(x[:5])
    
    ranking20 = rank20(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank20"] = ranking20
        
    Ranking_holder.loc[i, "Max Rank20"] = max(ranking20)
    Ranking_holder.loc[i, "Sum Ranks20"] = np.sum(ranking20)
            
    x = -np.sort(-ranking20)
    Ranking_holder.loc[i, "Top 5 Avg20"] = np.mean(x[:5])
    
    ranking21 = rank21(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank21"] = ranking21
    
    Ranking_holder.loc[i, "Max Rank21"] = max(ranking21)
    Ranking_holder.loc[i, "Sum Ranks21"] = np.sum(ranking21)
            
    x = -np.sort(-ranking21)
    Ranking_holder.loc[i, "Top 5 Avg21"] = np.mean(x[:5])
    
    ranking22 = rank22(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank22"] = ranking22
        
    Ranking_holder.loc[i, "Max Rank22"] = max(ranking22)
    Ranking_holder.loc[i, "Sum Ranks22"] = np.sum(ranking22)
            
    x = -np.sort(-ranking22)
    Ranking_holder.loc[i, "Top 5 Avg22"] = np.mean(x[:5])
    
    ranking23 = rank23(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank23"] = ranking23
        
    Ranking_holder.loc[i, "Max Rank23"] = max(ranking23)
    Ranking_holder.loc[i, "Sum Ranks23"] = np.sum(ranking23)
            
    x = -np.sort(-ranking23)
    Ranking_holder.loc[i, "Top 5 Avg23"] = np.mean(x[:5])
    
    ranking24 = rank24(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank24"] = ranking24
        
    Ranking_holder.loc[i, "Max Rank24"] = max(ranking24)
    Ranking_holder.loc[i, "Sum Ranks24"] = np.sum(ranking24)
            
    x = -np.sort(-ranking24)
    Ranking_holder.loc[i, "Top 5 Avg24"] = np.mean(x[:5])
    
    ranking25 = rank25(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank25"] = ranking25
    
    Ranking_holder.loc[i, "Max Rank25"] = max(ranking25)
    Ranking_holder.loc[i, "Sum Ranks25"] = np.sum(ranking25)
            
    x = -np.sort(-ranking25)
    Ranking_holder.loc[i, "Top 5 Avg25"] = np.mean(x[:5])
    
    ranking26 = rank26(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank26"] = ranking26
    
    Ranking_holder.loc[i, "Max Rank26"] = max(ranking26)
    Ranking_holder.loc[i, "Sum Ranks26"] = np.sum(ranking26)
            
    x = -np.sort(-ranking26)
    Ranking_holder.loc[i, "Top 5 Avg26"] = np.mean(x[:5])
    
    ranking27 = rank27(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank27"] = ranking27
    
    Ranking_holder.loc[i, "Max Rank27"] = max(ranking27)
    Ranking_holder.loc[i, "Sum Ranks27"] = np.sum(ranking27)
            
    x = -np.sort(-ranking27)
    Ranking_holder.loc[i, "Top 5 Avg27"] = np.mean(x[:5])
    
    ranking28 = rank28(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank28"] = ranking28
    
    Ranking_holder.loc[i, "Max Rank28"] = max(ranking28)
    Ranking_holder.loc[i, "Sum Ranks28"] = np.sum(ranking28)
            
    x = -np.sort(-ranking28)
    Ranking_holder.loc[i, "Top 5 Avg28"] = np.mean(x[:5])
    
    ranking29 = rank29(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank29"] = ranking28
    
    Ranking_holder.loc[i, "Max Rank29"] = max(ranking29)
    Ranking_holder.loc[i, "Sum Ranks29"] = np.sum(ranking29)
            
    x = -np.sort(-ranking29)
    Ranking_holder.loc[i, "Top 5 Avg29"] = np.mean(x[:5])
    
    fin_ra = np.asarray(df_sliced[["RA"]].values.tolist()); fin_dec = np.asarray(df_sliced[["dec"]].values.tolist())
    ## Storing values and extending the reduced catalog
    
    df_sliced = (pd.DataFrame.sort_values(df_sliced, by = ["Rank"], ascending = False)) ## Orders resultant sliced array
    df_sliced2 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank2"], ascending = False))
    df_sliced3 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank3"], ascending = False))
    df_sliced4 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank4"], ascending = False))
    df_sliced5 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank5"], ascending = False))
    df_sliced6 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank6"], ascending = False))
    df_sliced7 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank7"], ascending = False))
    df_sliced8 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank8"], ascending = False)) ## Orders resultant sliced array
    
    df_sliced9 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank9"], ascending = False))
    
    df_sliced10 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank10"], ascending = False))
    df_sliced11 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank11"], ascending = False))
    
    df_sliced12 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank12"], ascending = False))
    
    df_sliced13 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank13"], ascending = False))
    df_sliced14 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank14"], ascending = False))
    df_sliced15 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank15"], ascending = False))
    df_sliced16 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank16"], ascending = False))
    df_sliced17 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank17"], ascending = False))
    
    df_sliced18 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank18"], ascending = False))
    df_sliced19 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank19"], ascending = False))
    df_sliced20 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank20"], ascending = False))
    df_sliced21 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank21"], ascending = False))
    df_sliced22 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank22"], ascending = False))
    df_sliced23 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank23"], ascending = False))
    
    df_sliced24 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank24"], ascending = False))
    df_sliced25 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank25"], ascending = False))
    df_sliced26 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank26"], ascending = False))
    df_sliced27 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank27"], ascending = False))
    df_sliced28 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank28"], ascending = False))
    df_sliced29 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank29"], ascending = False))
    
    idi = df_sliced[["Identifier"]].values.tolist() ##Mask handling to check for values
    id2 = df_sliced2[["Identifier"]].values.tolist()
    id3 = df_sliced3[["Identifier"]].values.tolist()
    id4 = df_sliced4[["Identifier"]].values.tolist()
    id5 = df_sliced5[["Identifier"]].values.tolist()
    id6 = df_sliced6[["Identifier"]].values.tolist()
    id7 = df_sliced7[["Identifier"]].values.tolist()
    id8 = df_sliced8[["Identifier"]].values.tolist() ##Mask handling to check for values
    
    id9 = df_sliced9[["Identifier"]].values.tolist()
    
    id10 = df_sliced10[["Identifier"]].values.tolist()
    id11 = df_sliced11[["Identifier"]].values.tolist()
    
    id12 = df_sliced12[["Identifier"]].values.tolist()
    
    id13 = df_sliced13[["Identifier"]].values.tolist()
    id14 = df_sliced14[["Identifier"]].values.tolist()
    id15 = df_sliced15[["Identifier"]].values.tolist()
    id16 = df_sliced16[["Identifier"]].values.tolist()
    id17 = df_sliced17[["Identifier"]].values.tolist()
    
    id18 = df_sliced18[["Identifier"]].values.tolist()
    id19 = df_sliced19[["Identifier"]].values.tolist()
    id20 = df_sliced20[["Identifier"]].values.tolist()
    id21 = df_sliced21[["Identifier"]].values.tolist()
    id22 = df_sliced22[["Identifier"]].values.tolist()
    id23 = df_sliced23[["Identifier"]].values.tolist()
    
    id24 = df_sliced24[["Identifier"]].values.tolist()
    id25 = df_sliced25[["Identifier"]].values.tolist() 
    id26 = df_sliced26[["Identifier"]].values.tolist() 
    id27 = df_sliced27[["Identifier"]].values.tolist() 
    id28 = df_sliced28[["Identifier"]].values.tolist()
    id29 = df_sliced29[["Identifier"]].values.tolist()
    
    
    mask_check = [i for i, val in enumerate(idi) if val == [1]]
    mask_check2 = [i for i, val in enumerate(id2) if val == [1]]
    mask_check3 = [i for i, val in enumerate(id3) if val == [1]]
    mask_check4 = [i for i, val in enumerate(id4) if val == [1]]
    mask_check5 = [i for i, val in enumerate(id5) if val == [1]]
    mask_check6 = [i for i, val in enumerate(id6) if val == [1]]
    mask_check7 = [i for i, val in enumerate(id7) if val == [1]]
    mask_check8 = [i for i, val in enumerate(id8) if val == [1]]
    
    mask_check9 = [i for i, val in enumerate(id9) if val == [1]]
    
    mask_check10 = [i for i, val in enumerate(id10) if val == [1]]
    mask_check11 = [i for i, val in enumerate(id11) if val == [1]]
    
    mask_check12 = [i for i, val in enumerate(id12) if val == [1]]
    
    mask_check13 = [i for i, val in enumerate(id13) if val == [1]]
    mask_check14 = [i for i, val in enumerate(id14) if val == [1]]
    mask_check15 = [i for i, val in enumerate(id15) if val == [1]]
    mask_check16 = [i for i, val in enumerate(id16) if val == [1]]
    mask_check17 = [i for i, val in enumerate(id17) if val == [1]]
    
    mask_check18 = [i for i, val in enumerate(id18) if val == [1]]
    mask_check19 = [i for i, val in enumerate(id19) if val == [1]]
    mask_check20 = [i for i, val in enumerate(id20) if val == [1]]
    mask_check21 = [i for i, val in enumerate(id21) if val == [1]]
    mask_check22 = [i for i, val in enumerate(id22) if val == [1]]
    mask_check23 = [i for i, val in enumerate(id23) if val == [1]]
    
    mask_check24 = [i for i, val in enumerate(id24) if val == [1]]
    mask_check25 = [i for i, val in enumerate(id25) if val == [1]]
    mask_check26 = [i for i, val in enumerate(id26) if val == [1]]
    mask_check27 = [i for i, val in enumerate(id27) if val == [1]]
    mask_check28 = [i for i, val in enumerate(id28) if val == [1]]
    mask_check29 = [i for i, val in enumerate(id29) if val == [1]]
    
    Luminosity =  np.asarray(Luminosity)
    
    if len(mask_check20) == 0:
        print("Did not place\n\n\n")
        next
    else:
        length = len(id20) + 1
        
        placement[i] = mask_check[0] + 1; length = len(idi) + 1
        placement2[i] = mask_check2[0] + 1
        placement3[i] = mask_check3[0] + 1
        placement4[i] = mask_check4[0] + 1
        placement5[i] = mask_check5[0] + 1
        placement6[i] = mask_check6[0] + 1
        placement7[i] = mask_check7[0] + 1
        placement8[i] = mask_check8[0] + 1
        
        placement9[i] = mask_check9[0] + 1
        
        placement10[i] = mask_check10[0] + 1
        placement11[i] = mask_check11[0] + 1
        
        placement12[i] = mask_check12[0] + 1
        
        placement13[i] = mask_check13[0] + 1
        placement14[i] = mask_check14[0] + 1
        placement15[i] = mask_check15[0] + 1
        placement16[i] = mask_check16[0] + 1
        placement17[i] = mask_check17[0] + 1
        
        placement18[i] = mask_check18[0] + 1
        placement19[i] = mask_check19[0] + 1
        placement20[i] = mask_check20[0] + 1
        placement21[i] = mask_check21[0] + 1
        placement22[i] = mask_check22[0] + 1
        placement23[i] = mask_check23[0] + 1
        
        placement24[i] = mask_check24[0] + 1
        placement25[i] = mask_check25[0] + 1
        placement26[i] = mask_check26[0] + 1
        placement27[i] = mask_check27[0] + 1
        placement28[i] = mask_check28[0] + 1
        placement29[i] = mask_check29[0] + 1
        
        #display(Markdown("The keplerian orbit appears to be happening at r ={0:.2f} km"  .format(float(kepler(M_kep, w))/1000)))
        print("Galaxy data: \nDistance is {0:.2f} Mpc\nLuminosity is {1:.3e}\nra and dec [{2:.2f}, {3:.2f}] compared to reported ra and dec [{4:.2f}, {5:.2f}] \nTrue luminosity {6:.3e} W" .format(dl[int(placement[i] - 1)][0], Luminosity[int(placement[i] - 1)][0], fin_ra[int(placement[i] - 1)][0], fin_dec[int(placement[i] - 1)][0], testr[0][0], testd[0][0], b[i]))
        
        print("Galaxy placed", int(placement[i]), "out of", str(length), "with statistic 1\n\n\n")
        print("Galaxy placed", int(placement2[i]), "out of", str(length), "with statistic 2\n\n\n")
        print("Galaxy placed", int(placement3[i]), "out of", str(length), "with statistic 3\n\n\n")        
        print("Galaxy placed", int(placement4[i]), "out of", str(length), "with statistic 4\n\n\n")
        print("Galaxy placed", int(placement5[i]), "out of", str(length), "with statistic 5\n\n\n")
        print("Galaxy placed", int(placement6[i]), "out of", str(length), "with statistic 6\n\n\n")
        print("Galaxy placed", int(placement7[i]), "out of", str(length), "with statistic 7\n\n\n")
        print("Galaxy placed", int(placement8[i]), "out of", str(length), "with statistic 8\n\n\n")
        
        print("Galaxy placed", int(placement9[i]), "out of", str(length), "with statistic 9\n\n\n")
        
        print("Galaxy placed", int(placement10[i]), "out of", str(length), "with statistic 10\n\n\n")
        print("Galaxy placed", int(placement11[i]), "out of", str(length), "with statistic 11\n\n\n")
        
        print("Galaxy placed", int(placement12[i]), "out of", str(length), "with statistic 12\n\n\n")
        
        print("Galaxy placed", int(placement13[i]), "out of", str(length), "with statistic 13\n\n\n")
        print("Galaxy placed", int(placement14[i]), "out of", str(length), "with statistic 14\n\n\n")
        print("Galaxy placed", int(placement15[i]), "out of", str(length), "with statistic 15\n\n\n")
        print("Galaxy placed", int(placement16[i]), "out of", str(length), "with statistic 16\n\n\n")
        print("Galaxy placed", int(placement17[i]), "out of", str(length), "with statistic 17\n\n\n")
        
        print("Galaxy placed", int(placement18[i]), "out of", str(length), "with statistic 18\n\n\n")
        print("Galaxy placed", int(placement19[i]), "out of", str(length), "with statistic 19\n\n\n")
        print("Galaxy placed", int(placement20[i]), "out of", str(length), "with statistic 20\n\n\n")
        print("Galaxy placed", int(placement21[i]), "out of", str(length), "with statistic 21\n\n\n")
        print("Galaxy placed", int(placement22[i]), "out of", str(length), "with statistic 22\n\n\n")
        print("Galaxy placed", int(placement23[i]), "out of", str(length), "with statistic 23\n\n\n")
        
        print("Galaxy placed", int(placement24[i]), "out of", str(length), "with statistic 24\n\n\n")
        print("Galaxy placed", int(placement25[i]), "out of", str(length), "with statistic 25\n\n\n")
        print("Galaxy placed", int(placement26[i]), "out of", str(length), "with statistic 26\n\n\n")
        print("Galaxy placed", int(placement27[i]), "out of", str(length), "with statistic 27\n\n\n")
        print("Galaxy placed", int(placement28[i]), "out of", str(length), "with statistic 28\n\n\n")
        print("Galaxy placed", int(placement29[i]), "out of", str(length), "with statistic 29\n\n\n")
        
        percentages[i] = placement[i]/length 
        percentages2[i] = placement2[i]/length 
        percentages3[i] = placement3[i]/length 
        percentages4[i] = placement4[i]/length 
        percentages5[i] = placement5[i]/length 
        percentages6[i] = placement6[i]/length 
        percentages7[i] = placement7[i]/length 
        percentages8[i] = placement8[i]/length 
        
        percentages9[i] = placement9[i]/length 
        
        percentages10[i] = placement10[i]/length 
        percentages11[i] = placement11[i]/length 
        
        percentages12[i] = placement12[i]/length 
        
        percentages13[i] = placement13[i]/length 
        percentages14[i] = placement14[i]/length 
        percentages15[i] = placement15[i]/length 
        percentages16[i] = placement16[i]/length
        percentages17[i] = placement17[i]/length 
        
        percentages18[i] = placement18[i]/length 
        percentages19[i] = placement19[i]/length 
        percentages20[i] = placement20[i]/length
        percentages21[i] = placement21[i]/length            
        percentages22[i] = placement22[i]/length            
        percentages23[i] = placement23[i]/length            
        
        percentages24[i] = placement24[i]/length
        percentages25[i] = placement25[i]/length
        percentages26[i] = placement26[i]/length
        percentages27[i] = placement27[i]/length
        percentages28[i] = placement28[i]/length
        percentages29[i] = placement29[i]/length
        
        distances[i] = int(dl[int(placement[i]) - 1][0]); luminosity_i[i] = int(Luminosity[int(placement[i]) - 1][0])
        ras_dex[i, 0] = fin_ra[int(placement[i] - 1)]; ras_dex[i, 1] = fin_dec[int(placement[i] - 1)]; test_case[i, 0] = testr[0][0]; test_case[i, 1] = testd[0][0]
        
        #rank_host[i] = df_sliced20[["Rank20"]].values.tolist()[id20.index(max(id20))][0]

    faulty[i, 0] = df_master[["RA"]].values.tolist()[current_i][0]      #ra of galaxy
    faulty[i, 1] = ra_prime[i]                                          #ra of grb
    faulty[i, 2] = df_master[["dec"]].values.tolist()[current_i][0]     #dec of galaxy
    faulty[i, 3] = dec_prime[i]                                         #dec of grb
    
    if math.isnan(rank_host[i]) == True:    
        faulty[i, 4] = 1 #Mask
        no_se_func.append(i)
        #break
    
    else:
        faulty[i, 4] = 0 #Mask
        next
  


#Saving the ranking number data to a csv file

Ranking_holder.to_csv("Max and Sum of rank values within 250 Mpc.csv", header = True, index = False)


"""
f_v = [i for i, val in enumerate(faulty[:, 4]) if val == 0]
f_1v = [i for i, val in enumerate(faulty[:, 4]) if val == 1]

sets = set(np.arange(0, len(faulty), 1)) - set(f_v)
ft = pd.DataFrame(faulty)
faulty_cols = ["Galaxy RA", "GRB RA", "Galaxy dec", "GRB dec", "Mask"]
ft.columns = faulty_cols

ab_fault = ft.take(list(sets), axis = 0)
ab_vals = ab_fault.values.tolist()[0]
"""
place_array = np.zeros(shape = (N, 29))
place_array[:, 0] = percentages
place_array[:, 1] = percentages2
place_array[:, 2] = percentages3
place_array[:, 3] = percentages4
place_array[:, 4] = percentages5
place_array[:, 5] = percentages6
place_array[:, 6] = percentages7
place_array[:, 7] = percentages8
place_array[:, 8] = percentages9
place_array[:, 9] = percentages10
place_array[:, 10] = percentages11
place_array[:, 11] = percentages12
place_array[:, 12] = percentages13
place_array[:, 13] = percentages14
place_array[:, 14] = percentages15
place_array[:, 15] = percentages16
place_array[:, 16] = percentages17
place_array[:, 17] = percentages18
place_array[:, 18] = percentages19
place_array[:, 19] = percentages20
place_array[:, 20] = percentages21
place_array[:, 21] = percentages22
place_array[:, 22] = percentages23

place_array[:, 23] = percentages24
place_array[:, 24] = percentages25
place_array[:, 25] = percentages26
place_array[:, 26] = percentages27
place_array[:, 27] = percentages28
place_array[:, 28] = percentages29

zeros = [i for i, val in enumerate(place_array[:, 28]) if val == 0]
df_place_array = pd.DataFrame(place_array)

plus_one = [i for i, val in enumerate(place_array[:, 28]) if val > 0.9]

indices_to_keep = set(range(df_place_array.shape[0])) - set(zeros) - set(plus_one) #- set(no_se_func)
df_place_array = np.asarray(pd.DataFrame.take(df_place_array, list(indices_to_keep), axis = 0).values.tolist())

df_dist = pd.DataFrame(distances)
df_distance = np.asarray(pd.DataFrame.take(df_dist, list(indices_to_keep), axis = 0).values.tolist())

df_ang = pd.DataFrame(rotation_angle)
df_ang = np.asarray(pd.DataFrame.take(df_ang, list(indices_to_keep), axis = 0).values.tolist())

df_lumin = pd.DataFrame(b)
df_lumin = np.asarray(pd.DataFrame.take(df_lumin, list(indices_to_keep), axis = 0).values.tolist())
"""
plt.figure(3)
for p in range(20):
    plt.plot(df_place_array[:, p], np.log10(df_distance), "x", alpha = 2/(p/2 + 1), label = "Statistic" + str(p))
    
plt.title("Distance vs. percentage performance")
plt.ylabel("Log$_{10}$ Distance /Mpc"); plt.xlabel("Percentage placement"); plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("Statistic_Comparison.png")
"""



### The following can be used to investigate any values that flag up as false

f_v = [i for i, val in enumerate(faulty[:, 4]) if val == 0]
f_1v = [i for i, val in enumerate(faulty[:, 4]) if val == 1]

sets = set(np.arange(0, len(faulty), 1)) - set(f_v)
ft = pd.DataFrame(faulty)
faulty_cols = ["Galaxy RA", "GRB RA", "Galaxy dec", "GRB dec", "Mask"]
ft.columns = faulty_cols

place_array = np.zeros(shape = (N, 29))
place_array[:, 0] = percentages
place_array[:, 1] = percentages2
place_array[:, 2] = percentages3
place_array[:, 3] = percentages4
place_array[:, 4] = percentages5
place_array[:, 5] = percentages6
place_array[:, 6] = percentages7
place_array[:, 7] = percentages8
place_array[:, 8] = percentages9
place_array[:, 9] = percentages10
place_array[:, 10] = percentages11
place_array[:, 11] = percentages12
place_array[:, 12] = percentages13
place_array[:, 13] = percentages14
place_array[:, 14] = percentages15
place_array[:, 15] = percentages16
place_array[:, 16] = percentages17
place_array[:, 17] = percentages18
place_array[:, 18] = percentages19
place_array[:, 19] = percentages20
place_array[:, 20] = percentages21
place_array[:, 21] = percentages22
place_array[:, 22] = percentages23
place_array[:, 23] = percentages24
place_array[:, 24] = percentages25
place_array[:, 25] = percentages26
place_array[:, 26] = percentages27
place_array[:, 27] = percentages28
place_array[:, 28] = percentages29

zeros = [i for i, val in enumerate(place_array[:, 19]) if val == 0]
df_place_array = pd.DataFrame(place_array)

plus_one = [i for i, val in enumerate(place_array[:, 19]) if val > 0.9]

indices_to_keep = set(range(df_place_array.shape[0])) - set(zeros) - set(plus_one) #- set(no_se_func)
df_place_array = np.asarray(pd.DataFrame.take(df_place_array, list(indices_to_keep), axis = 0).values.tolist())

df_dist = pd.DataFrame(distances)
df_distance = np.asarray(pd.DataFrame.take(df_dist, list(indices_to_keep), axis = 0).values.tolist())

df_ang = pd.DataFrame(rotation_angle)
df_ang = np.asarray(pd.DataFrame.take(df_ang, list(indices_to_keep), axis = 0).values.tolist())

df_lumin = pd.DataFrame(b)
df_lumin = np.asarray(pd.DataFrame.take(df_lumin, list(indices_to_keep), axis = 0).values.tolist())

rankN = np.zeros(shape = (len(df_place_array), 29))
for i in range(len(df_place_array)):
    
    df_array_init = pd.DataFrame(df_place_array[i, :]) ## Takes percentage placement for each run
    
    counting_mask = np.arange(df_array_init.shape[0])
    df_array_init["Mask"] = counting_mask    ## Creates a matching mask for keeping track of where the entries end up
    
    df_array = (pd.DataFrame.sort_values(df_array_init, by = [0], ascending = True)) ## Orders resultant sliced array
    
    for k in range(df_array.shape[0]):  
        rankN[i, k] = [i for i, val in enumerate(df_array[["Mask"]].values.tolist()) if val == [k]][0] ## 
counter = 5

'''
for p in range(29):
    df_rank = pd.DataFrame(rankN[:, p])
    
    plt.figure(p + 4)
    
    val = df_rank[0].value_counts()
    vals = df_rank[0].value_counts().values.tolist()
    quantities = np.zeros(29)
    idq = val.index.values.tolist()

    for j in range(len(vals)):
        quantities[int(idq[j])] = vals[j]
    
    for o in range(29):
        plt.bar((o + 1), quantities[o], color = "black")
    
    plt.xlabel("Placement"); plt.ylabel("Frequency")
    plt.title("Statistic " + str(p + 1))
    plt.grid()
    plt.savefig("Statistic " + str(p + 1) + ".png")
    counter += 1
    
for i in range(29):
    plt.figure(counter)
    plt.plot(np.log10(df_distance), df_place_array[:, i], "kx", label = "Statistic " + str(i + 1))
    plt.ylabel("Percentage performance")
    plt.xlabel("Log$_{10}$ Distance /Mpc")
    plt.grid()
    plt.legend(loc = "best")
    plt.savefig("OmittedGalaxies_Statistic" + str(i + 1) + ".png")
    counter += 1
    
for j in range(29):
    plt.figure(counter)
    plt.plot(np.log10(df_lumin), df_place_array[:, j], "kx", label = "Statistic " + str(j + 1))
    plt.ylabel("Percentage performance")
    plt.xlabel("Log$_{10}$ Luminosity /W")
    plt.grid()
    plt.legend(loc = "best")
    plt.savefig("OmittedGalaxies_Lumin_Statistic" + str(j + 1) + ".png")
    counter += 1
    
for k in range(29):
    plt.figure(counter)
    plt.plot((df_ang), df_place_array[:, k], "kx", label = "Statistic " + str(k + 1))
    plt.ylabel("Percentage performance")
    plt.xlabel("Angular Offset /$^o$")
    plt.grid()
    plt.legend(loc = "best")
    plt.savefig("OmittedGalaxies_Ang_Statistic" + str(k + 1) + ".png")
    counter += 1
'''

elapsed_time = timer() - start # in seconds
print('The code took {:.5g} s to complete'.format(elapsed_time))



#%%

'''
This will do the same thing as above but for outside the 250Mpc sphere and so 
we don't have to select a galaxy we just have to generate a random GRB at some 
random location in the sky with a random luminosity
'''

start = timer()

df_master = pd.read_csv("Data Files/GLADE_Master_comma_100Mpc.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

L1 = np.linspace(56, 59, 101) #In J now
L2, c = L_func1(L1) #  ##Builds broken power law

cumuL = cumulative(L2) ##Luminosity Distribution
df_cumLum =   pd.read_csv("Data Files/Cumulative_Luminosity_100Mpc.csv") #calls the cumulative luminosity from the csv
df_cumLum.columns = ["NaN", "Cumulative Luminosity"] 
normal_c = df_cumLum[["Cumulative Luminosity"]].values[-1][0]
L_rank = df_cumLum[["Cumulative Luminosity"]].values * 1/normal_c
df_cumLum = df_cumLum[["Cumulative Luminosity"]].values# This is all to do with building a usable and callable power law
                                                                         

lum_N = np.linspace(0, df_cumLum.shape[0], df_cumLum.shape[0])
df_dL = df_master[["Luminosity Distance"]] #grabbing the luminosity distance from master file

#using HEALPix to split the sky into equal area sectors
hp = HEALPix(nside=16, order='ring', frame=ICRS())

#Creates from Numbers for dummies and Luminosity fraction between 0->1
dummies = random(N)
RandL = random(N)

gals = np.zeros(N)          ## Picks out a luminosity
gal_index = np.zeros(N)

#making empty arrays to store data for later
lum_atsource = np.zeros(N)
received_luminosity = np.zeros(N)
cumul_N = np.zeros(N)

lum_list = list(L_rank)

df_dL = df_dL.values.tolist()   ## Luminosity distance values for use
a = np.zeros(N)                 ## For storing temporary and unimportant values
b = np.zeros(N)                 ## For storing temporary and unimportant values

#grabs the positions of every galaxy in the catalogue
test_ra = df_master[["RA"]]
test_dec = df_master[["dec"]]

indices = list(np.arange(df_master.shape[0]))

error_radius = 2 * (2.62)  ## Change as necessary - this is an example value from HEARSCH 

distances = np.zeros(N)
luminosity_i = np.zeros(N)    
rank_host = np.zeros(N)

'''
This is to produce a random point on the sky
'''
angles = np.arccos(2 * random(N) - 1)

ra_rand = uniform(0, 360, size = N)
dec_rand = (np.pi/2 - angles) * (180/np.pi) ##Gives you random ra and dec

#makes a random distance from us, at least 250 Mpc away
r_min = 250 # Mpc
r = r_min / random(N)**(1/3)

#I want to try and keep the summation of the ranks and the largest rank so that I can work with them later

names_out = np.array(["Max Rank", "Sum Ranks", "Top 5 Avg"], dtype = str)
number_ranks_out = np.arange(2, 30, 1)

for i in number_ranks_out:
    #sets the names of the columns in the rank holder
    name_hold = np.array([names_out[0] + str(i), names_out[1] + str(i)])
    names_out = np.append(names_out, name_hold)
    
Ranking_holder_out = pd.DataFrame(columns = names_out, index = range(N))


for i in range(N):
    '''
    This makes N random luminosities from the luminosity power law for grbs, it
    then matches that to a random distance from here to '''       
    a[i] = find_nearest(cumuL, RandL[i])
    if a[i] == len(L1):
        a[i] = len(L1) - 1
    b[i] = 10**L1[int(a[i])] #int(a[i])
    received_luminosity[i] = Luminosity_for_convert(b[i], r[i])

error_radius = 2 * (2.62)

ranks = np.zeros(shape = (N, 5))

ra_prime = np.zeros(N); dec_prime = np.zeros(N)#1612646.0

#keep track of the progress
count = 0

for i in range(N):
    phi = 2 * np.pi * random()
    theta = np.arccos(2 * random() - 1)
    
    thph = spherical_convert(theta, phi)
    mod = modulus(thph)
    thph /= mod
        
    xyz = np.transpose(spherical_convert(float(ra_rand[i]), float(dec_rand[i])))
    m = np.zeros(shape = (N, 3))
    
    #for j in range(1):
        
    m[i, :] = np.transpose(np.cross(xyz, thph))#np.transpose(xyz[i, :] * thph[:, i])
    m_mod = modulus(m[i, :])
    m[i, :] /= m_mod
        
    rotation_angle = error_radius * normal(size = N)

    #for k in range(1):
        #rota = rotation(m[i, :], rotation_angle[i]) ###Ammend this!!
        
        #x_prime = mat_mul(rota, xyz) #rota * xyz[i, :]
        
    x_prime = axis_rotation(m[i, :], xyz, rotation_angle[i])
    
    xmod = modulus(x_prime)
    x_prime /= xmod
    ra_prime[i], dec_prime[i] = back_convert(x_prime)
    
    
    ''''My new function'''
    #selects the corresponding sectors to look through
    df_sliced = Sector_find_reduced(ra_prime[i], dec_prime[i], error_radius)
    df_sliced = df_sliced.rename(columns = {"Unnamed: 0.1": "Unnamed: 0"})
     
 
    
    ra = np.array(df_sliced[["RA"]].values.tolist())[:, 0]
    dec = np.array(df_sliced[["dec"]].values.tolist())[:, 0]
    
    Luminosity = np.array(df_sliced[["B Luminosity"]].values.tolist()) #Luminosity_Handling(np.array(df_sliced[["Absolute B Magnitude"]].values.tolist())) ## Converts A    
    dl = np.array(df_sliced[["Luminosity Distance"]].values.tolist())    
    
    lum_prob, SGR_test = L_func(received_luminosity[i], c, dl) ##Uses the luminosity function to calculate probabilities
    df_sliced["Luminosity Probability"] = lum_prob
    df_sliced["SGR flag"] = SGR_test
    
    angular_distaance = np.zeros(df_sliced.shape[0])
    
    
    for k in range(df_sliced.shape[0]):
        angular_distaance[k] = Ang_Dist(ra[k], ra_prime[i], dec[k], dec_prime[i])
    
    # Spit out comparison ra and dec
    # Sky position and true luminosity 
    # We might find that knowing the details might help better interpret the results
    # Test revisions
    
    df_sliced["Angular Distance"] = angular_distaance

    ranking = rank(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank"] = ranking

    Ranking_holder_out.loc[i, "Max Rank"] = max(ranking)
    Ranking_holder_out.loc[i, "Sum Ranks"] = np.sum(ranking)
    
    x = -np.sort(-ranking)
    Ranking_holder_out.loc[i, "Top 5 Avg"] = np.mean(x[:5])
    
    ranking2 = rank2(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank2"] = ranking2
    
    Ranking_holder_out.loc[i, "Max Rank2"] = max(ranking2)
    Ranking_holder_out.loc[i, "Sum Ranks2"] = np.sum(ranking2)
        
    x = -np.sort(-ranking2)
    Ranking_holder_out.loc[i, "Top 5 Avg2"] = np.mean(x[:5])
    
    ranking3 = rank3(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank3"] = ranking3
        
    Ranking_holder_out.loc[i, "Max Rank3"] = max(ranking3)
    Ranking_holder_out.loc[i, "Sum Ranks3"] = np.sum(ranking3)
        
    x = -np.sort(-ranking3)
    Ranking_holder_out.loc[i, "Top 5 Avg3"] = np.mean(x[:5])
    
    ranking4 = rank4(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank4"] = ranking4
        
    Ranking_holder_out.loc[i, "Max Rank4"] = max(ranking4)
    Ranking_holder_out.loc[i, "Sum Ranks4"] = np.sum(ranking4)
        
    x = -np.sort(-ranking4)
    Ranking_holder_out.loc[i, "Top 5 Avg4"] = np.mean(x[:5])
    
    ranking5 = rank5(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank5"] = ranking5
        
    Ranking_holder_out.loc[i, "Max Rank5"] = max(ranking5)
    Ranking_holder_out.loc[i, "Sum Ranks5"] = np.sum(ranking5)
        
    x = -np.sort(-ranking5)
    Ranking_holder_out.loc[i, "Top 5 Avg5"] = np.mean(x[:5])
    
    ranking6 = rank6(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank6"] = ranking6
        
    Ranking_holder_out.loc[i, "Max Rank6"] = max(ranking6)
    Ranking_holder_out.loc[i, "Sum Ranks6"] = np.sum(ranking6)
        
    x = -np.sort(-ranking6)
    Ranking_holder_out.loc[i, "Top 5 Avg6"] = np.mean(x[:5])
    
    ranking7 = rank7(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank7"] = ranking7
        
    Ranking_holder_out.loc[i, "Max Rank7"] = max(ranking7)
    Ranking_holder_out.loc[i, "Sum Ranks7"] = np.sum(ranking7)
        
    x = -np.sort(-ranking7)
    Ranking_holder_out.loc[i, "Top 5 Avg7"] = np.mean(x[:5])
    
    ranking8 = rank8(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank8"] = ranking8
    
    Ranking_holder_out.loc[i, "Max Rank8"] = max(ranking8)
    Ranking_holder_out.loc[i, "Sum Ranks8"] = np.sum(ranking8)
        
    x = -np.sort(-ranking8)
    Ranking_holder_out.loc[i, "Top 5 Avg8"] = np.mean(x[:5])
    
    ranking9 = rank9(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank9"] = ranking9
        
    Ranking_holder_out.loc[i, "Max Rank9"] = max(ranking9)
    Ranking_holder_out.loc[i, "Sum Ranks9"] = np.sum(ranking9)
        
    x = -np.sort(-ranking9)
    Ranking_holder_out.loc[i, "Top 5 Avg9"] = np.mean(x[:5])
    
    ranking10 = rank10(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank10"] = ranking10
        
    Ranking_holder_out.loc[i, "Max Rank10"] = max(ranking10)
    Ranking_holder_out.loc[i, "Sum Ranks10"] = np.sum(ranking10)
        
    x = -np.sort(-ranking10)
    Ranking_holder_out.loc[i, "Top 5 Avg10"] = np.mean(x[:5])
    
    ranking11 = rank11(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank11"] = ranking11
        
    Ranking_holder_out.loc[i, "Max Rank11"] = max(ranking11)
    Ranking_holder_out.loc[i, "Sum Ranks11"] = np.sum(ranking11)
        
    x = -np.sort(-ranking11)
    Ranking_holder_out.loc[i, "Top 5 Avg11"] = np.mean(x[:5])
    
    ranking12 = rank12(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank12"] = ranking12
        
    Ranking_holder_out.loc[i, "Max Rank12"] = max(ranking12)
    Ranking_holder_out.loc[i, "Sum Ranks12"] = np.sum(ranking12)
        
    x = -np.sort(-ranking12)
    Ranking_holder_out.loc[i, "Top 5 Avg12"] = np.mean(x[:5])
    
    ranking13 = rank13(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank13"] = ranking13
        
    Ranking_holder_out.loc[i, "Max Rank13"] = max(ranking13)
    Ranking_holder_out.loc[i, "Sum Ranks13"] = np.sum(ranking13)
        
    x = -np.sort(-ranking13)
    Ranking_holder_out.loc[i, "Top 5 Avg13"] = np.mean(x[:5])
    
    ranking14 = rank14(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank14"] = ranking14
    
    Ranking_holder_out.loc[i, "Max Rank14"] = max(ranking14)
    Ranking_holder_out.loc[i, "Sum Ranks14"] = np.sum(ranking14)
        
    x = -np.sort(-ranking14)
    Ranking_holder_out.loc[i, "Top 5 Avg14"] = np.mean(x[:5])
    
    ranking15 = rank15(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank15"] = ranking15
        
    Ranking_holder_out.loc[i, "Max Rank15"] = max(ranking15)
    Ranking_holder_out.loc[i, "Sum Ranks15"] = np.sum(ranking15)
        
    x = -np.sort(-ranking15)
    Ranking_holder_out.loc[i, "Top 5 Avg15"] = np.mean(x[:5])
    
    ranking16 = rank16(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank16"] = ranking16
        
    Ranking_holder_out.loc[i, "Max Rank16"] = max(ranking16)
    Ranking_holder_out.loc[i, "Sum Ranks16"] = np.sum(ranking16)
        
    x = -np.sort(-ranking16)
    Ranking_holder_out.loc[i, "Top 5 Avg16"] = np.mean(x[:5])
    
    ranking17 = rank17(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank17"] = ranking17
        
    Ranking_holder_out.loc[i, "Max Rank17"] = max(ranking17)
    Ranking_holder_out.loc[i, "Sum Ranks17"] = np.sum(ranking17)
        
    x = -np.sort(-ranking17)
    Ranking_holder_out.loc[i, "Top 5 Avg17"] = np.mean(x[:5])
    
    ranking18 = rank18(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank18"] = ranking18
        
    Ranking_holder_out.loc[i, "Max Rank18"] = max(ranking18)
    Ranking_holder_out.loc[i, "Sum Ranks18"] = np.sum(ranking18)
        
    x = -np.sort(-ranking18)
    Ranking_holder_out.loc[i, "Top 5 Avg18"] = np.mean(x[:5])
    
    ranking19 = rank19(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank19"] = ranking19
        
    Ranking_holder_out.loc[i, "Max Rank19"] = max(ranking19)
    Ranking_holder_out.loc[i, "Sum Ranks19"] = np.sum(ranking19)
        
    x = -np.sort(-ranking19)
    Ranking_holder_out.loc[i, "Top 5 Avg19"] = np.mean(x[:5])
    
    ranking20 = rank20(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank20"] = ranking20
        
    Ranking_holder_out.loc[i, "Max Rank20"] = max(ranking20)
    Ranking_holder_out.loc[i, "Sum Ranks20"] = np.sum(ranking20)
        
    x = -np.sort(-ranking20)
    Ranking_holder_out.loc[i, "Top 5 Avg20"] = np.mean(x[:5])
    
    ranking21 = rank21(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank21"] = ranking21
    
    Ranking_holder_out.loc[i, "Max Rank21"] = max(ranking21)
    Ranking_holder_out.loc[i, "Sum Ranks21"] = np.sum(ranking21)
        
    x = -np.sort(-ranking21)
    Ranking_holder_out.loc[i, "Top 5 Avg21"] = np.mean(x[:5])
    
    ranking22 = rank22(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank22"] = ranking22
        
    Ranking_holder_out.loc[i, "Max Rank22"] = max(ranking22)
    Ranking_holder_out.loc[i, "Sum Ranks22"] = np.sum(ranking22)
        
    x = -np.sort(-ranking22)
    Ranking_holder_out.loc[i, "Top 5 Avg22"] = np.mean(x[:5])
    
    ranking23 = rank23(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank23"] = ranking23
        
    Ranking_holder_out.loc[i, "Max Rank23"] = max(ranking23)
    Ranking_holder_out.loc[i, "Sum Ranks23"] = np.sum(ranking23)
        
    x = -np.sort(-ranking23)
    Ranking_holder_out.loc[i, "Top 5 Avg23"] = np.mean(x[:5])
    
    ranking24 = rank24(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank24"] = ranking24
        
    Ranking_holder_out.loc[i, "Max Rank24"] = max(ranking24)
    Ranking_holder_out.loc[i, "Sum Ranks24"] = np.sum(ranking24)
        
    x = -np.sort(-ranking24)
    Ranking_holder_out.loc[i, "Top 5 Avg24"] = np.mean(x[:5])
    
    ranking25 = rank25(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank25"] = ranking25
    
    Ranking_holder_out.loc[i, "Max Rank25"] = max(ranking25)
    Ranking_holder_out.loc[i, "Sum Ranks25"] = np.sum(ranking25)
        
    x = -np.sort(-ranking25)
    Ranking_holder_out.loc[i, "Top 5 Avg25"] = np.mean(x[:5])
    
    ranking26 = rank26(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank26"] = ranking26
    
    Ranking_holder_out.loc[i, "Max Rank26"] = max(ranking26)
    Ranking_holder_out.loc[i, "Sum Ranks26"] = np.sum(ranking26)
        
    x = -np.sort(-ranking26)
    Ranking_holder_out.loc[i, "Top 5 Avg26"] = np.mean(x[:5])
    
    ranking27 = rank27(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank27"] = ranking27
    
    Ranking_holder_out.loc[i, "Max Rank27"] = max(ranking27)
    Ranking_holder_out.loc[i, "Sum Ranks27"] = np.sum(ranking27)
        
    x = -np.sort(-ranking27)
    Ranking_holder_out.loc[i, "Top 5 Avg27"] = np.mean(x[:5])
    
    ranking28 = rank28(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank28"] = ranking28
    
    Ranking_holder_out.loc[i, "Max Rank28"] = max(ranking28)
    Ranking_holder_out.loc[i, "Sum Ranks28"] = np.sum(ranking28)
        
    x = -np.sort(-ranking28)
    Ranking_holder_out.loc[i, "Top 5 Avg28"] = np.mean(x[:5])
    
    ranking29 = rank29(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank29"] = ranking29
    
    Ranking_holder_out.loc[i, "Max Rank29"] = max(ranking29)
    Ranking_holder_out.loc[i, "Sum Ranks29"] = np.sum(ranking29)
        
    x = -np.sort(-ranking29)
    Ranking_holder_out.loc[i, "Top 5 Avg29"] = np.mean(x[:5])
    
    fin_ra = np.asarray(df_sliced[["RA"]].values.tolist()); fin_dec = np.asarray(df_sliced[["dec"]].values.tolist())
    ## Storing values and extending the reduced catalog
    
    df_sliced = (pd.DataFrame.sort_values(df_sliced, by = ["Rank"], ascending = False)) ## Orders resultant sliced array
    df_sliced2 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank2"], ascending = False))
    df_sliced3 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank3"], ascending = False))
    df_sliced4 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank4"], ascending = False))
    df_sliced5 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank5"], ascending = False))
    df_sliced6 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank6"], ascending = False))
    df_sliced7 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank7"], ascending = False))
    df_sliced8 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank8"], ascending = False)) ## Orders resultant sliced array
    
    df_sliced9 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank9"], ascending = False))
    
    df_sliced10 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank10"], ascending = False))
    df_sliced11 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank11"], ascending = False))
    
    df_sliced12 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank12"], ascending = False))
    
    df_sliced13 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank13"], ascending = False))
    df_sliced14 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank14"], ascending = False))
    df_sliced15 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank15"], ascending = False))
    df_sliced16 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank16"], ascending = False))
    df_sliced17 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank17"], ascending = False))
    
    df_sliced18 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank18"], ascending = False))
    df_sliced19 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank19"], ascending = False))
    df_sliced20 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank20"], ascending = False))
    df_sliced21 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank21"], ascending = False))
    df_sliced22 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank22"], ascending = False))
    df_sliced23 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank23"], ascending = False))
    
    df_sliced24 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank24"], ascending = False))
    df_sliced25 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank25"], ascending = False))
    df_sliced26 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank26"], ascending = False))
    df_sliced27 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank27"], ascending = False))
    df_sliced28 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank28"], ascending = False))
    df_sliced29 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank29"], ascending = False))
    
    count += 1
    if count % 50 == 0:
        print(count)

elapsed_time = timer() - start # in seconds
print('The code took {:.5g} s to complete'.format(elapsed_time))

Ranking_holder_out.to_csv("Max and Sum of rank values outside.csv", header = True, index = False)

#this sets the number of bins and there edges so that the graphs can be adequately compared
#bin_sides = np.linspace(0, 1, 75)


#%%
"""
This section will look into the histograms and find possibile rank values which 
provide a significant distinction between inside the outside sphere GRBs, it also
compares the relative performance for each statistic, finding the relative position
that certain statistics place X % of the time  
"""

#A distinction value at which we accept and the corresponding cutoff value
distinct = 0.1*N
cutoff = 1 - distinct/N

#opening a txt file to store the info
histtxt = open("Histogram Distinction for cutoff at {:.3g}.txt".format(cutoff), "w")

for j in names:
    #finding the range (upper)
    if max(Ranking_holder[j]) >= max(Ranking_holder_out[j]):
        upper = max(Ranking_holder[j])
    elif max(Ranking_holder[j]) < max(Ranking_holder_out[j]):
        upper = max(Ranking_holder_out[j])
    
    #finding the range (lower)
    if min(Ranking_holder[j]) <= min(Ranking_holder_out[j]):
        lower = min(Ranking_holder[j])
    elif min(Ranking_holder[j]) > min(Ranking_holder_out[j]):
        lower = min(Ranking_holder_out[j])
    
        '''
    #plots Histograms
    plt.figure()

    plt.hist(Ranking_holder[j], bins = 75, range = (lower, upper), color ="blue", alpha =0.75, label = "Within")
    
    plt.hist(Ranking_holder_out[j], bins = 75, range = (lower, upper), color = "red", alpha =0.75, label = "Outside")
    plt.title(j)
    plt.xlabel(j)
    plt.legend(loc = "best")
    
    plt.savefig("{}.png".format(j))
        '''
        
    #finds values for each histogram
    values, bin_edge = np.histogram(Ranking_holder_out[j], bins = 75, range = (lower, upper))
    valuesin, bin_edgein = np.histogram(Ranking_holder[j], bins = 75, range = (lower, upper))
    
    #setps up the initial count of GRB's inside vs outside
    num = 0 
    numin = 0
    for i in range(len(values)):  
        #adds up the number of GRBs in this bin which are inside and outside of 200Mpc
        numin += valuesin[i]
        num += values[i]
        
        #how many outside are left
        left = 1000 - num
        
        #finding the relative difference to determine if this is significant
        diff = (numin - num)
        
        if abs(diff) >= distinct:
            #sets off when the difference in nunmber exceeds our distinction parameter
            
            #finds the value of the centre of the bin
            AVG = (bin_edge[i] + bin_edge[i+1])/2
            if diff*-1 > 0:
                #prints the information and writes to the txt file then breaks this loop 
               # print("Less than {} in 500 GRB's outside 200MPc produce ".format(left) + j + " values greater than: {:.3g}\n".format(AVG))
                #print("For this " + j + " value, {} out of 500 GRBs inside 200Mpc will be missed\n\n".format(numin))
                
                histtxt.write("lOOK INTO " + j + " THIS COULD BE USEFUL\n")
                histtxt.write("Less than {} in 500 GRB's outside 200MPc produce ".format(left) + j + " values greater than: {:.3g}\n".format(AVG))
                histtxt.write("For this " + j + " value, {} out of 500 GRBs inside 200Mpc will be missed\n".format(numin))
                histtxt.write("##############################################################################\n\n")
            
            elif diff*-1 < 0:
                #prints the information and writes to the txt file then breaks this loop 
                #print("Less than {} in 500 GRB's outside 200MPc produce ".format(num) + j + " values less than: {:.3g}\n".format(AVG))
                #print("For this " + j + " value, {} out of 500 GRBs inside 200Mpc will be missed\n\n".format(sum(valuesin) - numin))
                
                histtxt.write("lOOK INTO " + j + " THIS COULD BE USEFUL\n")
                histtxt.write("Less than {} in 500 GRB's outside 200MPc produce ".format(num) + j + " values less than: {:.3g}\n".format(AVG))
                histtxt.write("For this " + j + " value, {} out of 500 GRBs inside 200Mpc will be missed\n".format(sum(valuesin) - numin))
                histtxt.write("##############################################################################\n\n")
            
            break
    
        elif ((numin == N) & (num != N)):
            #print("less than {} in 500 GRBs outside produce".format(num) + j + " values less than: {:.3g}\n".format(max(Ranking_holder[j])))
            #print("For this "  + j + " all 500 inside GRBs are contained\n") 
            
            histtxt.write("less than {} in 500 GRBs outside produce ".format(num) + j + " values less than: {:.3g}\n\n".format(max(Ranking_holder[j])))
            histtxt.write("For this "  + j + " all 500 inside GRBs are contained\n\n")
            break
        
    #what to print if the difference between the number inside and out never exceeds our distinction parameter
    if num == N:
       # print("For Hist " + j + " there never reaches a good enough distinction of {} between the number of GRBs inside vs outside\n\n".format(distinct))
        histtxt.write("For Hist " + j + " there never reaches a good enough distinction of {} between the number of GRBs inside vs outside\n\n".format(distinct))     

histtxt.close()

"""
#setting an acceptance level for how well a stat performs compared to the others 
acceptance = 0.75

#making an array to see which stat is best
best = np.array([])

#making a txt file to store the information
perform = open("Relative placement of stats for {:.3g} % of the time.txt".format(acceptance*100), "w")
for p in range(29):
    '''
    This bit is the same as the code above (inside) it is basically finding the
    number of times each stat ranked a certain place out of the number of stats
    and then adds one to that location, this is what produces the black histograms
    '''
    df_rank = pd.DataFrame(rankN[:, p])
    
    val = df_rank[0].value_counts()
    vals = df_rank[0].value_counts().values.tolist()
    quantities = np.zeros(29)
    idq = val.index.values.tolist()

    for j in range(len(vals)):
        quantities[int(idq[j])] = vals[j]
    
    #this is my added bit which makes a count and then counts the number of times the stat scores a certain place
    count = 0
    for k in range(29):
        count += quantities[k]
        
        #once the count exceeds the desired percentage of the total number of events this writes the information to the txt and breaks the loop
        if count >= acceptance*sum(quantities):         
            print("{0:.3g} % of the time statistic {1:.3g} was at least {2:.3g} best compared to the other stats\n".format(acceptance*100, p+1, k+1))
            perform.write("{0:.3g} % of the time statistic {1:.3g} was at least {2:.3g} best compared to the other stats\n\n".format(acceptance*100, p+1, k+1))
            break
    
    #keeping track of the counts for each stat
    best = np.append(best, k+1)

stats_best = np.asarray(np.where(best == min(best))).flatten() + 1

if len(stats_best) == 1:
    perform.write("The statistic which performs best {0:.3g} % of the time is stat {1:.3g}".format(acceptance*100, stats_best[0]))
    print("The statistic which performs best {0:.3g} % of the time is stat {1:.3g}".format(acceptance*100, stats_best[0]))

else:
    perform.write("The statistics which perform best {0:.3g} % of the time are stats".format(acceptance*100) + str(stats_best))
    print("The statistics which perform best {0:.3g} % of the time are stats".format(acceptance*100) + str(stats_best))

#closing the txt
perform.close()
"""