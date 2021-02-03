# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:55:25 2021
@author: Jordan
Automating the rank production
"""

import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d; from astropy_healpix import HEALPix; 
from astropy.coordinates import ICRS, SkyCoord; from astropy import units as u;
from timeit import default_timer as timer

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

def Ang_Dist(ra1, ra2, dec1, dec2):## Calculates the angular distance between apparent position and galaxy
    
    ra1 *= (np.pi/180); ra2 *= (np.pi/180)
    dec1 *= (np.pi/180); dec2 *= (np.pi/180)
    
    return (180/np.pi) * np.arccos(np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
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
def Luminosity_back_convert(L_given, d_L): #  ##Converts luminosity to luminosity at source
    #L = L0/4 *np.pi * d_l**2
    return (L_given) * (4 * np.pi * (3.086e22 * d_L)**2)
#################################################################
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

def rank(theta, sigma, d_lum, luminosity, luminosity_probability, a, b, c, d): ## Normal
    ## Implements a ranking statistic with arbitrary powers
    return np.exp(-((a*theta**2)/(2 * (sigma)**2))) *(1/d_lum**b * luminosity**c)[:, 0] * luminosity_probability**d

#%%
start = timer()

df_master = pd.read_csv("Data Files/GLADE_Master_comma_100Mpc.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

#setting the luminosity function of the GRBs and the cumulative distribution
L1 = np.linspace(56, 59, 101) #In J now
L2, c = L_func1(L1) #  ##Builds broken power law
cumuL = cumulative(L2) ##Luminosity Distribution

