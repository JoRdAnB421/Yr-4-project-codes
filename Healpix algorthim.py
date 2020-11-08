# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 16:05:44 2020
@author: Jordan
Using a type ofHealpix alogrithm to separate the sky into equal areas to section 
the galaxies off into individual sectors
"""

import numpy as np
import pandas as pd


#importing the master GLADE file
#df_master = pd.read_csv("Data Files/GLADE_Master.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

def rankw(theta, sigma, d_lum, luminosity, luminosity_probability): ## Normal
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

    