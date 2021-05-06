<<<<<<< HEAD
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

=======
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 16:05:44 2020
@author: Jordan
Using Healpix alogrithm to separate the sky into equal areas to section the 
galaxies off into individual sectors
"""

import numpy as np
import pandas as pd

#importing the master GLADE file
df_master = pd.read_csv("Data Files/GLADE_Master.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

#RA around the celestial sphere varies from 0 -> 360, made an array in steps of 10 deg
RA_full =  np.arange(0, 370, 20) #degrees

#Dec varies from 90 -> -90, made an array in steps of 20 deg
Dec_full = np.arange(-90, 100, 20) #degrees


count = 0

for i in range(len(RA_full)-1):
    '''
    This for loop will loop over the RA sections of 20 degrees, at each section
    a holding array is made to select the galaxies in the main dataframe where the
    RA is above the lower value. A second holding array is then made which selects
    from those the ones below the upper limit. This holding array is then fed into
    the next loop to check the declination under the same method.
    '''
    RA_hold = df_master.iloc[np.where(RA_full[i] <= df_master["RA"])]
    RA_hold = RA_hold.iloc[np.where(RA_hold["RA"] <= RA_full[i+1])]
    
    for j in range(len(Dec_full)-1):
        '''
        Declination sectors are found in the same way as the RA before and then
        the Galaxies that have survived the tests will go into their own data frame
        for that sector
        '''
        
        Dec_hold = RA_hold.iloc[np.where(Dec_full[j] <= RA_hold["dec"])]
        Dec_hold = Dec_hold.iloc[np.where(Dec_hold["dec"] <= Dec_full[j+1])]
        
        name = str('Sector_{0},{1}'.format(RA_full[i], Dec_full[j]))
        
        Dec_hold.to_csv(r'Data Files\GLADE_Sectioned\{}.csv'.format(name), header = True)
    
    count += 1
    print(count)
    
>>>>>>> 1c9aa4c7d0c5291d1fcd3a218fcf5e23efec9c37
    