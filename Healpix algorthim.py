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
df_master = pd.read_csv("Data Files/GLADE_Master.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

err_rad = 5

#RA around the celestial sphere varies from 0 -> 360, made an array in steps of 10 deg
RA_full =  np.arange(0, 370, err_rad) #degrees

#Dec varies from 90 -> -90, made an array in steps of 20 deg
Dec_full = np.arange(-90, 100, err_rad) #degrees

#this count is just to check how the loop is running, it will print the number of cycles each cycle
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
        
        #creating the name of the sector
        name = str('Sector_{0},{1}'.format(RA_full[i], Dec_full[j]))
        
        #Saving all the Galaxies into a csv file which is labelled after the lower RA and Dec
        Dec_hold.to_csv(r'Data Files\GLADE_Sectioned\{}.csv'.format(name), header = True)
    
    count += 1
    print(count)
    
    