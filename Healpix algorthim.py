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

err_rad = 10

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
"""

##############################################################################
'''
Trying to define a function which will check the location of the GRB and
place it into one of the sectors, from there it will check if it is close to 
the edges of the sector and if so grab the galaxies from that sector also
'''

def Sect_sel(RA_grb, Dec_grb, err = 5):
    '''
    Take the simulated GRB position and find what sector it belongs to and where 
    that position lies in relation to the centre of the sector. As the size
    of the sectors is dictated by the error radius of the GRB position, only 
    when the grb is directly in the centre would we only use that one sector
    '''
    
    #RA around the celestial sphere varies from 0 -> 360, made an array in steps of 10 deg
    ra =  np.arange(0, 370, err) #degrees

    #Dec varies from 90 -> -90, made an array in steps of 20 deg
    dec = np.arange(-90, 100, err) #degrees
    
    for i in range(len(ra)-1):
        if ra[i] <= RA_grb <= ra[i+1]:    
            rahold = ra[i]
            for j in range(len(dec)-1):
                if dec[j] <= Dec_grb <= dec[j+1]:
                    dechold = dec[j]
                    break
                else:
                    continue
            break
        else:
            continue
    
    
    Sec_contain = str('Sector_{0},{1}'.format(rahold, dechold))
    
    #grabbing the corresponding sector containing the grb
    df_contain = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(Sec_contain)\
                             , delimiter = ",")     
            
    #finding which side of the ra & dec centre lines the grb lies
    
    #sector to the left
    if RA_grb < ((ra[i]+ra[i+1])/2):
        lindex = ra[i-1]
        left_sec = str('Sector_{0},{1}'.format(lindex, dechold))
        
        #calling the corresponding sector csv file and appending it to other neccesary dataframes
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(left_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
    
    #sector to the right
    elif RA_grb > ((ra[i]+ra[i+1])/2):
        rindex = ra[i+1]
        right_sec = str('Sector_{0},{1}'.format(rindex, dechold))
          
        #calling the corresponding sector csv file and appending it to other neccesary dataframes      
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(right_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
        
    #sector below
    elif Dec_grb < ((dec[j]+dec[j+1])/2):
        dindex = dec[j-1]
        down_sec = str('Sector_{0},{1}'.format(rahold, dindex))
                
        #calling the corresponding sector csv file and appending it to other neccesary dataframes
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(down_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
            
    #sector above
    elif RA_grb > ((ra[i]+ra[i+1])/2):
        uindex = dec[j+1]
        right_sec = str('Sector_{0},{1}'.format(rahold, uindex))
                
        #calling the corresponding sector csv file and appending it to other neccesary dataframes
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(right_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
            
    #sector left and below
    elif RA_grb < ((ra[i]+ra[i+1])/2) and Dec_grb < ((dec[j]+dec[j+1])/2):
        dindex = dec[j-1]
        lindex = ra[i-1]
        leftdown_sec = str('Sector_{0},{1}'.format(lindex, dindex))      
                
        #calling the corresponding sector csv file and appending it to other neccesary dataframes
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(leftdown_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
            
    #sector left and above
    elif RA_grb < ((ra[i]+ra[i+1])/2) and Dec_grb > ((dec[j]+dec[j+1])/2):
        uindex = dec[j+1]
        lindex = ra[i-1]
        leftup_sec = str('Sector_{0},{1}'.format(lindex, uindex))
                
        #calling the corresponding sector csv file and appending it to other neccesary dataframes
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(leftup_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
            
    #sector right and below
    elif RA_grb > ((ra[i]+ra[i+1])/2) and Dec_grb < ((dec[j]+dec[j+1])/2):
        dindex = dec[j-1]
        rindex = ra[i+1]
        rightdown_sec = str('Sector_{0},{1}'.format(rindex, dindex))
                
        #calling the corresponding sector csv file and appending it to other neccesary dataframes
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(rightdown_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
            
    #sector right and above
    elif RA_grb > ((ra[i]+ra[i+1])/2) and Dec_grb > ((dec[j]+dec[j+1])/2):
        uindex = dec[j+1]
        rindex = ra[i+1]
        rightup_sec = str('Sector_{0},{1}'.format(rindex, uindex))
                
        #calling the corresponding sector csv file and appending it to other neccesary dataframes
        holder = pd.read_csv("Data Files/GLADE_Sectioned/{}.csv".format(rightup_sec)\
                             , delimiter = ",")
        
        df_contain = df_contain.append(holder)
        
    return df_contain
"""
    