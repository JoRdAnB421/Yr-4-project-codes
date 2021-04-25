# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 12:59:00 2020
@author: Jordan
testing astropy.healpix
"""

import numpy as np; import pandas as pd; from astropy_healpix import HEALPix; 
from astropy.coordinates import ICRS, SkyCoord; from astropy import units as u;

#importing the master GLADE file
df_master = pd.read_csv("Data Files/GLADE_Master_comma_100Mpc.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

#grabbing just the coordinates of each galaxy
ra, dec = df_master["RA"], df_master["dec"]

#using HEALPix to split the celestial sky into equal area sectors
hp = HEALPix(nside=16, order='ring', frame=ICRS())

#making a holding array to hold the sectors for each galaxy
hold = np.array([])    


for i in range(len(ra)):
    '''
    This loops over all the galaxies, at each one it takes the galaxies coordinates
    and uses skycoord_to_healpix to find what HEALPix sector that galaxies belongs to,
    the result is then appended to the hold array
    '''
    coords = SkyCoord(ra[i], dec[i], unit = "deg")
    sector = hp.skycoord_to_healpix(coords)
    hold = np.append(hold, sector)


#adding the sector identifyier data to the master file
df_master["Sector"] = hold

for j in range(hp.npix):
    '''
    This loops over the sectors defined by HEALPix, grabbing the corresponding 
    galaxies from the main file and creating a new csv file for that sector
    '''
    #isolating the galaxies in the sector 
    sec = df_master.iloc[np.where(df_master["Sector"] == j)]
    
    #making a name for the current sector
    name = str("Sector_{}".format(j))
    
    #turning dataframe to csv file
    sec.to_csv(r'Data Files\GLADE_Sectioned_Reduced\{}.csv'.format(name), header = True)

def Sector_find(RA_grb, Dec_grb, err_radius):
    '''
    Give coordinates of the grb location and an error in the position, this function
    will use cone_search to find all sky sectors that the cone intersects and 
    will read the corresponding csv files and compile them into one dataframe
    '''
    
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