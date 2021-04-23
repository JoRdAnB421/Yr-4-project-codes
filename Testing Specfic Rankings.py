# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 17:40:56 2021
@author: Jordan
Running specific Ranking statistics for various reasons
"""

import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d; from astropy_healpix import HEALPix; 
from astropy.coordinates import ICRS, SkyCoord; from astropy import units as u;
from timeit import default_timer as timer

#number of GRBs to be tested
N = 1000

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

def rank(theta, sigma, d_lum, luminosity, luminosity_probability, n1, n2, n3, n4): ## Normal
    ## Implements a ranking statistic with arbitrary powers
    return np.exp(-((n1*theta**2)/(2 * (sigma)**2))) *(1/d_lum**n2 * luminosity**n3)[:, 0] * luminosity_probability**n4

#%%

start = timer()

df_master = pd.read_csv("Data Files/GLADE_Master_comma_100Mpc.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

#setting the luminosity function of the GRBs and the cumulative distribution
L1 = np.linspace(56, 59, 101) #In J now
L2, c = L_func1(L1) #  ##Builds broken power law
cumuL = cumulative(L2) ##Luminosity Distribution

#This is grabbing the galaxy luminosity information and formatting it for use
cumuL = cumulative(L2) ##Luminosity Distribution
df_cumLum =   pd.read_csv("Data Files/Cumulative_Luminosity_100Mpc.csv")
df_cumLum.columns = ["NaN", "Cumulative Luminosity"]
normal_c = df_cumLum[["Cumulative Luminosity"]].values[-1][0]
L_rank = df_cumLum[["Cumulative Luminosity"]].values * 1/normal_c
df_cumLum = df_cumLum[["Cumulative Luminosity"]].values#    

#grabbing the luminsotity distance of the galaxies
df_dL = df_master[["Luminosity Distance"]]

#using HEALPix to split the sky into equal area sectors
hp = HEALPix(nside=16, order='ring', frame=ICRS())

#setting up to choose random GRB values and associate to a galaxy within the range
dummies = random(N)
RandL = random(N)

#setting arrays for the host galaxies and the Luminosity recieved
gals = np.zeros(N)          ## Picks out a luminosity
received_luminosity = np.zeros(N)

#creating empty arrays to store information
df_dL = df_dL.values.tolist()   ## Luminosity distance values for use
a = np.zeros(N)                 ## For storing temporary and unimportant values
b = np.zeros(N)                 ## For storing temporary and unimportant values

#taking the sky positions of all the galaxies
test_ra = df_master[["RA"]]
test_dec = df_master[["dec"]]

#used to create an index of the galaxies (for keeping track of them)
indices = list(np.arange(df_master.shape[0]))

error_radius = 2 * (2.62)  ## Change as necessary - this is an example value from HEARSCH 

#sets arrays which will check for if the    
rank_host = np.zeros(N)
faulty = np.zeros(shape = (N, 5))       ## All of this used to store values

#setting random rotations
phi = 2 * np.pi * random(N) * (180/np.pi)               ## Random positions for rotations
theta = np.arccos(2 * random(N) - 1) * (180/np.pi)      

#converting to spherical coordinates
thph = spherical_convert(theta, phi)
mod = np.zeros(N)

for i in range(N):
    mod[i] = modulus(thph[:, i])
    thph[:, i] /= mod[i]

#setting empty cartesian arrays for the loop
xyz = np.zeros(shape = (N, 3))
m = np.zeros(shape = (N, 3))

#empty arrays for storing the GRB coordinate
ra_prime = np.zeros(N); dec_prime = np.zeros(N)
rotation_angle = error_radius * normal(size = N) * (np.pi/180) 

###################################################################################################
'''
This is producing N GRBs associated with galaxies inside the sphere
'''
count = 0
for i in range(N):
    '''
    This loop will simulate N GRBs originating from galaxies with the 200Mpc
    observable sphere and then return the coordinates of the GRBs along with the 
    identifier of the host galaxies.
    '''
    
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
####################################################################################################    
    ##Rotation of test ra and dec
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
    
    count += 1
    if count % (0.1*N) == 0:
        print(count)
###################################################################################################
'''
This next section is to produce N GRBs simulated to exist outside of the sphere
'''



a = np.zeros(N) ## For storing temporary and unimportant values
b = np.zeros(N) ## For storing temporary and unimportant values

#Creates from Numbers for dummies and Luminosity fraction between 0->1
dummies_out = random(N)
RandL_out = random(N)

#empty array for the recieved luminosity outside
received_luminosity_out = np.zeros(N)

#This is to produce a random point on the sky
angles = np.arccos(2 * random(N) - 1)

ra_rand = uniform(0, 360, size = N)
dec_rand = (np.pi/2 - angles) * (180/np.pi) ##Gives you random ra and dec

#makes a random distance from us, at least 250 Mpc away
r_min = 250 # Mpc
r = r_min / random(N)**(1/3)

#empty arrays for the outside GRB locations
ra_prime_out = np.zeros(N); dec_prime_out = np.zeros(N)

#setting random rotations of the position of the GRBs
rotation_angle_out = error_radius * normal(size = N)

count = 0
for i in range(N):
    '''
    This makes N random luminosities from the luminosity power law for grbs, it
    then matches that to a random distance from here to '''       
    a[i] = find_nearest(cumuL, RandL[i])
    if a[i] == len(L1):
        a[i] = len(L1) - 1
    b[i] = 10**L1[int(a[i])] #int(a[i])
    received_luminosity_out[i] = Luminosity_for_convert(b[i], r[i])
    
    #applying the randomised error
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
        
    x_prime = axis_rotation(m[i, :], xyz, rotation_angle_out[i])
    
    xmod = modulus(x_prime)
    x_prime /= xmod
    ra_prime_out[i], dec_prime_out[i] = back_convert(x_prime)
    
    count += 1
    if count % (0.1*N) == 0:
        print(count)

elapsed_time = timer() - start # in seconds
print('Simulating the inside and outside GRBs took {:.4g} minutes'.format(elapsed_time/60))

#%%

'''
This section will start with some intial conditions for the powers of the terms
and then vary them all a little bit to see if that produces a better Ranking 
Statistic.
'''

start = timer()

#sets the number of different combinations of powers to test
Num = 50

#sets what powers you want to use
powers = [0.0942490673932074, 0.342760123472853, -0.0258831489973406, -4.6668065925182]

#powers = [1, 1, 0, 1]

#I want to measure for various different styles of ranking and so this defines those types
names = np.array(["Max Rank", "Sum Ranks", "Top 5 Avg"], dtype = str)    
Ranking_holder = pd.DataFrame(columns = names, index = range(N))

#this is the names for the outside GRBs
Ranking_holder_out = pd.DataFrame(columns = names, index = range(N))

#setting up a dataframe to store the Best ranks
column_names = ["Rank Type", "Distinction", "Cutoff Value", "a", "b", "c", "d", "Mean Placement of Host"]
Iterative_rank = pd.DataFrame(index = range(1), columns = column_names)

#a holder for any faults
no_se_func = []

#keeps track of how well placed the host galaxies were
placement = np.zeros(N)
percentages = np.zeros(N)

count = 0
for i in range(N):
    
    '''
    The first part of this loop will work on the given ranking statistic
    for the 1000 GRBs inside the sphere.
    '''
    current_i = indices.index(gals[i])
    
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
    
    #####################################################################
    #taking the angular position of the galaxies
    ra = np.array(df_sliced[["RA"]].values.tolist())[:, 0]
    dec = np.array(df_sliced[["dec"]].values.tolist())[:, 0]
    
    #finding the luminosity of the galaxies 
    Luminosity = np.array(df_sliced[["B Luminosity"]].values.tolist()) #Luminosity_Handling(np.array(df_sliced[["Absolute B Magnitude"]].values.tolist())) ## Converts A    
    dl = np.array(df_sliced[["Luminosity Distance"]].values.tolist())    
    
    #finding the luminosity probability of this GRB 
    lum_prob, SGR_test = L_func(received_luminosity[i], c, dl) ##Uses the luminosity function to calculate probabilities
    df_sliced["Luminosity Probability"] = lum_prob
    df_sliced["SGR flag"] = SGR_test
    #####################################################################
    
    #some work on finding the angular distance of this GRB to the galaxies which are in the sectors
    angular_distaance = np.zeros(df_sliced.shape[0])
    
    for k in range(df_sliced.shape[0]):
        angular_distaance[k] = Ang_Dist(ra[k], ra_prime[i], dec[k], dec_prime[i])
  
    id_check = [i for i, val in enumerate(angular_distaance) if math.isnan(val) == True]
    
    for k in range(len(id_check)):
        angular_distaance[int(id_check[k])] = Ang_Dist(ra_prime[i], testr, dec_prime[i], testd)
    
    #keeping track of the angular distance values
    df_sliced["Angular Distance"] = angular_distaance
    
    #####################################################################
    #using the proposed ranking statistic and storing it 
    ranking = rank(angular_distaance, error_radius, dl, Luminosity, lum_prob, *powers)  ## Uses defined ranking statistic
    df_sliced["Rank"] = ranking

    #finding values for maximum rank and rank summation etc
    Ranking_holder.loc[i, "Max Rank"] = max(ranking)
    Ranking_holder.loc[i, "Sum Ranks"] = np.sum(ranking)
    
    x = -np.sort(-ranking)
    Ranking_holder.loc[i, "Top 5 Avg"] = np.mean(x[:5])
    ######################################################################
    
    fin_ra = np.asarray(df_sliced[["RA"]].values.tolist()); fin_dec = np.asarray(df_sliced[["dec"]].values.tolist())
    ## Storing values and extending the reduced catalog
    
    df_sliced = (pd.DataFrame.sort_values(df_sliced, by = ["Rank"], ascending = False)) ## Orders resultant sliced array
    
    idi = df_sliced[["Identifier"]].values.tolist() ##Mask handling to check for values
    mask_check = [i for i, val in enumerate(idi) if val == [1]]
    
    Luminosity =  np.asarray(Luminosity)

    if len(mask_check) == 0:
        print("Did not place\n\n\n")
        next
    else:
        length = len(idi) + 1
        
        placement[i] = mask_check[0] + 1; length = len(idi) + 1
    
        percentages[i] = placement[i]/length
    
    #####################################################################
    #checking for faults
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
        
    
    #####################################################################
    '''
    This is the second part of the loop which will run for the 1000 outside 
    GRBs
    '''
    #selects the corresponding sectors to look through
    df_sliced_out = Sector_find_reduced(ra_prime_out[i], dec_prime_out[i], error_radius)
    df_sliced_out = df_sliced_out.rename(columns = {"Unnamed: 0.1": "Unnamed: 0"})
    
    if df_sliced_out.empty == True:
        print("No Galaxies within the error radius of {} found assume it must have come from outside the sphere")
        Ranking_holder_out.loc[i] = [0, 0, 0]
        count += 1
        continue
 
    #selecting the galaxy RA and Dec
    ra_out = np.array(df_sliced_out[["RA"]].values.tolist())[:, 0]
    dec_out = np.array(df_sliced_out[["dec"]].values.tolist())[:, 0]
    
    #collecting the luminosities and lum probabilities
    Luminosity_out = np.array(df_sliced_out[["B Luminosity"]].values.tolist()) #Luminosity_Handling(np.array(df_sliced[["Absolute B Magnitude"]].values.tolist())) ## Converts A    
    dl_out = np.array(df_sliced_out[["Luminosity Distance"]].values.tolist())    
    
    lum_prob_out, SGR_test_out = L_func(received_luminosity_out[i], c, dl_out) ##Uses the luminosity function to calculate probabilities
    df_sliced_out["Luminosity Probability"] = lum_prob_out
    df_sliced_out["SGR flag"] = SGR_test_out
    
    #sorting out the angular distances 
    angular_distaance_out = np.zeros(df_sliced_out.shape[0])
    
    
    for k in range(df_sliced_out.shape[0]):
        angular_distaance_out[k] = Ang_Dist(ra_out[k], ra_prime_out[i], dec_out[k], dec_prime_out[i])
    
    
    df_sliced_out["Angular Distance"] = angular_distaance_out

    ####################################################################
    #defining the ranking statistic and finding the various values
    ranking_out = rank(angular_distaance_out, error_radius, dl_out, Luminosity_out, lum_prob_out, *powers)  ## Uses defined ranking statistic
    df_sliced_out["Rank"] = ranking_out

    Ranking_holder_out.loc[i, "Max Rank"] = max(ranking_out)
    Ranking_holder_out.loc[i, "Sum Ranks"] = np.sum(ranking_out)
    
    x = -np.sort(-ranking_out)
    Ranking_holder_out.loc[i, "Top 5 Avg"] = np.mean(x[:5])
    #####################################################################

    df_sliced_out = (pd.DataFrame.sort_values(df_sliced_out, by = ["Rank"], ascending = False)) ## Orders resultant sliced array
    
    count += 1
    if count % (0.1*N) == 0:
        print(count)
    
print("The Ranking Statistic has been completed on all inside and outside GRBs.")




'''
This last part of the loop is going to check if the distinction level between
the inside GRBs and outside GRBs is bettered for these new proposed ranking
when considering the three ranking values Max Rank, Sum Ranks and top 5 avg Rank values.

If the distinction level is bettered outside to inside GRBs is reached then, that ranking 
stat will be saved and the powers on its components will be become the new intial 
conditions for the next change
'''
#%%
#distinction set as 0 to start with to see how well the base statistic works
distinction_level = 0

for k in names:
    #finding the range (upper)
    if max(Ranking_holder[k]) >= max(Ranking_holder_out[k]):
        upper = max(Ranking_holder[k])
    elif max(Ranking_holder[k]) < max(Ranking_holder_out[k]):
        upper = max(Ranking_holder_out[k])
    
    #finding the range (lower)
    if min(Ranking_holder[k]) <= min(Ranking_holder_out[k]):
        lower = min(Ranking_holder[k])
    elif min(Ranking_holder[k]) > min(Ranking_holder_out[k]):
        lower = min(Ranking_holder_out[k])
    #finds values for each histogram
    values, bin_edge = np.histogram(Ranking_holder_out[k], bins = 75, range = (lower, upper))
    valuesin, bin_edgein = np.histogram(Ranking_holder[k], bins = 75, range = (lower, upper))
    #sets up the initial count of GRB's inside vs outside
    num = 0 
    numin = 0
    
    
    for t in range(len(values)):  
        #adds up the number of GRBs in this bin which are inside and outside of 200Mpc
        numin += valuesin[t]
        num += values[t]
        
        #how many outside are left
        left = 1000 - num
        
        #finding the relative difference to determine if this is significant
        diff = (numin - num)
        
        if abs(diff)/N > distinction_level:
            #sets off when the difference in number exceeds our distinction parameter
            #i.e., it finds a better distinction value
            x = numin
            y = num
            
            #storing the values and bin edges to plot later
            goodvals = values
            goodvalsin = valuesin
            goodedge = bin_edge
            goodedgein = bin_edgein
            
            
            New_distinct = diff/N
            distinction_level = abs(New_distinct)
            
            #finds the value of the centre of the bin
            AVG = (bin_edge[t] + bin_edge[t+1])/2
            
            #sets the current best Ranking type
            best_type = k

mean_host_placement = np.mean(percentages)

initial = [best_type, distinction_level, AVG, *powers, mean_host_placement]

Iterative_rank.loc[0] = initial

elapsed_time = timer() - start # in seconds
print('The code took {:.3g} minutes to complete'.format(elapsed_time/60))
            
#%%
plt.close()

plt.figure()
plt.bar(goodedge[:-1], goodvals, width= (goodedge[1]-goodedge[0]), color='red', alpha=0.5, label = "outside")
plt.bar(goodedgein[:-1], goodvalsin, width= (goodedgein[1]-goodedgein[0]), color='blue', alpha=0.5, label = "inside")
plt.axvline(AVG, ls = "--", color = "black", label = "Best distinction point")


plt.legend(loc = "best")
plt.xlabel("Stat values")
plt.ylabel("Number")
plt.title("Histogram of {0} producing the best distinction \n value for {1} simulated SGRBs".format(best_type, N))
plt.savefig("Histogram.png")