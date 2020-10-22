
import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d

N = 50       ##Change to alter the number of loops the code runs for

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


no_se_func = []
ras_dex = np.zeros(shape = (N, 2))
test_case = np.zeros(shape = (N, 2))


def Ang_Dist(ra1, ra2, dec1, dec2):## Calculates the angular distance between apparent position and galaxy
    
    ra1 *= (np.pi/180); ra2 *= (np.pi/180)
    dec1 *= (np.pi/180); dec2 *= (np.pi/180)
    
    return (180/np.pi) * np.arccos(np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
#################################################################
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
#########################################################################################
#########################################################################################
df_master = pd.read_csv("Data Files/GLADE_Master.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined

L1 = np.linspace(56, 59, 101) #In J now
L2, c = L_func1(L1) #  ##Builds broken power law

cumuL = cumulative(L2) ##Luminosity Distribution
df_cumLum =   pd.read_csv("Data Files/Cumulative Luminosity.csv")
df_cumLum.columns = ["NaN", "Cumulative Luminosity"]
normal_c = df_cumLum[["Cumulative Luminosity"]].values[-1][0]
L_rank = df_cumLum[["Cumulative Luminosity"]].values * 1/normal_c
df_cumLum = df_cumLum[["Cumulative Luminosity"]].values#                ## This is all to do with building a usable and callable power law
                                                                         

lum_N = np.linspace(0, df_cumLum.shape[0], df_cumLum.shape[0])
df_dL = df_master[["Luminosity Distance"]]


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
"""
placement18 = np.zeros(N)
percentages18 = np.zeros(N)
for i in range(N):
    
    current_i = indices.index(gals[i])
    
    testr = np.array(test_ra.iloc[[current_i]].values.tolist())
    testd = np.array(test_dec.iloc[[current_i]].values.tolist())

    ident = np.zeros(df_master.shape[0])
      
    print(str(i + 1), "out of " + str(N))
    print("Test galaxy: ", str(gals[i]))
    
    ident[current_i] = 1
    df_master["Identifier"] = ident  ## Creates a mask for identifying the host galaxy
    
    
    q, t, df_sliced = reduction(abs(ra_prime[i]), dec_prime[i], df_master) ## Reduces the catalog by RA and dec
    ra = np.array(df_sliced[["RA"]].values.tolist())[:, 0]
    dec = np.array(df_sliced[["dec"]].values.tolist())[:, 0]
    
    Luminosity = np.array(df_sliced[["B Luminosity"]].values.tolist()) #Luminosity_Handling(np.array(df_sliced[["Absolute B Magnitude"]].values.tolist())) ## Converts A    
    dl = np.array(df_sliced[["Luminosity Distance"]].values.tolist())    
    
    lum_prob, SGR_test = L_func(received_luminosity[i], c, dl) ##Uses the luminosity function to calculate probabilities
    df_sliced["Luminosity Probability"] = lum_prob
    df_sliced["SGR flag"] = SGR_test
    
    angular_distaance = np.zeros(df_sliced.shape[0])
    
    for k in range(df_sliced.shape[0]):
        angular_distaance[k] = Ang_Dist(ra[k], testr[0][0], dec[k], testd[0][0])
    
    df_sliced["Angular Distance"] = angular_distaance

    ranking23 = rank23(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank23"] = ranking23
    
    df_sliced18 = (pd.DataFrame.sort_values(df_sliced, by = ["Rank23"], ascending = False))

    id18 = df_sliced18[["Identifier"]].values.tolist()
    mask_check18 = [i for i, val in enumerate(id18) if val == [1]]
    
    Luminosity =  np.asarray(Luminosity)
    
    if len(mask_check18) == 0:
        print("Did not place\n\n\n")
        next
    else:
        length = len(id18) + 1
        
        placement18[i] = mask_check18[0] + 1

        #display(Markdown("The keplerian orbit appears to be happening at r ={0:.2f} km"  .format(float(kepler(M_kep, w))/1000)))
        #print("Galaxy data: \nDistance is {0:.2f} Mpc\nLuminosity is {1:.3e}\nra and dec [{2:.2f}, {3:.2f}] compared to reported ra and dec [{4:.2f}, {5:.2f}] \nTrue luminosity {6:.3e} W" .format(dl[int(placement18[i] - 1)][0], Luminosity[int(placement18[i] - 1)][0], fin_ra[int(placement18[i] - 1)][0], fin_dec[int(placement18[i] - 1)][0], testr[0][0], testd[0][0], b[i]))
        
        print("Galaxy placed", int(placement18[i]), "out of", str(length), "with statistic 18\n\n\n")
        
        percentages18[i] = placement18[i]/length 


"""
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

    ident = np.zeros(df_master.shape[0])
      
    print(str(i + 1), "out of " + str(N))
    print("Test galaxy: ", str(gals[i]))
    
    ident[current_i] = 1
    df_master["Identifier"] = ident  ## Creates a mask for identifying the host galaxy
    
    
    q, t, df_sliced = reduction(abs(ra_prime[i]), dec_prime[i], df_master) ## Reduces the catalog by RA and dec
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
    
    ranking2 = rank2(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank2"] = ranking2
    
    ranking3 = rank3(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank3"] = ranking3
    
    ranking4 = rank4(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank4"] = ranking4
    
    ranking5 = rank5(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank5"] = ranking5
    
    ranking6 = rank6(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank6"] = ranking6
    
    ranking7 = rank7(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank7"] = ranking7
    
    ranking8 = rank8(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank8"] = ranking8

    ranking9 = rank9(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank9"] = ranking9
    
    ranking10 = rank10(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank10"] = ranking10
    
    ranking11 = rank11(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank11"] = ranking11
    
    ranking12 = rank12(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank12"] = ranking12
    
    ranking13 = rank13(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank13"] = ranking13
    
    ranking14 = rank14(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank14"] = ranking14

    ranking15 = rank15(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank15"] = ranking15
    
    ranking16 = rank16(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank16"] = ranking16
    
    ranking17 = rank17(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank17"] = ranking17
    
    ranking18 = rank18(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank18"] = ranking18
    
    ranking19 = rank19(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank19"] = ranking19
    
    ranking20 = rank20(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank20"] = ranking20
    
    ranking21 = rank21(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank21"] = ranking21

    ranking22 = rank22(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank22"] = ranking22
    
    ranking23 = rank23(angular_distaance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank23"] = ranking23
    
    
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
    """
    for k in range(5):        
        aa[i][k] = np.exp(-(df_sliced[["Angular Distance"]].head(5).values.tolist()[k][0])/error_radius)
        ab[i][k] = df_sliced[["Luminosity Distance"]].head(5).values.tolist()[k][0]
        ac[i][k] = df_sliced[["B Luminosity"]].head(5).values.tolist()[k][0]
        ad[i][k] = df_sliced[["Luminosity Probability"]].head(5).values.tolist()[k][0]
    """     
"""
plt.figure(0)
plt.plot(percentages19, np.log10(distances), "kx")
#plt.title("Distance vs. percentage performance")
plt.ylabel("Log$_{10}$ Distance /Mpc"); plt.xlabel("Percentage placement"); plt.grid()
#plt.xlim(1e-27, 1)
plt.savefig("Distances vs. percentage.png")

plt.figure(1)
plt.plot(percentages19, np.log10(b), "kx")
#plt.title("Intrinsic Luminosity vs. percentage performance")
plt.ylabel("Log$_{10}$ Luminosity /W"); plt.xlabel("Percentage placement"); plt.grid()
#plt.xlim(1e-27, 1)
plt.savefig("Luminosity vs. percentage.png")

plt.figure(2)
plt.plot(percentages19, rotation_angle, "kx")
plt.ylabel("Angular offset /$^o$"); plt.xlabel("Percentage performance")
plt.grid()
plt.savefig("Angular offset vs. percentage.png")
### The following can be used to investigate any values that flag up as false
"""
f_v = [i for i, val in enumerate(faulty[:, 4]) if val == 0]
f_1v = [i for i, val in enumerate(faulty[:, 4]) if val == 1]

sets = set(np.arange(0, len(faulty), 1)) - set(f_v)
ft = pd.DataFrame(faulty)
faulty_cols = ["Galaxy RA", "GRB RA", "Galaxy dec", "GRB dec", "Mask"]
ft.columns = faulty_cols
"""
ab_fault = ft.take(list(sets), axis = 0)
ab_vals = ab_fault.values.tolist()[0]
"""
place_array = np.zeros(shape = (N, 23))
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
"""
plt.figure(3)
for p in range(20):
    plt.plot(df_place_array[:, p], np.log10(df_distance), "x", alpha = 2/(p/2 + 1), label = "Statistic" + str(p))
    
plt.title("Distance vs. percentage performance")
plt.ylabel("Log$_{10}$ Distance /Mpc"); plt.xlabel("Percentage placement"); plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("Statistic_Comparison.png")
"""

rankN = np.zeros(shape = (len(df_place_array), 23))


for i in range(len(df_place_array)):
    
    df_array_init = pd.DataFrame(df_place_array[i, :]) ## Takes percentage placement for each run
    
    counting_mask = np.arange(df_array_init.shape[0])
    df_array_init["Mask"] = counting_mask    ## Creates a matching mask for keeping track of where the entries end up
    
    df_array = (pd.DataFrame.sort_values(df_array_init, by = [0], ascending = True)) ## Orders resultant sliced array
    
    for k in range(df_array.shape[0]):  
        rankN[i, k] = [i for i, val in enumerate(df_array[["Mask"]].values.tolist()) if val == [k]][0] ## 
counter = 5
for p in range(23):
    df_rank = pd.DataFrame(rankN[:, p])
    
    plt.figure(p + 4)
    
    val = df_rank[0].value_counts()
    vals = df_rank[0].value_counts().values.tolist()
    quantities = np.zeros(23)
    idq = val.index.values.tolist()

    for j in range(len(vals)):
        quantities[int(idq[j])] = vals[j]

    for o in range(23):
        plt.bar((o + 1), quantities[o], color = "black")
    
    plt.xlabel("Placement"); plt.ylabel("Frequency")
    plt.title("Statistic " + str(p + 1))
    plt.grid()
    plt.savefig("Statistic " + str(p + 1) + ".png")
    counter += 1
    
for i in range(23):
    plt.figure(counter)
    plt.plot(np.log10(df_distance), df_place_array[:, i], "kx", label = "Statistic " + str(i + 1))
    plt.ylabel("Percentage performance")
    plt.xlabel("Log$_{10}$ Distance /Mpc")
    plt.grid()
    plt.legend(loc = "best")
    plt.savefig("OmittedGalaxies_Statistic" + str(i + 1) + ".png")
    counter += 1
    
for j in range(23):
    plt.figure(counter)
    plt.plot(np.log10(df_lumin), df_place_array[:, j], "kx", label = "Statistic " + str(j + 1))
    plt.ylabel("Percentage performance")
    plt.xlabel("Log$_{10}$ Luminosity /W")
    plt.grid()
    plt.legend(loc = "best")
    plt.savefig("OmittedGalaxies_Lumin_Statistic" + str(j + 1) + ".png")
    counter += 1
    
for k in range(23):
    plt.figure(counter)
    plt.plot((df_ang), df_place_array[:, k], "kx", label = "Statistic " + str(k + 1))
    plt.ylabel("Percentage performance")
    plt.xlabel("Angular Offset /$^o$")
    plt.grid()
    plt.legend(loc = "best")
    plt.savefig("OmittedGalaxies_Ang_Statistic" + str(k + 1) + ".png")
    counter += 1

"""
plt.plot(L1, L2, "k", label = "Broken Power Law")
plt.xlabel("Log$_{10}$ Luminosity /W")
plt.ylabel("Luminosity Probability, $\phi$(L)")
plt.grid()

#df_distance = np.log10(df_distance)
plt.figure(0)
plt.plot(np.log10(df_distance), df_place_array[:, 17], "kx", label = "Statistic 18")
plt.ylabel("Percentage performance")
plt.xlabel("Log$_{10}$ Distance /Mpc")
plt.grid()
plt.legend(loc = "best")
plt.savefig("OmittedGalaxies_Statistic18.png")

plt.figure(1)
plt.plot(df_ang, df_place_array[:, 17], "kx", label = "Statistic 18")
plt.ylabel("Percentage performance")
plt.xlabel("Angular offset /$^o$")
plt.grid()
plt.legend(loc = "best")
plt.savefig("OmittedGalaxies_ang_Statistic18.png")

plt.figure(2)
plt.plot(np.log10(df_lumin), df_place_array[:, 17], "kx", label = "Statistic 18")
plt.ylabel("Percentage performance")
plt.xlabel("Intrinsic Luminosity /W")
plt.grid()
plt.legend(loc = "best")
plt.savefig("OmittedGalaxies_lumin_Statistic18.png")

"""