
import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d

N = 5      ##Change to alter the number of loops the code runs for

placement = np.zeros(shape = (N,5))
placement2 = np.zeros(shape = (N,5))
placement3 = np.zeros(shape = (N,5))
placement4 = np.zeros(shape = (N,5))
placement5 = np.zeros(shape = (N,5))
placement6 = np.zeros(shape = (N,5))
placement7 = np.zeros(shape = (N,5))
placement8 = np.zeros(shape = (N,5))
placement9 = np.zeros(shape = (N,5))
placement10 = np.zeros(shape = (N,5))
placement11 = np.zeros(shape = (N,5))
placement12 = np.zeros(shape = (N,5))
placement13 = np.zeros(shape = (N,5))
placement14 = np.zeros(shape = (N,5))
placement15 = np.zeros(shape = (N,5))
placement16 = np.zeros(shape = (N,5))
placement17 = np.zeros(shape = (N,5))
placement18 = np.zeros(shape = (N,5))
placement19 = np.zeros(shape = (N,5))
placement20 = np.zeros(shape = (N,5))
placement21 = np.zeros(shape = (N,5))
placement22 = np.zeros(shape = (N,5))
placement23 = np.zeros(shape = (N,5))


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
    return np.exp(-(theta**2/((sigma)**2))) * (1/d_lum**0 * luminosity)[:, 0] * luminosity_probability

def rank19(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Luminosity
    
    return np.exp(-(theta**2/(2 * (sigma)**2))) * (1/d_lum * luminosity**0)[:, 0] * luminosity_probability

def rank20(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Lum_Prob
    return np.exp(-(theta**2/(2 * (sigma)))) * (1/d_lum * luminosity)[:, 0] * luminosity_probability**0

def rank21(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Lum_Prob
    return np.exp(-(theta**2/(2 * (sigma)**2)))**(4) * (1/d_lum**8 * luminosity)[:, 0] * luminosity_probability**2

def rank22(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Lum_Prob
    return np.exp(-(theta**2/(2 * (sigma))))**(sigma**4) * (1/d_lum**8 * luminosity)[:, 0] * luminosity_probability**2

def rank23(theta, sigma, d_lum, luminosity, luminosity_probability): ## No Lum_Prob
    return np.exp(-((theta**2)**100/(2 * (sigma)))) * (1/d_lum**8 * luminosity)[:, 0] * luminosity_probability**2
#################################################################
def convert(h, m, s): #Degrees minutes seconds to degrees (More for applied code than here)
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

## You have an array and a value 
    #interp1D
    
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

L1 = np.linspace(56, 59, 100) #In J now
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

lum_atsource = np.zeros(N)
received_luminosity = np.zeros(N)
cumul_N = np.zeros(N)

lum_list = list(L_rank)


#error_radius = 10 * random(N)

angles = np.arccos(2 * random(N) - 1)

ra_rand = uniform(0, 360, size = N)
dec_rand = (np.pi/2 - angles) * (180/np.pi) ##Gives you random ra and dec

r_max = 5000 # 5 Gpc
r = r_max * random(N)**(1/3)

RandL = random(N)
received_luminosity = np.zeros(N)  ## Used to find dummy luminosities
cumul_N = np.zeros(N)


df_dL = df_dL.values.tolist()
a = np.zeros(N)
b = np.zeros(N)

for i in range(N):
       
    a[i] = find_nearest(cumuL, RandL[i])
    if a[i] == len(L1):
        a[i] = len(L1) - 1
    b[i] = L1[int(a[i])] #int(a[i])
    received_luminosity[i] = Luminosity_for_convert(b[i], r[i])

error_radius = 2 * (2.62)

ranks = np.zeros(shape = (N, 5))

for i in range(0, N):
    print(i + 1)
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

    ra_prime = np.zeros(N); dec_prime = np.zeros(N)#1612646.0
    
    #for k in range(1):
        #rota = rotation(m[i, :], rotation_angle[i]) ###Ammend this!!
        
        #x_prime = mat_mul(rota, xyz) #rota * xyz[i, :]
        
    x_prime = axis_rotation(m[i, :], xyz, rotation_angle[i])
    
    xmod = modulus(x_prime)
    x_prime /= xmod
    ra_prime[i], dec_prime[i] = back_convert(x_prime)
        
    ## These will be used as test 
    
    q, t, df_sliced = reduction(abs(ra_prime[i]), dec_prime[i], df_master) ## Reduces the catalog by RA and dec
    
    if df_sliced.shape[0] == 0:
        print("Did not place")
        next
    else:
        print()
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
        
    #angular_distance = Ang_Dist(ra, ra_prime, dec, dec_prime)
    
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
    
    for k in range(5):
        
        placement[i][k] = df_sliced[["Rank"]].head(5).values.tolist()[k][0]
        placement2[i][k] = df_sliced2[["Rank2"]].head(5).values.tolist()[k][0]
        placement3[i][k] = df_sliced3[["Rank3"]].head(5).values.tolist()[k][0]
        placement4[i][k] = df_sliced4[["Rank4"]].head(5).values.tolist()[k][0]
        placement5[i][k] = df_sliced5[["Rank5"]].head(5).values.tolist()[k][0]
        placement6[i][k] = df_sliced6[["Rank6"]].head(5).values.tolist()[k][0]
        placement7[i][k] = df_sliced7[["Rank7"]].head(5).values.tolist()[k][0]
        placement8[i][k] = df_sliced8[["Rank8"]].head(5).values.tolist()[k][0]
        placement9[i][k] = df_sliced9[["Rank9"]].head(5).values.tolist()[k][0]
        placement10[i][k] = df_sliced10[["Rank10"]].head(5).values.tolist()[k][0]
        placement11[i][k] = df_sliced11[["Rank11"]].head(5).values.tolist()[k][0]
        placement12[i][k] = df_sliced12[["Rank12"]].head(5).values.tolist()[k][0]
        placement13[i][k] = df_sliced13[["Rank13"]].head(5).values.tolist()[k][0]
        placement14[i][k] = df_sliced14[["Rank14"]].head(5).values.tolist()[k][0]
        placement15[i][k] = df_sliced15[["Rank15"]].head(5).values.tolist()[k][0]
        placement16[i][k] = df_sliced16[["Rank16"]].head(5).values.tolist()[k][0]
        placement17[i][k] = df_sliced17[["Rank17"]].head(5).values.tolist()[k][0]
        placement18[i][k] = df_sliced18[["Rank18"]].head(5).values.tolist()[k][0]
        placement19[i][k] = df_sliced19[["Rank19"]].head(5).values.tolist()[k][0]
        placement20[i][k] = df_sliced20[["Rank20"]].head(5).values.tolist()[k][0]
        placement21[i][k] = df_sliced21[["Rank21"]].head(5).values.tolist()[k][0]
        placement22[i][k] = df_sliced22[["Rank22"]].head(5).values.tolist()[k][0]
        placement23[i][k] = df_sliced22[["Rank23"]].head(5).values.tolist()[k][0]
    
    
