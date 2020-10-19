
import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d

N = 10      ##Change to alter the number of loops the code runs for

## GLADE Master clarify columns

def Ang_Dist(ra1, ra2, dec1, dec2):## Calculates the angular distance between apparent position and galaxy
    
    ra1 *= (np.pi/180); ra2 *= (np.pi/180)
    dec1 *= (np.pi/180); dec2 *= (np.pi/180)
        
    return (180/np.pi) * np.arccos(np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))

def rank(theta, sigma, d_lum, luminosity, luminosity_probability):
    ## Implements a ranking statistic defined in report
    return np.exp(-(theta**2/2*(sigma)**2)) * (1/d_lum * luminosity)[:, 0] * luminosity_probability #* Colour_factor

def convert(h, m, s): #Degrees minutes seconds to degrees (More for applied code than here)
    return h + (m/60) + (s/3600)

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

def modulus(array):  ##Test  ##Finds the modulus of a matrix/array
    
    return np.sqrt(array[0]**2 + array[1]**2 + array[2]**2)

def find_nearest(array, value): #Kind of a hash and not exactly interpolation, but for this point, should be okay
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return idx

## cumulative luminosity greater than random number and take the first

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

def Luminosity_back_convert(L_given, d_L): #  ##Converts luminosity to luminosity at source
    #L = L0/4 *np.pi * d_l**2
    return (L_given) * (4 * np.pi * (3.086e22 * d_L)**2)

def Luminosity_for_convert(L_given, d_L): #  ##Converts luminosity at source to apparent luminosity
    return(L_given)/(4 * np.pi * (3.086e22 * d_L)**2)

def L_func(L_test, c, d_L): ##   ##Takes an input and returns a probability based on the broken power law

    #L_test = received_luminosity[-1]
    #d_L = dl
    
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

def cumulative(array): ###   #Builds cumulative distributions
    
    N = array.shape[0]
    summing = np.zeros(N)
    array = L2
    
    for i in range(1, N):
        df = pd.DataFrame(array[:i])
        summing[i] = df.sum().values.tolist()[0]
        
    return summing
############################################################################################
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

df_master = pd.read_csv("Data Files/GLADE_Master_LGRB.csv", delimiter = ",", low_memory = False) ##GLADE_Master.csv previously defined
df_master.drop("Unnamed: 0", axis = 1, inplace = True)

L1 = np.linspace(55, 59, 100) #In J now
L2, c = L_func1(L1) #  ##Builds broken power law

cumuL = cumulative(L2) ##Luminosity Distribution
df_cumLum =   pd.read_csv("Data Files/Cumulative Luminosity.csv")
df_cumLum.columns = ["NaN", "Cumulative Luminosity"]
normal_c = df_cumLum[["Cumulative Luminosity"]].values[-1][0]
L_rank = df_cumLum[["Cumulative Luminosity"]].values * 1/normal_c
df_cumLum = df_cumLum[["Cumulative Luminosity"]].values#                ## This is all to do with building a usable and callable power law
                                                                        ## 

lum_N = np.linspace(0, df_cumLum.shape[0], df_cumLum.shape[0])

df_dL = df_master[["Luminosity Distance"]]


tests = randint(0, 2, size = N) ## If tests[i] = 0, use test galaxy, or if = 1, choose random point beyond the catalog
dummies = random(N)
RandL = random(N)

gals = np.zeros(N)          ## Picks out a luminosity
gal_index = np.zeros(N)

aa = np.zeros(shape = (N, 5))  # Storing Angular distance
ab = np.zeros(shape = (N, 5))  # Storing Luminosity Distance
ac = np.zeros(shape = (N, 5))  # Storing B Luminosity
ad = np.zeros(shape = (N, 5))  # Storing Luminosity Probability

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

for i in range(N):
    
    gals[i] = find_nearest(L_rank, dummies[i])  ## Picks out galaxies from the cumulative luminosity distribution
       
    a[i] = (find_nearest(cumuL, (RandL[i])))
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
    
    
    ra, dec, df_sliced = reduction(abs(ra_prime[i]), dec_prime[i], df_master) ## Reduces the catalog by RA and dec

    
    Luminosity = df_sliced[["B Luminosity"]].values.tolist()#Luminosity_Handling(np.array(df_sliced[["Absolute B Magnitude"]].values.tolist())) ## Converts A    
    dl = np.array(df_sliced[["Luminosity Distance"]].values.tolist())    
    
    lum_prob, SGR_test = L_func(received_luminosity[i], c, dl) ##Uses the luminosity function to calculate probabilities
    df_sliced["Luminosity Probability"] = lum_prob
    df_sliced["SGR flag"] = SGR_test
    
    angular_distance = Ang_Dist(ra, testr[0], dec, testd[0])
    df_sliced["Angular Distance"] = angular_distance
        
    ranking = rank(angular_distance, error_radius, dl, Luminosity, lum_prob)  ## Uses defined ranking statistic
    df_sliced["Rank"] = ranking
    
    ## Storing values and extending the reduced catalog
    
    df_sliced = (pd.DataFrame.sort_values(df_sliced, by = ["Rank"], ascending = False)) ## Orders resultant sliced array
    
    idi = df_sliced[["Identifier"]].values.tolist() ##Mask handling to check for values
    mask_check = [i for i, val in enumerate(idi) if val == [1]]
    
    Luminosity =  np.array(Luminosity)
    
    if len(mask_check) == 0:
        print("Did not place\n\n\n")
        next
    else:
        placement = mask_check[0] + 1; length = len(idi) + 1
        
        print("Galaxy placed", placement, "out of", str(length), "With distance", int(dl[placement - 1]), "\n\n\n")
        
        percentages[i] = placement/length; distances[i] = int(dl[placement - 1][0]); luminosity_i[i] = int(Luminosity[placement - 1][0])
        rank_host[i] = df_sliced[["Rank"]].values.tolist()[idi.index(max(idi))][0]

    faulty[i, 0] = df_master[["RA"]].values.tolist()[current_i][0]      #ra of galaxy
    faulty[i, 1] = ra_prime[i]                                          #ra of grb
    faulty[i, 2] = df_master[["dec"]].values.tolist()[current_i][0]     #dec of galaxy
    faulty[i, 3] = dec_prime[i]                                         #dec of grb
    
    if math.isnan(rank_host[i]) == True:    
        faulty[i, 4] = 1 #Mask
        break
    
    else:
        faulty[i, 4] = 0 #Mask
        next
    
    for k in range(5):        
        aa[i][k] = np.exp(-(df_sliced[["Angular Distance"]].head(5).values.tolist()[k][0])/error_radius)
        ab[i][k] = df_sliced[["Luminosity Distance"]].head(5).values.tolist()[k][0]
        ac[i][k] = df_sliced[["B Luminosity"]].head(5).values.tolist()[k][0]
        ad[i][k] = df_sliced[["Luminosity Probability"]].head(5).values.tolist()[k][0]

plt.figure(0)
plt.semilogy(percentages[:292], distances[:292], "kx")
plt.title("Distance vs. percentage performance")
plt.ylabel("Log$_{10}$ Distance /Mpc"); plt.xlabel("Percentage placement"); plt.grid()
#plt.xlim(1e-27, 1)
plt.savefig("Distances vs. percentage.png")

plt.figure(1)
plt.semilogy(percentages, b, "kx")
plt.title("Luminosity vs. percentage performance")
plt.ylabel("Log$_{10}$ Luminosity /W"); plt.xlabel("Percentage placement"); plt.grid()
#plt.xlim(1e-27, 1)
plt.savefig("Luminosity vs. percentage.png")


### The following can be used to investigate any values that flag up as false

f_v = [i for i, val in enumerate(faulty[:, 4]) if val == 0]

sets = set(np.arange(0, len(faulty), 1)) - set(f_v)
ft = pd.DataFrame(faulty)
faulty_cols = ["Galaxy RA", "GRB RA", "Galaxy dec", "GRB dec", "Mask"]
ft.columns = faulty_cols

ab_fault = ft.take(list(sets), axis = 0)
ab_vals = ab_fault.values.tolist()[0]
