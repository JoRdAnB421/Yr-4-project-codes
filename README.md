# Yr-4-project-codes
Coallating the various different codes developed off of David Evans' initial framework. Their are different files relating to different jobs, however the main codes used are the Automating Ranking Stats.py which implements the MCMC algorithm for finding the optimal ranking statistic, and the Testing functions from Ex_within.py which allows specific testing of individual statistics. 

# Dependencies 
pylab
numpy
pandas
math
json
scipy
astropy
HEALPix
timeit

# For Sam
I have kept the two important scripts that you will need for your research on the home directory (not in any folders).

# How Testing Specfic Rankings.py works

I will briefly explain what the main script does: 

"Testing Specfic Rankings.py" is the main script which simulates N SGRBs inside and outside of 100 Mpc, for each of these simulated SGRB, we find the relevant galaxies which lie within the error radius (I had it set to the error radius of the Fermi telescope but this could be changed), and apply the choosen statistic to these galaxies given a rank to each of the galaxies which is higher for those more likely to be the host galaxy. We then take measures of these rank values (I took the highest rank, the sum of all ranks and then average of the top 5 ranks) and store the values.

Once this runs over every single SGRB, the script will create histograms for each of the three measures (highest rank, the sum of all ranks and then average of the top 5 ranks), though you won't see these plotted, for each of these it will essential find the one which is able to best separate the two groups (inside and outside SGRBs), this best histogram is then plotted for you to see. Some properties of this histogram are saved into a .csv file which should be self-explanatory given their headers (but contact me for further information). 

The next part of the script takes this best histogram and creates a cumulative plot of the values which can help to better see the distinction between the two groups. Finally the script will calculate the efficiency (the fraction of "inside" SGRBs that will be caught by this statistic) and the false alarm rate (the fraction of "outside" SGRBs that will be incorrectly identified as being "inside").

# How to use/edit the code

To use "Testing Specfic Rankings.py" you will need 4 data files, these were too big to upload to github on their own and so I have put them in a zipped folder called "Data files.rar" in the folder Data files, you will need to unzip these files and have them in the same directory as where you run the code (these are information from the GLADE galaxy catalogue and some cumulative luminosity files).

Also to set the range for the inside SGRBs you will have to remove the galaxies that are past a certain luminosity distance in the csv file before you use that data file in either the "HealPix" or the "Testing Speciifc Rankings.py". Unfortunetly I don't not have a code that does this, however it should not be a had code to write.

in the script itself there are a few things that you can alter to change the statistic your testing and number of simulated SGRBs:

On line 17 you can change the number of simulated SGRBs

On line 347 you can change the error radius (useful for if you're looking with a different satellite

On line 505 you can change the powers on the statistic (this is synoymus with changing the statistic as they will always have the same base 4 terms)

It would be possible to add extra measures of the rank values if you want however you would have to change quite a few parts of the script and so I can lend a hand if it comes to it.

# Using a new catalogue

Patrick mentioned that you will be looking at other catalogues (like SDSS) and to do this you must first put the new catalogue into the "proper HealPix.py" script as you have to split up the catalogue into smaller data files which represent different equal area sectors of the sky. You would also need to make a separate file which is "reduced" that is to say has all of the galaxies outside of your observing limit removed (I had this at 100 Mpc due to the completeness of the Glade galaxy catalogue). 

When using this script you will just have to change the name of the file you're reducing on line 12 and what the sectors will be called on line 50.
