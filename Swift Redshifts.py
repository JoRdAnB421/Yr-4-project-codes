# -*- coding: utf-8 -*-
"""
Created on Sun May  2 18:39:47 2021
@author: Jordan
Swift SGRBs 
"""

import numpy  as np
import pandas as pd

#import data
file = pd.read_csv("Swift SGRBs.txt", delimiter = "\t")
redshift_file = file.Redshift.str.extract(r'(\d+\.?\d*)')

#pull out only redshifts
redshifts = np.array([], dtype = float)

for i in redshift_file[0]:
    if type(i) == str:
        redshifts = np.append(redshifts, i)