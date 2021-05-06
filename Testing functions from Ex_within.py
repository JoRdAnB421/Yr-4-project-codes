<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:50:35 2020
@author: Jordan
Testing functions made by David Evans
"""

import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d

def rotation(x, angle):##Test  #Rotation about the z axis
    #need angle in radians
    
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                          [np.sin(angle), np.cos(angle), 0],
                          [0, 0, 1]])
    
    return np.dot(rotation, x)


vec = np.array([1,1,1])


vec2 = rotation(vec, np.pi/2)

=======
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:50:35 2020
@author: Jordan
Testing functions made by David Evans
"""

import pylab as plt; import numpy as np; import pandas as pd
import math; import json; from numpy.random import random, normal, uniform, randint
from scipy.interpolate import interp1d

def rotation(x, angle):##Test  #Rotation about the z axis
    #need angle in radians
    
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                          [np.sin(angle), np.cos(angle), 0],
                          [0, 0, 1]])
    
    return np.dot(rotation, x)


vec = np.array([1,1,1])


vec2 = rotation(vec, np.pi/2)

>>>>>>> 1c9aa4c7d0c5291d1fcd3a218fcf5e23efec9c37
print(vec2)