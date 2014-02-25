'''
Created on 30 janv. 2014

@author: colley
'''

import tools as tl
import matplotlib.pyplot as pl
import burkeAtmModel as bam
import numpy as np
import scipy.optimize as spo
import lmfit as lm



    
#
# MAIN
#





try:
    pl.show()
except AttributeError:
    pass