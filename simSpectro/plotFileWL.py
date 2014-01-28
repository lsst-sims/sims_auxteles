# -*- coding: utf-8 -*-
'''
Created on 31 août 2012

@author: colley
'''

import sys

import numpy as np
import matplotlib.pyplot as pl


# en m/s
G_SpeedLight = 3e8

def freq2wavelength(a, unit=1e-9):
    """
    a en Hertz
    defaut unit nm
    """
    return (G_SpeedLight/a)/unit


def read3column(fName):
    c1 = []
    c2 = []
    c3 = []    
    file = open(fName, 'r')
    for line in file.readlines():
        if line[0] == "#":
            continue
        line = line.split()
        c1.append(line[0])
        c2.append(line[1])
        try:
            c3.append(line[2])
        except:
             c3.append(0)
    file.close()
    c1 = np.array(c1, dtype=np.float64)
    c2 = np.array(c2, dtype=np.float64)
    c3 = np.array(c3, dtype=np.float64)
    return c1, c2, c3



    
    
def plot_syscaldEdl(fName):
    c1, c2, c3 = read3column(fName)
    pl.figure()
    pl.xlabel("wavelength nm")
    pl.ylabel("?")
    pl.grid()
    pl.title("Spectrum file %s"%fName)
    pl.plot(c1, c2)




MyArgs = sys.argv

plot_syscaldEdl(MyArgs[1])
pl.show()
