# -*- coding: utf-8 -*-
'''
Created on 31 août 2012

@author: colley
'''

import sys

import numpy as np
import matplotlib.pyplot as pl


# en m/s
S_SpeedLight = 3e8

def freq2wavelength(a, unit=1e-9):
    """
    a en Hertz
    defaut unit nm
    """
    return (S_SpeedLight/a)/unit


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
        c3.append(line[2])
    file.close()
    c1 = np.array(c1, dtype=np.float64)
    c2 = np.array(c2, dtype=np.float64)
    c3 = np.array(c3, dtype=np.float64)
    return c1, c2, c3


def plot_syscaldEdn(fName):
    c1, c2, c3 = read3column(fName)
    pl.figure()
    pl.xlabel("wavelength nm")
    pl.ylabel("?")
    pl.grid()
    pl.title("Spectrum file %s"%fName)
    pl.plot(freq2wavelength(c1), c2)
    
    
def plot_syscaldEdl(fName):
    c1, c2, c3 = read3column(fName)
    pl.figure()
    pl.xlabel("wavelength nm")
    pl.ylabel("?")
    pl.grid()
    pl.title("Spectrum file %s"%fName)
    pl.plot(c1, c2)

#
#
#

fName1 = 'alpha_lyr_stis_004cut2.ascii_Modtran_10.dat__Modtran_10.dat_syscaldEdn.dat'
fName2 = 'alpha_lyr_stis_004cut2.ascii_Modtran_10.dat__Modtran_10.dat_syscaldEdl.dat'
fName1 = 'alpha_lyr_stis_004cut2.ascii_Modtran_10.dat_dEdn.dat'
fName2 = 'alpha_lyr_stis_004cut2.ascii_Modtran_10.dat__Modtran_10.dat_syscaldEdl.dat'

#MyArgs = sys.argv

plot_syscaldEdn(fName1)
plot_syscaldEdl(fName2)
pl.show()

