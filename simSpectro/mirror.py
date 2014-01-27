# -*- coding: utf-8 -*-

#############################################################################
#          mirror - a class defining the properties of the mirror           #
#                                                                           #
#        author:          Gurvan BAZIN                                      #
#                         Universitäts-Sternwarte-München                   #
#        latest revision: June 2010                                         #
#############################################################################

# include some useful libraries
# and class definitions
from scipy import *

from constants import *


class mirror:

    def __init__(self):
        # wavelength and response arrays
        self.wl = []
        self.re = []

    def read(self, filename):
        # read the file
        file = open(filename, 'r')
        for line in file.readlines():
            line = line.split()
            if line[0] == '#':
                continue
            self.wl.append(float(line[0]))
            self.re.append(float(line[1])/100.)
        file.close()
        self.wl = array(self.wl)
        self.re = array(self.re)

    def getwl(self):
        # return the wavelength array
        return self.wl

    def getnu(self):
        # return the frequency array
        return flipud(clight/self.wl)
    
    def getre(self):
        # return the response array in increasing wavelength
        return self.re
    
    def getrenu(self):
        # return the response array in increasing frequency
        return flipud(self.re)
