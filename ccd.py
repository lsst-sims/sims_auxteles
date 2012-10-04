# -*- coding: utf-8 -*-
# modif jmc
#############################################################################
#       ccd - a class defining the ccd of the LSST Auxiliary Telescope      #
#                                                                           #
#        author:          Gurvan BAZIN                                      #
#                         Universitäts-Sternwarte-München                   #
#        latest revision: June 2010                                         #
#############################################################################

# include some useful libraries
# and class definitions
from constants import *
from scipy import *
import pylab as pl
class ccd:

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
            wl = float(line[0])
            self.wl.append(wl)
            re = float(line[1])/100.
            # basic check use of ccd transmission 
            #if wl > 900: re = 0
            self.re.append(re)
        file.close()
        self.wl = array(self.wl)
        self.re = array(self.re)
        #self.plotResponse()
        
    def plotResponse(self):
        pl.figure()
        pl.xlim(200, 1200)
        pl.plot(self.wl, self.re)
        pl.xlabel("longueur d'onde nm")
        pl.title("Rendement quantique du CCD suivant la longueur d'onde")
        pl.grid()
       
        
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
