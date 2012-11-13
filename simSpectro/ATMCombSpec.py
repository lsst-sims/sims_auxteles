from scipy import *
from constants import *
import pylab as pl

class ATMCombSpec:

    def __init__(self):
        self.wl = []
        self.nu = []
        self.tr = []
        self.ab = []

    def read(self, filename):
        file = open(filename, 'r')
        for line in file.readlines():
            line = line.split()
            freq = float(line[0])*3.e10
            trans = float(line[1])
            self.nu.append(freq)
            self.wl.append(clight/freq)
            self.tr.append(trans)
            self.ab.append(1.-trans)
        file.close()
        self.nu = array(self.nu)
        self.wl = array(self.wl)
        self.tr = array(self.tr)
        self.ab = array(self.ab)
        #self.plot(filename)
        
    def gettrans(self):
        return array([self.nu, self.wl, self.tr])

    def getabs(self):
        return array([self.nu, self.wl, self.ab])
    
    def plot(self,fName):
        pl.figure()        
        pl.xlabel("wavelength nm")
        pl.ylabel("%")
        pl.grid()
        pl.title("Atmosphere file %s"%fName)
        pl.plot(self.wl, self.tr)