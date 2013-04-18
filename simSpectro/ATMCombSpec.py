from scipy import *
from constants import *
import pylab as pl
import numpy as np 

class ATMCombSpec:

    def __init__(self):
        self._initList()  
        
    def _initList(self):
        self.wl = []
        self.nu = []
        self.tr = []
        self.ab = []

    def setWithTrans(self, wl, tr):
        self.wl = np.flipud(wl)
        self.tr = np.flipud(tr)
        self.ab = 1.0 - self.tr
        self.nu = clight/self.wl
        
    def read(self, filename):
        self._initList()
        myfile = open(filename, 'r')
        for line in myfile.readlines():
            line = line.split()
            freq = float(line[0])*3.e10
            trans = float(line[1])
            self.nu.append(freq)
            self.wl.append(clight/freq)
            self.tr.append(trans)
            self.ab.append(1.-trans)
        myfile.close()
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
        print len(self.wl)
