'''
Created on 14 nov. 2012

@author: colley
'''


import numpy as np
import templatMODTRAN as tmod
import pylab as pl


class BurkeAtmModelv1(object):
    '''
    Burke and All model from paper 'precision  determination of atmosphere extinction at optical and NIR wavelengths' 2010
    '''

    def __init__(self, ModFileTempl, arrayWL, pressure=782):
        '''        
        '''
        self._Tpl = tmod.TemplateMODTRAN(ModFileTempl)
        self._PressureRef = pressure
        self._NbPar = 18
        self._Par  = None
        self.setSlopeAirMassComponents(0.98, 0.71, 0.97, 0.59)
        self._aWL = arrayWL
        self._WL0 = 6750
        
        
    def setSlopeAirMassComponents(self, Raygleigth, h2o, o3, o2):
        self._SlopeAMassRay = Raygleigth
        self._SlopeAMassh2o = h2o 
        self._SlopeAMasso3 = o3
        self._SlopeAMasso2 = o2
    
    
    def setBurkeParamModel(self, par):
        assert (len(par) == self._NbPar)
        self._Par = par
        
                
    def printBurkeModel(self):
        if self._Par == None:
            print "No parameter yet defined !"
            return 
        print "Constants"
        print "  airmass slope:"
        print "    H20 :   ", self._SlopeAMassh2o
        print "    Rayleigth :   ", self._SlopeAMassRay
        print "    O2 :   ", self._SlopeAMasso2
        print "    O3 :   ", self._SlopeAMasso3
        print "  ref. pressure:  %d  mb"% self._PressureRef
        print "  ref. wavelength:  %g  A"% self._WL0
        print "Parameters"
        print "  Tgray :   ", self._Par[0]
        print "  Tau0 :   ", self._Par[1]
        print "  Tau1 :   ", self._Par[2]
        print "  Tau2 :   ", self._Par[3]
        print "  alpha :   ", self._Par[4]
        print "  Cmol :   ", self._Par[5]
        print "  C_O3 :   ", self._Par[6]
        print "  C_H2O :   ", self._Par[7:18]
        pl.figure()
        xt = np.linspace(0,3600*6, 256)
        pl.plot(xt, self.getC_H20(xt))
        pl.xlabel("s")
        pl.grid()
        pl.title("C_H20")
        self._Tpl.plotTemplate()
        
        
    def setAerosolGrayModel(self, par, Tgray, tau0, tau1, tau2, alpha):
        par[0] = Tgray
        par[1] = tau0
        par[2] = tau1
        par[3] = tau2
        par[4] = alpha
        return par
    
    
    def setMolModel(self,par, Cmol):
        par[5] = Cmol
        return par
    
    
    def setOzoneModel(self, par,O3):
        par[6] = O3
        return par
        
    
    
    def ComputeAtmTransmission(self, alt, az, t, presure):
        self._AirMass = 1/np.cos(np.pi/2 - alt)
        self._Alt = alt
        self._Az = az
        self._EW = np.cos(alt)*np.sin(az)
        self._NS = np.cos(alt)*np.cos(az)
        self._Time = t
        self._PresRat = presure/self._PressureRef
        ret = 1
        ret *= self._transGrayAero()
        ret *= self._transMols()
        ret *= self._transMola()
        ret *= self._transO3()
        ret *= self._transH2O()
        return ret
    
    
    def getC_H20(self, t): 
        return self._Par[7] + t*(self._Par[7] + t*(self._Par[8] + t*self._Par[9]))
    
        
    def _transGrayAero(self):
        tau  = (self._Par[1] + self._EW*self._Par[2] + self._NS*self._Par[3])
        tau *= (self._aWL/self._WL0)**self._Par[4]
        return self._Par[0]*np.exp(-self._AirMass*tau)
    
    
    def _transMols(self):
        return (1 - self._Par[5]*self._PresRat*self._Tpl._Trmols*self._AirMass) 
    
    
    def _transMola(self):
        return (1 - np.sqrt(self._Par[5]*self._PresRat)*self._Tpl._Trmola*self._AirMass) 
    
    
    def _transO3(self):
        return (1 - self._Par[6]*self._Tpl._Tr03*self._AirMass)
    
    
    def _transH2O(self):
        return (1 - self.getC_H20(self._Time)*self._Tpl._TrH2O*self._AirMass)
    
        
        