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
    
#
# PRIVATE
#

    def __init__(self, ModFileTempl, timeObs, pressure=782):
        '''        
        '''
        self._Tpl = tmod.TemplateMODTRAN(ModFileTempl)
        # convert in angstrom
        self._Tpl.convertWaveLength(1.e-10)
        self._PressureRef = pressure*1.0
        self._NbObs = len(timeObs)
        self._TimeObs = timeObs
        self._NbParNoH20 = 7
        # nb paral H20 is egal to nbObs
        self._NbPar = self._NbParNoH20 + self._NbObs 
        self._Par  = None
        self.setModelSlopeAirMassComponents(0.98, 0.71, 0.97, 0.59)
        self._aWL = self._Tpl._wl
        self._WL0 = 6750.0

    def _transGrayAero(self):
        tau  = (self._Par[1] + self._EW*self._Par[2] + self._NS*self._Par[3])
        tau *= (self._aWL/self._WL0)**self._Par[4]
        return self._Par[0]*np.exp(-self._AirMass*tau)
        
    def _transMols(self):
        return (1 - self._Par[5]*self._PresRat*self._Tpl._Amols*self._AirMass) 
    
    def _transMola(self):
        return (1 - np.sqrt(self._Par[5]*self._PresRat)*self._Tpl._Amola*self._AirMass) 
    
    def _transO3(self):
        return (1 - self._Par[6]*self._Tpl._A03*self._AirMass)
    
    def _transH2O(self):
        return (1 - self.getC_H20(self._idxTime)*self._Tpl._AH2O*self._AirMass)

#       
# PUBLIC
#

# setter

    def setParamAerosolGray(self, Tgray, tau0, tau1, tau2, alpha):
        self._Par[0] = Tgray
        self._Par[1] = tau0
        self._Par[2] = tau1
        self._Par[3] = tau2
        self._Par[4] = alpha
                
    def setParamMol(self, Cmol):
        self._Par[5] = Cmol
    
    def setParamOzone(self, O3):
        self._Par[6] = O3
        
    def setParamH20(self, aIn):        
        assert (len(aIn) == self._NbObs)
        self._Par[self._NbParNoH20:self._NbPar] = aIn
        
    def setParamExample1(self):
        self._Par = np.zeros(self._NbPar, dtype=np.float32)
        # from [1], table 3, 2007, 2nov
        self.setParamAerosolGray(0.98, 3.9/100, 0.02/100, -0.03/100, -1.70)
        # from [1], 4.results , first line
        self.setParamMol(0.91)
        # from [1], table 3, 2007, 2nov
        self.setParamOzone(0.8)
        self.setParamH20(np.ones(self._NbObs))

    def setParamNoEffect(self):
        """
        to have template MODTRAN with any modification
        """
        self._Par = np.zeros(self._NbPar, dtype=np.float32)
        # from [1], table 3, 2007, 2nov
        self.setParamAerosolGray(1.0, 0, 0, 0, 0)
        # from [1], 4.results , first line
        self.setParamMol(1.0)
        # from [1], table 3, 2007, 2nov
        self.setParamOzone(1.0)
        self.setParamH20(np.ones(self._NbObs))
        
    def setModelSlopeAirMassComponents(self, Raygleigth, h2o, o3, o2):
        self._SlopeAMassRay = Raygleigth
        self._SlopeAMassh2o = h2o 
        self._SlopeAMasso3 = o3
        self._SlopeAMasso2 = o2
    
    def setParam(self, par):
        assert (len(par) == self._NbPar)
        self._Par = par
        
    def setVarObs(self, alt, az, t, presure):
        self._AirMass = np.fabs(1/np.cos(np.pi/2 - alt))
        self._Alt = alt
        self._Az = az
        self._EW = np.cos(alt)*np.sin(az)
        self._NS = np.cos(alt)*np.cos(az)
        self._Time = t
        aIdx = np.where(np.fabs(t-self._TimeObs) < 0.01)[0]
        if len(aIdx) == 0:
            print "No observation at this time ", t
            print self._TimeObs
            return None
        self._idxTime = aIdx[0]
        self._PresRat = presure*1.0/self._PressureRef        
        
# getter

    def computeAtmTrans(self):   
        ret = 1
        ret *= self._transGrayAero()
        ret *= self._transMols()
        ret *= self._transMola()
        ret *= self._transO3()
        ret *= self._transH2O()
        self._CurtTrans = ret
        return ret
    
    def computeAtmTransAt(self, alt, az, t, presure):
        self.setVarObs(alt, az, t, presure)  
        ret = self.computeAtmTrans()
        return ret
    
    def getC_H20(self,IdxObs): 
        return self._Par[self._NbParNoH20+IdxObs]

    def computeDeltaDataModel(self,alt, az, t, presure, param, data ):
        self.setParam(param)
        return data - self.computeAtmTrans(alt, az, t, presure)

# pretty print, plot

    def printAndPlotBurkeModel(self):
        if self.printBurkeModel() == None: return None
        pl.figure()
        pl.plot(self._TimeObs, self._Par[self._NbParNoH20:self._NbPar],"*")
        pl.xlabel("s")
        pl.grid()
        pl.title("C_H20")
        self._Tpl.plotTemplate()
                
    def printBurkeModel(self):
        if self._Par == None:
            print "No parameter yet defined !"
            return None
        self.printBurkeModelConst()
        self.printBurkeModelParam()
        
    def printBurkeModelConst(self):
        print "BurkeAndAll atmosphere model:"
        print "============================"
        print "Constants"
        print "  airmass slope:"
        print "    H20 :   ", self._SlopeAMassh2o
        print "    Rayleigth :   ", self._SlopeAMassRay
        print "    O2 :   ", self._SlopeAMasso2
        print "    O3 :   ", self._SlopeAMasso3
        print "  ref. pressure:  %d  mb"% self._PressureRef
        print "  ref. wavelength:  %g  A"% self._WL0
        
    def printBurkeModelParam(self):
        print "Parameters"
        print "  Tgray :   ", self._Par[0]
        print "  Tau0 :   ", self._Par[1]
        print "  Tau1 :   ", self._Par[2]
        print "  Tau2 :   ", self._Par[3]
        print "  alpha :   ", self._Par[4]
        print "  Cmol :   ", self._Par[5]
        print "  C_O3 :   ", self._Par[6]
        print "  C_H2O :   ", self._Par[7:self._NbPar]
                   
    def plotCurrentTrans(self):
        if self._CurtTrans == None:
            print "No current transmission !"
            return None    
        pl.figure()
        pl.plot(self._Tpl._wl, self._CurtTrans)
        pl.xlabel("Angstrom")
        pl.ylabel("%")
        pl.grid()
        pl.title("atm trans at: [alt:%.1f,az:%.1f], pres ratio %.2f, time %.1f"%(np.rad2deg(self._Alt),np.rad2deg(self._Az), self._PresRat, self._Time ))
        
        