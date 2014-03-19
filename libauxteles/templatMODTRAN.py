'''
Created on 6 nov. 2012

@author: colley
'''

import tools as tl 
import matplotlib.pyplot as pl
import numpy as np
import scipy.interpolate as sci


class TemplateMODTRAN(object):
    '''
    read output MODTRAN file to extract template absorption of
    * 03
    * molecular scattering
    * molecular absorption
    * H20
    '''
    def __init__(self, fileName, unit = 1e-9):
        out = tl.readTextFileColumn(fileName)
        self._wl = out[:,0]
        self._Unit = unit
        # _Axx absorption 
        self._AO3   = 1 - out[:,3]
        # keyword MODTRAN: MOL_SCAT
        self._Amols = 1 - out[:,5]
        # keyword MODTRAN: O2
        self._Amola = 1 - out[:,7]
        self._AH2O  = 1 - out[:,2]*out[:,4]
        perm = np.argsort(self._wl)
        self._AO3 = self._AO3[perm]
        self._Amols = self._Amols[perm]
        self._Amola = self._Amola[perm]
        self._AH2O = self._AH2O [perm]
        self._wl = self._wl[perm]
        self._ref = 650.0 
        self._res = 2600.0 # resolution read in file
        self._AM = 1.0 
    
    def setAirMass(self, pAM):
        self._AM = pAM

    
    def restrictWL(self, pWLmin, pWlmax):
        """
        restrict wavelength and all MODTRAN template to pWLmin, pWlmax 
        """
        iMin, iMax = tl.indexInIntervalCheck(self._wl, pWLmin, pWlmax)
        print  iMin, iMax
        # suppress wl
        self._wl = self._wl[iMin:iMax]
        self._AO3 = self._AO3[iMin:iMax]
        self._AH2O = self._AH2O[iMin:iMax]
        self._Amola = self._Amola[iMin:iMax]
        self._Amols = self._Amols[iMin:iMax]
        
        
    def resampleBetween(self, pWLmin, pWLmax, pNb):
        newWL = np.linspace(pWLmin, pWLmax, pNb, True)
        self.resample(newWL)
    
    
    def _setZeroNeg(self, a):
        idx = np.where(a <0) [0]
        if len(idx) > 0:
            a[idx] = 0.0
            print a[idx]
            raise
        return a
    
    
    def _threasoldOne(self, a):
        idx = np.where(a > 1.0) [0]
        if len(idx) > 0:
            a[idx] = 1.0
            print "threasoldOne detection"
        return a
    
    
    def downgradeTemplate(self, resOut):
        self._AO3 = self._setZeroNeg(tl.downgradeResol(self._wl, self._AO3, self._res, resOut, self._ref))
        self._AH2O = self._setZeroNeg(tl.downgradeResol(self._wl, self._AH2O, self._res, resOut, self._ref))
        self._Amola = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amola, self._res, resOut, self._ref))
        self._Amols = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amols, self._res, resOut, self._ref))        
        self._res = resOut
        
        
    def downgradeTemplateAndResample(self, resOut, pWL):
        self._AO3 = self._setZeroNeg(tl.downgradeResol(self._wl, self._AO3, self._res, resOut, self._ref, pWL))
        self._AH2O = self._setZeroNeg(tl.downgradeResol(self._wl, self._AH2O, self._res, resOut, self._ref, pWL))
        self._Amola = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amola, self._res, resOut, self._ref, pWL))
        self._Amols = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amols, self._res, resOut, self._ref, pWL))
        self._res = resOut       
        self._wl = pWL
        
        
    def resample(self, pWL):
        """
        resample template for a given wavelength vector
        """
        self._AO3 = tl.interpolLinear(self._wl, self._AO3, pWL)        
        self._AH2O = tl.interpolLinear(self._wl, self._AH2O, pWL)        
        self._Amola = tl.interpolLinear(self._wl, self._Amola, pWL)        
        self._Amols= tl.interpolLinear(self._wl, self._Amols, pWL)        
        self._wl = pWL
                
                
    def convertWaveLength(self, unit):
        fact = (self._Unit/unit)
        self._wl *= fact
        self._ref *= fact
        self._Unit = unit
                
                
    def getTrO3(self):
        return 1.0 - self._AO3
        
    def getTrmols(self):
        return 1.0 - self._Amols    
    
    def getTrmola(self):
        return 1.0 - self._Amola    
    
    def getTrH2O (self):
        return 1.0 - self._AH2O     
    
    def getTrAll(self):
        return self.getTrO3()*self.getTrmols()*self.getTrmola()*self.getTrH2O()
    
    def plotTemplate(self, WLmin=None, WLmax=None, pTitle=""):
        if WLmin ==None:
            idxMin = 0
        else:
            idx = np.where(self._wl >= WLmin)[0]
            idxMin = idx[0]
            
        if WLmax ==None:
            idxMax = len(self._wl)
        else:
            idx = np.where(self._wl <= WLmax)[0]
            idxMax = idx[-1]        
        pl.figure()
        pl.plot(self._wl[idxMin:idxMax], self.getTrO3()[idxMin:idxMax])
        pl.plot(self._wl[idxMin:idxMax], self.getTrmols()[idxMin:idxMax])
        pl.plot(self._wl[idxMin:idxMax], self.getTrmola()[idxMin:idxMax])
        pl.plot(self._wl[idxMin:idxMax], self.getTrH2O()[idxMin:idxMax])
        TrAll = self.getTrAll()
        pl.plot(self._wl[idxMin:idxMax], TrAll[idxMin:idxMax],'y')
        pl.title("Template MODTRAN transmission "+pTitle)
        pl.legend(["03","mols","mola/02","H2O","all"],loc=4)
        pl.grid()
        pl.xlabel("%gm"%self._Unit )
        pl.ylabel("%")
        pl.ylim([0,1.02])


###############################################################################
class TemplateMODTRANAirMass(TemplateMODTRAN):
###############################################################################
    """
    MODTRAN template with air mass linear interpolation
    """
    def __init__(self):
        rootPackage = tl.getRootPackage()
        self._aTplt = np.load(rootPackage+"/data/modtran/tpltMODTRAN.npy")
        self._aTplt = self._aTplt.astype(np.float64)
        self._aTplt = np.rollaxis(self._aTplt,1)
        assert(self._aTplt.shape[0] == 4)
        assert(self._aTplt.shape[2] > 100)
        for aTplt in self._aTplt:
            for tplt in aTplt:
                tplt[:] = self._threasoldOne(tplt)        
        # now shape is (4, ~12,  ~5212)
        self._Lcomp = ["O3", "MolS", "MolA", "H2O"]        
        self._wl = np.load(rootPackage+"/data/modtran/tpltMODTRAN_wlnm.npy")
        self._Aam = np.load(rootPackage+"/data/modtran/tpltMODTRAN_AM.npy")
        self._Unit = 1e-9
        self._ref = 650.0 
        self._res = 2600.0 # resolution read in file
        self._AM = 1.0
        self._initInterpolation()
        
    def setAirMass(self, pAM):
        print pAM
        assert(self._Aam[0] <= pAM)
        assert(pAM <= self._Aam[-1])
        self._AM = np.array([pAM])

       
    def restrictWL(self, pWLmin, pWlmax):
        """
        restrict wavelength and all MODTRAN template to pWLmin, pWlmax 
        """
        iMin, iMax = tl.indexInIntervalCheck(self._wl, pWLmin, pWlmax)
        self._wl = self._wl[iMin:iMax]
        self._aTplt = self._aTplt[...,iMin:iMax]
        self._initInterpolation()
    
    
    def downgradeTemplate(self, resOut):         
        for aTplt in self._aTplt:
            for tplt in aTplt:
                tplt[:] = self._threasoldOne(tl.downgradeResol(self._wl, tplt,\
                                                          self._res, resOut,\
                                                          self._ref))              
        self._res = resOut
        self._initInterpolation()
    
    
    def downgradeTemplateAndResample(self, resOut, pWL):
        out = np.empty(( 4,self._aTplt.shape[1], len(pWL)), dtype=np.float64)
        for aTplt, aTout in zip(self._aTplt, out):
            for tplt, tout in zip(aTplt, aTout):
                tout[:] = self._threasoldOne(tl.downgradeResol(self._wl, tplt,\
                                                          self._res, resOut,\
                                                          self._ref, pWL))        
        self._res = resOut
        del self._aTplt
        self._aTplt = out
        self._wl = pWL
        self._initInterpolation()
    
    
    def _initInterpolation(self):
        am  = self._Aam
        self._DoInter = {}
        for comp, tplt in zip(self._Lcomp, self._aTplt):
            self._DoInter[comp] = sci.interp1d(am, tplt.T, 'linear')


    def resample(self, pWL):
        """
        resample template for a given wavelength vector
        """
        out = np.empty(( 4,self._aTplt.shape[1], len(pWL)), dtype=np.float64)
        for tplt, tout in zip(self._aTplt, out):
            tout[...] = sci.interp1d(self._wl, tplt,'linear',\
                                     bounds_error=True)(pWL)
        del self._aTplt
        self._aTplt = out
        self._wl = pWL
        self._initInterpolation()
        
    
    def getTrO3(self):        
        return self._DoInter["O3"](self._AM).ravel()
    
    
    def getTrmols(self):
        return self._DoInter["MolS"](self._AM).ravel()
    
    
    def getTrmola(self):
        return self._DoInter["MolA"](self._AM).ravel()
    
    
    def getTrH2O(self):
        return self._DoInter["H2O"](self._AM).ravel()

