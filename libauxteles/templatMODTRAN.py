'''
Created on 6 nov. 2012

@author: colley
'''

import tools as tl 
import matplotlib.pyplot as pl
import numpy as np


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
        self._A03   = 1 - out[:,3]
        self._Amols = 1 - out[:,5]
        self._Amola = 1 - out[:,7]
        self._AH2O  = 1 - out[:,2]*out[:,4]
        perm = np.argsort(self._wl)
        self._A03 = self._A03[perm]
        self._Amols = self._Amols[perm]
        self._Amola = self._Amola[perm]
        self._AH2O = self._AH2O [perm]
        self._wl = self._wl[perm]
        self._ref = 650.0 
        self._res = 2600.0 # resolution read in file
    
    
    def restrictWL(self, pWLmin, pWlmax):
        """
        restrict wavelength and all MODTRAN template to pWLmin, pWlmax 
        """
        iMin, iMax = tl.indexInIntervalCheck(self._wl, pWLmin, pWlmax)
        print  iMin, iMax
        # suppress wl
        self._wl = self._wl[iMin:iMax]
        self._A03 = self._A03[iMin:iMax]
        self._AH2O = self._AH2O[iMin:iMax]
        self._Amola = self._Amola[iMin:iMax]
        self._Amols = self._Amols[iMin:iMax]
        
        
    def resampleBetween(self, pWLmin, pWLmax, pNb):
        newWL = np.linspace(pWLmin, pWLmax, pNb, True)
        self.resample(newWL)
    
    
    def _setZeroNeg(self, a):
        idx = np.where(a <0) [0]
        if len(a) > 0:
            a[idx] = 0.0
        return a
    
    def downgradeTemplate(self, resOut):
        self._A03 = self._setZeroNeg(tl.downgradeResol(self._wl, self._A03, self._res, resOut, self._ref))
        self._AH2O = self._setZeroNeg(tl.downgradeResol(self._wl, self._AH2O, self._res, resOut, self._ref))
        self._Amola = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amola, self._res, resOut, self._ref))
        self._Amols = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amols, self._res, resOut, self._ref))        
        self._res = resOut
        
        
    def downgradeTemplateAndResample(self, resOut, pWL):
        self._A03 = self._setZeroNeg(tl.downgradeResol(self._wl, self._A03, self._res, resOut, self._ref, pWL))
        self._AH2O = self._setZeroNeg(tl.downgradeResol(self._wl, self._AH2O, self._res, resOut, self._ref, pWL))
        self._Amola = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amola, self._res, resOut, self._ref, pWL))
        self._Amols = self._setZeroNeg(tl.downgradeResol(self._wl, self._Amols, self._res, resOut, self._ref, pWL))
        self._res = resOut       
        self._wl = pWL
        
        
    def resample(self, pWL):
        self._A03 = tl.interpolLinear(self._wl, self._A03, pWL)        
        self._AH2O = tl.interpolLinear(self._wl, self._AH2O, pWL)        
        self._Amola = tl.interpolLinear(self._wl, self._Amola, pWL)        
        self._Amols= tl.interpolLinear(self._wl, self._Amols, pWL)        
        self._wl = pWL
                
                
    def convertWaveLength(self, unit):
        fact = (self._Unit/unit)
        self._wl *= fact
        self._ref *= fact
        self._Unit = unit
                
                
    def getTr03(self):
        return 1.0 - self._A03
        
    def getTrmols(self):
        return 1.0 - self._Amols    
    
    def getTrmola(self):
        return 1.0 - self._Amola    
    
    def getTrH2O (self):
        return 1.0 - self._AH2O     
    
    def getTrAll(self):
        return self.getTr03()*self.getTrmols()*self.getTrmola()* self.getTrH2O()
    
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
        pl.plot(self._wl[idxMin:idxMax], self.getTr03()[idxMin:idxMax])
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