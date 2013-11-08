'''
Created on 4 dec. 2012

@author: colley
'''
import kurucz as kur
import numpy as np 
import tools as tl
import pylab as pl
        
        
        
class StarTargetSimuAll(object):
    def __init__(self, oKur, pNbpar=2):
        """
        oKur: kurucz model to compute star flux with param target stars
        pNbpar :  fix number of parameter other are considered as constant
        """
        self._oKur = oKur    
        assert isinstance(self._oKur , kur.Kurucz)
        self._NbPstar = pNbpar                  
#        self._aParam = np.array([ [0.0 , 7700, 1.74]], dtype=np.float64) 
        # metalicity, temperature, gravity in kurucz unit
        # F0I 
        # F5V
        # G5III
        # A0V
        self._aParam = np.array([ [0.0 , 7700.0, 1.74],
                                  [0.0 , 6440.0, 4.34],
                                  [0.0 , 5150.0, 2.54],
                                  [0.0 , 9520.0, 4.14] ], dtype=np.float64)
        # visual magnitude, parallaxe in mas 
        self._aMagParal = np.array([[ 6.0, 200.0],
                                    [ 7.5, 250.0],
                                    [ 7.0, 250.0],
                                    [ 6.4, 100.0] ], dtype=np.float64) 
        self._NbStar = self._aParam.shape[0]                          
        self.updateFlux()
        self._dictID = {}
            
           
    def _computeKuruczCoef(self):
        lCoef = []
        for iS in range(len(self._aParam)):
            coef = tl.coefKuruczEarth(self._aParam[iS, 1], self._aMagParal[iS, 0], self._aMagParal[iS, 1])
            lCoef.append(coef)
        self._Kcoef = np.array(lCoef)
    
    
    def getMeanTemp(self):
        return  self._aParam[:, 1].mean()
    
           
    def getAllTemperature(self):
        return self._aParam[:, 1]
    
    
    def addStar(self, pParm):
        self._NbStar += 1
        self._aParam = np.vstack((self._aParam, pParm))
        # pas optimal mais rapide a ecrire
        self._aFlux = np.zeros((self._NbStar, self._NbWL), dtype=np.float32)
        self._setFlux()
    
        
    def _setFlux(self):
        self._NbWL = len(self._oKur.getWL()) 
        self._aFlux = np.zeros((self._NbStar, self._NbWL), dtype=np.float32)
        for idx in range(self._NbStar):            
            self._aFlux[idx] = self._oKur.getFluxInterLin(self._aParam[idx])
    
    
    def addMetGra(self, idx, pTemp):
        par = np.copy(self._aParam[idx]).ravel()
        #print par, pTemp, pTemp[0]
        par[1] = pTemp[0]
        return par
    
    
    def addMet(self, idx, pTempGra):
        par = np.copy(self._aParam[idx]).ravel()       
        #print pTempGra, par[1:]
        par[1:] = pTempGra
        return par
  
  
    def getFluxIdx(self, pIdx):
        return  self._aFlux[pIdx]
  
    
    def IDstar2Idx(self, ID):
        try:
            idx = self._dictID[ID]            
        except KeyError:
            return None
        return idx
                         
    
    def getFluxIdxAtEarth(self, pIdx):
        return self._aFlux[pIdx] * self._Kcoef[pIdx]

    
    def getWL(self):
        return self._oKur.getWL()
    
    
    def plotAll(self):
        pl.figure()
        pl.title('star catalog flux density at earth (T:temp, G:gravity, m:magnitude, p[mas]:parallaxe)')
        lLeg = []
        for idx in range(self._NbStar):            
            pl.plot(self.getWL(), self.getFluxIdx(idx) * self._Kcoef[idx])
            sLeg = "T:%.0f G:%.1f m:%.1f p:%.0f" % (self._aParam[idx, 1], self._aParam[idx, 2], self._aMagParal[idx, 0], self._aMagParal[idx, 1])
            lLeg.append(sLeg)
        pl.ylabel(self._oKur.UnitFluxStr)
        pl.xlabel(self._oKur.UnitWLStr)
        pl.legend(lLeg)
        pl.grid()
        

    def setMetToZero(self):
        self._aParam[:, 0] = 0
        self.updateFlux() 
        
       
    def setAllStars(self, pTemp, pMet, pGra, pMag, pPara, pID):
        self._aParam = np.empty((len(pTemp), 3), dtype=np.float64)
        # kurucz parameters
        self._aParam[:, 0] = pMet
        self._aParam[:, 1] = pTemp
        self._aParam[:, 2] = pGra
        # visual magnitude, parallaxe in mas
        self._aMagParal = np.empty((len(pTemp), 2), dtype=np.float64)
        self._aMagParal[:, 0] = pMag
        self._aMagParal[:, 1] = pPara 
        # compute flux
        self._NbStar = self._aParam.shape[0]
        self._dictID = {}
        for idx, ID in enumerate(pID):
            self._dictID[ID] = idx
        self.updateFlux()                         


    def updateFlux(self):
        self._setFlux()
        self._computeKuruczCoef()        
    
