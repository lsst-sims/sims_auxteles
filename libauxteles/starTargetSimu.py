'''
Created on 4 dec. 2012

@author: colley
'''
import kurucz as kur
import numpy as np 
  
        
class StarTargetSimuAll(object):
    def __init__(self, oKur, pNbpar=2):
        """
        oKur: kurucz model to compute star flux with param target stars
        pNbpar :  fix number of parameter other are considered as constant
        """
        assert isinstance(oKur, kur.Kurucz)
        self._NbPstar = pNbpar 
        self._oKur = oKur             
#        self._aParam = np.array([ [0.0 , 7700, 1.74]], dtype=np.float64) 
        self._aParam = np.array([ [0.0 , 7700, 1.74],
                                  [0.0 , 6440, 4.34],
                                  [0.0 , 5150, 2.54],
                                  [0.0 , 9520, 4.14] ], dtype=np.float64) 
        self._NbStar = self._aParam.shape[0]      
        self._NbWL = len(oKur.getWL())       
        self._aFlux= np.zeros((self._NbStar, self._NbWL), dtype = np.float32)
        self._setFlux() 
    
    def getMeanTemp(self):
        return  self._aParam[:,1].mean()
           
    def getAllTemperature(self):
        return self._aParam[:,1]
    
    def addStar(self, pParm):
        self._NbStar += 1
        self._aParam = np.vstack((self._aParam, pParm))
        # pas optimal mais rapide a ecrire
        self._aFlux= np.zeros((self._NbStar, self._NbWL), dtype = np.float32)
        self._setFlux()
        
    def _setFlux(self):    
        for idx in range( self._NbStar):            
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
    
    def getWL(self):
        return self._oKur.getWL()
