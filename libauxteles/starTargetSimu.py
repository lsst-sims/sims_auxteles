'''
Created on 4 dec. 2012

@author: colley
'''
import kurucz as kur
import numpy as np 



class StarTargetSimu(object):

    def __init__(self, oKur):
        self.ParAll = np.zeros(3,  dtype = np.float32)
        assert isinstance(oKur, kur.Kurucz)
        self._oKur = oKur
        self._NbStar = 4
        self._NbWL = len(oKur.getWL())
        self._NbPstar = 2
        self._aParam = np.zeros((self._NbStar,2), dtype = np.float32)
        self._aFlux= np.zeros((self._NbStar, self._NbWL), dtype = np.float32)
        # F
        self._aParam[0,0] = 7200
        self._aParam[0,1] = 4.34
        # F
        self._aParam[1,0] = 6440
        self._aParam[1,1] = 4.2
        # G 
        self._aParam[2,0] = 5770
        self._aParam[2,1] = 4.49
        # A
        self._aParam[3,0] = 9520    
        self._aParam[3,1] = 4.14
        # precompute flux
        self.ParAll[0] = -3
        self._setaFlux()    
    
    
    def _setaFlux(self):    
        for idx in range( self._NbStar):
            self.getAllParam(self._aParam[idx])
            self._aFlux[idx] = self._oKur.getFluxInterLin(self.ParAll)
    
                
    def getFluxIdx(self, pIdx):
        return  self._aFlux[pIdx]
    
    def getAllParam(self, pAr):
        self.ParAll[1:] = pAr
        return self.ParAll
        
        
class StarTargetSimuAll(object):
    def __init__(self, oKur, pNbpar=2):
        """
        oKur: kurucz model to compute star flux with param target stars
        pNbpar :  fix number of parameter other are considered as constant
        """
        assert isinstance(oKur, kur.Kurucz)   
        self._NbPstar = pNbpar
        self._oKur = oKur     
        self._NbStar = 4
        self._aParam = np.zeros((self._NbStar, 3), dtype = np.float32)
        self._aParam = np.array([[0.0 , 7700, 1.74],
                                 [0.0 , 6440, 4.34],
                                 [0.0 , 5150, 2.54],
                                 [0.0 , 9520, 4.14 ]], dtype=np.float32)        
        self._NbWL = len(oKur.getWL())       
        self._aFlux= np.zeros((self._NbStar, self._NbWL), dtype = np.float32)
        self._setFlux() 
               
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
        par = np.copy(self._aParam[idx])
        par[1:] = pTempGra
        return par
  
    def getFluxIdx(self, pIdx):
        return  self._aFlux[pIdx]
        
        