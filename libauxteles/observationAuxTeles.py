'''
Created on 30 nov. 2012

@author: colley
'''
import numpy as np
import burkeAtmModel as atm
import starTargetSimu as star
import pylab as pl


class ObsSurvey(object):
    def __init__(self):
        self.aFlux = None
        self.aConst = None
        self.aIdxParAtm = None
        self.aIdxParStar= None
        self.NbFlux = 0
        self._NbNight = 0
        self._NbStar = 0
        self._NbConst = 0
        self._NbPatm = 0
        
    def readObsNight(self, pRep):
        """
        fill aFlux,aConst,aIdxParAtm,aIdxParStar
        """
        pass
    
    def getGuessParam(self):
        pass
    
    def plotObs(self,idx, pTitle=""):        
        pl.figure()
        pl.plot(self.aFlux[idx])
        pl.title('Obs flux '+pTitle)
        pl.grid()
        

class ObsSurveySimu01(ObsSurvey):
    """
    simu with random position, atm param
    """
    def __init__(self, pNbNight=1, pNbObsNight=1):
        ObsSurvey.__class__(self)        
        self._NbNight = pNbNight
        self._NbObsNight = pNbObsNight
        self._NbConst = 4
        self._Oatm = atm.BurkeAtmModelv2("")
        self._TrueParam = None
                         
    def setStarTarget(self, pOstar):
        assert isinstance(pOstar, star.StarTargetSimu)    
        self._Ostar = pOstar
        self._NbStar = self._Ostar._NbStar
        self._NbWL = self._Ostar._NbWL
        self._NbPstar  = self._Ostar._NbPstar
        self._NbSimul = self._NbNight*self._NbObsNight
        self._NbFlux  = self._NbSimul*self._NbStar
                
    def setAtmModel(self, pOatm):
        assert isinstance(pOatm, atm.BurkeAtmModelv2)                       
        self._Oatm = pOatm
        self._NbPatm = self._Oatm._NbPar
                
    def setAuxTeles(self, pOauxTeles):
        pass
    
    def _createTrueParAndConst(self):
        # create true param
        self._SimuParNight  = np.zeros((self._NbNight, 8), dtype=np.float32)
        self._SimuParTgray    = np.zeros(self._NbFlux , dtype=np.float32)
        self._SimuParC_H2O  = np.zeros(self._NbSimul , dtype=np.float32)
        self._SimuParObsIdx = np.zeros((self._NbFlux, 3) , dtype=np.int32)
        self._SimuConstObs  = np.zeros((self._NbFlux, self._NbConst) , dtype=np.float32)        
        
    def _randomTrueParAndConst(self):
        self._createTrueParAndConst()
        #######:  _SimuParNight       
        # tau 0
        self._SimuParNight[:,0] = np.random.uniform(0.5, 1.5, self._NbNight)
        # tau 1
        self._SimuParNight[:,1] = np.random.uniform(-0.05, 0.05, self._NbNight)
        # tau 2
        self._SimuParNight[:,2] = np.random.uniform(-0.05, 0.05, self._NbNight)
        # alpha
        self._SimuParNight[:,3] = np.random.uniform(-2, -0.05, self._NbNight)
        # Cmol
        self._SimuParNight[:,4] = np.random.uniform(0.5, 1.2, self._NbNight)
        # C_O3
        self._SimuParNight[:,5] = np.random.uniform(0.5, 1.2, self._NbNight)
        # dC_H2O/dWE
        self._SimuParNight[:,6] = np.random.uniform(-0.05, 0.05, self._NbNight)
        # dC_H2O/dNS
        self._SimuParNight[:,7] = np.random.uniform(-0.05, 0.05, self._NbNight)
        ######## _SimuParObs
        # Tgray
        self._SimuParTgray = np.random.uniform(0.9, 0.99, self._NbFlux)
        # C_H2O
        self._SimuParC_H2O = np.random.uniform(0.2, 1.2, self._NbSimul)
        ######## _SimuParObsIdx
        # idx night
        idx = np.arange(self._NbNight)
        self._SimuParObsIdx[:,0] = np.outer(idx, np.ones(self._NbObsNight*self._NbStar, dtype=np.int32)).ravel()
        # idx star
        idx = np.arange(self._NbStar)
        self._SimuParObsIdx[:,1] = np.outer( np.ones(self._NbSimul, dtype=np.int32), idx ).ravel()   
        # idx simul     
        idx = np.arange(self._NbSimul)
        self._SimuParObsIdx[:,2] = np.outer(idx,  np.ones(self._NbStar, dtype=np.int32) ).ravel()        
        ######## _SimuConstObs
        # time
        self._SimuConstObs[:,2] = np.arange(self._NbFlux)
        # ptg, alt , a
        self._SimuConstObs[:,0] = np.random.uniform(np.deg2rad(45), np.deg2rad(90), self._NbFlux)
        self._SimuConstObs[:,1] = np.random.uniform(0, 2*np.pi, self._NbFlux)
        # pressure
        self._SimuConstObs[:,3] = np.random.normal(782, 20, self._NbFlux)
    
    def getTrueParamAtmIdx(self, pIdx):
        par = np.zeros(10, dtype=np.float64)
        par[0] = self._SimuParTgray[pIdx]
        par[1:9] = self._SimuParNight[self._SimuParObsIdx[pIdx,0]]
        par[9] =  par[8]
        par[8] =  par[7]
        par[7] = self._SimuParC_H2O[self._SimuParObsIdx[pIdx,2]]
        return par
        
    def getConst(self, pIdx):
        return self._SimuConstObs[pIdx,:]
                
    def readObsNight(self, pRep):
        assert isinstance(self._Oatm, atm.BurkeAtmModelv1)
        #
        idxFlux = 0
        # index vector parameters : star param, atm param
        idxPar = self._NbStar*self._NbPstar
        self.aFlux = np.zeros((self._NbFlux, self._NbWL), dtype=np.float64)       
        self.aConst = np.zeros((self._NbFlux, self._NbConst), dtype=np.float64)
        self.aIdxParAtm = np.zeros((self._NbFlux, self._NbPatm), dtype=np.int16)
        self.aIdxParStar = np.zeros((self._NbFlux, self._NbPstar), dtype=np.int16)
        # random parameters
        ###
        self._randomTrueParAndConst()
        idxParNight = idxPar 
        idxParCH20  = idxParNight + 8 
        idxParTgray = idxParNight + 9        
        for idxN in range(self._NbNight): 
            for idxO in range(self._NbObsNight):                
                for idxS in range(self._NbStar):
                    print idxN,idxO,idxS
                    # compute flux 
                    self.aFlux[idxFlux,:] = self._Ostar.getFluxIdx(self._SimuParObsIdx[idxFlux,1])
                    # compute transmission atm 
                    self._Oatm.setConstObs(self.getConst(idxFlux))
                    self._Oatm.setParam(self.getTrueParamAtmIdx(idxFlux))                   
                    # multiply by transmission atm 
                    self.aFlux[idxFlux,:] *= self._Oatm.computeAtmTrans()
                    #self._Oatm.printBurkeModel()
                    #self._Oatm.plotCurrentTrans()
                    #self._Oatm._Tpl.plotTemplate()                    
                    # 
                    self.aConst[idxFlux,:] = self.getConst(idxFlux)
                    # 
                    self.aIdxParAtm[idxFlux, 0]    = idxParTgray                
                    self.aIdxParAtm[idxFlux, 1:7]  = np.arange(idxParNight, idxParNight+6)                    
                    self.aIdxParAtm[idxFlux, 7]    = idxParCH20
                    self.aIdxParAtm[idxFlux, 8:10] = np.arange(idxParNight+6, idxParNight+8)
                    #
                    self.aIdxParStar[idxFlux, 0] = idxS*2
                    self.aIdxParStar[idxFlux, 1] = idxS*2 +1                    
                    idxFlux +=1
                    # add 2 for Tgray, C_H20
                    idxParTgray += 1
                idxParCH20 = idxParTgray
                idxParTgray = idxParCH20 + 1
            # add 8 for tau 0,1,2, alpha, Cmol, CO3, dW_C_H20, dN_C_H20   
            idxParNight = idxParCH20
            idxParCH20 =   idxParNight + 8 
            idxParTgray =  idxParNight + 9                
        # with last observation    
        self._NbParam = self._NbStar*self._NbPstar + self._NbNight*8 + self._NbSimul*(1+self._NbStar)
                
    def getGuessDefault(self):
        guess = np.zeros(self._NbParam, dtype=np.float64)        
        self._Oatm.setParamNoEffect()
        g0 = self._Oatm._Par
        print "g0:",g0
        for idxS in range(self._NbStar):
            guess[idxS*2] = 6100
            guess[idxS*2+1] = 3
        for idx in range(self._NbFlux):
            print self.aIdxParAtm[idx,:]            
            guess[self.aIdxParAtm[idx,:]] = g0
        return guess
    
    def getFirstIdxAtm(self):
        return  2*self._NbStar    
    
    def getTrueParam(self):
        """
        return global true vector parameter
        """
        if self._TrueParam == None:
            par = np.zeros(self._NbParam)
            print "_NbParam",self._NbParam
            for idxS in range(self._NbStar):
                par[idxS*2] = self._Ostar._aParam[idxS, 0]
                par[idxS*2+1] = self._Ostar._aParam[idxS, 1]
            for idx in range(self._NbFlux):
                print self.aIdxParAtm[idx,:]            
                par[self.aIdxParAtm[idx]] = self.getTrueParamAtmIdx(idx)                
#            idx = 2*self._NbStar
#            idxObs = 0
#            for idxN in range(self._NbNight):
#                print idx,  self._SimuParNight[idxN]
#                par[idx:idx+8] = self._SimuParNight[idxN]
#                idx += 8
#                for idxO in range(self._NbObsNight):
#                    print idx
#                    par[idx:idx+2] = self._SimuParObs[idxObs]
#                    idxObs += 1
#                    idx += 2
            self._TrueParam = par
        return np.copy(self._TrueParam)
    
    def plotFlux(self,idx):
        pl.figure()
        pl.plot(self._Oatm._aWL, self.aFlux[idx,:])
        pl.title("Flux simu star/atm %d"%idx)
        pl.grid()
            
    def plotFluxAll(self):
        for idx in range(self._NbFlux):
            pl.figure()
            pl.plot(self._Oatm._aWL, self.aFlux[idx,:])
            pl.title("Flux simu star/atm %d"%idx)
            pl.grid()