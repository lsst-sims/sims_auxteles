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
        # flux measured
        self.aFlux = None
        # constant model vector
        self.aConst = None
        # index relation with a flux measured and atmosphere parameters
        self.aIdxParAtm = None
        # index relation with a flux measured and star parameters
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
    
    Before used called:
     * setAtmModel()
     * setStarTarget()
     
    principal method:
     * readObsNight():  simulated all data observation
    
    """
    def __init__(self, pNbNight=1, pNbObsNight=1):
        ObsSurvey.__class__(self)        
        self._NbNight = pNbNight
        self._NbObsNight = pNbObsNight
        self._NbConst = 4
        self._Oatm = atm.BurkeAtmModel("")
        self._TrueParam = None
        
# PRIVATE
#########
    def _createTrueParAndConst(self):
        # create true param
        self._SimuParNight  = np.zeros((self._NbNight, 8), dtype=np.float32)
        self._SimuParTgray  = np.zeros(self._NbFlux , dtype=np.float32)
        self._SimuParC_H2O  = np.zeros(self._NbSimul , dtype=np.float32)
        self._SimuParObsIdx = np.zeros((self._NbFlux, 3) , dtype=np.int32)
        self._SimuConstObs  = np.zeros((self._NbFlux, self._NbConst) , dtype=np.float32)        
        
    def _randomTrueParAndConst(self):
        self._createTrueParAndConst()
        #######:  _SimuParNight       
        # tau 0
        self._SimuParNight[:,0] = np.random.uniform(0.0, 0.015, self._NbNight)
        # tau 1
        self._SimuParNight[:,1] = np.random.uniform(-0.0005, 0.0005, self._NbNight)
        # tau 2
        self._SimuParNight[:,2] = np.random.uniform(-0.0005, 0.0005, self._NbNight)
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
        self._SimuParTgray = np.random.uniform(0.7, 0.8, self._NbFlux)
        # C_H2O
        self._SimuParC_H2O = np.random.uniform(0.5, 1.2, self._NbSimul)
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

# PUBLIC        
########    


# SETTER              
    def setStarTarget(self, pOstar):
        assert isinstance(pOstar, star.StarTargetSimu)    
        self._Ostar = pOstar
        self._NbStar = self._Ostar._NbStar
        self._NbWL = self._Ostar._NbWL
        self._NbPstar  = self._Ostar._NbPstar
        self._NbSimul = self._NbNight*self._NbObsNight
        self._NbFlux  = self._NbSimul*self._NbStar
                
    def setAtmModel(self, pOatm):
        assert isinstance(pOatm, atm.BurkeAtmModel)                       
        self._Oatm = pOatm
        self._NbPatm = self._Oatm._NbPar
                
    def setAuxTeles(self, pOauxTeles):
        pass
    
# GETTER
        
    def getTrueParamAtmIdx(self, pIdx):
        """
        return true parameter in Burke model format (10 values) for flux pIdx
        """
        par = np.zeros(10, dtype=np.float64)
        par[0] = self._SimuParTgray[pIdx]
        par[1:9] = self._SimuParNight[self._SimuParObsIdx[pIdx,0]]
        par[9] =  par[8]
        par[8] =  par[7]
        par[7] = self._SimuParC_H2O[self._SimuParObsIdx[pIdx,2]]
        return par
        
    def getConst(self, pIdx):
        """
        return observation constante for Burke model format for flux pIdx
        """
        return self._SimuConstObs[pIdx,:]
    
    def getGuessDefault(self):
        """
        return default parameter for all flux
        """
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
        """
        return index of first parameter relative to atmopshere in parameter vector 
        """
        return  2*self._NbStar    
    
    def getTrueParam(self):
        """
        return true parameter vector for all flux
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

        
    def readObsNight(self, pRep):
        """
        fill self.aFlux (array flux)  , here by simulation atmosphere and star
        """
        assert isinstance(self._Oatm, atm.BurkeAtmModel)
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
                    for idxPs in range(self._NbPstar):
                        self.aIdxParStar[idxFlux, idxPs] = idxS*self._NbPstar + idxPs
                        #self.aIdxParStar[idxFlux, 1] = idxS*2 +1                    
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
    
    def addNoisebySNR(self, snr=100, doPlot=False):
        mean = self.aFlux.ravel().mean()
        sigma =  mean/ snr 
        print "mean signal :", mean
        print "sigma for SNR %f: ", sigma
        size = self._NbFlux*self._NbWL
        noise = np.random.normal(0, sigma, size)
        if doPlot:
            pl.figure()
            pl.title("Add noise to spectrum")
            pl.plot(self._Oatm._aWL, self.aFlux[1])
        self.aFlux += noise.reshape(self._NbFlux,self._NbWL)
        if doPlot:           
            pl.plot(self._Oatm._aWL,self.aFlux[1])
            pl.xlabel("Angstrom")
            pl.legend(["no noise","with noise" ])
            pl.grid()
    
    def computeTransTheoIdx(self, idx, param):
        """
        idx   : index flux
        param : global parameter vector
        """        
        self._Oatm.setConstObs(self.aConst[idx])
        parBurke = param[self.aIdxParAtm[idx]]
        self._Oatm.setParam(parBurke)
        trans = self._Oatm.computeAtmTrans()
        return trans
    
    def computeFluxTheoIdx(self, idx, param):
        """
        idx   : index flux
        param : global parameter vector
        """
        parStar = np.array([-3, 0.0, 0.0])               
        # compute model transmission                
        self._Oatm.setConstObs(self.aConst[idx])
        par = param[self.aIdxParAtm[idx]]
        self._Oatm.setParam(par)
        trans = self._Oatm.computeAtmTrans()
        # compute model flux
        parStar[1:] = param[self.aIdxParStar[idx]]
        flux = self._Ostar._oKur.getFluxInterLin(parStar)            
        # model 
        atmStarMod = trans*flux            
        return atmStarMod
             
                
    def computeResidu(self, param):
        """
        return residu vector for parameter vector param
        """
        parStar = np.array([-3, 0.0, 0.0])
        residu = np.copy(self.aFlux)         
        for idx in range(self._NbFlux):           
            # compute model transmission                
            self._Oatm.setConstObs(self.aConst[idx])
            par = param[self.aIdxParAtm[idx]]
            self._Oatm.setParam(par)
            trans = self._Oatm.computeAtmTrans()
            # compute model flux
            parStar[1:] = param[self.aIdxParStar[idx]]
            flux = self._Ostar._oKur.getFluxInterLin(parStar)            
            # model 
            atmStarMod = trans*flux
            residu[idx,:] -= atmStarMod           
        return residu.ravel()
    
# PLOT FUNCTION
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


class ObsSurveySimu02(ObsSurveySimu01):
    """
    simu with random position, atm param, 
    star : only one parameter temperature
    """
    def __init__(self, pNbNight=1, pNbObsNight=1):
        ObsSurveySimu01.__init__(self, pNbNight, pNbObsNight)
        
    def getFirstIdxAtm(self):
        """
        return index of first parameter relative to atmopshere in parameter vector 
        """
        return  self._NbStar
    
    def getGuessDefault(self):
        """
        return default parameter for all flux
        """        
        guess = np.zeros(self._NbParam, dtype=np.float64)        
        self._Oatm.setParamNoEffect()
        g0 = self._Oatm._Par
        print "g0:",g0
        for idxS in range(self._NbStar):
            guess[idxS] = 8000            
        for idx in range(self._NbFlux):
            print self.aIdxParAtm[idx,:]            
            guess[self.aIdxParAtm[idx,:]] = g0
        return guess
    
    def getTrueParam(self):
        """
        return global true vector parameter
        """
        if self._TrueParam == None:
            par = np.zeros(self._NbParam)
            print "_NbParam",self._NbParam
            for idxS in range(self._NbStar):                
                par[idxS] = self._Ostar._aParam[idxS, 1]
            for idx in range(self._NbFlux):
                print self.aIdxParAtm[idx,:]            
                par[self.aIdxParAtm[idx]] = self.getTrueParamAtmIdx(idx)                
            self._TrueParam = par
        return np.copy(self._TrueParam)  
    
              
    def computeFluxTheoIdx(self, idxFlux, param):
        """
        idxFlux : index flux measure
        param   : parameter vector 
        """
        self._Oatm.setConstObs(self.aConst[idxFlux])
        par = param[self.aIdxParAtm[idxFlux]]
        #print "atm par:", par
        self._Oatm.setParam(par)
        trans = self._Oatm.computeAtmTrans()
        # compute model flux
        temp = param[self.aIdxParStar[idxFlux]]
        parStar =  self._Ostar.addMetGra(self.aIdxParStar[idxFlux], temp)
        flux = self._Ostar._oKur.getFluxInterLin(parStar)        
        atmStarMod = trans*flux
        return atmStarMod
   
    def computeResidu(self, param):        
        residu = np.copy(self.aFlux)  
        print "================================"  
        #print param
        for idx in range(self._NbFlux):
            #print idx
            # compute model transmission 
            self._Oatm.setConstObs(self.aConst[idx])
            par = param[self.aIdxParAtm[idx]]
            #print self.aIdxParAtm[idx]
            #print "atm par:", par
            self._Oatm.setParam(par)
            trans = self._Oatm.computeAtmTrans()
            # compute model flux
            temp = param[self.aIdxParStar[idx]]
            parStar =  self._Ostar.addMetGra(self.aIdxParStar[idx], temp)
            flux = self._Ostar._oKur.getFluxInterLin(parStar)
            #if idx < self._NbStar: print "star par:", parStar
            # model 
            atmStarMod = trans*flux
            residu[idx,:] -= atmStarMod        
        return residu.ravel()
    
    def setStarTarget(self, pOstar):
        assert isinstance(pOstar, star.StarTargetSimuAll)    
        self._Ostar = pOstar
        self._NbStar = self._Ostar._NbStar
        self._NbWL = self._Ostar._NbWL
        self._NbPstar  = self._Ostar._NbPstar
        self._NbSimul = self._NbNight*self._NbObsNight
        self._NbFlux  = self._NbSimul*self._NbStar
    
    