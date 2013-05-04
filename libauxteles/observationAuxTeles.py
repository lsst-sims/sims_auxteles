'''
Created on 30 nov. 2012

@author: colley
'''

import burkeAtmModel as atm
import starTargetSimu as star
import pylab as pl
import kurucz as kur
import numpy as np
import AuxSpecGen as aux
import copy as cp
import tools as tl



class ObsSurvey(object):
    def __init__(self):
        # array flux measured
        self.aFlux = None
        # constant model vector for each flux
        self.aConst = None
        # aIdxParAtm : indirection table for atmospheric parameters
        # return for idx flux , idx in global vector parameter of atmospheric associated to input idx flux
        # the order is coherent with convention used by BurkeAtmModel class ie:
        #   Tgray, Tau0, Tau1 , Tau2, alpha, Cmol, C_O3, C_H2O, dC_{H2O}/dEW, dC_{H2O}/dNS 
        self.aIdxParAtm = None
        # aIdxParStar : indirection table for star Kurucz parameters
        # return for idx flux , idx in global vector parameter Kurucz associated to input idx flux
        # the order is coherent with convention used by Kurucz class ie:
        #   metalicity, temperature, gravity        
        self.aIdxParStar= None
        self.NbFlux   = 0
        self._NbNight = 0
        self._NbStar  = 0
        self._NbConst = 0
        self._NbPatm  = 0
        self._NameNightAtm = [r'$\tau_0$',r'$\tau_1$', r'$\tau_2$', r'$\alpha$', '$C_{mol}$', '$C_{O3}$',r'$\partial_NC_{H_2O}$',r'$\partial_EC_{H_2O}$']
        self._NameTgray = '$T_{gray}$'
        self._NameWater = '$C_{H2O}$'
        self._NameStar = []
        
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
        ObsSurvey.__init__(self)        
        self._NbNight = pNbNight
        self._NbPeriodObsByNight = pNbObsNight
        self._NbConst = 4
        self._Oatm = atm.BurkeAtmModel("")
        self._TrueParam = None
       
        
# PRIVATE
#########
    def _createTrueParAndConst(self):
        # create true param
        # parameter night variation[idx flux]: tau0, tau1, tau2, alpha, Cmol,C_O3,dC_H2O/dWE,dC_H2O/dNS
        self._SimuParNight  = np.zeros((self._NbNight, 8), dtype=np.float32)
        # Tgray one by spectrum acquisition
        self._SimuParTgray  = np.zeros(self._NbFlux , dtype=np.float32)
        # C_H2O one by period of observation
        self._SimuParC_H2O  = np.zeros(self._NbPeriodObs , dtype=np.float32)
        # SimuParObsIdx [idx flux] = idx night, idx star, idx period observation
        self._SimuParObsIdx = np.zeros((self._NbFlux, 3) , dtype=np.int32)
        # observation constant [idx flux]= time, alt, az, pressure
        self._SimuConstObs  = np.zeros((self._NbFlux, self._NbConst) , dtype=np.float32)        
        
    def _randomTrueParAndConst(self):
        self._createTrueParAndConst()
        #######:  _SimuParNight       
        # tau 0
        self._SimuParNight[:,0] = np.random.uniform(0.008, 0.06, self._NbNight)
        #self._SimuParNight[:,0] = np.random.uniform(np.sqrt(0.01), np.sqrt(0.05), self._NbNight)
        # tau 1
        self._SimuParNight[:,1] = np.random.uniform(-0.005, 0.005, self._NbNight)
        # tau 2
        self._SimuParNight[:,2] = np.random.uniform(-0.005, 0.05, self._NbNight)
        # alpha
        self._SimuParNight[:,3] = np.random.uniform(-2, -0.05, self._NbNight)
        # Cmol
        self._SimuParNight[:,4] = np.random.uniform(0.5, 1.2, self._NbNight)
        # C_O3
        self._SimuParNight[:,5] = np.random.uniform(0.5, 1.2, self._NbNight)
        # dC_H2O/dWE
        self._SimuParNight[:,6] = np.random.uniform(-0.01, 0.01, self._NbNight)
        # dC_H2O/dNS
        self._SimuParNight[:,7] = np.random.uniform(-0.01, 0.01, self._NbNight)
        ######## _SimuParObs
        # Tgray
        self._SimuParTgray = np.random.uniform(0.6, 1.0, self._NbFlux)
        # C_H2O
        self._SimuParC_H2O = np.random.uniform(0.5, 1.2, self._NbPeriodObs)
        ######## _SimuParObsIdx
        # idx night
        idx = np.arange(self._NbNight)
        self._SimuParObsIdx[:,0] = np.outer(idx, np.ones(self._NbPeriodObsByNight*self._NbStar, dtype=np.int32)).ravel()
        # idx star
        idx = np.arange(self._NbStar)
        self._SimuParObsIdx[:,1] = np.outer( np.ones(self._NbPeriodObs, dtype=np.int32), idx ).ravel()   
        # idx  period observation  
        idx = np.arange(self._NbPeriodObs)
        self._SimuParObsIdx[:,2] = np.outer(idx,  np.ones(self._NbStar, dtype=np.int32) ).ravel()        
        ######## _SimuConstObs
        # ptg, alt , a
        self._SimuConstObs[:,0] = np.random.uniform(np.deg2rad(45), np.deg2rad(90), self._NbFlux)
        self._SimuConstObs[:,1] = np.random.uniform(0, 2*np.pi, self._NbFlux)
        # pressure
        self._SimuConstObs[:,2] = np.random.uniform(750, 820, self._NbFlux)
        #print self._SimuConstObs[:,3]
        # time
        self._SimuConstObs[:,3] = np.arange(self._NbFlux)

# PUBLIC        
########    


# SETTER              
    def setStarTarget(self, pOstar):
        assert isinstance(pOstar, star.StarTargetSimuAll)    
        self._oStarCat = pOstar
        self._NbStar = self._oStarCat._NbStar
        self._NbWL = self._oStarCat._NbWL
        self._NbPstar  = self._oStarCat._NbPstar
        self._NbPeriodObs = self._NbNight*self._NbPeriodObsByNight
        self._NbFlux  = self._NbPeriodObs*self._NbStar
        for idx in range(self._NbStar):
            self._NameStar.append('$t\degree$')   
            self._NameStar.append('$gra$')
                
    def setAtmModel(self, pOatm):
        assert isinstance(pOatm, atm.BurkeAtmModel)                       
        self._Oatm = pOatm
        self._NbPatm = self._Oatm._NbPar
                
    def setAuxTeles(self, pOauxTeles):
        pass
    
# GETTER
    def getWL(self):
        return self._Oatm._aWL
   
    def getNameVecParam(self):
        return self._NameStar+self._NameAtm    
    
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
        return observation constants for Burke model format for flux pIdx
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
        meanTemp = self._oStarCat.getMeanTemp()
        for idxS in range(self._NbStar):
            guess[idxS*2] = (1 + 0.1*(-1)**idxS)*self._oStarCat._aParam[idxS,1]          
            guess[idxS*2+1] = 2.5
        for idx in range(self._NbFlux):
            print self.aIdxParAtm[idx,:]            
            guess[self.aIdxParAtm[idx,:]] = g0
        print guess
        return guess
    
    def getFirstIdxAtm(self):
        """
        return index of first parameter relative to atmopshere in parameter vector 
        """
        return  2*self._NbStar
    
    def getIdxGravity(self):
        return [2*i+1 for i in range(self._NbStar)]
    
    def getIdxTgray(self):
        aIdx = []
        for idx, name in enumerate(self.getNameVecParam()):
            if name == self._NameTgray:
                aIdx.append(idx)
        return aIdx
    
    def getIdxAlpha(self):
        aIdx = []
        for idx, name in enumerate(self.getNameVecParam()):
            if name == r'$\alpha$':
                aIdx.append(idx)
        return aIdx
    
    def getIdxTau0(self):
        aIdx = []
        for idx, name in enumerate(self.getNameVecParam()):
            if name == r'$\tau_0$':
                aIdx.append(idx)
        return aIdx
    
    def getAtmPart(self, pParam):
        """
        get atmosphere part of parameter pParam
        """
        return pParam[self.getFirstIdxAtm():]
    
    def getStarPart(self, pParam):
        """
        get star part of parameter pParam
        """
        return pParam[:self.getFirstIdxAtm()]
    
    def _getxNightParam(self, offset, pParam):
        ret = np.zeros(self._NbNight, dtype=pParam.dtype)
        start = self.getFirstIdxAtm()       
        offsetNight = self.getNbParamAtmByNight()
        for idx in range(self._NbNight):
            ret[idx] = pParam[start + offset + offsetNight*idx]
        return ret
    
    def getTau0Param(self,pParam):
        """
        extract Tau0 parameters from pParam
        """
        return self._getxNightParam(0, pParam)
    
    def getAlphaParam(self,pParam):
        """
        extract Alpha parameters from pParam
        """
        return self._getxNightParam(3, pParam)
    
    def getTrueParam(self):
        """
        return true parameter vector for all flux
        """
        if self._TrueParam == None:
            par = np.zeros(self._NbParam, dtype=np.float64)
            print "_NbParam",self._NbParam
            for idxS in range(self._NbStar):
                par[idxS*2] = self._oStarCat._aParam[idxS, 1]
                par[idxS*2+1] = self._oStarCat._aParam[idxS, 2]
            for idx in range(self._NbFlux):
                print self.aIdxParAtm[idx,:]            
                par[self.aIdxParAtm[idx]] = self.getTrueParamAtmIdx(idx)
            self._TrueParam = par
        return np.copy(self._TrueParam)
        
    def readObsNight(self, pRep):
        """
        fill self.aFlux (array flux)  , here by simulation atmosphere and star
        """
        #assert isinstance(self._Oatm, atm.BurkeAtmModel)
        #
        idxFlux = 0
        # index vector parameters : star param, atm param
        idxPar = self._NbStar*self._NbPstar
        self.aFlux = np.zeros((self._NbFlux, self._NbWL), dtype=np.float64) 
        #  aConst =  _SimuConstObs    
        self.aConst = np.zeros((self._NbFlux, self._NbConst), dtype=np.float32)
        self.aIdxParAtm = np.zeros((self._NbFlux, self._NbPatm), dtype=np.int16)
        self.aIdxParStar = np.zeros((self._NbFlux, self._NbPstar), dtype=np.int16)       
        # random parameters
        ###
        self._randomTrueParAndConst()
        idxParNight = idxPar 
        idxParCH20  = idxParNight + 8 
        idxParTgray = idxParNight + 9  
        self._NameAtm = []        
        for idxN in range(self._NbNight): 
            [self._NameAtm.append(x) for x in self._NameNightAtm]
            for idxO in range(self._NbPeriodObsByNight):
                self._NameAtm.append(self._NameWater)            
                for idxS in range(self._NbStar):
                    self._NameAtm.append(self._NameTgray) 
                    print idxN,idxO,idxS
                    #
                    # compute star flux
                    # 
                    self.aFlux[idxFlux,:] = self._oStarCat.getFluxIdx(self._SimuParObsIdx[idxFlux,1])
                    # compute transmission atm 
                    self._Oatm.setConstObs(self.getConst(idxFlux))
                    self._Oatm.setParam(self.getTrueParamAtmIdx(idxFlux))
                    #                   
                    # * transmission atm
                    #  
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
        self._NbParam = self._NbStar*self._NbPstar + self._NbNight*(self.getNbParamAtmByNight())
    
    def getNbParamAtmByNight(self):
        return 8 + self._NbPeriodObsByNight*(1+self._NbStar)
    
    def addNoisebySNRglobal(self, snr=100, doPlot=False):
        mean = self.aFlux.ravel().mean()
        sigma =  mean/ snr 
        print "mean signal :", mean
        print "sigma for SNR : ", sigma
        size = self._NbFlux*self._NbWL
        noise = np.random.normal(0, sigma, size)
        if doPlot:
            pl.figure()
            pl.title("Add noise to spectrum")
            pl.plot(self.getWL(), self.aFlux[1])
            pl.plot(self.getWL(), self.aFlux[2])
        self.aFlux += noise.reshape(self._NbFlux,self._NbWL)
        if doPlot:           
            pl.plot(self.getWL(),self.aFlux[1])
            pl.plot(self.getWL(),self.aFlux[2])
            pl.xlabel("Angstrom")
            pl.legend(["no noise","with noise","no noise","with noise" ])
            pl.grid()
            
    def addNoisebySNR(self, snr=100, doPlot=[]):
        aFlux = self.aFlux.copy()
        for idxF in range(self._NbFlux):
            flux = self.aFlux[idxF]
            mean = flux.mean()
            sigma =  mean/ snr 
            print "mean signal :", mean
            print "sigma for SNR : ", sigma
            noise = np.random.normal(0, sigma, self._NbWL)
            self.aFlux[idxF] += noise
        if doPlot != []:
            pl.figure()
            pl.title("Add noise to flux SNR %f"%snr)
            for idx in doPlot:
                pl.plot(self.getWL(), self.aFlux[idx],'.')
                pl.plot(self.getWL(), aFlux[idx])
            pl.xlabel("Angstrom")            
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
               
        # compute model transmission                
        self._Oatm.setConstObs(self.aConst[idx])
        par = param[self.aIdxParAtm[idx]]
        self._Oatm.setParam(par)
        trans = self._Oatm.computeAtmTrans()
        # compute model flux
        tempGra = param[self.aIdxParStar[idx]]
        parStar = self._oStarCat.addMet(self._SimuParObsIdx[idx, 1], tempGra)
        flux = self._oStarCat._oKur.getFluxInterLin(parStar)            
        # model 
        atmStarMod = trans*flux            
        return atmStarMod
                                 
# PLOT FUNCTION
    def plotFlux(self,idx):
        pl.figure()
        pl.plot(self.getWL(), self.aFlux[idx,:])
        pl.title("Flux simu star/atm %d"%idx)
        pl.grid()
            
            
    def plotFluxAll(self):
        for idx in range(self._NbFlux):
            pl.figure()
            pl.plot(self.getWL(), self.aFlux[idx,:])
            pl.title("Flux simu star/atm %d"%idx)
            pl.grid()


class ObsSurveySimuTemp(ObsSurveySimu01):
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
        meanTemp = self._oStarCat.getMeanTemp()
        for idxS in range(self._NbStar):
            guess[idxS] = meanTemp
        for idx in range(self._NbFlux):
            print self.aIdxParAtm[idx,:]            
            guess[self.aIdxParAtm[idx,:]] = g0
        return guess
    
    def getTrueParam(self):
        """
        return global true vector parameter
        """
        if self._TrueParam == None:
            par = np.zeros(self._NbParam, dtype=np.float64)
            print "_NbParam",self._NbParam
            for idxS in range(self._NbStar):                
                par[idxS] = self._oStarCat._aParam[idxS, 1]
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
        parStar =  self._oStarCat.addMetGra(self.aIdxParStar[idxFlux], temp)
        flux = self._oStarCat._oKur.getFluxInterLin(parStar)        
        atmStarMod = trans*flux
        return atmStarMod
       
    def setStarTarget(self, pOstar):
        assert isinstance(pOstar, star.StarTargetSimuAll)    
        self._oStarCat = pOstar
        self._NbStar = self._oStarCat._NbStar
        self._NbWL = self._oStarCat._NbWL
        self._NbPstar  = self._oStarCat._NbPstar
        self._NbPeriodObs = self._NbNight*self._NbPeriodObsByNight
        self._NbFlux  = self._NbPeriodObs*self._NbStar
        for idx in range(self._NbStar):
            self._NameStar.append('$t\degree$')   
    

class ObsSurveySimuV2_1(ObsSurveySimu01):
    """
    Simu with :
     * star flux : catalog flux done by Kurucz model
     * star mvt  : random position 
     * atm       : Burke with constraints random parameters 
     * spectro   : used simulator designed by G. Bazin
     
    Principal method:
     * readObsNight() : simulated all data observation    
    """
    
    
    def __init__(self, pNbNight=1, pNbObsNight=1):
        ObsSurveySimu01.__init__(self, pNbNight, pNbObsNight)        
        self._AuxTeles = aux.AuxTeles()
        # 
        self._NbWL = -1       
        self._wlMin = 4500  # [angs]
        self._wlMax = 9500 # [angs]
        self.outpuRes = -1
        
    def setWLInterval(self):
        """
        [angstrom] wave length min max used to fit 
        """
        self._wlMin = wlMin
        self._wlMax = wlMax
        
        
    def getWL(self):
        return self._WL


    def setStarTarget(self, pOstar):  
        ObsSurveySimu01.setStarTarget(self, pOstar)
    
    
    def setAuxTeles(self, pAuxteles):
        """
        """        
        self._AuxTeles = pAuxteles
        assert isinstance(self._AuxTeles, aux.AuxTeles)
                
        
    def readObsNight(self, pRep):
        """
        fill self.aFlux (array flux)  , here by simulation atmosphere and star
        """
        # initialise AuxTeles
        # en nm
        wl = self._oStarCat.getWL()/10.
        # defined nb pixel to satisfy output resolution
        self._AuxTeles.ajustNBPixel(wl)
        self._NbWL = self._AuxTeles.nbpixel
        print "nbpixel def", self._AuxTeles.nbpixel
        # assert isinstance(self._Oatm, atm.BurkeAtmModel)      
        idxFlux = 0
        # index vector parameters : star param, atm param
        idxPar = self._NbStar*self._NbPstar
        self.aFlux = np.zeros((self._NbFlux, self._NbWL), dtype=np.float64) 
        #  aConst =  _SimuConstObs    
        self.aConst = np.zeros((self._NbFlux, self._NbConst), dtype=np.float32)
        self.aIdxParAtm = np.zeros((self._NbFlux, self._NbPatm), dtype=np.int16)
        self.aIdxParStar = np.zeros((self._NbFlux, self._NbPstar), dtype=np.int16)       
        # random parameters
        self._randomTrueParAndConst()
        idxParNight = idxPar 
        idxParCH20  = idxParNight + 8 
        idxParTgray = idxParNight + 9  
        self._NameAtm = []
        listAuxteles = []
        for idx in range(self._NbStar):
            # deepcopy to save all option choiced
            simSpec = cp.deepcopy(self._AuxTeles)
            print "nbpixel ", simSpec.nbpixel
            spectrum = self._oStarCat.getFluxIdx(idx)           
            simSpec.setStarSpectrum(wl, spectrum) 
            listAuxteles.append(simSpec)
            #simSpec.star.plotWL("star %d"%idx)   
#        for idx in range(self._NbStar):
#            listAuxteles[idx].star.plotWL("star %d"%idx)                   
        for idxN in range(self._NbNight): 
            [self._NameAtm.append(x) for x in self._NameNightAtm]
            for idxO in range(self._NbPeriodObsByNight):
                self._NameAtm.append(self._NameWater)            
                for idxS in range(self._NbStar):
                    self._NameAtm.append(self._NameTgray) 
                    print idxN,idxO,idxS
                    # Compute transmission atm 
                    self._Oatm.setConstObs(self.getConst(idxFlux))
                    self._Oatm.setParam(self.getTrueParamAtmIdx(idxFlux))
                    wlAtm = self._Oatm.getWL()
                    trans = self._Oatm.computeAtmTrans()
                    # calibrated star spectrum estimation
                    simSpec = listAuxteles[self._SimuParObsIdx[idxFlux,1]]
                    print "nbpixel loop", simSpec.nbpixel
                    # en nm
                    simSpec.setAtmTrans(wlAtm/10, trans)
                    #simSpec.atms.plot("")
                    wl, flux  = simSpec.getCalibFlux(240)
                    self.aFlux[idxFlux,:] = np.flipud(flux)
                    #self._WL = wl*10                                       
                    # defined const                         
                    self.aConst[idxFlux,:] = self.getConst(idxFlux)
                    # fill table atm parameters
                    self.aIdxParAtm[idxFlux, 0]    = idxParTgray                
                    self.aIdxParAtm[idxFlux, 1:7]  = np.arange(idxParNight, idxParNight+6)                    
                    self.aIdxParAtm[idxFlux, 7]    = idxParCH20
                    self.aIdxParAtm[idxFlux, 8:10] = np.arange(idxParNight+6, idxParNight+8) 
                    # fill table star parameters                   
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
        self._WL = np.flipud(wl) *10   
        self._NbParam = self._NbStar*self._NbPstar + self._NbNight*(self.getNbParamAtmByNight())
        # reduce flux to interval wlMin wlMax choiced
        iMin, iMax = tl.indexInIntervalCheck(self._WL, self._wlMin, self._wlMax)
        self._WL = self._WL[iMin: iMax]
        self._NbWL = len( self._WL )
        self.aFlux = self.aFlux[:,iMin: iMax]
        assert (self._NbWL == self.aFlux.shape[1])
        
