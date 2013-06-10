'''
Created on 5 dec. 2012

@author: colley
'''
import observationAuxTeles as obsAT
import scipy.optimize as spo
import numpy as np
import kurucz as kur
import burkeAtmModel as atm
import pylab as pl
import lmfit as lm
import tools as tl 



def getCorrelMatrixFromlmfit(pPars):
    sizeMat = len(pPars)
    corMat = np.ones((sizeMat, sizeMat), dtype=np.float32)
    for line, par in enumerate(pPars):
        print par
        for col in range(line):
            try:
                corMat[line, col] = corMat[ col, line] = pPars[par].correl['p%d'%col]
            except:
                print "getCorrelMatrixFromlmfit error with ", line
    return corMat
    

class StarFluxArray(object):
    """
    memory cache star flux 
    """
    def __init__(self, nbStar, oKurucz):
        assert isinstance(oKurucz, kur.Kurucz)
        self._oKur = oKurucz                
        self._nbStar = nbStar        
        self._aParKur = np.zeros((nbStar, 3), dtype=np.float64)
        sizeWl = len(self._oKur.getWL())
        print "sizeWl:", sizeWl
        self._aFlux = np.zeros((nbStar, sizeWl), dtype=np.float64)
        print self._aFlux.shape
        self._aCoef = np.ones(nbStar, dtype=np.float64)
        
        
    def setCoefEarth(self, aCoef):
        """
        define coefficient to compute flux observed at earth
        """
        self._aCoef = aCoef
        
        
    def getFlux(self, idx, parKur):
        if not np.array_equal(parKur, self._aParKur[idx]) :
            #print "New temp %d"%idx,parKur
            # comput new flux
            oKur = self._oKur
            assert isinstance(oKur, kur.Kurucz)
            flux = oKur.getFluxInterLin(parKur)
            self._aParKur[idx] = parKur
            #print flux.shape
            #print self._aFlux.shape
            self._aFlux[idx,:] = flux*self._aCoef[idx]
            #print flux
        else:
            #print "StarFluxArray: return memo flux"
            pass
        return self._aFlux[idx]


class AtmTransArray(object):
    """
    memory cache atmosphere transmission
    """    
    def __init__(self, nbTrans, oBurkeMod):
        assert isinstance(oBurkeMod, atm.BurkeAtmModel)
        self._oBurke = oBurkeMod               
        self._nbTrans = nbTrans      
        self._aParBurke = np.zeros((nbTrans, 10), dtype=np.float64)
        sizeWl = len(oBurkeMod._aWL)
        print sizeWl
        self._aTrans = np.zeros((nbTrans, sizeWl), dtype=np.float64)
    
    
    def _compteTrans(self):
        return self._oBurke.computeAtmTrans()
        
    
    def getTrans(self, idx, constObs, par):
        if not np.array_equal(par, self._aParBurke[idx]):
            diff =par- self._aParBurke[idx]
            #print "Atm change :", np.where(diff != 0.0)[0],' flux idx', idx
            # comput new flux
            oBurkeMod = self._oBurke
            assert isinstance(oBurkeMod, atm.BurkeAtmModel)
            #print oBurkeMod
            oBurkeMod.setConstObs(constObs)
            oBurkeMod.setParam(par)
            trans = self._compteTrans()
            self._aParBurke[idx] = par
            self._aTrans[idx,:]  = trans
        else:
            #print "AtmTransArray: return memo flux"
            pass
        return self._aTrans[idx]
        
        
class AtmTransArraySimuSpectro(AtmTransArray):
    """
    memory cache atmosphere transmission with simple simulation spectro (downgrade resolution)
    """    
    def __init__(self, nbTrans, oBurkeMod, res,  pWL):
        AtmTransArray.__init__(self, nbTrans, oBurkeMod)
        sizeWl = len(pWL)
        self._wl = pWL
        self._res  = res 
        self._aTrans = np.zeros((nbTrans, sizeWl), dtype=np.float64)
    
    def _compteTrans(self):
        self._oBurke.computeAtmTrans()
        trans = self._oBurke.downgradeTransAndResample(self._res, self._wl)
        return trans


class AtmStarSolver(object):
    """
    with lmfit optimize package
    """
    def __init__(self):
        # pseudo init pour completion EDI eclipse 
        self._oObs = obsAT.ObsSurveySimu01()
        #self._oObs = obsAT.ObsSurvey()
        #self._oStar = kur.Kurucz("")
        self._oAtm = atm.BurkeAtmModel("")
        self._CostFunc = []
        self._parEst = 0 # parameter estimation
        self._memoStarFlux = None
        self._memoAtmTrans = None
        
    def init(self, oObs, oKur, oAtm):
        """        
        """
        self._setObservation(oObs)
        self._setStarModel(oKur)
        self._setAtmModel(oAtm)            
        
    def _setObservation(self, oObs):
        self._oObs = oObs
        assert isinstance(self._oObs, obsAT.ObsSurvey)
        
    def _setStarModel(self, oKur):
        self._oStar = oKur
        assert isinstance(self._oStar, kur.Kurucz)
        
    def _setAtmModel(self, oAtm):
        self._oAtm = oAtm
        assert isinstance(self._oAtm, atm.BurkeAtmModel)  


#
#
#
    def _covar2Correl(self, mat):
        matTp = np.copy(mat)
        #print matTp
        diag = np.sqrt(matTp.diagonal())
        A = np.outer(np.ones(mat.shape[0]), diag)
        #print A
        res = matTp/(A*A.transpose())
        #print res
        return res
    
    def noAeroErrRel(self, estSol):
        guessTrue = self._oObs.getTrueParam()
        idxStart = self._oObs.getFirstIdxAtm()
        errRel = 100*(estSol[idxStart+4:] - guessTrue[idxStart+4:])/guessTrue[idxStart+4:]
        return errRel

    def aeroErrRel(self, estSol):
        guessTrue = self._oObs.getTrueParam()
        idxStart = self._oObs.getFirstIdxAtm()
        errRel = 100*(estSol[idxStart:idxStart+4] - guessTrue[idxStart:idxStart+4])/guessTrue[idxStart:idxStart+4]
        return errRel
    
    def transmisErrRelAtmAll(self):
        """
        return vector true relative error transmission function estimated
        """
        TruePar = self._oObs.getTrueParam()
        for idxFlux in range(self._oObs._NbFlux):  
            TransTrue = self._oObs.computeTransTheoIdx(idxFlux, TruePar)
            TransEst  = self._oObs.computeTransTheoIdx(idxFlux, self._parEst)
            errRel = 100*(TransEst - TransTrue)/TransTrue
            if  idxFlux == 0:       
                errRelTot = errRel
            else: 
                errRelTot = np.concatenate((errRelTot, errRel))
        return errRelTot
    
    def plotCostFuncHistory(self, pTitle=""):
        pl.figure()
        pl.title("Cost function history %s"%pTitle)
        #print self._CostFunc
        #print np.array(self._CostFunc)
        pl.plot(np.array(self._CostFunc))
        pl.yscale('log')
        pl.xlabel('iteration')
        pl.grid()
                
    def getChi2(self, param):
        return (self.getResidus(param)**2).sum()
    
#
# Plot
#

    
    def plotCorrelMatFromCovMat(self, mat, namePar, pTitle=''):
        self.plotCorrelMat(self._covar2Correl(mat), namePar, pTitle)

    def plotCorrelMatFromlmfit(self, params, namePar, pTitle=''):
        self.plotCorrelMat(getCorrelMatrixFromlmfit(params) , namePar, pTitle)
        
    def plotCorrelMat(self, mat, namePar, pTitle=''):  
        im = pl.matshow(mat,vmin=-1, vmax=1)  
        pl.title(pTitle)  
        aIdx = np.arange(len(namePar))
        pl.xticks(aIdx,  namePar)
        for label in im.axes.xaxis.get_ticklabels():
            label.set_rotation(90)
        pl.yticks(aIdx, namePar)
        pl.colorbar()
        
    def plotErrRel(self,estSol,pTitle=""):
        guessTrue = self._oObs.getTrueParam()
        errRel = 100*(estSol - guessTrue)/guessTrue
        pl.figure()
       
        aIdx = np.arange(len(errRel))
        pl.xticks(aIdx, self._oObs.getNameVecParam(), rotation=90)        
        pl.plot(errRel,"*")      
        pl.grid()
        pl.title('Relative error on parameters '+pTitle)
        pl.ylabel("%")
    
    def plotTrueEstimatePar(self, pLog=False):
        pl.figure()
        pl.title('Parameters true and estimated')
        if pLog:
            pl.semilogy(self._oObs.getTrueParam(),"*")          
            pl.semilogy(self._parEst,".")
        else:
            pl.plot(self._oObs.getTrueParam(),"*")          
            pl.plot(self._parEst,".")                 
        pl.legend(["True ","Estimate"])
        pl.xticks(np.arange(len(self._parEst)), self._oObs.getNameVecParam(), rotation=90)      
        pl.grid()
              
    def plotErrRelAtm(self, estSol, pTitle =""):
        guessTrue = self._oObs.getTrueParam()
        errRel = 100*(estSol - guessTrue)/guessTrue
        pl.figure()
        errRelAtm = self._oObs.getAtmPart(errRel)
        aIdx = np.arange(len(errRelAtm))
        pl.xticks(aIdx, self._oObs._NameAtm, rotation=90)        
        pl.plot(errRelAtm,"*")      
        pl.grid()
        pl.title('Relative error on parameters '+pTitle)
        pl.ylabel("%")
        
    def plotFluxRawTheo(self, idx, param, pTitle=""):
        pl.figure()
        pl.title("Flux index %d %s"% (idx, pTitle))
        pl.plot(self._oObs.getWL(),  self._oObs.aFlux[idx])
        pl.plot(self._oObs.getWL(),  self._oObs.aFlux[idx],"*")
        fluxTheo = self._oObs.computeFluxTheoIdx(idx, param)
        #print len (self._oStar.getWL()), len(fluxTheo)
        pl.plot(self._oObs._Oatm.getWL(),  fluxTheo)
        pl.xlabel("Angstrom")
        pl.legend(["raw flux ","raw flux ","theo. flux"])
        pl.grid()
        
    def plotTransTrueEst(self, idxFlux, pTitle=""):
        pl.figure()
        pl.title("Transmission associated flux %d %s"% (idxFlux, pTitle))
        TransTrue = 100*self._oObs.computeTransTheoIdx(idxFlux, self._oObs.getTrueParam())
        TransEst = 100* self._oObs.computeTransTheoIdx(idxFlux, self._parEst)                
        pl.plot(self._oObs._Oatm.getWL(), TransTrue )       
        pl.plot(self._oObs._Oatm.getWL(), TransEst )
        pl.xlabel("Angstrom")
        pl.ylabel("%")
        pl.legend(["True","Estimated"],loc=2)
        pl.grid()
        
    def plotFluxTheo(self, idx, param):
        pl.figure()
        pl.title("Theoric Flux %d"% idx)       
        fluxTheo = self._oObs.computeFluxTheoIdx(idx, param)
        pl.plot(fluxTheo)       
        pl.grid()
        
    def plotDistribErrRelAtm(self, idxFlux, pTitle=""):
        pl.figure()
        pl.title("distribution relative error transmission for flux %d. %s"% (idxFlux, pTitle))        
        TransTrue = self._oObs.computeTransTheoIdx(idxFlux, self._oObs.getTrueParam())
        TransEst = self._oObs.computeTransTheoIdx(idxFlux, self._parEst)        
        errRel = 100*(TransEst - TransTrue)/TransTrue
        pl.hist(errRel, 50, facecolor='green', log=True)
        pl.xlabel("%")       
        pl.grid()
            
    def plotDistribErrRelAtmAll(self, pTitle=""):
        pl.figure()
        pl.title("distribution relative error transmission for all flux. %s"% (pTitle))
        TruePar = self._oObs.getTrueParam()
        for idxFlux in range(self._oObs._NbFlux):  
            TransTrue = self._oObs.computeTransTheoIdx(idxFlux, TruePar)
            TransEst  = self._oObs.computeTransTheoIdx(idxFlux, self._parEst)
            errRel = 100.0*(TransEst - TransTrue)/TransTrue
            if  idxFlux == 0:       
                errRelTot = errRel
            else: 
                errRelTot = np.concatenate((errRelTot, errRel))
        pl.hist(errRelTot, 50, facecolor='green', log=True)
        pl.xlabel("%")       
        pl.grid()


#
# Slover method
#     
               
    def solveAtmStarTempGra(self, guess = None):
        """
        fit atmosphere parameters, temperature and gravity for star
        """
        if guess == None:
            guess = self._oObs.getGuessDefault()
        self._cptIte =0 
        def getResidu(param):
            print "================== START %d =============="%self._cptIte
            #print "getResidu ", param
            # add temp  
            #print "param: ", param            
            residu = self.getResidusTempGra(param)
            chi2 = (residu.ravel()**2).sum()
            print "chi2: ", chi2
            self._CostFunc.append(chi2)
            self._cptIte += 1
            #if self._cptIte == 20: raise
            return residu        
        res = self._FitRes = spo.leastsq(getResidu, guess, full_output=True)
        self._parEst = self._FitRes[0].copy()
        if res[4] >= 0 and res[4]<=4:
            print "FIT OK :",  res[3]
            return True
        else:
            print "FIT NOK : ",  res[3]
            return False
    
    def solveAtmStarTempGraWithBounds(self, guess = None):
        """
        fit atmosphere parameters, temperature and gravity for stars
        """
        if guess == None:
            guess = self._oObs.getGuessDefault()
        self._cptIte =0 
        def getResidu(params):
            print "================== START %d =============="%self._cptIte
            #print "getResidu ", param
            # add temp  
            #print "param: ", param            
#            v=[]
#            [v.append(params[a].value) for a in params]
            #print   params
            par = np.array([params[a].value for a in params], dtype=np.float64)
            #print par            
            #residu = self.getResidusTempGra(par)
            residu = self._ResidusFunc(par)
            chi2 = (residu**2).sum()
            Tgray = par[aIdxTgray]
            print "Tgray: min %f max %f"%(Tgray.min(), Tgray.max())
            print "chi2: ", chi2
            self._CostFunc.append(chi2)
            self._cptIte += 1
            #if self._cptIte == 20: raise
            return residu
        
        aIdxTgray = self._oObs.getIdxTgray()
        print "Tgray init", guess[aIdxTgray]    
        aIdxGra = self._oObs.getIdxGravity()
        print "Gravity init", guess[aIdxGra]
        aIdxalpha = self._oObs.getIdxAlpha()
        print "alpha init", guess[aIdxalpha]
        aIdxTau0 = self._oObs.getIdxTau0()
        print "Tau0 init", guess[aIdxTau0]
        # Create objets Parameters
        parGuess = lm.Parameters()
        aIdxGra = self._oObs.getIdxGravity()
        print aIdxGra      
        for idx, val in enumerate(guess):
            if idx in aIdxGra:
                parGuess.add("p%d"%idx, val, min= 0.5, max=4.5)
            elif idx in aIdxTgray:
                parGuess.add("p%d"%idx, val, min= 0.4, max=1.00001)
            elif idx in aIdxalpha:
                parGuess.add("p%d"%idx, val, min= -2.1, max=0)
            elif idx in aIdxTau0:
                parGuess.add("p%d"%idx, val, min= 0.005, max=0.1)
            else:
                parGuess.add("p%d"%idx, val)
        print parGuess 
        res = self._FitRes = lm.minimize(getResidu, parGuess, ftol=1e-5, xtol= 1e-5)
        if res.success:
            print "FIT OK :",  res.message
            self._parEst = np.array([res.params[a].value for a in res.params], dtype=np.float64)
            return True
        else:
            print "FIT NOK : ",  res.message, res.lmdif_message
            return False
                
    def _solveStarParam(self, guess):
        """
        guess: global guess
        """
        print "=========  solveStarParam"
        def residus(parTempGra):
            obs = self._oObs
            # fix metallicity
            parStar = np.array([-3, 0.0, 0.0])
            residu = np.copy(obs.aFlux) 
            print parTempGra            
            for idx in range(obs._NbFlux):                
                # compute model transmission                
                self._oAtm.setConstObs(obs.aConst[idx])
                par = guess[obs.aIdxParAtm[idx]]
                #print "atm par:", par
                self._oAtm.setParam(par)
                trans = self._oAtm.computeAtmTrans()
                # compute model flux
                parStar[1:] = parTempGra[obs.aIdxParStar[idx]]
                flux = self._oStar.getFluxInterLin(parStar)
                #print "star par:", parStar
                # model 
                atmStarMod = trans*flux
                residu[idx,:] -= atmStarMod
            print (residu.ravel()**2).sum()
            return residu.ravel()
        guessStar = np.copy(guess[:self._oObs.getFirstIdxAtm()]) 
        fit = spo.leastsq(residus,guessStar, full_output=True)
        print fit[0]
        parEst = np.copy(guess)
        parEst[:self._oObs.getFirstIdxAtm()] = fit[0]
        return parEst
               
    def _solveAtmParam(self, guess):
        """
        guess: global guess
        """
        print "=========  solveAtmParam"    
        def residus(parAtm):
            obs = self._oObs
            # fix metallicity
            parStar = np.array([-3, 0.0, 0.0])
            parGlo = np.copy(guess)
            parGlo[self._oObs.getFirstIdxAtm():] = parAtm
            residu = np.copy(obs.aFlux)             
            for idx in range(obs._NbFlux):
                print idx
                # compute model transmission                
                self._oAtm.setConstObs(obs.aConst[idx])
                par = parGlo[obs.aIdxParAtm[idx]]
                print "atm par:", par
                self._oAtm.setParam(par)
                trans = self._oAtm.computeAtmTrans()
                # compute model flux
                parStar[1:] = parGlo[obs.aIdxParStar[idx]]
                flux = self._oStar.getFluxInterLin(parStar)
                print "star par:", parStar
                # model 
                atmStarMod = trans*flux
                residu[idx,:] -= atmStarMod
            print atmStarMod.sum(),residu.ravel().sum()
            return residu.ravel()
        guessAtm = np.copy(guess[self._oObs.getFirstIdxAtm():]) 
        fit = spo.leastsq(residus, guessAtm, full_output=True)
        print fit[0]
        parEst = np.copy(guess)
        parEst[self._oObs.getFirstIdxAtm():] = fit[0]
        return parEst
        
    def solve2Step(self):
        """
        step one: temp, gra
        step two: atm. parameters
        """       
        guessTrue = self._oObs.getTrueParam()
        guess = np.copy(guessTrue)
        #guess += np.random.uniform(0, 0.00001, len(guess))*guess
        #guess = self._oObs.getGuessDefault()
        est1 = np.copy(guess)
        pl.figure()
        aChi2 = [self.getChi2(est1)]
        print self.getChi2(guessTrue) 
        print self.getChi2(guess)                    
        errRel = 100*(est1-guessTrue)/guessTrue
        pl.plot(errRel)
        leg=["0"]
        for Ite in range(5):            
            est2 = self._solveStarParam(est1)
            print est2
            aChi2.append(self.getChi2(est2))     
            est1 = self._solveAtmParam(est2)
            aChi2.append(self.getChi2(est1))   
            print Ite, (self.getResidus(est1)**2).sum()
            errRel = 100*(est1 - guessTrue)/guessTrue
            pl.plot(errRel,linewidth=2) 
            leg.append("%d"%(Ite+1))
        print est1, guess
        pl.legend(leg)
        print leg
        pl.grid()
        pl.figure()
        pl.semilogy(aChi2)
        pl.grid()
        pl.title("Chi2")
        print self.getChi2(guessTrue)  
        
    def computeResiduOld(self, param):
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
            flux = self._oStarCat._oKur.getFluxInterLin(parStar)            
            # model 
            atmStarMod = trans*flux
            residu[idx,:] -= atmStarMod           
        return residu.ravel()
       
    def solveAtmAndTemperature(self, guess=None):
        """
        used scale factor for star temperature
        """
        if guess == None:
            guess = self._oObs.getGuessDefault()
        self._cptIte = 0        
        def getResidu(param):
            print "================== START %d =============="%self._cptIte
            #print "getResidu ", param
            # add temp  
            #print "param: ", param            
            residu = self.getResidus(param)
            chi2 = (residu.ravel()**2).sum()
            print "chi2: ", chi2
            self._CostFunc.append(chi2)
            self._cptIte += 1
            #if self._cptIte == 20: raise
            return residu        
        res = self._FitRes = spo.leastsq(getResidu, guess, full_output=True)
        self._parEst = self._FitRes[0].copy()
        if res[4] >= 0 and res[4]<=4:
            print "FIT OK :",  res[3]
            return True
        else:
            print "FIT NOK : ",  res[3]
            return False


    def computeResiduSlow(self, param):
        if self._memoStarFlux == None:
            self._memoStarFlux =  StarFluxArray(self._NbStar, self._oStarCat._oKur) 
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
            parStar =  self._oStarCat.addMetGra(self.aIdxParStar[idx], temp)
            flux = self._memoStarFlux.getFlux(self._SimuParObsIdx[idx, 1], parStar)
            #flux = self._oStarCat._oKur.getFluxInterLin(parStar)
            #if idx < self._NbStar: print "star par:", parStar
            # model 
            atmStarMod = trans*flux
            residu[idx,:] -= atmStarMod        
        return residu.ravel()
    
    
    def getResidusTempGra(self, param):
        """
        with temperature star and gravity
        """
        if self._memoStarFlux == None:
            self._memoStarFlux =  StarFluxArray(self._oObs._NbStar, self._oStar) 
        if self._memoAtmTrans == None:
            self._memoAtmTrans =  AtmTransArray(self._oObs._NbFlux, self._oAtm) 
        residu = np.copy(self._oObs.aFlux)  
        #print "================================"  
        #print param
        for idx in range(self._oObs._NbFlux):            
            # compute model transmission 
            constObs = self._oObs.aConst[idx]            
            par = param[self._oObs.aIdxParAtm[idx]]
            #print par          
            trans = self._memoAtmTrans.getTrans(idx, constObs, par)
            # compute model flux
            tempGra = param[self._oObs.aIdxParStar[idx]]            
            parStar = self._oObs._oStarCat.addMet(self._oObs._SimuParObsIdx[idx, 1], tempGra)            
            flux = self._memoStarFlux.getFlux(self._oObs._SimuParObsIdx[idx, 1], parStar)           
            residu[idx,:] -= trans*flux
            if False:
                pl.figure()
                pl.plot(self._oAtm.getWL(), trans*flux) 
                pl.plot(self._oAtm.getWL(), self._oObs.aFlux[idx] )
                pl.legend(["theo","meas"])
                pl.grid() 
                  
        return residu.ravel()

        
    def getResidusTempGraPond(self, param):
        """
        with temperature star and gravity
        """
        if self._memoStarFlux == None:
            self._memoStarFlux =  StarFluxArray(self._oObs._NbStar, self._oStar)
            self._memoStarFlux.setCoefEarth(self._aEarthCoefForFlux)
        if self._memoAtmTrans == None:
            self._memoAtmTrans =  AtmTransArray(self._oObs._NbFlux, self._oAtm) 
        residu = np.copy(self._oObs.aFlux)  
        #print "================================"  
        #print param
        for idx in range(self._oObs._NbFlux):            
            # compute model transmission 
            constObs = self._oObs.aConst[idx]            
            par = param[self._oObs.aIdxParAtm[idx]]
            #print par          
            trans = self._memoAtmTrans.getTrans(idx, constObs, par)
            # compute model flux
            tempGra = param[self._oObs.aIdxParStar[idx]]            
            parStar = self._oObs._oStarCat.addMet(self._oObs._SimuParObsIdx[idx, 1], tempGra)            
            flux = self._memoStarFlux.getFlux(self._oObs._SimuParObsIdx[idx, 1], parStar)   
#            print   "trans", trans.shape      
#            print   "flux",flux.shape      
#            print   "efficiency", self._oObs.getInstruEfficiency().shape   
#            print   "residu",residu[idx,:].shape
            fluxTheo = trans * flux  * self._oObs.getInstruEfficiency()
            residu[idx,:] -= fluxTheo
            #residu[idx,:] /= self._oObs.aFlux.sum()
            print "sum", self._oObs.aFlux.sum()
            if True:
                pl.figure()
                pl.plot(self._oAtm.getWL(), fluxTheo) 
                pl.plot(self._oAtm.getWL(), self._oObs.aFlux[idx] )
                pl.legend(["theo","meas"])
                pl.grid()
            if idx == 3 and True:
                pl.show()
                raise        
            residu = residu/self._aSigma
        return residu.ravel()


    def getResidusWithSpecroSimu(self, param):
        """
        with temperature star and gravity
        """
        if self._memoStarFlux == None:
            self._memoStarFlux =  StarFluxArray(self._oObs._NbStar, self._oStar) 
        if self._memoAtmTrans == None:
            self._memoAtmTrans =  AtmTransArraySimuSpectro(self._oObs._NbFlux, self._oAtm, self._oObs.outpuRes, self._oObs.getWL()) 
        residu = np.copy(self._oObs.aFlux)  
        #print "================================"  
        #print param
        for idx in range(self._oObs._NbFlux):            
            # compute model transmission 
            constObs = self._oObs.aConst[idx]            
            par = param[self._oObs.aIdxParAtm[idx]]
            #print par          
            trans = self._memoAtmTrans.getTrans(idx, constObs, par)            
            # compute model flux
            tempGra = param[self._oObs.aIdxParStar[idx]]            
            parStar = self._oObs._oStarCat.addMet(self._oObs._SimuParObsIdx[idx, 1], tempGra)            
            flux = self._memoStarFlux.getFlux(self._oObs._SimuParObsIdx[idx, 1], parStar)           
            residu[idx,:] -= trans*flux
            if True:
                pl.figure()
                pl.plot(self._oAtm.getWL(), trans*flux) 
                pl.plot(self._oAtm.getWL(), self._oObs.aFlux[idx] )
                pl.legend(["theo","meas"])
                pl.grid() 
                  
        return residu.ravel()

    def setResidusFunc(self, fResidus):
        self._ResidusFunc = fResidus
        
        
    def setEarthCoefForFlux(self, aCoef):
        self._aEarthCoefForFlux =   aCoef
        
    def setSigmaMeasure(self, aSigma):
        self._aSigma =  aSigma/ aSigma.mean()    
        
        
    def getResidus(self, param):
        """
        with only temperature star
        """
        if self._memoStarFlux == None:
            self._memoStarFlux =  StarFluxArray(self._oObs._NbStar, self._oStar) 
        if self._memoAtmTrans == None:
            self._memoAtmTrans =  AtmTransArray(self._oObs._NbFlux, self._oAtm) 
        residu = np.copy(self._oObs.aFlux)  
        #print "================================"  
        #print param
        for idx in range(self._oObs._NbFlux):
            #print idx
            # compute model transmission 
            constObs = self._oObs.aConst[idx]            
            par = param[self._oObs.aIdxParAtm[idx]]            
            trans = self._memoAtmTrans.getTrans(idx, constObs, par)
            # compute model flux
            temp = param[self._oObs.aIdxParStar[idx]]            
            parStar = self._oObs._oStarCat.addMetGra(self._oObs._SimuParObsIdx[idx, 1], temp)
            #print 'temperature:', temp, parStar[1]
            flux = self._memoStarFlux.getFlux(self._oObs._SimuParObsIdx[idx, 1], parStar)
            #flux = self._oStarCat._oKur.getFluxInterLin(parStar)
            #if idx < self._NbStar: print "star par:", parStar
            # model 
            atmStarMod = trans*flux
            residu[idx,:] -= atmStarMod        
        return residu.ravel()

    def solveOnlyAtm(self, guess=None):
        if guess == None:
            guess = self._oObs.getGuessDefault()
        print "guess:",guess 
        tempStar = self._oObs._oStarCat.getAllTemperature()
        def getResidu(param):
            #print "getResidu ", param
            # add temp
            TempParam = np.concatenate((tempStar, param))
            residu = self.getResidus(TempParam)
            chi2 = (residu.ravel()**2).sum()
            print "chi2: ",chi2
            self._CostFunc.append(chi2)
            return residu
        self._FitRes = spo.leastsq(getResidu, guess, full_output=True)
        res = self._FitRes
        #self._FitRes = spo.minimize(getResiduself._oObs.computeResidu, guess, method='powell')
        print self._FitRes[0]
        self._parEst = self._FitRes[0].copy()
        if res[4] >= 0 and res[4]<=4:
            print "FIT OK :",  res[3]
            return True
        else:
            print "FIT NOK : ",  res[3]
            return False
              
        
        