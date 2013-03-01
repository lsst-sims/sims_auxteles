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

#from matplotlib import rc

#rc('text', usetex=True)
#rc('font', family='serif')



class StarFluxArray(object):
    """
    memory cache star flux 
    """
    def __init__(self, nbStar, oKurucz):
        assert isinstance(oKurucz, kur.Kurucz)
        self._oKur = oKurucz                
        self._nbStar = nbStar        
        self._aParKur = np.zeros((nbStar, 3), dtype=np.float32)
        sizeWl = len(self._oKur.getWL())
        print "sizeWl:", sizeWl
        self._aFlux = np.zeros((nbStar, sizeWl), dtype=np.float32)
        print self._aFlux.shape
                
    def getFlux(self, idx, parKur):
        if not np.array_equal(parKur, self._aParKur[idx]) :
            print "New parameter %d"%idx,parKur
            # comput new flux
            oKur = self._oKur
            assert isinstance(oKur, kur.Kurucz)
            flux = oKur.getFluxInterLin(parKur)
            self._aParKur[idx] = parKur
            print flux.shape
            print self._aFlux.shape
            self._aFlux[idx,:] = flux
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
        print "sizeWl:", sizeWl
        self._aTrans = np.zeros((nbTrans, sizeWl), dtype=np.float64)
        
    def getTrans(self, idx, constObs, par):
        if not np.array_equal(par, self._aParBurke[idx]):
            #print "New parameter %d"%idx,par
            # comput new flux
            oBurkeMod = self._oBurke
            assert isinstance(oBurkeMod, atm.BurkeAtmModel)
            oBurkeMod.setConstObs(constObs)
            oBurkeMod.setParam(par)
            trans = oBurkeMod.computeAtmTrans()
            self._aParBurke[idx] = par
            self._aTrans[idx,:]  = trans
        else:
            #print "AtmTransArray: return memo flux"
            pass
        return self._aTrans[idx]
        

class AtmSolver(object):
    """
    with leastsq() from scipy optimize package
    """
    def __init__(self):
        # pseudo init pour completion EDI eclipse 
        self._oObs = obsAT.ObsSurveySimu02()
        #self._oObs = obsAT.ObsSurvey()
        self._oStar = kur.Kurucz("")
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

    
    def plotCorrelMatrix(self, mat, namePar, pTitle=''):
        #pl.figure()        
        #pl.pcolor(mat)
        #pl.matshow(mat, cmap=pl.cm.gray)        
        im = pl.matshow(self._covar2Correl(mat),vmin=-1, vmax=1)  
        pl.title(pTitle)  
        aIdx = np.arange(len(namePar))
        pl.xticks(aIdx,  namePar)
        for label in im.axes.xaxis.get_ticklabels():
            label.set_rotation(90)
        pl.yticks(aIdx, namePar)
        pl.colorbar()

    def plotErrRel(self, estSol, pTitle =""):
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
        pl.plot(self._oStar.getWL(),  self._oObs.aFlux[idx])
        fluxTheo = self._oObs.computeFluxTheoIdx(idx, param)
        pl.plot(self._oStar.getWL(),  fluxTheo)
        pl.xlabel("Angstrom")
        pl.legend(["raw flux ","theo. flux"])
        pl.grid()
        
    def plotTransTrueEst(self, idxFlux, pTitle=""):
        pl.figure()
        pl.title("Transmission associated flux %d %s"% (idxFlux, pTitle))
        TransTrue = 100*self._oObs.computeTransTheoIdx(idxFlux, self._oObs.getTrueParam())
        TransEst =100* self._oObs.computeTransTheoIdx(idxFlux, self._parEst)                
        pl.plot(self._oAtm.getWL(), TransTrue )       
        pl.plot(self._oAtm.getWL(), TransEst )
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
            errRel = 100*(TransEst - TransTrue)/TransTrue
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
               
    def solve(self, guess = None):
        if guess == None:
            guess = self._oObs.getGuessDefault()        
        #guess += np.random.uniform(0, 0.05, len(guess))*guess
        print guess
        pl.figure()
        leg=[]
        def residus(param):
            obs = self._oObs
            parStar = np.array([-3, 0.0, 0.0])
            residu = np.copy(obs.aFlux)  
            print "================================"  
            for idx in range(obs._NbFlux):
                #print idx
                # compute model transmission                
                self._oAtm.setConstObs(obs.aConst[idx])
                par = param[obs.aIdxParAtm[idx]]
                #print "atm par:", par
                self._oAtm.setParam(par)
                trans = self._oAtm.computeAtmTrans()
                # compute model flux
                parStar[1:] = param[obs.aIdxParStar[idx]]
                flux = self._oStar.getFluxInterLin(parStar)
                if idx <= self._oObs._NbStar: print "star par:", parStar
                # model 
                atmStarMod = trans*flux
                residu[idx,:] -= atmStarMod
                if idx %4 ==1 and False:                    
                    #pl.plot(self._oAtm._aWL,atmStarMod)
                    pl.plot(self._oAtm._aWL,obs.aFlux[idx])                    
                    leg.append("Tgray %.2f ,C_H2O %.2f ,C_mol %.2f"%(par[0],par[7],par[5]))
            #pl.legend(leg)
            #pl.grid()
            #pl.xlabel("Angstrom")
            #pl.ylabel("flux in arbitrary unit")
            #pl.title("Same star spectrum for different atmospheric conditions")
            #print self._oAtm._aWL
            #print atmStarMod.sum(),residu.ravel().sum()
            #pl.show()
            return residu.ravel()
        self._FitRes = spo.leastsq(residus, guess, full_output=True)
        print self._FitRes[0]
        return self._FitRes[0]
    
                
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
       
    def solve3(self, guess=None):
        if guess == None:
            guess = self._oObs.getGuessDefault()
        self._FitRes = spo.leastsq(self._oObs.computeResidu, guess, full_output=True)
        #self._FitRes = spo.minimize(self._oObs.computeResidu, guess, method='powell')
        print self._FitRes[0]
        return self._FitRes[0]


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
    
    def getResidus(self, param):
        if self._memoStarFlux == None:
            self._memoStarFlux =  StarFluxArray(self._oObs._NbStar, self._oStar) 
        if self._memoAtmTrans == None:
            self._memoAtmTrans =  AtmTransArray(self._oObs._NbFlux, self._oAtm) 
        residu = np.copy(self._oObs.aFlux)  
        print "================================"  
        #print param
        for idx in range(self._oObs._NbFlux):
            #print idx
            # compute model transmission 
            constObs = self._oObs.aConst[idx]            
            par = param[self._oObs.aIdxParAtm[idx]]
            trans = self._memoAtmTrans.getTrans(idx, constObs, par)
            # compute model flux
            temp = param[self._oObs.aIdxParStar[idx]]
            parStar = self._oObs._oStarCat.addMetGra(self._oObs.aIdxParStar[idx], temp)
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
        #self._FitRes = spo.minimize(getResiduself._oObs.computeResidu, guess, method='powell')
        print self._FitRes[0]
        self._parEst = self._FitRes[0].copy()
        return self._parEst           
       
        
class AtmStarSolverv(AtmSolver):
    def __init__(self):
        pass
    
  
        