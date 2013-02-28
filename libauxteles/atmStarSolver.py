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
        

class AtmStarSolverv1(object):
    """
    with leastsq() from scipy optimize package
    """
    def __init__(self):
        # pseudo init pour completion EDI eclipse 
        self._Obs = obsAT.ObsSurveySimu02()
        #self._Obs = obsAT.ObsSurvey()
        self._StarMod = kur.Kurucz("")
        self._AtmMod = atm.BurkeAtmModel("")
        self._CostFunc = []
        
        
    def init(self, oObs, oKur, oAtm):
        """        
        """
        self._setObservation(oObs)
        self._setStarModel(oKur)
        self._setAtmModel(oAtm)
        
    def _setObservation(self, oObs):
        self._Obs = oObs
        assert isinstance(self._Obs, obsAT.ObsSurvey)
        
    def _setStarModel(self, oKur):
        self._StarMod = oKur
        assert isinstance(self._StarMod, kur.Kurucz)
        
    def _setAtmModel(self, oAtm):
        self._AtmMod = oAtm
        assert isinstance(self._AtmMod, atm.BurkeAtmModel)  
                    
    def solve(self, guess = None):
        if guess == None:
            guess = self._Obs.getGuessDefault()        
        #guess += np.random.uniform(0, 0.05, len(guess))*guess
        print guess
        pl.figure()
        leg=[]
        def residus(param):
            obs = self._Obs
            parStar = np.array([-3, 0.0, 0.0])
            residu = np.copy(obs.aFlux)  
            print "================================"  
            for idx in range(obs._NbFlux):
                #print idx
                # compute model transmission                
                self._AtmMod.setConstObs(obs.aConst[idx])
                par = param[obs.aIdxParAtm[idx]]
                #print "atm par:", par
                self._AtmMod.setParam(par)
                trans = self._AtmMod.computeAtmTrans()
                # compute model flux
                parStar[1:] = param[obs.aIdxParStar[idx]]
                flux = self._StarMod.getFluxInterLin(parStar)
                if idx <= self._Obs._NbStar: print "star par:", parStar
                # model 
                atmStarMod = trans*flux
                residu[idx,:] -= atmStarMod
                if idx %4 ==1 and False:                    
                    #pl.plot(self._AtmMod._aWL,atmStarMod)
                    pl.plot(self._AtmMod._aWL,obs.aFlux[idx])                    
                    leg.append("Tgray %.2f ,C_H2O %.2f ,C_mol %.2f"%(par[0],par[7],par[5]))
            #pl.legend(leg)
            #pl.grid()
            #pl.xlabel("Angstrom")
            #pl.ylabel("flux in arbitrary unit")
            #pl.title("Same star spectrum for different atmospheric conditions")
            #print self._AtmMod._aWL
            #print atmStarMod.sum(),residu.ravel().sum()
            #pl.show()
            return residu.ravel()
        self._FitRes = spo.leastsq(residus, guess, full_output=True)
        print self._FitRes[0]
        return self._FitRes[0]
    
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
    
    def getResidus(self, param):
        return self._Obs.computeResidu(param)     
        
    def _solveStarParam(self, guess):
        """
        guess: global guess
        """
        print "=========  solveStarParam"
        def residus(parTempGra):
            obs = self._Obs
            # fix metallicity
            parStar = np.array([-3, 0.0, 0.0])
            residu = np.copy(obs.aFlux) 
            print parTempGra            
            for idx in range(obs._NbFlux):                
                # compute model transmission                
                self._AtmMod.setConstObs(obs.aConst[idx])
                par = guess[obs.aIdxParAtm[idx]]
                #print "atm par:", par
                self._AtmMod.setParam(par)
                trans = self._AtmMod.computeAtmTrans()
                # compute model flux
                parStar[1:] = parTempGra[obs.aIdxParStar[idx]]
                flux = self._StarMod.getFluxInterLin(parStar)
                #print "star par:", parStar
                # model 
                atmStarMod = trans*flux
                residu[idx,:] -= atmStarMod
            print (residu.ravel()**2).sum()
            return residu.ravel()
        guessStar = np.copy(guess[:self._Obs.getFirstIdxAtm()]) 
        fit = spo.leastsq(residus,guessStar, full_output=True)
        print fit[0]
        parEst = np.copy(guess)
        parEst[:self._Obs.getFirstIdxAtm()] = fit[0]
        return parEst
               
    def _solveAtmParam(self, guess):
        """
        guess: global guess
        """
        print "=========  solveAtmParam"    
        def residus(parAtm):
            obs = self._Obs
            # fix metallicity
            parStar = np.array([-3, 0.0, 0.0])
            parGlo = np.copy(guess)
            parGlo[self._Obs.getFirstIdxAtm():] = parAtm
            residu = np.copy(obs.aFlux)             
            for idx in range(obs._NbFlux):
                print idx
                # compute model transmission                
                self._AtmMod.setConstObs(obs.aConst[idx])
                par = parGlo[obs.aIdxParAtm[idx]]
                print "atm par:", par
                self._AtmMod.setParam(par)
                trans = self._AtmMod.computeAtmTrans()
                # compute model flux
                parStar[1:] = parGlo[obs.aIdxParStar[idx]]
                flux = self._StarMod.getFluxInterLin(parStar)
                print "star par:", parStar
                # model 
                atmStarMod = trans*flux
                residu[idx,:] -= atmStarMod
            print atmStarMod.sum(),residu.ravel().sum()
            return residu.ravel()
        guessAtm = np.copy(guess[self._Obs.getFirstIdxAtm():]) 
        fit = spo.leastsq(residus, guessAtm, full_output=True)
        print fit[0]
        parEst = np.copy(guess)
        parEst[self._Obs.getFirstIdxAtm():] = fit[0]
        return parEst
        
    def solve2Step(self):
        """
        step one: temp, gra
        step two: atm. parameters
        """       
        guessTrue = self._Obs.getTrueParam()
        guess = np.copy(guessTrue)
        #guess += np.random.uniform(0, 0.00001, len(guess))*guess
        #guess = self._Obs.getGuessDefault()
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
        
        
    def solve3(self, guess=None):
        if guess == None:
            guess = self._Obs.getGuessDefault()
        self._FitRes = spo.leastsq(self._Obs.computeResidu, guess, full_output=True)
        #self._FitRes = spo.minimize(self._Obs.computeResidu, guess, method='powell')
        print self._FitRes[0]
        return self._FitRes[0]


    def solveOnlyAtm(self, guess=None):
        if guess == None:
            guess = self._Obs.getGuessDefault()
        print "guess:",guess 
        tempStar = self._Obs._Ostar._aParam[:, 1]   
        def getResidu(param):
            #print "getResidu ", param
            # add temp
            TempParam = np.concatenate((tempStar, param))
            residu = self._Obs.computeResidu(TempParam)
            chi2 = (residu.ravel()**2).sum()
            print "chi2: ",chi2
            self._CostFunc.append(chi2)
            return residu
        self._FitRes = spo.leastsq(getResidu, guess, full_output=True)
        #self._FitRes = spo.minimize(getResiduself._Obs.computeResidu, guess, method='powell')
        print self._FitRes[0]
        return self._FitRes[0]            
       
        
    
        