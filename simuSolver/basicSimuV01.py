'''
Created on 6 dec. 2012

@author: colley
'''
FileKuruczPic = '/home/colley/projet/lsst/stellar_spectra/k93.pic'
fileModtran = '/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt'

import atmStarSolver as asSol
import starTargetSimu as star
import observationAuxTeles as obsAT
import kurucz as kur
import pylab as pl
from burkeAtmModel import BurkeAtmModel
import numpy as np
import pylab as pl


class SimuAtmStarSolve():
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.resampleBetween(4000, 10000, 1000)
        self.oStar = star.StarTargetSimu(self.oKur)
        self.oAtm = BurkeAtmModel(fileModtran)
        self.oAtm.resample( self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimu01(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStar)
        self.oObs.readObsNight("")
        self.oSol = asSol.AtmStarSolverv1()
        self.oSol.init(self.oObs, self.oKur, self.oAtm)
        self._parEst = 0 # parameter estimation

    def plotErrRel(self, estSol, pTitle =""):
        guessTrue = self.oObs.getTrueParam()
        errRel = 100*(estSol - guessTrue)/guessTrue
        pl.figure()
        pl.plot(errRel,"*")
        pl.grid()
        pl.title('Relative error on parameters '+pTitle)
        pl.ylabel("%")
    
    def plotFluxRawTheo(self, idx, param, pTitle=""):
        pl.figure()
        pl.title("Flux index %d %s"% (idx, pTitle))
        pl.plot(self.oKur.getWL(),  self.oObs.aFlux[idx])
        fluxTheo = self.oObs.computeFluxTheoIdx(idx, param)
        pl.plot(self.oKur.getWL(),  fluxTheo)
        pl.xlabel("Angstrom")
        pl.legend(["raw flux ","theo. flux"])
        pl.grid()
        
    def plotTransTrueEst(self, idxFlux, pTitle=""):
        pl.figure()
        pl.title("Transmission associated flux %d %s"% (idxFlux, pTitle))
        TransTrue = 100*self.oObs.computeTransTheoIdx(idxFlux, self.oObs.getTrueParam())
        TransEst =100* self.oObs.computeTransTheoIdx(idxFlux, self._parEst)                
        pl.plot(self.oKur.getWL(), TransTrue )       
        pl.plot(self.oKur.getWL(), TransEst )
        pl.xlabel("Angstrom")
        pl.ylabel("%")
        pl.legend(["True","Estimated"])
        pl.grid()
        
    def plotFluxTheo(self, idx, param):
        pl.figure()
        pl.title("Theoric Flux %d"% idx)       
        fluxTheo = self.oObs.computeFluxTheoIdx(idx, param)
        pl.plot(fluxTheo)       
        pl.grid()
        
    def plotDistribErrRelAtm(self, idxFlux, pTitle=""):
        pl.figure()
        pl.title("distribution relative error transmission for flux %d. %s"% (idxFlux, pTitle))        
        TransTrue = self.oObs.computeTransTheoIdx(idxFlux, self.oObs.getTrueParam())
        TransEst = self.oObs.computeTransTheoIdx(idxFlux, self._parEst)        
        errRel = 100*(TransEst - TransTrue)/TransTrue
        pl.hist(errRel, 50, facecolor='green', log=True)
        pl.xlabel("%")       
        pl.grid()
            
    def plotDistribErrRelAtmAll(self, pTitle=""):
        pl.figure()
        pl.title("distribution relative error transmission for all flux. %s"% (pTitle))
        TruePar = self.oObs.getTrueParam()
        for idxFlux in range(self.oObs._NbFlux):  
            TransTrue = self.oObs.computeTransTheoIdx(idxFlux, TruePar)
            TransEst  = self.oObs.computeTransTheoIdx(idxFlux, self._parEst)
            errRel = 100*(TransEst - TransTrue)/TransTrue
            if  idxFlux == 0:       
                errRelTot = errRel
            else: 
                errRelTot = np.concatenate((errRelTot, errRel))
        pl.hist(errRelTot, 50, facecolor='green', log=True)
        pl.xlabel("%")       
        pl.grid()
       
       
class SimuAtmStarSolve2(SimuAtmStarSolve):
    """
    used ObsSurveySimu02 and StarTargetSimuAll class
    """
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.resampleBetween(4000, 10000, 1000)
        self.oStar = star.StarTargetSimuAll(self.oKur,1)
        self.oAtm = BurkeAtmModel(fileModtran)
        self.oAtm.resample( self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimu02(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStar)
        self.oObs.readObsNight("")
        self.oSol = asSol.AtmStarSolverv1()
        self.oSol.init(self.oObs, self.oKur, self.oAtm)



def simu01():    
    oSim = SimuAtmStarSolve(2,6)
    guess = oSim.oObs.getTrueParam() 
    print oSim.oSol.getChi2(guess)
#    oSim.plotFluxRawTheo(0, guess)
#    oSim.plotFluxRawTheo(1, guess)
#    oSim.plotFluxRawTheo(2, guess)
#    oSim.plotFluxRawTheo(3, guess)   
    guess = oSim.oObs.getTrueParam() 
    guess += np.random.normal(0,0.1,len(guess))*guess    
#    oSim.plotFluxRawTheo(0, guess)
#    oSim.plotFluxRawTheo(1, guess)
#    oSim.plotFluxRawTheo(2, guess)
#    oSim.plotFluxRawTheo(3, guess)    
#    #guess = oSim.oObs.getGuessDefault()    
    oSim.plotErrRel(guess, "guess")
    estSol = oSim.oSol.solve(guess)
    oSim.plotErrRel(estSol, 'estimated')
    
    
def simu02():   
    # create simulation observation
    oKur = kur.Kurucz(FileKuruczPic)
    oKur.resampleBetween(4000, 10000, 1000)
    oStar = star.StarTargetSimu(oKur)
    oAtm = BurkeAtmModel(fileModtran)
    oAtm.resample( oKur.getWL())
    oObs = obsAT.ObsSurveySimu01(2,4)
    oObs.setAtmModel(oAtm)
    oObs.setStarTarget(oStar)
    oObs.readObsNight("")
    # create object solver
    oSol = asSol.AtmStarSolverv1()
    oSol.init(oObs, oKur, oAtm)
    oSol.solve2Step()
    print oObs.getTrueParam()

def simuOnlyTemp():   
    # create simulation observation
    np.random.seed(100)
    oSim = SimuAtmStarSolve2(2,4)  
    guessTrue = oSim.oObs.getTrueParam() 
    guess = guessTrue*(1+np.random.normal(0, 0.05, len(guessTrue)))
    #guess =  oSim.oObs.getGuessDefault()
    oSim.plotFluxRawTheo(0, guess, 'with guess')
    oSim.plotFluxRawTheo(6, guess, 'with guess')
    oSim.plotErrRel(guess, "guess")
    print oSim.oSol.getChi2(guess)
    estSol = oSim.oSol.solve3(guess)
    print guessTrue
    oSim.plotErrRel(estSol, 'estimated')
    oSim.plotFluxRawTheo(0, estSol, 'with param estimated')
    oSim.plotFluxRawTheo(6, estSol, 'with param estimated')
    
    
def simuOnlyAtm():   
    # create simulation observation
    np.random.seed(103)
    oSim = SimuAtmStarSolve2(1,24)
    oSim.oObs.addNoisebySNR (500, True)
    guessTrue = oSim.oObs.getTrueParam() 
    guess = guessTrue*(1+np.random.normal(0,0.05,len(guessTrue)))
    guess = oSim.oObs.getGuessDefault()
    #guess =  oSim.oObs.getGuessDefault()
    oSim.plotFluxRawTheo(0, guess, 'with guess')
    oSim.plotFluxRawTheo(2, guess, 'with guess')
    oSim.plotErrRel(guess, "guess")
    print oSim.oSol.getChi2(guess)
    oSim._parEst = oSim.oSol.solveOnlyAtm(guess[4:])
    tempStar = np.array([ 7700. , 6440. , 5150. , 9520.])
    oSim._parEst = np.concatenate((tempStar, oSim._parEst))
    oSim.plotErrRel(oSim._parEst, 'estimated')
    oSim.oSol.plotCostFuncHistory()
    oSim.plotFluxRawTheo(0, oSim._parEst, 'with param estimated')
    oSim.plotFluxRawTheo(2, oSim._parEst, 'with param estimated')
    oSim.plotDistribErrRelAtm(1)
    oSim.plotDistribErrRelAtm(0)
    oSim.plotDistribErrRelAtmAll()
    oSim.plotTransTrueEst(1)
    oSim.plotTransTrueEst(2)
    oSim.plotTransTrueEst(3)
    
#simuOnlyTemp()
#simu01()
simuOnlyAtm()

try:
    #pl.show()
    pass

except AttributeError:
    pass

