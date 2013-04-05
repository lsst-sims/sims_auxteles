'''
Created on 6 dec. 2012

@author: colley
'''
FileKuruczPic = '/home/colley/projet/lsst/stellar_spectra/k93.pic'
fileModtran = '/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt'

import atmStarSolver as sol
import starTargetSimu as star
import observationAuxTeles as obsAT
import kurucz as kur
import pylab as pl
from burkeAtmModel import BurkeAtmModel, BurkeAtmModelTauPos
import numpy as np


class SimuAtmStarSolve():
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.resampleBetween(4000, 10000, 500)
        self.oStarCat = star.StarTargetSimuAll(self.oKur, 2)
        self.oAtm = BurkeAtmModel(fileModtran)
        self.oAtm.resample( self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimu01(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        self.oObs.readObsNight("")
        self.oSol = sol.AtmStarSolver()
        self.oSol.init(self.oObs, self.oKur, self.oAtm)
       
      
       
class SimuAtmStarSolve2(SimuAtmStarSolve):
    """
    used ObsSurveySimu02 and StarTargetSimuAll class
    """
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.resampleBetween(4000, 10000, 1000)
        self.oStarCat = star.StarTargetSimuAll(self.oKur,1)
        self.oAtm = BurkeAtmModel(fileModtran)
        self.oAtm.resample( self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimu02(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        self.oObs.readObsNight("")
        self.oSol = sol.AtmStarSolver()
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
    oSim.plotErrRelAtm(guess, "guess")
    estSol = oSim.oSol.solve(guess)
    oSim.plotErrRelAtm(estSol, 'estimated')
    
    
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
    oSol = sol.AtmStarSolver()
    oSol.init(oObs, oKur, oAtm)
    oSol.solve2Step()
    print oObs.getTrueParam()

def simuOnlyTemp():   
    # create simulation observation
    np.random.seed(16)
    nbNight = 1
    nbPerioByNight = 16
    oSim = SimuAtmStarSolve2(nbNight, nbPerioByNight)  
    snr = 400
    oSim.oObs.addNoisebySNR (snr)   
    guessTrue = oSim.oObs.getTrueParam() 
    #guess = guessTrue*(1+np.random.normal(0, 0.05, len(guessTrue)))
    guess =  oSim.oObs.getGuessDefault()
    for idx in range(oSim.oSol._oObs._NbStar):
        guess[idx] =  guessTrue[idx]*1.2
    print oSim.oSol.getChi2(guess)
    oSim.oSol.plotErrRel(guess, "guess") 
    if oSim.oSol.solveAtmAndTemperature(guess):
        oSim.oSol.plotCostFuncHistory()
        oSim.oSol.plotCorrelMatrix(oSim.oSol._FitRes[1], oSim.oObs.getNameVecParam(), "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(nbNight, nbPerioByNight, snr) )
        oSim.oSol.plotErrRel(oSim.oSol._parEst,"Estimated")
        oSim.oSol.plotDistribErrRelAtmAll()
        oSim.oSol.plotTransTrueEst(0)
        oSim.oSol.plotTransTrueEst(20)
#        oSim.oSol.plotTransTrueEst(100)
#        oSim.oSol.plotTransTrueEst(140)

def simuAtmTempGra():   
    # create simulation observation
    np.random.seed(16)
    nbNight = 1
    nbPerioByNight = 4
    oSim = SimuAtmStarSolve(nbNight, nbPerioByNight)  
    snr = 400000000000000
    oSim.oObs.addNoisebySNR (snr)   
    guessTrue = oSim.oObs.getTrueParam() 
    #guess = guessTrue*(1+np.random.normal(0, 0.05, len(guessTrue)))
    guess =  oSim.oObs.getGuessDefault()
#    for idx in range(oSim.oSol._oObs._NbStar):
#        guess[2*idx] =  guessTrue[idx]*1.2
    print guess
    oSim.oSol.plotErrRel(guess, "guess") 
    if oSim.oSol.solveAtmStarTempGra(guess):
        oSim.oSol.plotCostFuncHistory()
        oSim.oSol.plotCorrelMatrix(oSim.oSol._FitRes[1], oSim.oObs.getNameVecParam(), "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(nbNight, nbPerioByNight, snr) )
        oSim.oSol.plotErrRel(oSim.oSol._parEst,"Estimated")
        oSim.oSol.plotDistribErrRelAtmAll()
        oSim.oSol.plotTransTrueEst(0)
        oSim.oSol.plotTransTrueEst(20)
#        oSim.oSol.plotTransTrueEst(100)
#        oSim.oSol.plotTransTrueEst(140)

def simuAtmTempGraLM():   
    # create simulation observation
    # problem with seed 162, night 1, 24, but ok with tau and alpha bounds !
    np.random.seed(2629)
    nbNight = 1
    nbPerioByNight = 24
    oSim = SimuAtmStarSolve(nbNight, nbPerioByNight)  
    snr = 200
    oSim.oObs.addNoisebySNR (snr)   
    guessTrue = oSim.oObs.getTrueParam() 
    #guess = guessTrue*(1+np.random.normal(0, 0.05, len(guessTrue)))
    guess =  oSim.oObs.getGuessDefault()
#    for idx in range(oSim.oSol._oObs._NbStar):
#        guess[2*idx] =  guessTrue[idx]*1.2
    print guess
    oSim.oSol.plotErrRel(guess, "guess") 
    if oSim.oSol.solveAtmStarTempGraWithBounds(guess):
        oSim.oSol.plotCostFuncHistory()
        oSim.oSol.plotCorrelMatFromlmfit(oSim.oSol._FitRes.params, oSim.oObs.getNameVecParam(), "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(nbNight, nbPerioByNight, snr) )
        oSim.oSol.plotErrRel(oSim.oSol._parEst,"Estimated")
        oSim.oSol.plotDistribErrRelAtmAll()
        oSim.oSol.plotTransTrueEst(0)
        oSim.oSol.plotTransTrueEst(1)
        print oSim.oObs.getTrueParam()
        print oSim.oSol._parEst
        oSim.oSol.plotTrueEstimatePar()
#        oSim.oSol.plotTransTrueEst(100)
#        oSim.oSol.plotTransTrueEst(140)
        oSim.oSol.plotFluxRawTheo(0, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(2, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(3, oSim.oSol._parEst)


def simuOnlyAtm():   
    # create simulation observation
    np.random.seed(104)
    nbNight = 1
    nbPerioByNight = 15
    oSim = SimuAtmStarSolve2(nbNight, nbPerioByNight)
    #oSim.oObs.addNoisebySNR (400, [0,1,2,3])
    snr = 400000000
    oSim.oObs.addNoisebySNR (snr)   
    #guess = guessTrue*(1+np.random.normal(0,0.05,len(guessTrue)))
    guess = oSim.oObs.getGuessDefault()
    #guess =  oSim.oObs.getGuessDefault()
    #oSim.oSol.plotFluxRawTheo(1, guess, 'with guess')
    #oSim.oSol.plotFluxRawTheo(2, guess, 'with guess')
    oSim.oSol.plotErrRelAtm(guess, "guess")
    print oSim.oSol.getChi2(guess)
    if oSim.oSol.solveOnlyAtm(guess[4:]):
        tempStar = oSim.oStarCat.getAllTemperature()    
        oSim.oSol._parEst = np.concatenate((tempStar, oSim.oSol._parEst))
        oSim.oSol.plotErrRelAtm(oSim.oSol._parEst, 'estimated')
        oSim.oSol.plotCostFuncHistory()
        #oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst, 'with param estimated')
        #oSim.oSol.plotFluxRawTheo(2, oSim.oSol._parEst, 'with param estimated')
    #    oSim.plotDistribErrRelAtm(1)
    #    oSim.plotDistribErrRelAtm(0)
        oSim.oSol.plotDistribErrRelAtmAll()
        #oSim.oSol.plotTransTrueEst(1)
        #oSim.oSol.plotTransTrueEst(2)
        #oSim.oSol.plotTransTrueEst(3)
        oSim.oSol.plotCorrelMatFromCovMat(oSim.oSol._FitRes[1], oSim.oObs._NameAtm, "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(nbNight, nbPerioByNight, snr) )
        #oSim.oAtm.computeAtmTrans(True)
        oSim.oAtm.printAndPlotBurkeModel()
        print oSim.oObs.getTau0Param(oSim.oSol._parEst)
        print oSim.oObs.getTau0Param(oSim.oObs.getTrueParam())
        print oSim.oObs.getAlphaParam(oSim.oSol._parEst)
        print oSim.oObs.getAlphaParam(oSim.oObs.getTrueParam())
#simuOnlyTemp()
#simu01()
if __name__ == '__main__':
    #simuOnlyTemp() 
    #simuOnlyAtm()
    #simuAtmTempGra()
    simuAtmTempGraLM()
    
    try:
        pl.show()
        pass
    
    except AttributeError:
        pass

