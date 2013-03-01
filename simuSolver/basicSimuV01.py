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


class SimuAtmStarSolve():
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.resampleBetween(4000, 10000, 1000)
        self.oStarCat = star.StarTargetSimu(self.oKur)
        self.oAtm = BurkeAtmModel(fileModtran)
        self.oAtm.resample( self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimu01(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        self.oObs.readObsNight("")
        self.oSol = asSol.AtmStarSolverv()
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
        self.oSol = asSol.AtmSolver()
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
    np.random.seed(1)
    oSim = SimuAtmStarSolve2(2,24)
    #oSim.oObs.addNoisebySNR (400, [0,1,2,3])
    oSim.oObs.addNoisebySNR (400)   
    #guess = guessTrue*(1+np.random.normal(0,0.05,len(guessTrue)))
    guess = oSim.oObs.getGuessDefault()
    #guess =  oSim.oObs.getGuessDefault()
    #oSim.oSol.plotFluxRawTheo(1, guess, 'with guess')
    #oSim.oSol.plotFluxRawTheo(2, guess, 'with guess')
    oSim.oSol.plotErrRel(guess, "guess")
    print oSim.oSol.getChi2(guess)
    oSim.oSol._parEst = oSim.oSol.solveOnlyAtm(guess[4:])
    tempStar = oSim.oStarCat.getAllTemperature()    
    oSim.oSol._parEst = np.concatenate((tempStar, oSim.oSol._parEst))
    oSim.oSol.plotErrRel(oSim.oSol._parEst, 'estimated')
    oSim.oSol.plotCostFuncHistory()
    #oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst, 'with param estimated')
    #oSim.oSol.plotFluxRawTheo(2, oSim.oSol._parEst, 'with param estimated')
#    oSim.plotDistribErrRelAtm(1)
#    oSim.plotDistribErrRelAtm(0)
    #oSim.oSol.plotDistribErrRelAtmAll()
    #oSim.oSol.plotTransTrueEst(1)
    #oSim.oSol.plotTransTrueEst(2)
    #oSim.oSol.plotTransTrueEst(3)
    oSim.oSol.plotCorrelMatrix(oSim.oSol._FitRes[1], oSim.oObs._NameAtm )
    #oSim.oAtm.computeAtmTrans(True)
    oSim.oAtm.printAndPlotBurkeModel()
#simuOnlyTemp()
#simu01()
if __name__ == '__main__':    
    simuOnlyAtm()
    
    try:
        pl.show()
        pass
    
    except AttributeError:
        pass

