'''
Created on 16 avr. 2013

@author: colley
'''

FileKuruczPic = '/home/colley/projet/lsst/stellar_spectra/k93.pic'
fileModtran = '/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt'

import atmStarSolver as sol
import starTargetSimu as star
import observationAuxTeles as obsAT
import kurucz as kur
import pylab as pl
import burkeAtmModel as atm
import numpy as np
import AuxSpecGen as aux


class SimuVersion2_1():
    """
    used ObsSurveySimuV2_1 class to simulate spectrum ie: 
     Simu with :
     * star flux : catalog flux done by Kurucz model
     * star mvt  : random position 
     * atm       : Burke with constraints random parameters 
     * spectro   : used simulator designed by G. Bazin    
     
    order call :
        * setAuxTeles()
        * doSimu()
        * initSolver()
    """
    def __init__(self, night, obsByNight):
        self.oAtm = atm.BurkeAtmModel(fileModtran)
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.resample(self.oAtm.getWL())
        self.oStarCat = star.StarTargetSimuAll(self.oKur, 2)
        print "size Kurucz ", len(self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimuV2_1(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        
        
    def setAuxTeles(self, pAuxteles=None):
        if pAuxteles==None: 
            pAuxteles = aux.AuxTeles()
        self.oObs.setAuxTeles(pAuxteles)
        
        
    def doSimu(self):        
        self.oObs.readObsNight("")
    
    
    def initSolver(self):
        self.oSol = sol.AtmStarSolver()
        self.oKurSol = kur.Kurucz(FileKuruczPic)
        self.oKurSol.setCoefUnit(1e-8)
        # resample in ccd wl      
        self.oKurSol.resample(self.oObs._WL)
        self.oAtmSol = atm.BurkeAtmModel(fileModtran)
        self.oAtmSol.downgradeTemplateAndResample(self.oObs.outpuRes, self.oObs._WL)      
        self.oSol.init(self.oObs, self.oKurSol, self.oAtmSol)


    def getDefaultAuxTeles(self):
        return aux.AuxTeles()
    

def test_doSimu01():
    np.random.seed(1000)
    oSolSim = SimuVersion2_1(1,1)
    at = oSolSim.getDefaultAuxTeles()
    #at.doPlot = True
    inputRes  = 3200
    outputRes = 300
    at.setResolution(inputRes, outputRes)    
    oSolSim.setAuxTeles(at)
    oSolSim.doSimu()
    oSolSim.initSolver()
    parTrue = oSolSim.oObs.getTrueParam() 
    #oSolSim.oObs.plotFluxAll()
    oSolSim.oSol.plotFluxRawTheo(0, parTrue)
    oSolSim.oSol.plotFluxRawTheo(1, parTrue)
    oSolSim.oSol.plotFluxRawTheo(2, parTrue)
    oSolSim.oSol.plotFluxRawTheo(3, parTrue)
#    oSolSim.oObs.plotFlux(1)
#    oSolSim.oObs.plotFlux(2)
#    oSolSim.oObs.plotFlux(3)
    

def test_residu():
    nbNight = 1
    nbPerioByNight = 2
    oSim = SimuVersion2_1(nbNight, nbPerioByNight )
    at = oSim.getDefaultAuxTeles()
    #at.doPlot = True
    oSim.oObs.outpuRes =  600
    at.setResolution(oSim.oAtm.getResolution(), oSim.oObs.outpuRes)    
    oSim.setAuxTeles(at)
    oSim.doSimu()
    oSim.initSolver()
    guess = oSim.oObs.getTrueParam()
    oSim.oSol.getResidusTempGra(guess)
   
   
def test_SimuAndSolve01():
    nbNight = 1
    nbPerioByNight = 2
    oSim = SimuVersion2_1(nbNight, nbPerioByNight )
    at = oSim.getDefaultAuxTeles()
    #at.doPlot = True
    outputRes = 600
    at.setResolution(oSim.oAtm.getResolution(), outputRes)    
    oSim.setAuxTeles(at)
    oSim.doSimu()
    oSim.initSolver()
    #parTrue = oSim.oObs.getTrueParam() 
    #oSim.oObs.plotFluxAll()
#    oSim.oSol.plotFluxRawTheo(0, parTrue)
#    oSim.oSol.plotFluxRawTheo(1, parTrue)
#    oSim.oSol.plotFluxRawTheo(2, parTrue)
#    oSim.oSol.plotFluxRawTheo(3, parTrue)
    guess = oSim.oObs.getGuessDefault()
    guess = oSim.oObs.getTrueParam()
    if oSim.oSol.solveAtmStarTempGraWithBounds(guess):
        oSim.oSol.plotCostFuncHistory()
        oSim.oSol.plotFluxRawTheo(0, guess)
        oSim.oSol.plotFluxRawTheo(0, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(1, guess)
        oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(2, guess)
        oSim.oSol.plotFluxRawTheo(2, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(3, guess)
        oSim.oSol.plotFluxRawTheo(3, oSim.oSol._parEst)
        #oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst)
        oSim.oSol.plotCorrelMatFromlmfit(oSim.oSol._FitRes.params, oSim.oObs.getNameVecParam(), "Correlation matrix %d night(s) %d obs. period with 4 stars"%(nbNight, nbPerioByNight) )
        oSim.oSol.plotErrRel(oSim.oSol._parEst,"Estimated")
        oSim.oSol.plotDistribErrRelAtmAll()
        oSim.oSol.plotTransTrueEst(0)
        oSim.oSol.plotTransTrueEst(1)
        oSim.oSol.plotTransTrueEst(2)
        oSim.oSol.plotTransTrueEst(3)
    

#
# MAIN 
#
#test_doSimu01()
#test_SimuAndSolve01()
test_residu()

try:
    pl.show()
except AttributeError:
    pass
