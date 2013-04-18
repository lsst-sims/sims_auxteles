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
    def __init__(self, night, obsByNight, wlMin = 4400, wlMax=10100):
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.restrictToWLinterval(wlMin , wlMax)
        self.oStarCat = star.StarTargetSimuAll(self.oKur, 2)
        self.oAtm = atm.BurkeAtmModel(fileModtran)
        self.oAtm.resample( self.oKur.getWL())
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
        self.oSol.init(self.oObs, self.oKur, self.oAtm)

    def getDefaultAuxTeles(self):
        return aux.AuxTeles()
    

def test_doSimu01():
    oSolSim = SimuVersion2_1(2,2)
    at = oSolSim.getDefaultAuxTeles()
    #at.doPlot = True
    oSolSim.setAuxTeles(at)
    oSolSim.doSimu()
    oSolSim.oObs.plotFluxAll()  
#    oSolSim.oObs.plotFlux(1)
#    oSolSim.oObs.plotFlux(2)
#    oSolSim.oObs.plotFlux(3)
    


#
# MAIN 
#
test_doSimu01()

try:
    pl.show()
    pass

except AttributeError:
    pass
