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


class SimuAtmStarSolve():
    """
    used ObsSurveySimuV2_1 class to simulate spectrum ie: 
     Simu with :
     * star flux : catalog flux done by Kurucz model
     * star mvt  : random position 
     * atm       : Burke with constraints random parameters 
     * spectro   : used simulator designed by G. Bazin    
    """
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.resampleBetween(4000, 10000, 500)
        self.oStarCat = star.StarTargetSimuAll(self.oKur, 2)
        self.oAtm = atm.BurkeAtmModel(fileModtran)
        self.oAtm.resample( self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimu01(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        self.oObs.readObsNight("")
        self.oSol = sol.AtmStarSolver()
        self.oSol.init(self.oObs, self.oKur, self.oAtm)
       
