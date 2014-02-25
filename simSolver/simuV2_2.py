'''
Created on 30 oct. 2013

@author: colley
'''

import os
import sys

import AuxSpecGen as aux
import atmStarSolver as sol
import burkeAtmModel as atm
import kurucz as kur
import numpy as np
import obsStrategy as obS
import observationAuxTeles as obsAT
import pickle as pk
import matplotlib.pyplot as pl
import simuV2_1 as simu
import starTargetSimu as star
import tools as tl


sys.path.append('../libauxteles')
sys.path.append('../simSpectro')


#S_StarPos = "../data/simuRef/starPos2013-12-21.pkl"
S_StarPos = "../data/simuRef/select4Star_20.pkl"




###############################################################################
class SimuVersion2_2(simu.SimuVersion2_1):
###############################################################################    
    """
    simulation with ObsSurveySimuV2_2 class to simulate spectrum ie: 
     Simu with :
     * star flux : catalog flux done by Kurucz model observed at earth
     * star mvt  : real position with Hipp 
     * atm       : Burke with constraints random parameters 
     * data flux : with spectro simulator designed by G. Bazin
     * solver    : Levenberg-Marquardt with bounds 
     
    used, call with your own parameters method with this order
     * computeExposureTime()
     * doSimuWithAuxTeles()
     * initSimpleModelSpectro()
     * initSolver()
     
    """
        
    def __init__(self):
        super(SimuVersion2_2, self).__init__()
        self.dataStar = obS.DataPositionStarTarget()
            
        
    def setFileStarPos(self, pFileStarPos):
        self.FileStarPos = pFileStarPos
        tmp = obS.DataPositionStarTarget()
        self.dataStar = tmp.read(pFileStarPos)
        oKur = kur.Kurucz(simu.G_FileKuruczPic)       
        self.oStarCat = star.StarTargetSimuAll(oKur)
        #
        # Use star from observation strategy with Hipparcos catalog
        #
        d = self.dataStar
        self.oStarCat.setAllStars(d.Kur[:,0], d.Kur[:,2], d.Kur[:,1], d.Mag, d.Par, d.IdCst)
        self.oStarCat.setMetToZero()
        print self.oStarCat._aParam
     

        
    def computeExposureTime(self, snr, resol,wlmin=None, wlmax =None):
        """
        snr : mean snr expected
        
        description:
          compute for each star in catalog the exposure time to attain input snr
        """        
        self._SNR = snr
        self._Resol = resol*1.0
        self.oStarCat._oKur.restrictToWLinterval(self.oAtm.getWL().min()-100, self.oAtm.getWL().max()+100)        
        self.oStarCat._oKur.setCoefUnitFlux(1e-4, r"$J.m^{-2}.s^{-1}.nm^{-1}$")
        self.oStarCat.updateFlux()
        # define default atmospheric transmission
        self.oAtm.setDefaultConstObs()        
        self.oAtm.setDefaultParam()
        self.oAtm.setParamAerosolGray(1.0, 0.02, 0.0, 0.0, -1)
#        self.oAtm.setTgray(1.0)        
        transAtm = self.oAtm.computeAtmTrans(False)
        print "mean trans amt ", transAtm.mean()
        # init aux teles
        self.oAuxTel.setAtmTrans(self.oAtm.getWL()/10.0, transAtm)
        self.oAuxTel.setStarSpectrum(self.oStarCat._oKur.getWL()/10.0, self.oStarCat.getFluxIdxAtEarth(0))
        self.oAuxTel.setResolutionAndAjustNBPixel(self.oAtm.getResolution(), self._Resol)
        print "nb pixel ", self.oAuxTel.nbpixel
        #self.oStarCat.plotAll()
        lTimeInt = []
        #for idx in range(1):
        for idx in range(self.oStarCat._NbStar):
            #self.oAuxTel.doPlot = True        
            self.oAuxTel.setStarSpectrum(self.oStarCat._oKur.getWL()/10.0, self.oStarCat.getFluxIdxAtEarth(idx))
            self.oAuxTel.getCalibFlux(addNoise=False)
            #self.oAuxTel.plotRawFlux_photonElectron()
            lTimeInt.append( self.oAuxTel.getMeanTimeExposForSNR(self._SNR, wlmin, wlmax) )
            print "time to SNR %f: %f"%(self._SNR, lTimeInt[idx])
        self.TimeInt = np.array(lTimeInt)
        
        
    def doSimuWithAuxTeles(self, wlMin=350, wlMax=1000):
        """
        pFileStarPos : night observation position and star target simulation
        wlMin: [nm]  cut spectro wavelength under wlMin
        wlMax: [nm]  cut spectro wavelength upper wlMax
        """
        self.oObs = obsAT.ObsSurveySimuV2_2(self.dataStar)        
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        # set resolution in and out
        self.oObs.setAuxTeles(self.oAuxTel)
        self.oObs.setWLInterval(wlMin, wlMax)
        self.oObs.setExposureTime(self.TimeInt)
        self.oObs.readObsNight("")


    def doSimuWithFastModel(self, wlMin=350, wlMax=1000):
        # define restricted wl
        self.oObs = obsAT.ObsSurveySimuV2_2(self.dataStar)        
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        # set resolution in and out
        self.oObs.setAuxTeles(self.oAuxTel)
        self.oObs.setWLInterval(wlMin, wlMax)
        self.oObs.setExposureTime(self.TimeInt)
        self.oObs.readObsNightFast(self._Resol)
        
        
################################## END CLASS ##################################        

#
# tests
#

def test_SimuwithFastModelSpectro():
    """
    
    """
    np.random.seed(29)
    oSim = SimuVersion2_2()
    oSim.setFileStarPos(S_StarPos)
    snr = 2000000000000.0
    resolution = 500
    wlMin = 400
    wlMax = 950
    oSim.computeExposureTime( snr, resolution, wlMin, wlMax)  
    oSim.doSimuWithFastModel(wlMin, wlMax)
    oSim.initSolver()
    guess = oSim.oObs.getGuessDefault()    
    if oSim.oSol.solveAtmStarTempGraWithBounds(guess):
        titre = "resolution=%d snr=%d"%(resolution, snr) 
        oSim.oSol.plotCostFuncHistory()
        #oSim.oSol.plotCorrelMatFromlmfit(oSim.oSol._FitRes.params, oSim.oObs.getNameVecParam(), "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(night, obsByNight, snr) )
        oSim.oSol.plotErrRel(guess,"guess")
        oSim.oSol.plotErrRel(oSim.oSol._parEst,"Estimated")
        oSim.oSol.plotDistribErrRelAtmAll(titre)
        oSim.oSol.plotTransTrueEst(0,titre)
        oSim.oSol.plotTransTrueEst(1,titre)        
        oSim.oSol.plotFluxRawTheo(0, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(2, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(3, oSim.oSol._parEst)
        #oSim.oSol.plotFluxRawTheo(36, oSim.oSol._parEst)
        

def simuWithDifferentResolution():
    obsByNight = 20
    wlMin = 400
    wlMax = 950
    snr = 200
    lRes = [50, 100, 200, 300, 400, 500, 600, 700,800]    
    lRes = [ 200, 500]
    lEC = []
    lMean = []
    nbLoop = 7
    nbRes = len(lRes)
    aSTDTol = np.zeros((nbLoop, nbRes, obsByNight*4), dtype=np.float64)
    for Li in range(nbLoop):
        for iRes,resol in enumerate(lRes):
            np.random.seed(260+Li)
            oSim = SimuVersion2_2()
            oSim.setFileStarPos(S_StarPos)
            oSim.computeExposureTime( snr, resol, wlMin, wlMax)  
            oSim.doSimuWithFastModel(wlMin, wlMax)
            oSim.initSolver()
            guess = oSim.oObs.getGuessDefault()
            oSim.oSol.solveAtmStarTempGraWithBounds(guess)
            errRelTot, aSTD = oSim.oSol.transmisErrRelAtmAllAndSTD()
            aSTDTol[Li, iRes, :] = aSTD
            lMean.append( errRelTot.mean() )
            lEC.append( errRelTot.std() )
    aSTDres = np.array([[a.mean(),a.std()] for a in [aSTDTol[:,i,:] for i in range(nbRes)]])
    aMean = np.array(lMean).reshape(nbLoop, len(lRes))
    aEC = np.array(lEC).reshape(nbLoop, len(lRes))
    fOut = "stdERvResSNR%dspectro2-%d.pkl"%(snr,nbLoop)
    f=open(fOut, "wb")
    pk.dump(aEC, f, pk.HIGHEST_PROTOCOL)
    f.close()
    fOut = "meanERvResSNR%dspectro2-%d.pkl"%(snr,nbLoop)
    f=open(fOut, "wb")
    pk.dump(aMean, f, pk.HIGHEST_PROTOCOL)
    f.close()
    pl.figure()
    pl.title("Standard deviation of true relative error atmospheric transmission, SNR ~%.1f"%oSim.oObs.snrMean)
    print aEC.mean(axis=0), aEC.std(axis=0)
    pl.errorbar(lRes, aEC.mean(axis=0), yerr=aEC.std(axis=0), fmt='ro')
    pl.errorbar(lRes, aSTDres[:,0], yerr=aSTDres[:,1], fmt='bo')
    pl.xlim((20, 900))
    pl.ylim((0, 0.4))
    pl.grid()
    pl.ylabel("standard deviation in %")
    pl.xlabel("spectro resolution")   
    pl.figure()
    pl.title("Standard deviation of true relative error atmospheric transmission, SNR ~%.1f"%oSim.oObs.snrMean)
    print aEC.mean(axis=0), aEC.std(axis=0)
    pl.errorbar(lRes, aSTDres[:,0], yerr=aSTDres[:,1], fmt='bo')
    pl.xlim((20, 900))
    pl.ylim((0, 0.4))
    pl.grid()
    pl.ylabel("standard deviation in %")
    pl.xlabel("spectro resolution")   
    pl.figure()
    pl.title("Mean of true relative error atmospheric transmission, SNR ~%.1f"%oSim.oObs.snrMean)
    print aMean.mean(axis=0), aMean.std(axis=0)
    pl.errorbar(lRes, aMean.mean(axis=0), yerr=aMean.std(axis=0), fmt='ro')
    pl.xlim((20, 900))
    pl.grid()
    pl.ylabel("mean relative error in %")
    pl.xlabel("spectro resolution")     
    pl.figure()
    pl.title("true relative error atmospheric transmission, SNR ~%.1f"%oSim.oObs.snrMean)
    pl.errorbar(lRes, aMean.mean(axis=0), yerr=aMean.std(axis=0), fmt='bo',label="mean")
    pl.errorbar(lRes, aSTDres[:,0], yerr=aSTDres[:,1], fmt='ro',label="std dev")
    pl.xlim((20, 900))
    pl.legend(numpoints=1)
    pl.grid()
    pl.ylabel(" %")
    pl.xlabel("spectro resolution")   
    

        
#
# MAIN 
#
        
if __name__ == "__main__":
    #test_SimuwithFastModelSpectro()
    simuWithDifferentResolution()

try:
    pl.show()
except AttributeError:
    pass