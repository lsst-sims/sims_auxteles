'''
Created on 16 avr. 2013

@author: colley
'''
import sys
sys.path.append('../libauxteles')
sys.path.append('../simSpectro')


G_FileKuruczPic = '../data/kurucz93/k93.pic'
G_fileModtran = '../data/modtran/TemplateT04.01_1.txt'


import atmStarSolver as sol
import starTargetSimu as star
import observationAuxTeles as obsAT
import kurucz as kur
import pylab as pl
import burkeAtmModel as atm
import numpy as np
import AuxSpecGen as aux
import pickle as pk


###############################################################################
class SimuVersion2_1(object):
###############################################################################    
    """
    simulation with ObsSurveySimuV2_1 class to simulate spectrum ie: 
     Simu with :
     * star flux : catalog flux done by Kurucz model observed at earth
     * star mvt  : random position 
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
        self.oAuxTel = aux.AuxTeles()
        self.oAtm = atm.BurkeAtmModel(G_fileModtran)
        
        
    def computeExposureTime(self, snr, resol,wlmin=None, wlmax =None):
        """
        snr : mean snr expected
        
        description:
          compute for each star in catalog the exposure time to attain input snr
        """        
        self._SNR = snr
        self._Resol = resol*1.0
        oKur = kur.Kurucz(G_FileKuruczPic)
        oKur.restrictToWLinterval(self.oAtm.getWL().min()-100, self.oAtm.getWL().max()+100)
        #print self.oAtm.getWL()
        oKur.setCoefUnitFlux(1e-4, r"$J.m^{-2}.s^{-1}.nm^{-1}$")
        self.oStarCat = star.StarTargetSimuAll(oKur)
        # define default atmospheric transmission
        self.oAtm.setDefaultConstObs()        
        self.oAtm.setDefaultParam()
        self.oAtm.setParamAerosolGray(1.0, 0.02, 0.0, 0.0, -1)
#        self.oAtm.setTgray(1.0)
        
        transAtm = self.oAtm.computeAtmTrans(False)
        print "mean trans amt ", transAtm.mean()
        # init aux teles
        self.oAuxTel.setAtmTrans(self.oAtm.getWL()/10.0, transAtm)
        self.oAuxTel.setStarSpectrum(oKur.getWL()/10.0, self.oStarCat.getFluxIdxAtEarth(0))
        self.oAuxTel.setResolutionAndAjustNBPixel(self.oAtm.getResolution(), self._Resol)
        #self.oAuxTel.setResolutionAndAjustNBPixel(60, resol)
        print "nb pixel ", self.oAuxTel.nbpixel
        #self.oStarCat.plotAll()
        lTimeInt = []
        #for idx in range(1):
        for idx in range(self.oStarCat._NbStar):
            #self.oAuxTel.doPlot = True        
            self.oAuxTel.setStarSpectrum(oKur.getWL()/10.0, self.oStarCat.getFluxIdxAtEarth(idx))
            self.oAuxTel.getCalibFlux(addNoise=False)
            #self.oAuxTel.plotRawFlux_photonElectron()
            lTimeInt.append( self.oAuxTel.getMeanTimeExposForSNR(self._SNR, wlmin, wlmax) )
            print "time to SNR %f: %f"%(self._SNR, lTimeInt[idx])
        self.TimeInt = np.array(lTimeInt)


    def doSimuWithAuxTeles(self, night=1, obsByNight=24, wlMin=350, wlMax=1000):
        """
        night: number of night observation , 1 advise
        obsByNight : nunber of onservation period by night , max ~ 24 period every 20mn during 8 hours
        wlMin: [nm]  cut spectro wavelength under wlMin
        wlMax: [nm]  cut spectro wavelength upper wlMax
        """
        self.oObs = obsAT.ObsSurveySimuV2_1(night, obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        # set resolution in and out
        self.oObs.setAuxTeles(self.oAuxTel)
        self.oObs.setWLInterval(wlMin, wlMax)
        self.oObs.setExposureTime(self.TimeInt)
        self.oObs.readObsNight("")

    
    def doSimuWithFastModel(self,night=1, obsByNight=24, wlMin=350, wlMax=1000):
        # define restricted wl
        self.oObs = obsAT.ObsSurveySimuV2_1(night, obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        # set resolution in and out
        self.oObs.setAuxTeles(self.oAuxTel)
        self.oObs.setWLInterval(wlMin, wlMax)
        self.oObs.setExposureTime(self.TimeInt)
        self.oObs.readObsNightFast(self._Resol)
        

    def initSimpleModelSpectro(self):
        """
        init kurucz model and atmospheric Burke model
        """
        # define star catalog at wl ccd        
        self.oStarCat._oKur.resample(self.oObs.getWL()*10)
        # recompute flux catalog with the new wavelength
        self.oStarCat._setFlux()      
        # downgrade resolution MODTRAN template for atmospheric model and resample at wl ccd
        self.oAtm.downgradeTemplateAndResample( self._Resol, self.oObs.getWL()*10 )
        
        
    def initSolver(self):
        """
        Init solver and theoric spectro function model FluxTheo
        
        """
        self.oSol = sol.AtmStarSolver()
        self.oSol.init(self.oObs, self.oStarCat._oKur, self.oAtm)
        self.oSol.setEarthCoefForFlux(self.oStarCat._Kcoef)
        # define residus function where simple spectro model is defined, ie :
        #
        # FluxTheo =  FluxStar(k).Cearth.Burke@res=output_AuxTeles_resolution(b).TolalInstruEfficiency
        #
        # where :
        #   * k are Kurucz star parameter
        #   * b are Burke model parameters
        #
        self.oSol.setResidusFunc(self.oSol.getResidusTempGraPond)
        self.oSol.setSigmaMeasure(self.oObs.aSigma)
        
        
    
        
################################## END CLASS ##################################        

#
# tests
#
def test_computeExposureTime():
    oSim =  SimuVersion2_1()
    #oSim.oAuxTel.setDiamTelescope(4)
    snr = 200
    resolution = 500
    oSim.computeExposureTime(snr, resolution)  
          
     
def test_doSimuWithAuxTeles():
    oSim =  SimuVersion2_1()
    snr = 200
    resolution = 200
    oSim.computeExposureTime( snr, resolution)  
    night =  1
    obsByNight = 1
    wlMin=352
    wlMax=957
    oSim.doSimuWithAuxTeles( night, obsByNight, wlMin, wlMax)
    if True:
        pl.figure()
        pl.title("instrument efficiency")
        pl.ylabel("%")
        pl.plot(oSim.oObs.getWL(), oSim.oObs.getInstruEfficiency()*100)
        pl.grid()        
    oSim.oObs.plotWithErrorBar(0)
    oSim.oObs.plotWithErrorBar(1)
    oSim.oObs.plotWithErrorBar(2)
    oSim.oObs.plotWithErrorBar(3)


def test_initSimpleModelSpectro():
    oSim = SimuVersion2_1()
    snr = 200.0
    resolution = 500
    oSim.computeExposureTime( snr, resolution)  
    night =  1
    obsByNight = 1
    wlMin = 352
    wlMax = 957
    oSim.oAuxTel.verbose = 10
    oSim.doSimuWithFastModel( night, obsByNight, wlMin, wlMax)
    
    oSim.oAuxTel.getInstruEfficiency(True)
    oSim.oAtm.setDefaultParam()
    oSim.oAtm.setDefaultConstObs()    
    oSim.oAtm.computeAtmTrans(True)
    oSim.oStarCat.plotAll()
   
   
def test_SimuwithFastModelSpectro():
    np.random.seed(269)
    oSim = SimuVersion2_1()
    snr = 200.0
    resolution = 500
    wlMin = 400
    wlMax = 950
    oSim.computeExposureTime( snr, resolution, wlMin, wlMax)  
    night =  1
    obsByNight = 24
    oSim.doSimuWithFastModel( night, obsByNight, wlMin, wlMax)
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
   
    
def test_initSolver():
    np.random.seed(504)
    oSim = SimuVersion2_1()
    snr = 40.0
    resolution = 1000
    oSim.computeExposureTime( snr, resolution)  
    night =  1
    obsByNight = 16
    wlMin = 400
    wlMax = 950
    oSim.doSimuWithAuxTeles( night, obsByNight, wlMin, wlMax)
    oSim.initSimpleModelSpectro()
    oSim.initSolver() 
    oSim.oAtm.setDefaultParam()
    oSim.oAtm.setDefaultConstObs()    
    oSim.oAtm.computeAtmTrans(True)
    oSim.oStarCat.plotAll()   
    guess = oSim.oObs.getGuessDefault()
    if oSim.oSol.solveAtmStarTempGraWithBounds(guess):
        titre = "resolution=%d snr=%d"%(resolution, snr) 
        oSim.oSol.plotCostFuncHistory()
        oSim.oSol.plotCorrelMatFromlmfit(oSim.oSol._FitRes.params, oSim.oObs.getNameVecParam(), "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(night, obsByNight, snr) )
        oSim.oSol.plotErrRel(guess,"guess")
        oSim.oSol.plotErrRel(oSim.oSol._parEst,"Estimated")
        oSim.oSol.plotDistribErrRelAtmAll(titre)
        oSim.oSol.plotTransTrueEst(0,titre)
        oSim.oSol.plotTransTrueEst(1,titre)        
        oSim.oSol.plotFluxRawTheo(0, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(2, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawThele(3, oSim.oSol._parEst)
        oSim.oSol.plotFluxRawTheo(36, oSim.oSol._parEst)


def test_diffFluxRawTheo():
    np.random.seed(505)
    oSim = SimuVersion2_1()
    snr = 4000000000000000000.0
    resolution = 500
    oSim.computeExposureTime( snr, resolution)  
    night =  1
    obsByNight = 20
    wlMin = 400
    wlMax = 950
    oSim.doSimuWithAuxTeles( night, obsByNight, wlMin, wlMax)
    oSim.initSimpleModelSpectro()
    oSim.initSolver()
    #diff = oSim.oSol.getDiffStarFluxRawTheo()
    diff = oSim.oSol.getDiffSNRStarFluxRawTheo()
    pl.figure()
    pl.plot(oSim.oObs.getWL(), diff[0])
    pl.plot(oSim.oObs.getWL(),diff[1])
    pl.plot(oSim.oObs.getWL(),diff[2])
    pl.plot(oSim.oObs.getWL(),diff[3])
    pl.plot(oSim.oObs.getWL(),diff[4])
    pl.legend(["0","1","2","3","4",])
    pl.grid()
    pl.figure()
    pl.title("Mean")
    pl.plot(oSim.oObs.getWL(),np.mean(diff, axis=0))     
    pl.grid()
    pl.figure()
    pl.title("histogram snr model error, pixel number %d "% oSim.oObs._NbWL)
    pl.hist(np.mean(diff, axis=0), 40, (-500,500),facecolor='green')
    pl.xlabel('SNR')
    pl.figure()
    pl.title("histogram snr model error, pixel number %d "%oSim.oObs._NbWL)
    pl.hist(np.mean(diff, axis=0), 40, (-250,250),facecolor='green')
    pl.xlabel('SNR')
    diff = oSim.oSol.getDiffStarFluxRawTheo()
    pl.figure()
    pl.title("Mean error model for all observation")
    pl.plot(oSim.oObs.getWL(),np.mean(diff, axis=0))
    pl.xlabel('nm')  
    pl.grid()    
    pl.figure()
    pl.title("relative error model for all observation")
    er = diff/oSim.oObs.aFlux
    pl.plot(oSim.oObs.getWL(),100*np.mean(er, axis=0))
    pl.xlabel('nm')  
    pl.ylabel('%')  
    pl.grid()    
    oSim.oSol.plotFluxRawTheo(0, oSim.oObs.getTrueParam())
    oSim.oSol.plotFluxRawTheo(1, oSim.oObs.getTrueParam())
  
#
# mains 
#
def simuWithDifferentResolution():    
    night =  1
    obsByNight = 20
    wlMin = 400
    wlMax = 950
    snr = 200
    lRes = [50, 100, 200, 300, 400, 500, 600, 700,800]
    #lRes = [ 100]
    lEC = []
    nbLoop = 7
    for Li in range(nbLoop):
        for resol in lRes:
            np.random.seed(260+Li)
            oSim = SimuVersion2_1()
            oSim.computeExposureTime( snr, resol, wlMin, wlMax)  
            oSim.doSimuWithFastModel( night, obsByNight, wlMin, wlMax)
            oSim.initSolver()
            guess = oSim.oObs.getGuessDefault()
            oSim.oSol.solveAtmStarTempGraWithBounds(guess)
            errRelTot = oSim.oSol.transmisErrRelAtmAll()
            lEC.append(errRelTot.std())
    aEC = np.array(lEC).reshape(nbLoop, len(lRes))
    fOut = "stdERvResSNRspectro%d.pic"%snr
    f=open(fOut, "wb")
    pk.dump(aEC, f, pk.HIGHEST_PROTOCOL)
    f.close()
    pl.figure()
    pl.title("Standard deviation of true relative error atmospheric transmission, SNR ~%.1f"%oSim.oObs.snrMean)
    print aEC.mean(axis=0), aEC.std(axis=0)
    pl.errorbar(lRes, aEC.mean(axis=0), yerr=aEC.std(axis=0), fmt='ro')
    pl.xlim((20, 900))
    pl.ylim((0, 0.4))
    pl.grid()
    pl.ylabel("standard deviation in %")
    pl.xlabel("spectro resolution")   
    
        
############################################################################     
class SimuVersion2_1old():
############################################################################    
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
        self.oAtm = atm.BurkeAtmModel(G_fileModtran)
        self.oKur = kur.Kurucz(G_FileKuruczPic)
        self.oKur.setCoefUnitFlux(1e-8)
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
        self.oKurSol = kur.Kurucz(G_FileKuruczPic)
        self.oKurSol.setCoefUnitFlux(1e-8)
        # resample in ccd wl      
        self.oKurSol.resample(self.oObs._WL)
        self.oAtmSol = atm.BurkeAtmModel(G_fileModtran)
        self.oAtmSol.downgradeTemplateAndResample(self.oObs.outpuRes, self.oObs._WL)      
        self.oSol.init(self.oObs, self.oKurSol, self.oAtmSol)


    def getDefaultAuxTeles(self):
        return aux.AuxTeles()
    
################################## END CLASS ##################################    

#
# tests
#
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

#
# SimuVersion2_1
#test_computeExposureTime()
#test_doSimuWithAuxTeles()
#test_initSimpleModelSpectro()
#test_initSolver()
#test_diffFluxRawTheo()
test_SimuwithFastModelSpectro()

#simuWithDifferentResolution()

#
# SimuVersion2_1OLD
#
#test_doSimu01()
#test_SimuAndSolve01()
#test_residu()



try:
    pl.show()
except AttributeError:
    pass
