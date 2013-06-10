'''
Created on 6 dec. 2012

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
from burkeAtmModel import BurkeAtmModel
import numpy as np
import tools as tl 
import pickle as pk


###############################################################################
class SimuAtmStarSolve():
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(G_FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.resampleBetween(4000, 10000, 500)
        self.oStarCat = star.StarTargetSimuAll(self.oKur, 2)
        self.oAtm = BurkeAtmModel(G_fileModtran)
        self.oAtm.resample( self.oKur.getWL() )
        self.oObs = obsAT.ObsSurveySimu01(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        self.oObs.readObsNight("")
        self.oSol = sol.AtmStarSolver()
        self.oSol.init(self.oObs, self.oKur, self.oAtm)

#
# tests
#     

#
# mains
# 

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


def simuAtmTempGraLM():   
    # create simulation observation
    # problem with seed 162, night 1, 24, but ok with tau and alpha bounds !
    np.random.seed(2629)
    nbNight = 1
    nbPerioByNight = 24
    oSim = SimuAtmStarSolve(nbNight, nbPerioByNight)  
    snr = 400
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
#        oSim.oSol.plotFluxRawTheo(0, oSim.oSol._parEst)
#        oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst)
#        oSim.oSol.plotFluxRawTheo(2, oSim.oSol._parEst)
#        oSim.oSol.plotFluxRawTheo(3, oSim.oSol._parEst)


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
      
      
###############################################################################      
class SimuAtmStarSolveWithResol():
    """
    used data with given resolution
    """
    def __init__(self, night, obsByNight, resol):
        self.oAtm = BurkeAtmModel(G_fileModtran)
        self.oAtm.downgradeTemplate(resol)
        # create regular grid, number of node agree with resolution and Nyquist criterion  
        newWL = np.linspace(self.oAtm._aWL[0], self.oAtm._aWL[-1], self.oAtm.getNBins())
        iMin, iMax = tl.indexInIntervalCheck(newWL, 4000, 9500)
        newWL = newWL[iMin:iMax]
        print self.oAtm.getNBins()
        print len(newWL)
        self.oAtm.resample( newWL )
        self.oKur = kur.Kurucz(G_FileKuruczPic)
        self.oKur.setCoefUnitFlux(1e-8)
        self.oKur.resample(newWL)
        self.oStarCat = star.StarTargetSimuAll(self.oKur, 2)
        self.oObs = obsAT.ObsSurveySimu01(night, obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        self.oObs.readObsNight("")
        self.oSol = sol.AtmStarSolver()
        self.oSol.init(self.oObs, self.oKur, self.oAtm)
        self.oSol.setResidusFunc(self.oSol.getResidusTempGra)
        
#
# tests
#     

#
# mains
# 
def simuWithResolution():   
    # create simulation observation
    # problem with seed 162, night 1, 24, but ok with tau and alpha bounds !
    np.random.seed(268)
    nbNight = 1
    nbPerioByNight = 24
    resol = 500
    oSim = SimuAtmStarSolveWithResol(nbNight, nbPerioByNight, resol)  
    snr = 200
    oSim.oObs.addNoisebySNR (snr)   
    guessTrue = oSim.oObs.getTrueParam() 
    #guess = guessTrue*(1+np.random.normal(0, 0.05, lsimuWithResolutionen(guessTrue)))
    guess =  oSim.oObs.getGuessDefault()
#    for idx in range(oSim.oSol._oObs._NbStar):
#        guess[2*idx] =  guessTrue[idx]*1.2
    print guess
    oSim.oSol.plotErrRel(guess, "guess")    
    if oSim.oSol.solveAtmStarTempGraWithBounds(guess):
        titre = "resolution=%d"%resol
        oSim.oSol.plotCostFuncHistory()
        oSim.oSol.plotCorrelMatFromlmfit(oSim.oSol._FitRes.params, oSim.oObs.getNameVecParam(), "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(nbNight, nbPerioByNight, snr) )
        oSim.oSol.plotErrRel(oSim.oSol._parEst,"Estimated")
        oSim.oSol.plotDistribErrRelAtmAll(titre)
        oSim.oSol.plotTransTrueEst(0,titre)
        oSim.oSol.plotTransTrueEst(1,titre)
        print oSim.oObs.getTrueParam()
        print oSim.oSol._parEst
        oSim.oSol.plotTrueEstimatePar()
#        oSim.oSol.plotTransTrueEst(100)
#        oSim.oSol.plotTransTrueEst(140)
#        oSim.oSol.plotFluxRawTheo(0, oSim.oSol._parEst)        
        oSim.oSol.plotFluxRawTheo(0, oSim.oSol._parEst,titre)
        oSim.oSol.plotFluxRawTheo(1, oSim.oSol._parEst,titre)
        print oSim.oSol.transmisErrRelAtmAll().std()


def simuWithDifferentResolution():    
    nbNight = 1
    nbPerioByNight = 20
    snr = 300
    lRes = [50, 100, 200, 300, 400,500,600, 700, 800]
    #lRes = [50, 100]
    lEC = []
    nbLoop = 6
    for Li in range(nbLoop):
        for resol in lRes:
            np.random.seed(260+Li)
            oSim = SimuAtmStarSolveWithResol(nbNight, nbPerioByNight, resol)
            oSim.oObs.addNoisebySNR (snr)
            guess =  oSim.oObs.getGuessDefault()
            oSim.oSol.solveAtmStarTempGraWithBounds(guess)
            errRelTot = oSim.oSol.transmisErrRelAtmAll()
            lEC.append(errRelTot.std())
    aEC = np.array(lEC).reshape(nbLoop, len(lRes))
    fOut = "stdERvResSNR%d.pic"%snr
    f=open(fOut, "wb")
    pk.dump(aEC, f, pk.HIGHEST_PROTOCOL)
    f.close()
    pl.figure()
    pl.title("Standard deviation of true relative error atmospheric transmission, SNR %d"%snr)
    print aEC.mean(axis=0), aEC.std(axis=0)
    pl.errorbar(lRes, aEC.mean(axis=0), yerr=aEC.std(axis=0), fmt='ro')
    pl.xlim((20, 900))
    pl.ylim((0, 0.4))
    pl.grid()
    pl.ylabel("standard deviation in %")
    pl.xlabel("spectro resolution")


def test_errobar():
    lRes = [50, 100]
    a = np.random.normal(0,1, 6)
    aEC = a.reshape((3,2))
    pl.figure()
    pl.title("Standard deviation true relative error atmospheric transmission")
    print aEC.mean(axis=0), aEC.std(axis=0)
    pl.errorbar(lRes, aEC.mean(axis=0), yerr=aEC.std(axis=0), fmt='ro')
    pl.xlim((20, 900))
    pl.grid()
    pl.ylabel("standard deviation en %")
    pl.xlabel("Resolution")


     
###############################################################################       
class SimuAtmStarSolve2(SimuAtmStarSolve):
    """
    used ObsSurveySimu02 and StarTargetSimuAll class
    """
    def __init__(self, night, obsByNight):
        self.oKur = kur.Kurucz(G_FileKuruczPic)
        self.oKur.setCoefUnit(1e-8)
        self.oKur.resampleBetween(4000, 10000, 1000)
        self.oStarCat = star.StarTargetSimuAll(self.oKur,1)
        self.oAtm = BurkeAtmModel(G_fileModtran)
        self.oAtm.resample( self.oKur.getWL())
        self.oObs = obsAT.ObsSurveySimu02(night,obsByNight)
        self.oObs.setAtmModel(self.oAtm)
        self.oObs.setStarTarget(self.oStarCat)
        self.oObs.readObsNight("")
        self.oSol = sol.AtmStarSolver()
        self.oSol.init(self.oObs, self.oKur, self.oAtm)


#
# tests
#     

#
# mains
# 
 
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
        #oSim.oSol.plotTransTrueEst(2)simuWithDifferentResolution
        #oSim.oSol.plotTransTrueEst(3)
        oSim.oSol.plotCorrelMatFromCovMat(oSim.oSol._FitRes[1], oSim.oObs._NameAtm, "Correlation matrix %d night(s) %d obs. period with 4 stars, snr=%.1f"%(nbNight, nbPerioByNight, snr) )
        #oSim.oAtm .computeAtmTrans(True)
        oSim.oAtm.printAndPlotBurkeModel()
        print oSim.oObs.getTau0Param(oSim.oSol._parEst)
        print oSim.oObs.getTau0Param(oSim.oObs.getTrueParam())
        print oSim.oObs.getAlphaParam(oSim.oSol._parEst)
        print oSim.oObs.getAlphaParam(oSim.oObs.getTrueParam())

#
# MAIN
#

##################
# CLASS SimuAtmStarSolve
##################
#simuAtmTempGra()
#simuAtmTempGraLM()


##################
# CLASS SimuAtmStarSolve2
##################
#simuOnlyTemp() 
#simuOnlyAtm()


##################
# CLASS SimuAtmStarSolveWithResol
##################
simuWithResolution()
#simuWithDifferentResolution()
#test_errobar()



try:
    pl.show()
    pass

except AttributeError:
    pass

