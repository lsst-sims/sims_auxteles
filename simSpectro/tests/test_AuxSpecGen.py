'''
Created on 2 avr. 2013

@author: colley
'''
import unittest

from AuxSpecGen import *
import kurucz as kur


S_DirExple = getModuleDirectory()+"/example/"
S_StarFile = S_DirExple+'alpha_lyr_stis_004nm.txt'
S_AtmFile = S_DirExple+'Modtran_10.dat'
S_AtmFile2 = S_DirExple+'Modtran_50.dat'
S_AtmFile3 = S_DirExple+'Modtran_100.dat'
G_FileKuruczPic = '../../data/kurucz93/k93.pic'
S_fileModtran = '../../data/modtran/TemplateT04.01_1.txt'



def test_getRawFlux_Wmnm():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile)   
    simSpec.getCalibFlux(addNoise=False)
    if True:
        pl.figure()
        pl.title("RawFlux")
        pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')        
        pl.plot(simSpec.getWLccd(), simSpec.getRawFlux_Wmnm())       
        pl.grid()
    
    
def test_getInstruEfficiency():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile, 1e-6)   
    simSpec.getCalibFlux(addNoise=False)
    if True:
        pl.figure()
        pl.title("instrumental efficiency")
        pl.ylabel("%")        
        pl.plot(simSpec.getWLccd(), simSpec.getInstruEfficiency()*100)       
        pl.grid()
        
        
def test_FluxInstruCor():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile, 1e-6)   
    simSpec.getCalibFlux(addNoise=False)
    if True:
        pl.figure()
        pl.title("Flux instru correc")
        pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')        
        pl.plot(simSpec.getWLccd(), simSpec.getRawFlux_Wmnm()/simSpec.getInstruEfficiency())       
        pl.grid()
                
def test_FluxInstruAtmCor():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile, 1e-6)   
    simSpec.getCalibFlux(addNoise=False)
    atmTr1 = simSpec.atms.getTransAtwlDowngrade(simSpec.getWLccd(), simSpec.inputres, simSpec.outputres)
    atmTr2 = simSpec.atms.getTransAtwl(simSpec.getWLccd())
    if True:
        pl.figure()
        pl.title("Flux instru and atm correc")
        pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')        
        pl.plot(simSpec.getWLccd(), simSpec.getRawFlux_Wmnm()/simSpec.getInstruEfficiency()/atmTr1)       
        pl.plot(simSpec.getWLccd(), simSpec.getRawFlux_Wmnm()/simSpec.getInstruEfficiency()/atmTr2)     
        pl.legend(["smooth Atm"," raw Atm"])  
        pl.grid()

def test_comparFluxInstruAtmCor():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile)   
    simSpec.getCalibFlux(addNoise=False)
    atmTr1 = simSpec.atms.getTransAtwlDowngrade(simSpec.getWLccd(), simSpec.inputres, simSpec.outputres)
    star = tl.interpolBSpline(simSpec.star.wl[::-1], simSpec.star.dEdl[::-1], simSpec.getWLccd())
    starEst= simSpec.getRawFlux_Wmnm()
    simSpec.getCalibFlux()
    starEstNoise= simSpec.getRawFlux_Wmnm()
    if True:
        pl.figure()
        pl.title("Flux instru correc")
        pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')        
        pl.plot(simSpec.getWLccd(), starEst/simSpec.getInstruEfficiency())       
        pl.plot(simSpec.getWLccd(), starEstNoise/simSpec.getInstruEfficiency())       
        pl.plot(simSpec.getWLccd(), star*atmTr1,'.-')     
        pl.legend(["star*atm est ","star*atm est noise","star*atm raw"])  
        pl.grid()
        
def test_noisePhot():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile, 1e-4)   
    simSpec.computePhotonElec()
    if True:
        pl.figure()
        pl.title("Noise photon")        
        pl.plot(simSpec.getWLccd(),simSpec.getSigmaPhotonNoise())               
        pl.grid()
        pl.figure()
        pl.title("SNR")        
        pl.plot(simSpec.getWLccd()[50:-50],simSpec.getRawFlux_photonElectron()[50:-50]/simSpec.getSigmaPhotonNoise()[50:-50])               
        pl.grid()


def test_getMeanTimeExposForSNR():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile, 1e-4)
    simSpec.getCalibFlux(addNoise=False) 
    print "t snr=200: ",simSpec.getMeanTimeExposForSNR(200)
    print "t snr=300: ",simSpec.getMeanTimeExposForSNR(300)
    print "t snr=400: ",simSpec.getMeanTimeExposForSNR(400)


def test_getFluxWithStarLowResolution():
    resol = 500
    oAtm = atm.BurkeAtmModel(S_fileModtran)
    oKur = kur.Kurucz(G_FileKuruczPic)
    oKur.restrictToWLinterval(oAtm.getWL().min()-100, oAtm.getWL().max()+100) 
    oKur.setCoefUnitFlux(1e-4, r"$J.m^{-2}.s^{-1}.nm^{-1}$")
    oAtm.setDefaultParam()
    oAtm.setTgray(1.0)
    oAtm.setDefaultConstObs()
    transAtm = oAtm.computeAtmTrans(True)
    simSpec = AuxTeles()
    simSpec.setAtmTrans(oAtm.getWL()/10.0, transAtm)
    densityFlux = tl.coefKuruczEarth(7000, 7, 300)*oKur.getFluxInterLin(np.array([0, 7000, 2]))
    simSpec.setStarSpectrum(oKur.getWL()/10.0,densityFlux)
    simSpec.setResolutionAndAjustNBPixel(oAtm.getResolution(), resol)
    simSpec.doPlot = True
    simSpec.getFluxWithStarLowResolution(240, False)
      
      
def snr_estimation():
    simSpec = AuxTeles()
    #simSpec.doPlot = True
    simSpec.setAtmFile(S_AtmFile3)
    simSpec.setStarFile(S_StarFile, 1e-6)    
    #simSpec.ajustNBPixel()
    ret = simSpec.getCalibFlux()
    simSpec.star.plotNbPhotons("total photon-electron",  simSpec.star.elecccd)
    spectre =  ret[1].copy()[10:]
    print spectre.max()
    retnn = simSpec.getCalibFlux(addNoise=False)
    spectreNoNoise = ret[1].copy()[10:]
    print spectreNoNoise.max()
    diff = spectreNoNoise - spectre
    print "mean std min max:"
    print diff.mean(), diff.std(), diff.min(), diff.max()
    if True:
        pl.figure()
        pl.title("snr_estimation")
        pl.plot(ret[0][10:], spectre)
        pl.plot(retnn[0][10:], spectreNoNoise)
        pl.grid()

    
    
#snr_estimation()
#test_getRawFlux_Wmnm()
#test_getInstruEfficiency()
#test_FluxInstruCor()
#test_FluxInstruAtmCor()
#test_comparFluxInstruAtmCor()
#test_noisePhot()
#test_getMeanTimeExposForSNR()
test_getFluxWithStarLowResolution()
FlagPlot= True

class Test(unittest.TestCase):
    def plotSimProd(self,xNew, estInter,  prodFluxAtm):
        pl.figure()                    
        pl.plot(xNew, prodFluxAtm)
        pl.plot(xNew, estInter)
        pl.ylim([0, 6.0e-9])
        pl.grid()
        pl.title("comparaison flux*Atm and calib star Flux without atm")
        pl.legend(["flux*Atm ", "sim inter"])
            
            
    def notest_SimAcquSpectro01(self):
        simSpec = AuxTeles()
        simSpec.setAtmFile(S_AtmFile3)
        simSpec.setStarFile(S_StarFile)
        xNew = simSpec.atms.wl[::-1][:-250]
        fluxInter = tl.interpolLinear(simSpec.star.wl[::-1], simSpec.star.dEdl[::-1], xNew)
        prodFluxAtm = fluxInter*simSpec.atms.tr[::-1][:-250]
        simSpec.doPlot = False      
        ret = simSpec.getCalibFlux()
        estInter = tl.interpolLinear(ret[0][::-1], ret[1][::-1], xNew)
        if FlagPlot: 
            pl.figure()            
            pl.plot(ret[0], ret[1])
            pl.plot(xNew, prodFluxAtm)
            pl.plot(xNew, estInter)
            pl.ylim([0, 6.0e-9])
            pl.grid()
            pl.title("comparaison flux*Atm and calib star Flux without atm")
            pl.legend(["sim spectro ","flux*Atm ", "sim inter"])
        diff = prodFluxAtm-estInter
        meanSig = estInter.mean()
        diffmean = diff.mean()
        relErr = diffmean/meanSig
        relStd = diff.std()/meanSig
        self.assertTrue(np.fabs(relErr) < 1e-3)    
        self.assertTrue(relStd < 1)    
     
     
    def SimAcquSpectro02(self, simSpec, fileAtm): 
        simSpec.setAtmFile(fileAtm)
        xNew = simSpec.atms.wl[::-1][:-250]
        fluxInter = tl.interpolLinear(simSpec.star.wl[::-1], simSpec.star.dEdl[::-1], xNew)
        prodFluxAtm = fluxInter*simSpec.atms.tr[::-1][:-250]
        simSpec.doPlot = True      
        ret = simSpec.getCalibFlux()
        estInter = tl.interpolLinear(ret[0][::-1], ret[1][::-1], xNew)
        self.plotSimProd(xNew, estInter, prodFluxAtm)
        diff = prodFluxAtm-estInter
        meanSig = estInter.mean()
        diffmean = diff.mean()
        relErr = diffmean/meanSig
        relStd = diff.std()/meanSig
        self.assertTrue(np.fabs(relErr) < 1e-3)    
        self.assertTrue(relStd < 1)    


    def test_SimAcquSpectroMulti(self):
        simSpec = AuxTeles()
        simSpec.doPlot = True
        simSpec.setStarFile(S_StarFile)
        self.SimAcquSpectro02(simSpec, S_AtmFile)
        simSpec.atms.plot( S_AtmFile)
        #self.SimAcquSpectro02(simSpec, S_AtmFile2)
        #self.SimAcquSpectro02(simSpec, S_AtmFile3)
        
        
    def test_zz(self):                
        if FlagPlot:   
            try:
                pl.show()
            except AttributeError:
                pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main()
    pass

try:
    pl.show()
except AttributeError:
    pass
