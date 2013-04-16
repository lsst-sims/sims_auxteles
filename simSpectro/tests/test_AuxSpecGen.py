'''
Created on 2 avr. 2013

@author: colley
'''
import unittest

from AuxSpecGen import *
import tools as tl

S_DirExple = getModuleDirectory()+"/example/"
S_StarFile = S_DirExple+'alpha_lyr_stis_004nm.txt'
S_AtmFile = S_DirExple+'Modtran_10.dat'
S_AtmFile2 = S_DirExple+'Modtran_50.dat'
S_AtmFile3 = S_DirExple+'Modtran_100.dat'



def test_simuSpectroAuxTeles():
    simuSpectroAuxTeles(410, 1, 0.1, S_StarFile, S_AtmFile)
    pl.show()


#test_simuSpectroAuxTeles()
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
        simSpec.doPlot = False      
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
        simSpec.setStarFile(S_StarFile)
        self.SimAcquSpectro02(simSpec, S_AtmFile)
        self.SimAcquSpectro02(simSpec, S_AtmFile2)
        self.SimAcquSpectro02(simSpec, S_AtmFile3)
        
        
    def test_zz(self):                
        if FlagPlot:   
            try:
                pl.show()
            except AttributeError:
                pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

#try:
#    pl.show()
#except AttributeError:
#    pass
