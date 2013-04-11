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



def test_simuSpectroAuxTeles():
    simuSpectroAuxTeles(410, 1, 0.1, S_StarFile, S_AtmFile)
    pl.show()


#test_simuSpectroAuxTeles()
FlagPlot= True

class Test(unittest.TestCase):

    def test_SimAcquSpectro01(self):
        simSpec = AuxTeles()
        simSpec.setAtmFile(S_AtmFile)
        simSpec.setStarFile(S_StarFile)
#        print simSpec.star.wl[::-1][:10]
#        print simSpec.atms.wl[::-1][:10]
#        print simSpec.star.wl[::-1][-10:]
#        print simSpec.atms.wl[::-1][-300:]
        xNew = simSpec.atms.wl[::-1][:-250]
        fluxInter = tl.interpolLinear(simSpec.star.wl[::-1], simSpec.star.dEdl[::-1], xNew)
        prodFluxAtm = fluxInter*simSpec.atms.tr[::-1][:-250]
        simSpec.doPlot = FlagPlot       
        ret = simSpec.getCalibFlux()
        estInter = tl.interpolLinear(ret[0][::-1], ret[1][::-1], xNew)
        if FlagPlot: 
            pl.figure()
            pl.plot(ret[0], ret[1])
            pl.plot(xNew, prodFluxAtm)
            pl.plot(xNew, estInter)
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
