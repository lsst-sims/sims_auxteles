'''
Created on 14 nov. 2012

@author: colley
'''
import unittest
from kurucz import *
import scipy.optimize as spo


FileKuruczPic = '/home/colley/projet/lsst/stellar_spectra/k93.pic'

def test_plotFlux():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(1000, 12000)
    oKur.plotFlux(1540)
    oKur.plotMultiFluxesCont(1540,5)
    oKur.plotFlux(1544)
    
    
def test_setWLInterval():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(3000, 9000)
    print oKur._IdxMin,  oKur._Flux[oKur._IdxMin-1,0], oKur._Flux[oKur._IdxMin,0], oKur._Flux[oKur._IdxMin+1,0]
    print oKur._IdxMax,  oKur._Flux[oKur._IdxMax-1,0], oKur._Flux[oKur._IdxMax,0], oKur._Flux[oKur._IdxMax+1,0]    
    oKur.setWLInterval(3030, 9010)
    print oKur._IdxMin,  oKur._Flux[oKur._IdxMin-1,0], oKur._Flux[oKur._IdxMin,0], oKur._Flux[oKur._IdxMin+1,0]
    print oKur._IdxMax,  oKur._Flux[oKur._IdxMax-1,0], oKur._Flux[oKur._IdxMax,0], oKur._Flux[oKur._IdxMax+1,0]    
     
    
    
def test_getFluxInterLin():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(2000, 10000)
    pl.figure()
    fl = oKur.getFluxInterLin(-1.3, 6000, 2.7)
    pl.plot(oKur.getWL(), fl.T)
    fl = oKur.getFluxInterLin(-1.3, 6242, 2.7)
    pl.plot(oKur.getWL(), fl.T)
    fl = oKur.getFluxInterLin(-1.3, 6300, 2.7)
    pl.plot(oKur.getWL(), fl.T)
    pl.legend(['6000','6242','6300'])
    pl.grid()
    
         
class Test(unittest.TestCase):

    def notest_nearestFlux(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(1000, 12000)
        idx= oKur.getFluxIdxNPG(1.3, 6100, 2.7)
        oKur.plotFlux(idx)
        
    
    def notest_nearestFluxes(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)
        aIdx = []
        aIdx.append(oKur.getFluxIdxNPG(-1, 6100, 2.7))
        aIdx.append(oKur.getFluxIdxNPG(-2, 6100, 2.7))
        aIdx.append(oKur.getFluxIdxNPG(-3, 6100, 2.7))
        aIdx.append(oKur.getFluxIdxNPG(-4, 6100, 2.7))
        oKur.plotMultiFluxes(aIdx)
        aIdx = []
        aIdx.append(oKur.getFluxIdxNPG(-1, 6000, 2.7))
        aIdx.append(oKur.getFluxIdxNPG(-1, 7000, 2.7))
        aIdx.append(oKur.getFluxIdxNPG(-1, 8000, 2.7))
        aIdx.append(oKur.getFluxIdxNPG(-1, 9000, 2.7))
        oKur.plotMultiFluxes(aIdx)
        aIdx = []
        aIdx.append(oKur.getFluxIdxNPG(-1, 6100, 2))
        aIdx.append(oKur.getFluxIdxNPG(-1, 6100, 3))
        aIdx.append(oKur.getFluxIdxNPG(-1, 6100, 4))
        aIdx.append(oKur.getFluxIdxNPG(-1, 6100, 5))
        oKur.plotMultiFluxes(aIdx)
        pl.show()   
    
    def test_leastsq(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)        
        flux =oKur.getFluxIdxNPG(-2, 6100, 2.7)        
        def errorModel(par, kur, ydata):
            print par 
            res = ydata - kur.getFluxInterLin(par[0], par[1], par[2])
            return res.ravel()
        print "True parameter"
        print "-2, 6100, 2.7"
        print "\nguess parameter"
        p0 = np.array([-2.1, 6300., 2.5])
        print p0
        res = spo.leastsq(errorModel, p0, args=(oKur, flux), full_output=True)
        if res[4] >= 0 and res[4]<=4:
            print "FIT OK :",  res[3]
            #print res[0], res[2]["nfev"]
            print "\nSolution scipy.leastsq in %d calls "%res[2]["nfev"]
            print res[0]            
            if res[1] != None:                             
                print res[1]
            else: 
                print "no covariance estimated"
        else:
            print res
            print "FIT NOK : ",  res[3]       
#test_setWLInterval()
#test_plotFlux()
#test_getFluxInterLin()
unittest.main()
