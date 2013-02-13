'''
Created on 4 dec. 2012

@author: colley
'''
import unittest
from observationAuxTeles import *
import kurucz as kur
from burkeAtmModel import BurkeAtmModelv2
import starTargetSimu as star
import pylab as pl

FileKuruczPic = '/home/colley/projet/lsst/stellar_spectra/k93.pic'
fileModtran = '/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt'

def test_readObsNight():
    oKur = kur.Kurucz(FileKuruczPic)
    oKur.resampleBetween(5000, 10000, 1024)
    oStar = star.StarTargetSimu(oKur)
    oAtm = BurkeAtmModelv2(fileModtran)
    oAtm.resample( oKur.getWL())
    oObs = ObsSurveySimu01(2,2)
    oObs.setAtmModel(oAtm)
    oObs.setStarTarget(oStar)
    oObs.readObsNight("")
    print oObs.aConst
    print oObs.getTrueParam()
    print oObs.getGuessDefault()
    
test_readObsNight()

try:
    pl.show()
except AttributeError:
    pass

    
 
class Test(unittest.TestCase):


    def testName(self):
        pass

#unittest.main()