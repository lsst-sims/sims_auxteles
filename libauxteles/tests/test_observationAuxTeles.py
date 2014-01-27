'''
Created on 4 dec. 2012

@author: colley
'''
import unittest
from observationAuxTeles import *
import kurucz as kur
from burkeAtmModel import BurkeAtmModel
import starTargetSimu as star
import pylab as pl
import simuV2_2

S_FileKuruczPic = '../../data/kurucz93/k93.pic'
S_fileModtran = '../../data/modtran/TemplateT04.01_1.txt'
S_StarPos = "../../data/simuRef/starPos2013-12-21.pkl"


##################################################
# ObsSurveySimu01
##################################################
def test_readObsNight():
    oKur = kur.Kurucz(S_FileKuruczPic)
    oKur.resampleBetween(5000, 10000, 1024)
    oStar = star.StarTargetSimuAll(oKur)
    oAtm = BurkeAtmModel(S_fileModtran)
    oAtm.resample( oKur.getWL())
    oObs = ObsSurveySimu01(2,2)
    oObs.setAtmModel(oAtm)
    oObs.setStarTarget(oStar)
    oObs.readObsNight("")
    print oObs.aConst
    print oObs.getTrueParam()
    print oObs.getGuessDefault()

#
# TEST
#
#test_readObsNight()

##################################################
# ObsSurveySimuV2_2
##################################################

def test_ObsSurveySimuV2_2_init():
    print "============== test_ObsSurveySimuV2_2_init"
    tmp = obS.DataPositionStarTarget()
    dataStar = tmp.read(S_StarPos)    
    oObs = ObsSurveySimuV2_2(dataStar)
    print oObs.DataStar


def test_ObsSurveySimuV2_2_doSimu():
    print "============== test_ObsSurveySimuV2_2_doSimu"
    simu = simuV2_2.SimuVersion2_2()
    simu.setFileStarPos(S_StarPos)
    simu.computeExposureTime(200, 500)
    print simu.dataStar
    simu.doSimuWithFastModel()    
    simu.oObs.plotFlux(0)
    simu.oObs.plotFlux(1)
    simu.oObs.plotFlux(2)
    simu.oObs.plotFlux(3)
    simu.oObs.plotFlux(4)
    simu.oObs.plotFlux(5)
    

def test_ObsSurveySimuV2_2_posStar():
    print "============== test_ObsSurveySimuV2_2_doSimu"
    simu = simuV2_2.SimuVersion2_2()
    simu.setFileStarPos(S_StarPos)
    simu.computeExposureTime(200, 500)
    print simu.dataStar
    simu.doSimuWithFastModel()
    simu.oObs.plotAllPositionStar("", 'starpos') 
 
 
#
# TEST
#

#test_ObsSurveySimuV2_2_init()
#test_ObsSurveySimuV2_2_doSimu()
test_ObsSurveySimuV2_2_posStar()


try:
    pl.show()
except AttributeError:
    pass

    
 
class Test(unittest.TestCase):


    def testName(self):
        pass

#unittest.main()