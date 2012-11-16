'''
Created on 15 nov. 2012

@author: colley
'''
import unittest
from burkeAtmModel import *


fileModtran = '/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt'


def test_init():
    timeObs= np.linspace(0, 4*3600, 20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)
    par = np.zeros(oAtm._NbPar)
    oAtm.setModelParam(par)
    oAtm.printAndPlotBurkeModel()
        
def test_setDefaultModel():   
    timeObs= np.linspace(0, 4*3600,20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setDefaultParam()
    oAtm.printAndPlotBurkeModel()
    
def test_ComputeAtmTransmission():   
    timeObs= np.linspace(0, 4*3600,20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setDefaultParam()
    oAtm.computeAtmTransmission(np.pi/3, np.pi/2, timeObs[3], 800)    
    oAtm.plotCurrentTrans()
    oAtm.computeAtmTransmission(np.pi/4, np.pi/2, timeObs[3], 850)
    oAtm.plotCurrentTrans()
    
def test_AltVariation():
    timeObs= np.linspace(0, 4*3600,20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setDefaultParam()
    pl.figure()
    aAlt= np.deg2rad(np.linspace(30, 90, 4))
    aLgd = []
    for alt in aAlt:
        tr = oAtm.computeAtmTransmission(alt, 0, timeObs[3], 780)          
        pl.plot(oAtm._aWL, tr)
        aLgd.append("alt %.1f degree"%np.rad2deg(alt))
    pl.legend(aLgd, loc=4)
    pl.xlabel("Angstrom")
    pl.ylabel("%")
    pl.grid()    
    pl.title("Burke model atm. trans. var. with altitude ")
    
def test_PresVariation():
    timeObs= np.linspace(0, 4*3600,20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setDefaultParam()
    pl.figure()
    aPres= np.linspace(650, 950, 4)
    aLgd = []
    for pres in aPres:
        tr = oAtm.computeAtmTransmission(np.pi/2, 0, timeObs[3], pres)          
        pl.plot(oAtm._aWL, tr)
        aLgd.append("pres ratio %.2f"%oAtm._PresRat)
    pl.legend(aLgd, loc=4)
    pl.xlabel("Angstrom")
    pl.ylabel("%")
    pl.grid()    
    pl.title("Burke model atm. trans. varariation with pressure ")
 
def test_CH20Variation():
    timeObs= np.linspace(0, 4*3600,4)
    aVar = np.linspace(0.6, 1.4,4)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setDefaultParam()
    oAtm.setH20Param(aVar)
    pl.figure()    
    aLgd = []
    for idx in range(len(aVar)):
        tr = oAtm.computeAtmTransmission(np.pi/2, 0, timeObs[idx], 750)          
        pl.plot(oAtm._aWL, tr)
        aLgd.append("C_H2O %.2f"%aVar[idx])
    pl.legend(aLgd, loc=4)
    pl.xlabel("Angstrom")
    pl.ylabel("%")
    pl.grid()    
    pl.title("Burke model atm. trans. varariation with C_H2O ")
        
class Test(unittest.TestCase):
    def testName(self):
        pass

#test_init()
#test_setDefaultModel()
#test_ComputeAtmTransmission()
test_AltVariation()
test_PresVariation()
test_CH20Variation()


try:
    pl.show()
except AttributeError:
    pass

