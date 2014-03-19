'''
Created on 16 nov. 2012

@author: colley
'''
import unittest

from templatMODTRAN import *


FileMOD = '../../data/modtran/TemplateT04.01_1.txt'


def test_plotTemplate():
    oMOD = TemplateMODTRAN(FileMOD)
    oMOD.plotTemplate()
    
def test_resampleBetween():
    oMOD = TemplateMODTRAN(FileMOD)
    oMOD.plotTemplate(500, 1000, "raw")
    #print oMOD._wl[0], oMOD._wl[-1]
    #print oMOD._wl
    oMOD.resampleBetween(500, 1000, 512)
    oMOD.plotTemplate(500, 1000, "resample")
    
    
def test_restrictWL():
    oMOD = TemplateMODTRAN(FileMOD)
    oMOD.plotTemplate()
    oMOD.restrictWL(300, 1100)
    oMOD.plotTemplate()

class Test(unittest.TestCase):


    def testName(self):
        pass

#
# TEST, MAIN
#

#test_plotTemplate()
#test_resampleBetween()
#test_restrictWL()

#
# TemplateMODTRANAirMass
#

def test_TemplateMODTRANAirMass_init():
    oTpl = TemplateMODTRANAirMass()
    
    
def test_TemplateMODTRANAirMass_interpol():
    oTpl = TemplateMODTRANAirMass()
    plotInterpolation(oTpl)
    
    
def plotInterpolation(oTpl):
    oTpl.setAirMass(1)
    O3_1 = oTpl.getTrO3()
    oTpl.setAirMass(1.1)
    O3_11 = oTpl.getTrO3()
    oTpl.setAirMass(1.05)
    O3_105 = oTpl.getTrO3()
    oTpl.setAirMass(1.45)
    O3_245 = oTpl.getTrO3()
    # 
    pl.figure()
    pl.plot(oTpl._wl, O3_1 )
    pl.plot(oTpl._wl, O3_11 )
    pl.plot(oTpl._wl, O3_105,'-.' )
    pl.plot(oTpl._wl, O3_245)    
    pl.legend(["O3 AM=1","O3 AM=1.1","O3 AM=1.05","O3 AM=1.45"], loc="best")
    pl.grid()
    
    
def test_TemplateMODTRANAirMass_resample():
    oTpl = TemplateMODTRANAirMass()
    oTpl.resampleBetween(500, 1000, 512)
    plotInterpolation(oTpl)


def test_TemplateMODTRANAirMass_restrictWL():
    oTpl = TemplateMODTRANAirMass()
    plotInterpolation(oTpl)
    oTpl.restrictWL(500, 700)
    plotInterpolation(oTpl)
    
    
def test_TemplateMODTRANAirMass_downgradeTemplate():
    oTpl = TemplateMODTRANAirMass()
    oTpl.setAirMass(1.5)
    pl.figure()
    pl.plot(oTpl._wl, oTpl.getTrAll())
    oTpl.downgradeTemplate(600)
    pl.plot(oTpl._wl, oTpl.getTrAll())
    pl.grid()
    
def test_TemplateMODTRANAirMass_downgradeTemplateAndResample():
    oTpl = TemplateMODTRANAirMass()
    oTpl.setAirMass(1.5)
    pl.figure()
    pl.plot(oTpl._wl, oTpl.getTrAll())
    newWL = np.linspace(600, 900, 512, True)
    oTpl.downgradeTemplateAndResample(600,newWL )
    pl.plot(oTpl._wl, oTpl.getTrAll())
    pl.grid()
   

#
# TEST, MAIN
#
#test_TemplateMODTRANAirMass_init()
#test_TemplateMODTRANAirMass_interpol()
#test_TemplateMODTRANAirMass_resample()
#test_TemplateMODTRANAirMass_restrictWL()
#test_TemplateMODTRANAirMass_downgradeTemplate()
test_TemplateMODTRANAirMass_downgradeTemplateAndResample()
 
try:
    pl.show()
except AttributeError:
    pass

