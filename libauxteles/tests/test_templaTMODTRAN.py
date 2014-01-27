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
    
class Test(unittest.TestCase):


    def testName(self):
        pass


#test_plotTemplate()
test_resampleBetween()

try:
    pl.show()
except AttributeError:
    pass

