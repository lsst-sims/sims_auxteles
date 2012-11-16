'''
Created on 16 nov. 2012

@author: colley
'''
import unittest
from templatMODTRAN import *


FileMOD = '/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt'


def test_plotTemplate():
    oMOD = TemplateMODTRAN(FileMOD)
    oMOD.plotTemplate()
    
class Test(unittest.TestCase):


    def testName(self):
        pass


test_plotTemplate()

try:
    pl.show()
except AttributeError:
    pass

