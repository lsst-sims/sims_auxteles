'''
Created on 15 nov. 2012

@author: colley
'''
import unittest
from burkeAtmModel import *


def test_init():
    wl = np.linspace(2000, 10000, 10000, 256)
    oAtm= BurkeAtmModelv1('/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt', wl)
    par = np.zeros(oAtm._NbPar)
    oAtm.setBurkeParamModel(par)
    oAtm.printBurkeModel()
    
    
class Test(unittest.TestCase):
    def testName(self):
        pass

test_init()
pl.show()
