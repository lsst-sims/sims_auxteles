'''
Created on 14 nov. 2012

@author: colley
'''
import unittest

from kurucz import *



FileKuruczPic = '/home/colley/projet/lsst/stellar_spectra/k93.pic'

def test_plotFlux():
    oKur = Kurucz(FileKuruczPic)
    oKur.plotFlux(1540)
    oKur.plotMultiFlux(1540,5)
    oKur.plotFlux(1544)
    
class Test(unittest.TestCase):


    def testName(self):
        pass



test_plotFlux()
pl.show()