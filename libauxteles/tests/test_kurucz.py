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

    def test_nearestFlux(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.getNearestFlux(1.3, 8800, 2.7)
        

test_plotFlux()
pl.show()
unittest.main()
