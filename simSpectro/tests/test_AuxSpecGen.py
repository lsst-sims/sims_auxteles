'''
Created on 2 avr. 2013

@author: colley
'''
import unittest

from AuxSpecGen import *

S_DirExple = getModuleDirectory()+"/example/"
S_StarFile = S_DirExple+'alpha_lyr_stis_004nm.txt'
S_AtmFile = S_DirExple+'Modtran_10.dat'

def test_simuSpectroAuxTeles():
    simuSpectroAuxTeles(410, 1, 0.1, S_StarFile, S_AtmFile)
    pl.show()



test_simuSpectroAuxTeles()

class Test(unittest.TestCase):

    def test_SimAcquSpectro01(self):
        simSpec = SimAcquSpectro()
        simSpec.setAtm(S_AtmFile)
        simSpec.setStarFile(S_StarFile)
        ret = simSpec.getCalibFlux()
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    try:
        pl.show()
    except AttributeError:
        pass
