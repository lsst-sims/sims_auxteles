'''
Created on 6 nov. 2012

@author: colley
'''


import unittest
from tools import *


def test_deltaMagnitudeABConst():
    print "================test_deltaMagnitudeABConst"
    size = 2000
    wl_SFDstar = np.linspace(250, 1200, size)
    SFDstar = np.ones(size)    
    TransTrue = SFDstar
    TransEst = 0.95*TransTrue
    wl_Trans = wl_SFDstar
    delta = deltaMagnitudeAB(SFDstar, wl_SFDstar, TransTrue, TransEst, wl_Trans, 'U')
    print "delta:", delta
    delta = deltaMagnitudeAB(SFDstar, wl_SFDstar, TransTrue, TransEst, wl_Trans, 'I')
    print "delta:", delta
    delta = deltaMagnitudeAB(SFDstar, wl_SFDstar, TransTrue, TransEst, wl_Trans, 'Y')
    print "delta:", delta
    
    
def test_deltaMagnitudeABLin():
    print "================test_deltaMagnitudeABLin"
    size = 2000
    wl_SFDstar = np.linspace(250, 1200, size)
    SFDstar = np.ones(size)    
    TransTrue = SFDstar
    errTr = np.linspace(0.95, 1.3, size)
    TransEst = errTr*TransTrue
    wl_Trans = wl_SFDstar
    delta = deltaMagnitudeAB(SFDstar, wl_SFDstar, TransTrue, TransEst, wl_Trans, 'U')
    print "delta:", delta
    delta = deltaMagnitudeAB(SFDstar, wl_SFDstar, TransTrue, TransEst, wl_Trans, 'I')
    print "delta:", delta
    delta = deltaMagnitudeAB(SFDstar, wl_SFDstar, TransTrue, TransEst, wl_Trans, 'Y')
    print "delta:", delta
    
def test_vega(): 
    print coefKuruczEarth(9600, 0.03, 129)

def test_getOutputDir():
    print "test_getOutputDir: "
    print getOutputDir()
 
test_vega()    
test_deltaMagnitudeABConst()
test_deltaMagnitudeABLin()
test_getOutputDir()


class Test(unittest.TestCase):


    def test_01(self):
        ret = readTextFileColumn('file_tools_01.txt',';')       
        print ret
        
        
    def test_02(self):
        # mon comm
        ret = readTextFileColumn(getRootPackage()+'/data/modtran/TemplateT04.01_1.txt')        
        print ret


    def test_stellarRadius(self):
        """
        test with Sirius 
        ref : http://cas.sdss.org/dr6/en/proj/advanced/hr/radius1.asp        
        """
        siriusRadius  = stellarRadius(9500, -1.44,  379.21)
        rap = siriusRadius/696342
        self.assertTrue(np.fabs(rap - 1.76) < 0.02)  

    def test_coefKuruczEarth(self):
        """
        test with Sirius
          r = 1.76*696342
          d= 2.64*3.085e13
          In [8]: (r/d)**2
          Out[8]: 2.2643957244259747e-16
        """
        coef = coefKuruczEarth(9500, -1.44,  379.21) 
        tv = 2.2643957244259747e-16
        #err = np.fabs(coef-tv)   /tv     
        test = np.allclose(np.array([coef]) ,np.array([tv]), atol=0, rtol=1e-2)
        self.assertTrue(test)  
        
       
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()