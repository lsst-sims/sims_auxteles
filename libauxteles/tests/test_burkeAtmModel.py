'''
Created on 15 nov. 2012

@author: colley
'''
import unittest
from burkeAtmModel import *
import scipy.optimize as spo

fileModtran = '/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt'
CptCall = 0

def test_init():
    timeObs= np.linspace(0, 4*3600, 20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)
    par = np.zeros(oAtm._NbPar)
    oAtm.setParam(par)
    oAtm.printAndPlotBurkeModel()
        
def test_setDefaultModel():   
    timeObs= np.linspace(0, 4*3600,20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setParamExample1()
    oAtm.printAndPlotBurkeModel()
    
def test_ComputeAtmTransmission():   
    timeObs= np.linspace(0, 4*3600,20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setParamExample1()
    oAtm.computeAtmTrans(np.pi/3, np.pi/2, timeObs[3], 800)    
    oAtm.plotCurrentTrans()
    oAtm.computeAtmTrans(np.pi/4, np.pi/2, timeObs[3], 850)
    oAtm.plotCurrentTrans()
    
def test_AltVariation():
    timeObs= np.linspace(0, 4*3600,20)
    oAtm= BurkeAtmModelv1(fileModtran, timeObs)   
    oAtm.setParamExample1()
    pl.figure()
    aAlt= np.deg2rad(np.linspace(30, 90, 4))
    aLgd = []
    for alt in aAlt:
        tr = oAtm.computeAtmTransAt(alt, 0, timeObs[3], 780)          
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
    oAtm.setParamExample1()
    pl.figure()
    aPres= np.linspace(650, 950, 4)
    aLgd = []
    for pres in aPres:
        tr = oAtm.computeAtmTransAt(np.pi/2, 0, timeObs[3], pres)          
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
    oAtm.setParamExample1()
    oAtm.setParamH20(aVar)
    pl.figure()    
    aLgd = []
    for idx in range(len(aVar)):
        tr = oAtm.computeAtmTransAt(np.pi/2, 0, timeObs[idx], 750)          
        pl.plot(oAtm._aWL, tr)
        aLgd.append("C_H2O %.2f"%aVar[idx])
    pl.legend(aLgd, loc=4)
    pl.xlabel("Angstrom")
    pl.ylabel("%")
    pl.grid()    
    pl.title("Burke model atm. trans. varariation with C_H2O ")


def test_leastsq01():
    oAtm= BurkeAtmModelv1(fileModtran, np.array([0]))
    oAtm.setVarObs(np.pi/2, 0, np.array([0.]), 750)
    oAtm.setParamExample1()
    #oAtm.setParamNoEffect()
    sed = oAtm._Tpl.getTrAll()
    sed += np.random.normal(0, 1e-5, len(sed))  
    print "size sed", len(sed) 
    np.set_printoptions(precision=5) 
    def errorModel(param, atm, ydata):
        global CptCall
        CptCall += 1
        atm.setParam(param)
        print CptCall, param
        return ydata - atm.computeAtmTrans()
    oAtm._Par[0] = 0.8
    oAtm.setParamH20(np.array([1.4]))
    oAtm.printBurkeModel()
    p0 = oAtm._Par
    res = spo.leastsq(errorModel, p0, args=(oAtm, sed), full_output=True)
    print res
    if res[4] >= 0 and res[4]<=4:
        print "FIT OK :",  res[3]
        print res[0]
        oAtm.setParam(res[0])
        oAtm.printBurkeModelParam()
        if res[1] != None:             
            pl.figure()
            pl.pcolor(res[1])
            pl.colorbar()
    else:
        print "FIT NOK : ",  res[3]


def test_leastsq02():
    oAtm= BurkeAtmModelv1(fileModtran, np.array([0]))
    oAtm.setVarObs(np.pi/2, 0, np.array([0.]), 750)
    oAtm.setParamExample1()
    #oAtm.setParamNoEffect()
    sed = oAtm.computeAtmTrans()
    sed += np.random.normal(0, 1e-7, len(sed))  
    print "size sed", len(sed) 
    np.set_printoptions(precision=5) 
    def errorModel(param, atm, ydata):
        global CptCall
        CptCall += 1
        atm.setParam(param)
        #print CptCall, param
        return ydata - atm.computeAtmTrans()
    print "True parameter"
    oAtm.printBurkeModelParam()
    oAtm._Par += np.random.normal(0, 0.5, len(oAtm._Par))  
    print "\nguess parameter"
    oAtm.printBurkeModelParam()
    p0 = oAtm._Par
    res = spo.leastsq(errorModel, p0, args=(oAtm, sed), full_output=True)
    if res[4] >= 0 and res[4]<=4:
        print "FIT OK :",  res[3]
        #print res[0], res[2]["nfev"]
        print "\nSolution scipy.leastsq in %d calls "%res[2]["nfev"]
        oAtm.setParam(res[0])
        oAtm.printBurkeModelParam()
        if res[1] != None:             
            pl.figure()
            pl.pcolor(res[1])
            pl.colorbar()
        else: 
            print "no covariance estimated"
    else:
        print res
        print "FIT NOK : ",  res[3]
       
    

        
class Test(unittest.TestCase):
    def testName(self):
        pass

#test_init()
#test_setDefaultModel()
#test_ComputeAtmTransmission()
#test_AltVariation()
#test_PresVariation()
#test_CH20Variation()
test_leastsq02()

try:
    pl.show()
except AttributeError:
    pass

