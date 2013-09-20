import unittest
from burkeAtmModel import *
import scipy.optimize as spo


G_fileModtran = '../../data/modtran/TemplateT04.01_1.txt'
CptCall = 0


def plotErrRelAtm(pTrue, pEst, pTitle=""):    
    pl.figure()
    pl.title(pTitle)
    Err = 100*(pTrue- pEst)/pTrue
    pl.plot(Err)
    pl.ylabel("%")
    pl.grid()
    
    
def test_init():    
    oAtm= BurkeAtmModel(G_fileModtran)
    par = np.zeros(oAtm._NbPar)
    oAtm.setParam(par)
    oAtm.printBurkeModelParam()
        
        
def test_setDefaultModel():   
    oAtm= BurkeAtmModel(G_fileModtran)   
    oAtm.setParamExample1()
    oAtm.printBurkeModelParam()
    
    
def test_ComputeAtmTransmission():       
    oAtm= BurkeAtmModel(G_fileModtran)   
    oAtm.setParamExample1()
    oAtm.computeAtmTransAtConstObs(np.pi/3, np.pi/2, 800)    
    oAtm.plotCurrentTrans()
#    oAtm.computeAtmTransAtConstObs(np.pi/4, np.pi/2, 850)
#    oAtm.plotCurrentTrans()
    oAtm.printBurkeModel()
    
    
def test_AltVariation():
    oAtm= BurkeAtmModel(G_fileModtran)   
    oAtm.setDefaultParam()
    pl.figure()
    aAlt= np.deg2rad(np.linspace(30, 90, 4))
    aLgd = []
    for alt in aAlt:
        tr = oAtm.computeAtmTransAtConstObs(alt, 0, 780)          
        pl.plot(oAtm._aWL, tr)
        aLgd.append("alt %.1f degree"%np.rad2deg(alt))
    pl.legend(aLgd, loc=4)
    pl.xlabel("Angstrom")
    pl.ylabel("%")
    pl.grid()    
    pl.title("Burke model atm. trans. var. with altitude ")

    
def test_PresVariation():    
    oAtm= BurkeAtmModel(G_fileModtran)   
    oAtm.setParamExample1()
    pl.figure()
    aPres= np.linspace(650, 950, 4)
    aLgd = []
    for pres in aPres:
        tr = oAtm.computeAtmTransAtConstObs(np.pi/2, 0,  pres)          
        pl.plot(oAtm._aWL, tr)
        aLgd.append("pres ratio %.2f"%oAtm._PresRat)
    pl.legend(aLgd, loc=4)
    pl.xlabel("Angstrom")
    pl.ylabel("%")
    pl.grid()    
    pl.title("Burke model atm. trans. varariation with pressure ")
 
 
def test_CH20Variation():   
    aVar = np.linspace(0.6, 1.4,4)
    oAtm= BurkeAtmModel(G_fileModtran)   
    oAtm.setParamExample1()
    oAtm.setConstObsComp(np.pi/2, 0, 750)
    par= np.copy(oAtm._Par)
    pl.figure()    
    aLgd = []
    for idx in range(len(aVar)):
        par[7] = aVar[idx]
        oAtm.setParam(par)
        tr = oAtm.computeAtmTrans()       
        pl.plot(oAtm._aWL, tr)
        aLgd.append("C_H2O %.2f"%aVar[idx])
    pl.legend(aLgd, loc=4)
    pl.xlabel("Angstrom")
    pl.ylabel("%")
    pl.grid()    
    pl.title("Burke model atm. trans. varariation with C_H2O ")


def test_leastsq01():
    oAtm= BurkeAtmModelv1(G_fileModtran, np.array([0]))
    oAtm.setConstObsComp(np.pi/2, 0, np.array([0.]), 750)
    oAtm.setParamExample1()
    #oAtm.setDefaultParam()
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
            pl.pcolor(res[1][::-1])
            pl.colorbar()
    else:
        print "FIT NOK : ",  res[3]


def test_leastsq02():
    oAtm= BurkeAtmModel(G_fileModtran, np.array([0]))
    oAtm.setConstObsComp(np.pi/3, np.pi/4, np.array([0.]), 750)
    oAtm.setParamExample1()
    #oAtm.setDefaultParam()
    sed = oAtm.computeAtmTrans()
    sed += np.random.normal(0, 1e-5, len(sed))  
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
    oAtm._Par += np.random.normal(0, 0.2, len(oAtm._Par))
    oAtm._Par[2] = oAtm._Par[3] =oAtm._Par[1] =0
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
            oAtm.plotCovarMatrix(res[1], "Covariance matrix estimated with 1 transmission")
            print res[1]
        else: 
            print "no covariance estimated"
    else:
        print res
        print "FIT NOK : ",  res[3]
       
       
def test_leastsq03():    
    oAtm= BurkeAtmModel(G_fileModtran)
    matVarObs = np.array([[np.pi/2, 0,  750], 
                       [np.pi/3.5, np.pi/3, 750],
                       [np.pi/4, np.pi,  750]])    
    MyLegv=[]
    oAtm.setConstObs(matVarObs[0,:])
    oAtm.setParamExample1()
    parTrue= np.copy(oAtm._Par)
    sed = oAtm.computeAtmTrans()
    sed1 = np.copy(sed)
    MyLegv.append('%f'%oAtm._AirMass)
    #oAtm.plotCurrentTrans()
    oAtm.printBurkeModel()
    oAtm.setConstObs(matVarObs[1,:])
    MyLegv.append('%f'%oAtm._AirMass)
    #oAtm.plotCurrentTrans()
    sed2 = oAtm.computeAtmTrans().copy()
    oAtm.setConstObs(matVarObs[2,:])
    MyLegv.append('%f'%oAtm._AirMass)
    #oAtm.plotCurrentTrans()
    sed3 = oAtm.computeAtmTrans().copy()    
    sed = np.concatenate((sed,sed2))
    sed = np.concatenate((sed,sed3)) 
    sed += np.random.normal(0, 1e-4, len(sed))* sed
    pl.figure()
    pl.plot(oAtm._aWL, sed1)
    pl.plot(oAtm._aWL,sed2)
    pl.plot(oAtm._aWL,sed3)
    pl.plot(oAtm._aWL,sed[:len(sed1)])
    MyLegv.append('with Noise')
    pl.legend(MyLegv)
    pl.grid()      
    #sed += np.random.normal(0, 1e-3, len(sed))    
    def errorModel(param, atm, ydata, varObs):
        #print "Param: ", param        
        atm.setParam(param)
        atm.setConstObs(matVarObs[0,:])
        sed1 = np.copy(atm.computeAtmTrans())
        atm.setConstObs(matVarObs[1,:])
        sed = np.concatenate((sed1,atm.computeAtmTrans()))
        atm.setConstObs(matVarObs[2,:])
        sed = np.concatenate((sed,atm.computeAtmTrans()))
        print "chi2: ",((ydata - sed)**2).sum()      
        return ydata - sed
    print "True parameter"    
    oAtm.printBurkeModelParam()
    #oAtm._Par += np.random.normal(0, 0.05, len(oAtm._Par))*oAtm._Par
    #oAtm.setDefaultParam() 
    print "\nguess parameter"
    oAtm.printBurkeModelParam()
    p0 = oAtm._Par
    p0 += np.random.normal(0, 0.5, len(oAtm._Par))*p0
    oAtm.setDefaultParam()
    p0 = oAtm._Par
    plotErrRelAtm(parTrue, p0, "err ral guess")
    res = spo.leastsq(errorModel, p0, args=(oAtm, sed, matVarObs), full_output=True)
    print 'erreur relative', np.fabs((parTrue-res[0])/parTrue)
    if res[4] >= 0 and res[4]<=4:
        print "FIT OK :",  res[3]
        #print res[0], res[2]["nfev"]
        print "\nSolution scipy.leastsq in %d calls "%res[2]["nfev"]
        oAtm.setParam(res[0])
        print parTrue
        oAtm.printBurkeModelParam()
        pl.figure()        
        pl.plot(p0)
        pl.plot(res[0])
        pl.plot(parTrue)
        pl.legend(['guess','estimated','True'])
        pl.grid()
        oAtm.plotErrRelAtm(parTrue, res[0], "err rel. estimated")
        print "compare parameter"
        print p0
        print parTrue
        print res[0]
        if res[1] != None:             
            oAtm.plotCovarMatrix(res[1], "Covariance matrix estimated with 3 trans. obser")
            oAtm.plotCorrelMatFromCovMat(res[1], "Correlation matrix estimated with 3 trans. obser")
            #print res[1]
        else: 
            print "no covariance estimated"
    else:
        #print res
        print "FIT NOK : ",  res[3]   


def test_downgradeTemplate():
    oAtm= BurkeAtmModel(G_fileModtran)
    oAtm.setParamExample1()
    trFull = oAtm.computeAtmTransAtConstObs(np.pi/3, np.pi/2, 800)
    res = 200
    oAtm.downgradeTemplate(res)
    trDown = oAtm.computeAtmTransAtConstObs(np.pi/3, np.pi/2, 800)
    pl.figure()
    pl.title("downgradeTemplate resolution to %d"%res)
    pl.plot(oAtm.getWL(), trFull)
    pl.plot(oAtm.getWL(), trDown)
    pl.grid()
    
    
def test_downgradeTrans():
    oAtm= BurkeAtmModel(G_fileModtran)
    oAtm.setParamExample1()
    trFull = oAtm.computeAtmTransAtConstObs(np.pi/3, np.pi/2, 800)
    x = oAtm.getWL().copy()
    res = 400
    newx = np.linspace(3000, 9500, 800)
    trDownTrans = oAtm.downgradeTransAndResample(res, newx)
    oAtm.downgradeTemplateAndResample(res,newx)
    trDownTemplate = oAtm.computeAtmTransAtConstObs(np.pi/3, np.pi/2, 800)
    pl.figure()
    pl.title("downgradeTemplate and resample resolution to %d"%res)
    pl.plot(x, trFull)
    pl.plot(newx, trDownTrans)
    pl.plot(newx, trDownTemplate)
    pl.legend(["raw","down trans","down template"])
    pl.grid()
    
    
def test_downgradeTemplateAndResample(res = 600):
    newx = np.linspace(3000, 9000, 800)
    oAtm= BurkeAtmModel(G_fileModtran)
    oAtm.setParamExample1()
    trFull = oAtm.computeAtmTransAtConstObs(np.pi/3, np.pi/2, 800)
    x = oAtm.getWL().copy()    
    oAtm.downgradeTemplateAndResample(res,newx)
    trDown = oAtm.computeAtmTransAtConstObs(np.pi/3, np.pi/2, 800)
    oAtm.plotCurrentTrans()
    pl.figure()
    pl.title("downgradeTemplate and resample resolution to %d"%res)
    pl.plot(x, trFull)
    pl.plot(oAtm.getWL(), trDown)
    pl.grid()
       
        
class Test(unittest.TestCase):
    def testName(self):
        pass
    
    
    def test_covar2Correl(self):
        oAtm= BurkeAtmModel(G_fileModtran)
        mat = np.array([[2,4],[4,9]])
        oAtm._covar2Correl(mat)
        
        
#test_init()
#test_setDefaultModel()
#test_ComputeAtmTransmission()
#test_downgradeTemplate()
test_downgradeTemplateAndResample(500)
test_downgradeTemplateAndResample(200)
test_downgradeTemplateAndResample(100)
#test_downgradeTrans()

#test_comparedowngrade()
#test_AltVariation()
#test_PresVariation()
#test_CH20Variation()

#np.random.seed(150)
#test_leastsq03()
#test_leastsq03_v2()


try:
    pl.show()
except AttributeError:
    pass
