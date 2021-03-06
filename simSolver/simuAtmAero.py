'''
Created on 25 feb. 2014

@author: colley
'''
import tools as tl
import matplotlib.pyplot as pl
import burkeAtmModel as bam
import numpy as np
import scipy.optimize as spo
import lmfit as lm
import modtranTools as mdt

#
# Static Constant
#

S_DoPlot = 0

S_file = "../data/simuRef/atmAeroDiffZangle.txt"
S_file = "../data/simuRef/trans_final.plt"
S_fileModtran = '../data/modtran/TemplateT04.01_1.txt'
S_fileAero = "../data/simuRef/trans_aero.plt"


def coherenceTP7TransFinal(pFileTP7, pFileFinal, pFileAero):
    """
    """
    wlMin = 3000
    wlMax =  10000
    lSelect = ["COMBIN"]
    aTP7 = mdt.readTP7(pFileTP7, lSelect, wlMin, wlMax)
    aAtm = tl.readTextFileColumn(pFileFinal)
    pl.figure()
    pl.title("no correction aero")
    pl.xlabel("nm")
    pl.plot(aAtm[:,0], aAtm[:,1])
    pl.plot(aTP7[:,0]/10, aTP7[:,1])
    pl.grid()
    pl.legend(["output simu","tp7"], loc="best")
    aAero = tl.readTextFileColumn(pFileAero)
    aAtm[:,1] /=  aAero[:,1]
    pl.figure()
    pl.title("correction aero")    
    pl.xlabel("nm")
    pl.plot(aAtm[:,0], aAtm[:,1])
    pl.plot(aTP7[:,0]/10, aTP7[:,1])
    pl.grid()
    pl.legend(["output simu","tp7"], loc="best")
    # test freqency
    pl.figure()
    pl.title("freqency")    
    pl.xlabel("nm")
    pl.plot(aAtm[:,0],'*')
    pl.plot(aTP7[:,0]/10,'*')
    pl.grid()
    atmInter = tl.interpolLinear(aAtm[:,0], aAtm[:,1], aTP7[:,0]/10)
    diff = atmInter - aTP7[:,1]
    pl.figure()
    pl.title("difference final/aero - tp7")    
    pl.xlabel("nm")
    pl.plot(aTP7[:,0]/10, diff )
    pl.grid()
    

def check_selectedTP7(pFile):
    """ 
    find tp7 major composant to approximate " COMBIN" result
    """
    wlMin = 3000
    wlMax =  10000
    lSelect = ["COMBIN","H2O", 'O2','O3','MOLEC']
    aTP7 = mdt.readTP7(pFile, lSelect, wlMin, wlMax)
    totalSelect = aTP7[:,5]*aTP7[:,2]*aTP7[:,3]*aTP7[:,4]
    diff = totalSelect - aTP7[:,1]
    fig = pl.figure()
    pl.title("TP7: COMBIN - %s"%("*".join(lSelect[1:])))
    pl.plot(aTP7[:,0], diff)
    pl.grid()
    pl.xlabel('Angstrom')
    lSelect = ["COMBIN","H2O","H2Ob" ,'O2','O3','MOLEC']
    aTP7 = mdt.readTP7(pFile, lSelect, wlMin, wlMax)
    totalSelect = aTP7[:,2]*aTP7[:,3]*aTP7[:,4]*aTP7[:,5]*aTP7[:,6]
    diff = totalSelect - aTP7[:,1]
    pl.figure()
    pl.title("TP7:  COMBIN - %s"%("*".join(lSelect[1:])))
    pl.plot(aTP7[:,0], diff)
    pl.grid()
    pl.xlabel('Angstrom')
    lSelect = ["COMBIN","H2O","H2Ob" ,'O2','O3','MOLEC', 'NO2']
    aTP7 = mdt.readTP7(pFile, lSelect, wlMin, wlMax)
    totalSelect = aTP7[:,2]*aTP7[:,3]*aTP7[:,4]*aTP7[:,5]*aTP7[:,6]*aTP7[:,7]
    diff = totalSelect - aTP7[:,1]
    pl.figure()
    pl.title("TP7:  COMBIN - %s"%("*".join(lSelect[1:])))
    pl.plot(aTP7[:,0], diff)
    pl.grid()    
    pl.xlabel('Angstrom')
    lSelect = ['NO2']
    aTP7 = mdt.readTP7(pFile, lSelect, wlMin, wlMax)
    pl.figure()
    pl.title("TP7:  NO2")
    pl.plot(aTP7[:,0], aTP7[:,1])
    pl.grid()    
    pl.xlabel('Angstrom')
    

def atmAeroPlot(pFile):
    aAtm = tl.readTextFileColumn(pFile)
    print aAtm
    pl.figure()
    pl.plot(aAtm[:,0], aAtm[:,1])
    pl.grid()
    
    
def check_wl_atmAero(pFile):
    """
    check, compare wl definition atm and aero
    """
    aAtm = tl.readTextFileColumn(pFile)
    print aAtm    
    atmTheo = bam.BurkeAtmModel(S_fileModtran)
    wlTH = atmTheo.getWL()/10.
    wlMes = aAtm[:,0]    
    print wlTH, wlMes
    diff = wlTH - wlMes
    atmTheo.setDefaultConstObs()
    atmTheo.setDefaultParam()
    defAtm = atmTheo.computeAtmTrans(True)
    pl.figure()
    pl.title(pFile)
    pl.plot(wlMes, defAtm)
    pl.plot(wlMes, aAtm[:,1])


#
# Fit for Burke model and variation
#   

def fit_BurkeNoAero(aAtm, pConstTgray=False):
    """
    fit with Burke model but with only template and Tgray
    Parameters:
     * aAtm : array transmission, column 0 is frequency
    """
    # define position of template coefficient
    #idx2Full = [0,5,6,7]
    #p0 = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    idx2Full = [5,6,7]
    p0 = np.array([ 1.0, 1.0, 1.0], dtype=np.float64)
    myParm = np.zeros(10, dtype=np.float64)
    myParm[0] = 1.0
    myParm[1] = 0.0
    def getFullParam(smallParam):
        myParm[idx2Full] = smallParam
        print myParm
        return myParm    
    def errorModel(param, atm, ydata):
        atm.setParam(getFullParam(param))
        diff = ydata - atm.computeAtmTrans()
        print (diff*diff).sum()
        return diff        
    def errorModelSum(param, atm, ydata):        
        diff = errorModel(param, atm, ydata)
        return (diff*diff).sum()
    transMeas = aAtm[:,1]
    atmTheo = bam.BurkeAtmModel(S_fileModtran)
    wlMin = 3000
    wlMax = 9000
    rangeInterval = tl.getIndexInInterval(aAtm[:,0]*10, wlMin, wlMax)
    transMeas = aAtm[rangeInterval,1]
    atmTheo.restrictWL(wlMin, wlMax)
    print atmTheo.getWL()
    print aAtm[rangeInterval,0]
    assert np.allclose(atmTheo.getWL(), aAtm[rangeInterval,0]*10, 1e-3)
    atmTheo.setConstObsComp(np.pi/2, 0, 775)
    # parameter: tau0, alpha, Cmol, C_O3, C_H20    
    if pConstTgray:
        myBounds = ((None, 1), (None,None), (None, None), (None, None))
        res =  spo.minimize(errorModelSum, p0, args=(atmTheo, transMeas), bounds=myBounds, method='L-BFGS-B' )
        solEst = res.x
    else:
        res = spo.leastsq(errorModel, p0, ftol=1e-10, args=(atmTheo, transMeas), full_output=True)
        print res
        if not (res[4] >= 0 and res[4]<=4):
            print "FIT NOK : ",  res[3]
            print res[0]
            return
        if res[1] != None and False:             
            pl.figure()
            pl.pcolor(res[1][::-1])
            pl.colorbar()
        print "FIT OK :",  res[3]
        print res[0]
        solEst = res[0]
    wl = atmTheo.getWL()
    atmTheo.printBurkeModelParam()
    atmTheo.setParam(getFullParam(solEst))
    atmTheo.printBurkeModelParam()
    baseTilte = "Burke model without aerosol"
    if 0 in idx2Full : 
        baseTilte += " with Tgray"
    else: 
        baseTilte += " no Tgray"
    pl.figure()
    #
    pl.title(baseTilte+": fit")
    pl.plot(wl, transMeas)
    pl.plot(wl, atmTheo.computeAtmTrans())  
    pl.legend(["simu", "fit"], loc="best")
    pl.grid() 
    residu = (transMeas- atmTheo.computeAtmTrans())
    pl.figure()
    #
    pl.title(baseTilte + ': residus')
    pl.plot(wl, residu) 
    pl.xlabel('Angstrom')
    pl.grid()
    pl.figure()
    #
    pl.title(baseTilte + ": histogram residus")
    pl.hist(residu, 20)
    pl.xlabel('%')
    print "residu (m, std): ", residu.mean(),residu.std()
    
    
    
def fit_atmAeroOnly(pAero):
    """
    fit 3 parameters of aerosol composante
     * Tgray
     * aerosol 1,4
     
     Parameters:
      * pAero : column 0 frequency , column 1 aero composant
    """
    idx2Full = [0, 1, 4]
    myParm = np.zeros(10, dtype=np.float64)
    def getFullParam(smallParam):
        myParm[idx2Full] = smallParam
        print myParm
        return myParm    
    def errorModel(param, atm, ydata):
        atm.setParam(getFullParam(param))
        diff = ydata - atm.computeAtmTrans()
        print (diff*diff).sum()
        return diff    
    #print aAtm
    atmTheo = bam.BurkeAtmModel(S_fileModtran)
    wlMin = 3000
    wlMax = 10000
    rangeIdx = tl.getIndexInInterval(atmTheo._aWL, wlMin, wlMax)
    pAero = pAero[rangeIdx]        
    atmTheo.restrictWL(wlMin, wlMax)
    print pAero[:,0]
    print atmTheo.getWL()/10
    # check frequency identity
    assert np.allclose(pAero[:,0], atmTheo.getWL()/10, 1e-3)
    atmTheo.setDefaultParam()
    atmTheo.setConstObsComp(np.pi/2, 0, 775)
    # parameter: Tgray , Tau0, alpha
    p0 = np.array([1.0, 0.1, -1.0], dtype=np.float64)
    res = spo.leastsq(errorModel, p0, args=(atmTheo, pAero[:,1]), full_output=True)
    print res
    wl = atmTheo.getWL()
    atmTheo.printBurkeModelParam()
    if not (res[4] >= 0 and res[4]<=4):
        print "FIT NOK : ",  res[3]
        print res[0]
        return
    print "FIT OK :",  res[3]
    print res[0]
    atmTheo.setParam(getFullParam(res[0]))
    transEst = atmTheo.computeAtmTrans()
    atmTheo.printBurkeModelParam()
    if res[1] != None and False:      
        pl.figure()
        pl.pcolor(res[1][::-1])
        pl.colorbar() 
    pl.figure()
    #
    pl.title("aerosol transmission Burke model: fit")
    pl.plot(wl, pAero[:,1])
    pl.plot(wl, transEst)
    pl.legend(["measure", "fit"], loc="best")
    pl.grid()      
    pl.figure()
    #
    pl.title('aerosol transmission Burke model: residus') 
    residu = (pAero[:,1]- transEst)
    pl.plot(wl, residu)    
    pl.xlabel('Angstrom')
    pl.grid()
    pl.figure()
    #
    pl.title("aerosol transmission Burke model: histogram residu")
    pl.hist(residu, 40)
    pl.xlabel('%') 
    print "residu (m,std):", residu.mean(),residu.std() 
    

def fit_atmBurkeAll(aAtm):
    """
    fit 6 parameters of atmospherique transmission Burke model:
     * Tgray 0
     * aerosol 1,4
     * Cmol 5
     * C_O3 6
     * C_H2O 7
    no constraint
    """
    idx2Full = [0,1,4,5,6,7]
    myParm = np.zeros(10, dtype=np.float64)
    myParm[0] = 1.0
    def getFullParam(smallParam):
        myParm[idx2Full] = smallParam
        print myParm
        return myParm    
    def errorModelLM(param, atm, ydata):
        par = np.array([param[a].value for a in param], dtype=np.float64)
        atm.setParam(getFullParam(par))
        diff = ydata - atm.computeAtmTrans()
        print (diff*diff).sum()
        return diff
    def errorModel(param, atm, ydata):
        atm.setParam(getFullParam(param))
        diff = ydata - atm.computeAtmTrans()
        print (diff*diff).sum()
        return diff    
    #atmTheo = bam.BurkeAtmModel(S_fileModtran)
    atmTheo = bam.BurkeAtmModelAM()
    wlMin = 3000
    wlMax = 9000
    rangeInterval = tl.getIndexInInterval(aAtm[:,0]*10, wlMin, wlMax)
    transMeas = aAtm[rangeInterval,1]   
    atmTheo.restrictWL(wlMin, wlMax)
    wl = atmTheo.getWL()
    assert np.allclose(wl, aAtm[rangeInterval,0]*10, 1e-3)
    atmTheo.setConstObsComp(np.pi/2, 0, 775)
    # guess    :   Tgray, tau0, alpha, Cmol, C_O3, C_H20
    p0 = np.array([0.95,  0.05, -0.05, 1.0 , 1.0 , 1.0], dtype=np.float64)        
    atmTheo.setParam(getFullParam(p0))
    atmTheo.computeAtmTrans()
    pl.figure()
    pl.title("Transmission ")
    pl.plot(wl, transMeas)
    pl.plot(wl, atmTheo._CurtTrans)
    pl.legend(["measures","guess"],loc="best")
    pl.grid()
    if True:
        parGuess = lm.Parameters() 
        parGuess.add("Tgray", p0[0], max=1.001)  
        parGuess.add("tau0", p0[1], min=0.0)
        parGuess.add("alpha", p0[2])
        parGuess.add("Cmol", p0[3])
        parGuess.add("C_O3", p0[4])
        parGuess.add("C_H20", p0[5])
        res =  lm.minimize(errorModelLM, parGuess ,args=(atmTheo, transMeas))
        if not res.success:
            print "FIT NOK : ", res.message
            return
        parEst = np.array([res.params[a].value for a in res.params], dtype=np.float64)
    else:            
        res = spo.leastsq(errorModel, p0, args=(atmTheo, transMeas),maxfev= 20000, full_output=True)
        if not (res[4] >= 0 and res[4]<=4):
            print "FIT NOK : ",  res[3]
            print res[0]
            return
        if res[1] != None:             
            pl.figure()
            pl.pcolor(res[1][::-1])
            pl.colorbar() 
        parEst = res[0]
    print "FIT OK, solution: ",  parEst    
    atmTheo.printBurkeModelParam()
    atmTheo.setParam(getFullParam(parEst))
    transEst = atmTheo.computeAtmTrans(True)
    atmTheo.printBurkeModelParam()
    pl.figure()
    pl.title("atmospheric transmission with Burke model")
    pl.plot(wl, transMeas)
    pl.plot(wl, transEst)
    pl.legend(["measures", "fit"], loc="best")
    pl.grid()
    #       
    pl.figure()
    residu = (transMeas- transEst)
    pl.plot(wl, residu)
    pl.title('atmosphere fit with Burke model, residus') 
    #pl.ylabel('%')
    pl.xlabel('Angstrom')
    pl.grid()
    pl.figure()
    pl.hist(residu, 40)
    pl.xlabel('%')
    pl.title("atmosphere fit with Burke model, histogram residu")
    print residu.mean(),residu.std() 
    # plot only aerosol
    pl.figure()
    pl.title("aerosol with Tgray")
    aAtm = tl.readTextFileColumn("../data/simuRef/trans_aero.plt")
    pl.plot(aAtm[:,0]*10, aAtm[:,1])
    pl.plot(atmTheo.getWL(), atmTheo._transGrayAero())
    pl.legend(["measures","fit"], loc="best")
    pl.grid()
    return atmTheo, transMeas
    
    
#
# Fit from file
#
  
def fit_FileAtm(pFile):
    aAtm = tl.readTextFileColumn(pFile)
    fit_atmBurkeAll(aAtm)
    
    
def fit_FileAtmPlotComposant(pFile, pFileAero, pFileTP7):
    aAtm = tl.readTextFileColumn(pFile)
    atmTheo, transMeas = fit_atmBurkeAll(aAtm)
    assert isinstance(atmTheo, bam.BurkeAtmModel)
    wl = atmTheo.getWL()
    # aerosol
    aAero = tl.readTextFileColumn(pFileAero)
    aero = tl.interpolLinear(aAero[:,0], aAero[:,1], wl/10)
    residuAero = aero - atmTheo._transGrayAero()
    pl.figure()
    pl.title("aerosol transmission")
    pl.plot(atmTheo.getWL(), aero)
    pl.plot(atmTheo.getWL(), atmTheo._transGrayAero())
    pl.legend(["simu","Burke model"], loc="best")
    pl.grid()
    # MODTRAN composant
    lSelect = ["H2O", "H2Ob", 'O2','O3','MOLEC']
    aTP7 = mdt.readTP7(pFileTP7, lSelect, wl[0]*.95, wl[-1]*1.05)
    # H2O
    cH2O = tl.interpolLinear(aTP7[:,0], aTP7[:,1]*aTP7[:,2], wl)
    residuH2O = cH2O - atmTheo._transH2O()
    pl.figure()
    pl.title("H2O transmission")
    pl.plot(atmTheo.getWL(), cH2O)    
    pl.plot(atmTheo.getWL(), atmTheo._transH2O())
    pl.legend(["simu","Burke model"], loc="best")
    pl.grid()
    # O3
    cO3 = tl.interpolLinear(aTP7[:,0], aTP7[:,4], wl)
    residuO3 = cO3 - atmTheo._transO3()
    pl.figure()
    pl.title("O3 transmission")
    pl.plot(atmTheo.getWL(), cO3)    
    pl.plot(atmTheo.getWL(), atmTheo._transO3())
    pl.legend(["simu","Burke model"], loc="best")
    pl.grid()
    # mol scattering
    cMolS = tl.interpolLinear(aTP7[:,0], aTP7[:,5], wl)
    residucMolS = cMolS - atmTheo._transMols()
    pl.figure()
    pl.title("mol scattering transmission")
    pl.plot(atmTheo.getWL(), cMolS)    
    pl.plot(atmTheo.getWL(), atmTheo._transMols())
    pl.legend(["simu","Burke model"], loc="best")
    pl.grid()
    # mol absorption : O2
    cMolA = tl.interpolLinear(aTP7[:,0], aTP7[:,3], wl)
    residucMolA = cMolA - atmTheo._transMola()
    pl.figure()
    pl.title("mol absorption transmission")
    pl.plot(atmTheo.getWL(), cMolA)    
    pl.plot(atmTheo.getWL(), atmTheo._transMola())
    pl.legend(["simu","Burke model"], loc="best")
    pl.grid()
    pl.figure()
    pl.title('residus composant')
    pl.plot(atmTheo.getWL(), residuAero)
    pl.plot(atmTheo.getWL(), residuH2O)
    pl.plot(atmTheo.getWL(), residuO3)
    pl.plot(atmTheo.getWL(), residucMolA)
    pl.plot(atmTheo.getWL(), residucMolS)
    pl.xlabel("angstrom")
    pl.legend(["Aero","H2O", "O3", "MolA","MolS"], loc="best")
    pl.grid()
    
    
def fit_FileAero(pFile, nCol=1):
    aAtm = tl.readTextFileColumn(pFile)
    aAtm = aAtm.T[[0,nCol]]
    fit_atmAeroOnly(aAtm.T)
   
    
def fit_FileAtmCorAero(pFileAtm, pFileAero):
    aAtm = tl.readTextFileColumn(pFileAtm)    
    aAero = tl.readTextFileColumn(pFileAero)
    aAtm[:,1] /=  aAero[:,1] 
    fit_BurkeNoAero(aAtm) 
   
   
def fit_FileAtmCorAeroN02(pFileAtm, pFileAero, pFileTP7):
    print "="*80+"\n fit_FileAtmCorAeroN02\n"+"="*80
    aAtm = tl.readTextFileColumn(pFileAtm)    
    aAero = tl.readTextFileColumn(pFileAero)
    aAtm[:,1] /=  aAero[:,1]
    print "aAtm", aAtm[:,0]
    IdxRet = tl.indexInIntervalCheck(aAtm[:,0], 295, 1100)
    print aAtm.shape
    print IdxRet
    aAtm = aAtm[IdxRet[0]:IdxRet[1]]
    # read TP7 N02 composant
    aTP7 = mdt.readTP7(pFileTP7, ["NO2"])
    # interpole to frequence aATm
    print "aTP7:", aTP7[:,0]
    print "aAtm", aAtm[:,0]
    transNO2 = tl.interpolLinear(aTP7[:,0]/10, aTP7[:,1] , aAtm[:,0], True)
    # correction NO2 to  aAtm
    aAtm[:,1] /=  transNO2
    fit_BurkeNoAero(aAtm) 
         
         
def fit_FileAtmAeroBounds(pFile):
    idx2Full = [1,4,5,6,7]
    myParm = np.zeros(10, dtype=np.float64)
    myParm[0] = 1.0
    def getFullParam(smallParam):
        myParm[idx2Full] = smallParam
        print "Alpha=",myParm[4]
        return myParm    
    def errorModel(param, atm, ydata):
        atm.setParam(getFullParam(param))
        diff = ydata - atm.computeAtmTrans()
        residu2= (diff*diff).sum()
        return   residu2    
    aAtm = tl.readTextFileColumn(pFile)
    print aAtm
    transMeas = aAtm[:,1]
    atmTheo = bam.BurkeAtmModel(S_fileModtran)
    wlMin=3000
    wlMax = 11000
    transMeas = aAtm[:,1][tl.getIndexInInterval(atmTheo._aWL, wlMin, wlMax)]
    atmTheo.restrictWL(wlMin, wlMax)    
    atmTheo.setConstObsComp(np.pi/2, 0, 782)
    # parameter: tau0, alpha, Cmol, C_O3, C_H20
    p0 = np.array([1.0, -0.05, 1.0, 1.0, 1.0], dtype=np.float64)
    myBounds = ((0, None), (None,0), (None, None), (None, None), (None, None))
    #myBounds = ((0, None), (None,0), (None, None), (None, None), (None, None))
    #res =  spo.minimize(errorModel, p0, args=(atmTheo, transMeas), method='L-BFGS-B' , bounds= myBounds)
    res =  spo.minimize(errorModel, p0, args=(atmTheo, transMeas), bounds=myBounds, method='L-BFGS-B' )
    print res
    wl = atmTheo.getWL()
    atmTheo.printBurkeModelParam()
    print res.message
    print res.x
    if not res.success:
        print "FIT NOK : "
        return
    print "FIT OK :"    
    atmTheo.setParam(getFullParam(res.x))
    atmTheo.printBurkeModelParam()
    pl.figure()
    pl.plot(wl, transMeas)
    pl.plot(wl, atmTheo.computeAtmTrans(True))        
    pl.figure()
    residu = (transMeas- atmTheo.computeAtmTrans())*100
    pl.plot(wl, residu)
    pl.title('atmosphere fit with Burke model, residu') 
    pl.ylabel('%')
    pl.xlabel('Angstrom')
    pl.grid()
    pl.figure()
    pl.hist(residu, 20)
    pl.xlabel('%')
    pl.title("atmosphere fit with Burke model, histogram residu")
    print residu.mean(),residu.std()    
    
    
def fit_FileatmAeroBoundslmfit(pFile):
    """
    
    """
    idx2Full = [1,4,5,6,7]
    myParm = np.zeros(10, dtype=np.float64)
    myParm[0] = 1.0
    def getFullParam(smallParam):
        myParm[idx2Full] = smallParam
        print "Alpha=",myParm[4]
        return myParm    
    def errorModel(param, atm, ydata):
        par = np.array([param[a].value for a in param], dtype=np.float64)
        atm.setParam(getFullParam(par))
        diff = ydata - atm.computeAtmTrans()
        residu2 = (diff*diff).sum()
        print residu2
        return diff   
    aAtm = tl.readTextFileColumn(pFile)
    print aAtm
    transMeas = aAtm[:,1]
    atmTheo = bam.BurkeAtmModel(S_fileModtran)
    wlMin = 3000
    wlMax = 11000
    transMeas = aAtm[:,1][tl.getIndexInInterval(atmTheo._aWL, wlMin, wlMax)]
    atmTheo.restrictWL(wlMin, wlMax)    
    atmTheo.setConstObsComp(np.pi/2, 0, 775)
    # parameter: tau0, alpha, Cmol, C_O3, C_H20
    p0 = np.array([1.0, -0.05, 1.0, 1.0, 1.0], dtype=np.float64)
    parGuess = lm.Parameters()    
    parGuess.add("tau0", p0[0], min =0)
    parGuess.add("alpha", p0[1], max=0)
    parGuess.add("Cmol", p0[2])
    parGuess.add("C_O3", p0[3])
    parGuess.add("C_H20", p0[4])
    res =  lm.minimize(errorModel, parGuess ,args=(atmTheo, transMeas))
    print res
    wl = atmTheo.getWL()
    atmTheo.printBurkeModelParam()
    print res.message
    if not res.success:
        print "FIT NOK : "
        return
    print "FIT OK :" 
    parEst = np.array([res.params[a].value for a in res.params], dtype=np.float64) 
    print parEst
    atmTheo.setParam(getFullParam(parEst))
    atmTheo.printBurkeModelParam()
    pl.figure()
    pl.plot(wl, transMeas)
    pl.plot(wl, atmTheo.computeAtmTrans(True))        
    pl.figure()
    residu = (transMeas- atmTheo.computeAtmTrans())*100
    pl.plot(wl, residu)
    pl.title('atmosphere fit with Burke model, residu') 
    pl.ylabel('%')
    pl.xlabel('Angstrom')
    pl.grid()
    pl.figure()
    pl.hist(residu, 20)
    pl.xlabel('%')
    pl.title("atmosphere fit with Burke model, histogram residu")
    print residu.mean(),residu.std()


    
#
# MAIN
#

#fit_atmBurkeAll(S_file)
#fit_atmBurkeAll(S_file)
#check_wl_atmAero("../data/simuRef/trans_final.plt")
#check_wl_atmAero("../data/simuRef/trans_aero.plt")
#fit_atmBurkeAll(S_file)
#fit_FileatmAeroBoundslmfit(S_file)
#fit_FileAero(S_fileAero,1)

lSelect = ["COMBIN","H2O", 'O2','O3','MOLEC']
# lSelect = ["H2Ob","H2O"]
tp7File = "/home/colley/temp/tmp.tp7"
#readMODTRANtp7(tp7File, lSelect, 3000, 9000)

#check_selectedTP7(tp7File)
#coherenceTP7TransFinal(tp7File, S_file, S_fileAero)
#fit_FileAtmCorAero(S_file, S_fileAero)

#check_selectedTP7(tp7File)
#fit_FileAtmCorAero(S_file,S_fileAero )
#fit_FileAtmCorAeroN02(S_file, S_fileAero, tp7File)
#fit_FileAtm(S_file)
fit_FileAtmPlotComposant(S_file, S_fileAero, tp7File)


try:
    pl.show()
except AttributeError:
    pass
