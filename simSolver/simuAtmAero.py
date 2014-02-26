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
import copy 

#
# Static Constant
#

S_DoPlot = 0

S_file = "../data/simuRef/atmAeroDiffZangle.txt"
S_file = "../data/simuRef/trans_final.plt"
S_fileModtran = '../data/modtran/TemplateT04.01_1.txt'
S_fileAero = "../data/simuRef/trans_aero.plt"

_keyTP7 = "FREQ COMBIN    H2O   UMIX     O3  TRACE     N2    H2Ob MOLEC AER+CLD  HNO3 AER+CLD    -LOG    CO2     CO    CH4    N2O     O2    NH3     NO    NO2    SO2  CFC11  CFC12  CFC13  CFC14  CFC22 CFC113 CFC114 CFC115 CLONO2   HNO4 CHCL2F   CCL4   N2O5"
S_KeyTP7 = _keyTP7.split()




def readMODTRANtp7(pFile, pListCol, pWLMin, pWLmax):
    """
    read tp7 output modtran and select column
    WARNING : must add # to catch header !
    parameters: 
     * pFile :  string with path and name file
     * pListCol :  list of selected column
     * pWLMin, pWLmax : min max wavelength in Angstrom
     
    return:
     * numpy array where column 0 is  wavelength in Angstrom and other columns are selected by pListCol
    """  
    # reorder keyword like tp7 file sense, import for pyplot.legend()
    lSelectIdx = [S_KeyTP7.index(key) for key in pListCol]
    IdxSort = np.argsort(lSelectIdx)    
    ListColStr = [pListCol[idx] for idx in IdxSort]
    # read tp7 file
    aRet = tl.readTextFileColumn(pFile)
    if aRet is None: 
        print "ERROR: ", aRet
        return
    # selected columns
    remCol = set(range(aRet.shape[1]))
    ListCol = copy.copy(lSelectIdx)
    ListCol.append(0)
    remCol = list(remCol - set(ListCol))
    aRet = np.delete(aRet, remCol, 1)    
    # cm^-1 to angstrom
    wlAngs = 1.0e8/aRet[:,0]
    IdxRet = tl.indexInIntervalCheck(-wlAngs, -pWLmax, -pWLMin)
    if IdxRet[0] is None:
        return None
    # suppress line
    aRet = aRet[IdxRet[0]: IdxRet[1]]
    # replace wl in angstrom
    aRet[:,0] = 1.0e8/aRet[:,0]
    # growing sense
    aRet = np.flipud(aRet)
    # plot
    if S_DoPlot:
        pl.figure()
        pl.title('MODTRAN TP7 transmission selected')
        pl.xlabel("Angstrom")
        for col in aRet[:,1:].T:
            pl.plot(aRet[:,0], col)
        pl.grid()
        pl.legend(ListColStr,loc="best")
    return aRet


def coherenceTP7TransFinal(pFileTP7, pFileFinal, pFileAero):
    """
    """
    wlMin = 3000
    wlMax =  10000
    lSelect = ["COMBIN"]
    aTP7 = readMODTRANtp7(pFileTP7, lSelect, wlMin, wlMax)
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
    aTP7 = readMODTRANtp7(pFile, lSelect, wlMin, wlMax)
    totalSelect = aTP7[:,5]*aTP7[:,2]*aTP7[:,3]*aTP7[:,4]
    diff = totalSelect - aTP7[:,1]
    pl.figure()
    pl.title("TP7: COMBIN - %s"%("*".join(lSelect[1:])))
    pl.plot(aTP7[:,0], diff)
    pl.grid()
    pl.xlabel('Angstrom')
    lSelect = ["COMBIN","H2O","H2Ob" ,'O2','O3','MOLEC']
    aTP7 = readMODTRANtp7(pFile, lSelect, wlMin, wlMax)
    totalSelect = aTP7[:,2]*aTP7[:,3]*aTP7[:,4]*aTP7[:,5]*aTP7[:,6]
    diff = totalSelect - aTP7[:,1]
    pl.figure()
    pl.title("TP7:  COMBIN - %s"%("*".join(lSelect[1:])))
    pl.plot(aTP7[:,0], diff)
    pl.grid()
    pl.xlabel('Angstrom')
    lSelect = ["COMBIN","H2O","H2Ob" ,'O2','O3','MOLEC', 'NO2']
    aTP7 = readMODTRANtp7(pFile, lSelect, wlMin, wlMax)
    totalSelect = aTP7[:,2]*aTP7[:,3]*aTP7[:,4]*aTP7[:,5]*aTP7[:,6]*aTP7[:,7]
    diff = totalSelect - aTP7[:,1]
    pl.figure()
    pl.title("TP7:  COMBIN - %s"%("*".join(lSelect[1:])))
    pl.plot(aTP7[:,0], diff)
    pl.grid()    
    pl.xlabel('Angstrom')
    lSelect = ['NO2']
    aTP7 = readMODTRANtp7(pFile, lSelect, wlMin, wlMax)

    

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

def fit_BurkeNoAero(aAtm):
    """
    fit with Burke model but with only template and Tgray
    Parameters:
     * aAtm : array transmission, column 0 is frequency
    """
    # define position of template coefficient
    idx2Full = [0,5,6,7]
    p0 = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    #idx2Full = [5,6,7]
    #p0 = np.array([ 1.0, 1.0, 1.0], dtype=np.float64)
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
    transMeas = aAtm[:,1]
    atmTheo = bam.BurkeAtmModel(S_fileModtran)
    wlMin = 3000
    wlMax = 9000
    transMeas = aAtm[:,1][tl.getIndexInInterval(atmTheo._aWL, wlMin, wlMax)]
    atmTheo.restrictWL(wlMin, wlMax)    
    atmTheo.setConstObsComp(np.pi/2, 0, 775)
    # parameter: tau0, alpha, Cmol, C_O3, C_H20    
    res = spo.leastsq(errorModel, p0, ftol=1e-10, args=(atmTheo, transMeas), full_output=True)
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
    atmTheo.printBurkeModelParam()
    if res[1] != None and False:             
        pl.figure()
        pl.pcolor(res[1][::-1])
        pl.colorbar() 
    pl.figure()
    pl.plot(wl, transMeas)
    pl.plot(wl, atmTheo.computeAtmTrans())  
    pl.legend(["simu", "fit"], loc="best")
    pl.grid() 
    #
    pl.figure()
    residu = (transMeas- atmTheo.computeAtmTrans())
    pl.plot(wl, residu)
    mesTitle = 'residu fit with Burke model transmission, no aerosol'
    if 0 in idx2Full : 
        mesTitle += " with Tgray"
    else: 
        mesTitle += " no Tgray"
    pl.title(mesTitle) 
    pl.xlabel('Angstrom')
    pl.grid()
    pl.figure()
    pl.hist(residu, 20)
    pl.xlabel('%')
    pl.title("atmosphere fit with Burke model, histogram residu")
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
    

def fit_atmBurkeAll(pFile):
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
    def errorModel(param, atm, ydata):
        atm.setParam(getFullParam(param))
        diff = ydata - atm.computeAtmTrans()
        print (diff*diff).sum()
        return diff        
    aAtm = tl.readTextFileColumn(pFile)
    print aAtm
    transMeas = aAtm[:,1]
    atmTheo = bam.BurkeAtmModel(S_fileModtran)
    wlMin = 3000
    wlMax = 9000
    transMeas = aAtm[:,1][tl.getIndexInInterval(atmTheo._aWL, wlMin, wlMax)]
    atmTheo.restrictWL(wlMin, wlMax)    
    atmTheo.setConstObsComp(np.pi/2, 0, 775)
    #atmTheo.setConstObsComp(np.pi/2, 0,782)
    # parameter: tau0, alpha, Cmol, C_O3, C_H20
    p0 = np.array([1.0,1.0, -0.05, 1.0, 1.0, 1.0], dtype=np.float64)
    res = spo.leastsq(errorModel, p0, args=(atmTheo, transMeas), full_output=True)
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
    transEst = atmTheo.computeAtmTrans(True)
    atmTheo.printBurkeModelParam()
    if res[1] != None:             
        pl.figure()
        pl.pcolor(res[1][::-1])
        pl.colorbar() 
    pl.figure()
    pl.title("atmospheric transmission fit with Tgray")
    pl.plot(wl, transMeas)
    pl.plot(wl, transEst)
    pl.legend(["measure", "fit"], loc="best")
    pl.grid()
    #       
    pl.figure()
    residu = (transMeas- transEst)
    pl.plot(wl, residu)
    pl.title('atmosphere fit with Burke model with Tgray, residu') 
    #pl.ylabel('%')
    pl.xlabel('Angstrom')
    pl.grid()
    pl.figure()
    pl.hist(residu, 40)
    pl.xlabel('%')
    pl.title("atmosphere fit with Burke model with Tgray, histogram residu")
    print residu.mean(),residu.std() 
    # plot only aerosol
    pl.figure()
    pl.title("aerosol with Tgray")
    aAtm = tl.readTextFileColumn("../data/simuRef/trans_aero.plt")
    pl.plot(aAtm[:,0]*10, aAtm[:,1])
    pl.plot(atmTheo.getWL(), atmTheo._transGrayAero())
    pl.legend(["simu","fit"], loc="best")
    pl.grid()
    
    
#
# Fit from file
#
  
def fit_FileAtm(pFile):
    aAtm = tl.readTextFileColumn(pFile)
    fit_BurkeNoAero(aAtm)
    
def fit_FileAero(pFile):
    aAtm = tl.readTextFileColumn(pFile)
    fit_atmAeroOnly(aAtm)
   
    
def fit_FileAtmCorAero(pFileAtm, pFileAero):
    aAtm = tl.readTextFileColumn(pFileAtm)    
    aAero = tl.readTextFileColumn(pFileAero)
    aAtm[:,1] /=  aAero[:,1] 
    fit_BurkeNoAero(aAtm) 
   
   
def fit_FileAtmCorAeroN02(pFileAtm, pFileAero, pFileTP7):
    aAtm = tl.readTextFileColumn(pFileAtm)    
    aAero = tl.readTextFileColumn(pFileAero)
    aAtm[:,1] /=  aAero[:,1] 
    fit_BurkeNoAero(aAtm) 
    # to DO : 
    # read TP7 N02 composant
    # interpole to frequence aATm
    # correction NO2 to  aAtm
    
         
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
fit_FileAero(S_fileAero)

lSelect = ["COMBIN","H2O", 'O2','O3','MOLEC']
# lSelect = ["H2Ob","H2O"]
tp7File = "/home/colley/temp/tmp.tp7"
#readMODTRANtp7(tp7File, lSelect, 3000, 9000)

#check_selectedTP7(tp7File)
#coherenceTP7TransFinal(tp7File, S_file, S_fileAero)
#fit_FileAtmCorAero(S_file, S_fileAero)


try:
    pl.show()
except AttributeError:
    pass
