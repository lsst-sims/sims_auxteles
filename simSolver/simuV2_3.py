'''
Created on 30 janv. 2014

@author: colley
'''

import tools as tl
import matplotlib.pyplot as pl
import burkeAtmModel as bam
import numpy as np
import scipy.optimize as spo


S_file = "../data/simuRef/atmAeroDiffZangle.txt"
G_fileModtran = '../data/modtran/TemplateT04.01_1.txt'



def test_atmAeroPlot(pFile):
    aAtm = tl.readTextFileColumn(pFile)
    print aAtm
    pl.figure()
    pl.plot(aAtm[:,0], aAtm[:,1])
    pl.grid()
    
    
def test_wl_atmAero(pFile):
    aAtm = tl.readTextFileColumn(pFile)
    print aAtm    
    atmTheo = bam.BurkeAtmModel(G_fileModtran)
    wlTH = atmTheo.getWL()/10.
    wlMes = aAtm[:,0]    
    print wlTH, wlMes
    diff = wlTH - wlMes
    atmTheo.setDefaultConstObs()
    atmTheo.setDefaultParam()
    defAtm = atmTheo.computeAtmTrans(True)
    pl.figure()
    pl.plot(wlMes, defAtm)
    pl.plot(wlMes, aAtm[:,1])

    
def fit_atmAero(pFile):
    idx2Full = [1,4,5,6,7]
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
    wl = aAtm[:,0]
    transMeas = aAtm[:,1]
    atmTheo = bam.BurkeAtmModel(G_fileModtran)
    atmTheo.setConstObsComp(np.pi/2, 0, 782)
    # parameter: tau0, alpha, Cmol, C_O3, C_H20
    p0 = np.array([1.0, 0.05, 1.0, 1.0, 1.0], dtype=np.float64)
    res = spo.leastsq(errorModel, p0, args=(atmTheo, transMeas), full_output=True)
    print res
    atmTheo.printBurkeModelParam()
    if res[4] >= 0 and res[4]<=4:
        print "FIT OK :",  res[3]
        print res[0]
        atmTheo.setParam(getFullParam(res[0]))
        atmTheo.printBurkeModelParam()
        if res[1] != None:             
            pl.figure()
            pl.pcolor(res[1][::-1])
            pl.colorbar()
        pl.figure()
        pl.plot(wl, transMeas)
        pl.plot(wl, atmTheo.computeAtmTrans())        
    else:
        print "FIT NOK : ",  res[3]
        print res[0]
    
    
    
    
    
    
    
    
#
# MAIN
#

fit_atmAero(S_file)


try:
    pl.show()
except AttributeError:
    pass