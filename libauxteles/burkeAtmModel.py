'''
Created on 14 nov. 2012

@author: colley
'''

import numpy as np
import templatMODTRAN as tmod
import pylab as pl

    
class BurkeAtmModel(object):
    """
        Tgray :   self._Par[0]
        Tau0 :   self._Par[1]
        Tau1 :   self._Par[2]
        Tau2 :   self._Par[3]
        alpha :   self._Par[4]
        Cmol :   self._Par[5]
        C_O3 :   self._Par[6]
        C_H2O :   self._Par[7]
        dC_{H2O}/dEW :   self._Par[8]
        dC_{H2O}/dNS : self._Par[9]
    """
    
    def _razNeg(self, pAr):
        idx = np.where(pAr < 0.0)[0]
#        if len(idx)  > 0:
#            print "======================================================================Some negative value ", idx  
#            pl.figure()
#            pl.plot(pAr)
#            pAr[idx] = 0
#            pl.plot(pAr,'*')
#            pl.ploself._AbsH2OConst*t(pAr)
#            pl.grid()                    
        pAr[idx] = 0
        return pAr
        
    def __init__(self, ModFileTempl, pressure=782):
        '''        
        '''
        if ModFileTempl == "": return 
        self._Tpl = tmod.TemplateMODTRAN(ModFileTempl)
        # convert in angstrom
        self._Tpl.convertWaveLength(1.e-10)
        self._PressureRef = pressure*1.0
        # nb paral H20 is egal to nbObs
        self._NbPar = 10 
        self._Par  = None        
        self._aWL = self._Tpl._wl
        self._WL0 = 6750.0
        self._NameParam = ['$T_{gray}$',r'$\tau_0$',r'$\tau_1$', r'$\tau_2$', r'$\alpha$', '$C_{mol}$', '$C_{O3}$','$C_{H2O}$']
        self._NameParam.append('$dC_{H2O}/dEW$')
        self._NameParam.append('$dC_{H2O}/dNS$')
        self._AbsH2OConst = 10.0**(-0.4*0.01)
        self._AbsO2Const  = 10.0**(-0.4*0.005)
    
    def getWL(self):
        return self._aWL 
    
    def resample(self, pWL):
        self._Tpl.resample(pWL)
        self._aWL = self._Tpl._wl
          
    def _transGrayAero(self):
        tau  = self._Par[1] + self._EW*self._Par[2] + self._NS*self._Par[3]
        #print 'rap:', self._aWL[0]/self._WL0, self._aWL[-1]/self._WL0
        tau *= np.power(self._aWL/self._WL0, self._Par[4])        
        return self._razNeg(self._Par[0]*np.exp(-self._AirMass*tau))
        
    def _transMols(self):
        #return self._razNeg(1.0 - self._Par[5]*self._PresRat*np.power(self._Tpl._Amols, self._AirMass))
        return self._razNeg(1.0 - self._Par[5]*self._PresRat*(1.0 - np.power(self._Tpl.getTrmols(), self._AirMass)))
    
    def _transMola(self): 
        #return self._razNeg(1.0 - np.sqrt(self._Par[5]*self._PresRat)*self._AbsO2Const*np.power(self._Tpl._Amola, self._AirMass))
        return self._razNeg(1.0 - np.sqrt(self._Par[5]*self._PresRat)*(1.0 - self._AbsO2Const*np.power(self._Tpl.getTrmola(), self._AirMass)))
     
    def _transO3(self):
        #return self._razNeg(1.0 - self._Par[6]*np.power(self._Tpl._A03, self._AirMass))
        return self._razNeg(1.0 - self._Par[6]*(1.0 - np.power(self._Tpl.getTr03(), self._AirMass)))
        
    def _transH2O(self):
        C_H20 = self._Par[7] + self._EW*self._Par[8] + self._NS*self._Par[9]
        #return self._razNeg(1.0 - C_H20*self._AbsH2OConst*np.power(self._Tpl._AH2O, self._AirMass))
        return self._razNeg(1.0 - C_H20*(1.0 - self._AbsH2OConst*np.power(self._Tpl.getTrH2O(), self._AirMass)))
    
#       
# PUBLIC
#

# setter

    def setParamAerosolGray(self, Tgray, tau0, tau1, tau2, alpha):
        self._Par[0] = Tgray
        self._Par[1] = tau0
        self._Par[2] = tau1
        self._Par[3] = tau2
        self._Par[4] = alpha
               
    def setParamMol(self, Cmol):
        self._Par[5] = Cmol
    
    def setParamOzone(self, O3):
        self._Par[6] = O3
        
    def setParamH20(self, aIn):        
        assert (len(aIn) == 3)
        self._Par[7:10] = aIn
        
    def setParamExample1(self):
        self._Par = np.zeros(self._NbPar, dtype=np.float32)
        # from [1], table 3, 2007, 2nov
        self.setParamAerosolGray(0.98, 3.9/100, 0.02/100, -0.03/100, -1.70)
        # from [1], 4.results , first line
        self.setParamMol(0.91)
        # from [1], table 3, 2007, 2nov
        self.setParamOzone(0.8)
        # 
        self.setParamH20(np.array([1, 0.01, -0.05]))
        
    def setParamNoEffect(self):
        """
        to have template MODTRAN with any modification
        """
        self._Par = np.zeros(self._NbPar, dtype=np.float32)
        # from [1], table 3, 2007, 2nov
        self.setParamAerosolGray(1.0, 0.0, 0, 0, -1)
        # from [1], 4.results , first line
        self.setParamMol(1.0)
        # from [1], table 3, 2007, 2nov
        self.setParamOzone(1.0)
        #
        self.setParamH20(np.array([1, 0, 0]))
        
    def setParam(self, par):
        assert (len(par) == self._NbPar)
        self._Par = par
        
    def setConstObs(self, pVecObs):
        """        
        alt , az, pressure
        alt : 0 Hori, pi/2 zenith
        """
        return self.setConstObsComp(pVecObs[0], pVecObs[1], pVecObs[2])
    
    def setConstObsComp(self, alt, az, pressure):        
        self._AirMass = np.fabs(1/np.cos(np.pi/2 - alt))        
        self._Alt = alt
        self._Az = az
        self._EW = np.cos(alt)*np.sin(az)  
        self._NS = np.cos(alt)*np.cos(az)            
        self._PresRat = pressure*1.0/self._PressureRef
        #print "PresRat: ",self._PresRat, pressure
        
# getter
    def computeAtmTrans(self, flagPlot=False):
        #print "AirMass:",self._AirMass
        if flagPlot: 
            pl.figure()
            ret = self._transGrayAero()
            pl.plot(self._aWL,ret,'r--')
            ret1 = self._transMols()
            pl.plot(self._aWL,ret1,'y--.')
            ret2= self._transMola()
            pl.plot(self._aWL,ret2,'b:')
            ret3= self._transO3()
            pl.plot(self._aWL,ret3,'y--')
            ret4= self._transH2O()
            pl.plot(self._aWL,ret4,'r--.')
            pl.plot(self._aWL,ret*ret1*ret2*ret3*ret4)
            pl.legend(["gray-aero","Mols","Mola","O3","H2O","all"],loc=3)
            #pl.legend(["gray"])
            pl.title("composant modtran")
            pl.grid()
            self._CurtTrans = ret
        else:
            ret = self._transGrayAero()
            ret *= self._transMols()
            ret *= self._transMola()
            ret *= self._transO3()
            ret *= self._transH2O()
            self._CurtTrans = ret        
        return ret
    
    def computeAtmTransAtConstObs(self, alt, az, pressure):
        self.setConstObsComp(alt, az, pressure)  
        return self.computeAtmTrans()
    
    
# pretty print, plot

    def printAndPlotBurkeModel(self):
        if self.printBurkeModel() == None: return None
        pl.figure()
        pl.plot(self._TimeObs, self._Par[self._NbParNoH20:self._NbPar],"*")
        pl.xlabel("s")
        pl.grid()
        pl.title("C_H2O")
        self._Tpl.plotTemplate()
                
    def printBurkeModel(self):
        if self._Par == None:
            print "No parameter yet defined !"
            return None
        self.printBurkeModelConst()
        self.printBurkeModelParam()
        self.printBurkeObserConst()
        
    def printBurkeModelConst(self):
        print "Constants"
        print "  ref. pressure:  %d  mb"% self._PressureRef
        print "  ref. wavelength:  %g  A"% self._WL0
        
    def printBurkeObserConst(self):
        print "Constants observation"
        print "  EW       :   ", self._EW
        print "  SN       :   ", self._NS
        print "  pres. rat:   ", self._PresRat
               
    def printBurkeModelParam(self):
        print "Parameters"
        print "  Tgray :   ", self._Par[0]
        print "  Tau0 :   ", self._Par[1]
        print "  Tau1 :   ", self._Par[2]
        print "  Tau2 :   ", self._Par[3]
        print "  alpha :   ", self._Par[4]
        print "  Cmol :   ", self._Par[5]
        print "  C_O3 :   ", self._Par[6]
        print "  C_H2O :   ", self._Par[7:]
                   
    def plotCurrentTrans(self):
        if self._CurtTrans == None:
            print "No current transmission !"
            return None    
        pl.figure()
        pl.plot(self._Tpl._wl, self._CurtTrans)
        pl.xlabel("Angstrom")
        pl.ylabel("%")
        pl.grid()
        pl.title("atm trans at: [alt:%.1f,az:%.1f], pres ratio %.2f"%(np.rad2deg(self._Alt),np.rad2deg(self._Az), self._PresRat ))
    
    def plotCovarMatrix(self, mat, pTitle=''):
        #pl.figure()        
        #pl.pcolor(mat)
        #pl.matshow(mat, cmap=pl.cm.gray)  
        im = pl.matshow(mat)  
        pl.title(pTitle)  
        aIdx = np.arange(len(self._NameParam))
        pl.xticks(aIdx,  self._NameParam)
        for label in im.axes.xaxis.get_ticklabels():
            label.set_rotation(45)
        #pl.yticks(aIdx, self._NameParam)
        pl.colorbar()
    
    def _covar2Correl(self, mat):
        matTp = np.copy(mat)
        #print matTp
        diag = np.sqrt(matTp.diagonal())
        A = np.outer(np.ones(mat.shape[0]), diag)
        #print A
        res = matTp/(A*A.transpose())
        #print res
        return res
    
    def plotCorrelMatrix(self, mat, pTitle=''):
        #pl.figure()        
        #pl.pcolor(mat)
        #pl.matshow(mat, cmap=pl.cm.gray)        
        im = pl.matshow(self._covar2Correl(mat))  
        pl.title(pTitle)  
        aIdx = np.arange(len(self._NameParam))
        pl.xticks(aIdx,  self._NameParam)
        for label in im.axes.xaxis.get_ticklabels():
            label.set_rotation(45)
        pl.yticks(aIdx, self._NameParam)
        pl.colorbar()

    def plotErrRel(self, pTrue, pEst, pTitle=""):    
        pl.figure()       
        pl.title(pTitle)
        Err = 100*(pTrue- pEst)/pTrue
        pl.plot(Err)
        pl.ylabel("%")
        aIdx = np.arange(len(self._NameParam))
        pl.xticks(aIdx,  self._NameParam)
        pl.grid()



class BurkeAtmModelTauPos(BurkeAtmModel):
    """
    fit sqrt(Tau0) to assume Tau0 positive
    
        Tgray :   self._Par[0]
        sqrt(Tau0) :   self._Par[1]
        Tau1 :   self._Par[2]
        Tau2 :   self._Par[3]
        alpha :   self._Par[4]
        Cmol :   self._Par[5]
        C_O3 :   self._Par[6]
        C_H2O :   self._Par[7]
        dC_{H2O}/dEW :   self._Par[8]
        dC_{H2O}/dNS : self._Par[9]
    """
    
    def __init__(self, ModFileTempl, pressure=782):
        BurkeAtmModel.__init__(self, ModFileTempl, pressure)
        
    def _transGrayAero(self):
        tau  = self._Par[1]**2 + self._EW*self._Par[2] + self._NS*self._Par[3]
        tau *= np.power(self._aWL/self._WL0, self._Par[4])        
        return self._razNeg(self._Par[0]*np.exp(-self._AirMass*tau))
# setter

    def setParamAerosolGray(self, Tgray, tau0, tau1, tau2, alpha):
        BurkeAtmModel.setParamAerosolGray(Tgray, tau0, tau1, tau2, alpha)        
        self._Par[1] = np.sqrt(tau0)
    
    def printBurkeModelParam(self):
        print "Parameters"
        print "  Tgray :   ", self._Par[0]
        print "  Tau0 :   ", self._Par[1]**2
        print "  Tau1 :   ", self._Par[2]
        print "  Tau2 :   ", self._Par[3]
        print "  alpha :   ", self._Par[4]
        print "  Cmol :   ", self._Par[5]
        print "  C_O3 :   ", self._Par[6]
        print "  C_H2O :   ", self._Par[7:]
