'''
Created on 14 nov. 2012

@author: colley
'''

import numpy as np
import templatMODTRAN as tmod
import pylab as pl

class BurkeAtmModelv1(object):
    '''
    Burke and All model from paper 'precision  determination of atmosphere extinction at optical and NIR wavelengths' 2010
    '''   
#
# PRIVATE
#

    def __init__(self, ModFileTempl, timeObs, pressure=782):
        '''        
        '''
        if ModFileTempl == "": return 
        self._Tpl = tmod.TemplateMODTRAN(ModFileTempl)
        # convert in angstrom
        self._Tpl.convertWaveLength(1.e-10)
        self._PressureRef = pressure*1.0
        self._NbObs = len(timeObs)
        self._TimeObs = timeObs
        self._NbParNoH20 = 7
        # nb paral H20 is egal to nbObs
        self._NbPar = self._NbParNoH20 + self._NbObs 
        self._Par  = None        
        self._aWL = self._Tpl._wl
        self._WL0 = 6750.0
        self._NameParam = ['$T_{gray}$',r'$\tau_0$',r'$\tau_1$', r'$\tau_2$', r'$\alpha$', '$C_{mol}$', '$C_{O3}$','$C_{H2O}$']
        self._AbsH2OConst = 10.0**(-0.4*0.01)
        self._AbsO2Const  = 10.0**(-0.4*0.005)
        self.setConstObsPtgRefComp(0,0)
        
    def resample(self, pWL):
        self._Tpl.resample(pWL)
        self._aWL = self._Tpl._wl  
         
    def _transGrayAero(self):
        tau  = (self._Par[1] + self._EW*self._Par[2] + self._NS*self._Par[3])
        tau *= np.power(self._aWL/self._WL0, self._Par[4])
        return self._Par[0]*np.exp(-self._AirMass*tau)
        
    def _transMols(self):
        return (1 - self._Par[5]*self._PresRat*np.power(self._Tpl._Amols, self._AirMass)) 
    
    def _transMola(self): 
        return (1 - np.sqrt(self._Par[5]*self._PresRat)*self._AbsO2Const*np.power(self._Tpl._Amola, self._AirMass)) 
    
    def _transO3(self):
        return (1 - self._Par[6]*np.power(self._Tpl._A03, self._AirMass))
    
    def _transH2O(self):
        return (1 - self.getC_H20(self._idxTime)*self._AbsH2OConst*np.power(self._Tpl._AH2O, self._AirMass))

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
        assert (len(aIn) == self._NbObs)
        self._Par[self._NbParNoH20:self._NbPar] = aIn
        
    def setParamExample1(self):
        self._Par = np.zeros(self._NbPar, dtype=np.float32)
        # from [1], table 3, 2007, 2nov
        self.setParamAerosolGray(0.98, 3.9/100, 0.02/100, -0.03/100, -1.70)
        # from [1], 4.results , first line
        self.setParamMol(0.91)
        # from [1], table 3, 2007, 2nov
        self.setParamOzone(0.8)
        

    def setParamNoEffect(self):
        """
        to have template MODTRAN with any modification
        """
        self._Par = np.zeros(self._NbPar, dtype=np.float32)
        # from [1], table 3, 2007, 2nov
        self.setParamAerosolGray(0.9, 1.0, 0, 0, -1)
        # from [1], 4.results , first line
        self.setParamMol(1.0)
        # from [1], table 3, 2007, 2nov
        self.setParamOzone(1.0)
        self._Par[self._NbParNoH20:self._NbPar] = np.array([1, 0,0])
        
    
    def setParam(self, par):
        assert (len(par) == self._NbPar)
        self._Par = par
        
    def setConstObs(self, pVecObs):
        return self.setConstObsComp(pVecObs[0], pVecObs[1],pVecObs[2],pVecObs[3])
    
    def setConstObsComp(self, alt, az, t, pressure):
        self._AirMass = np.fabs(1/np.cos(np.pi/2 - alt))
        self._Alt = alt
        self._Az = az
        self._EW = np.cos(alt)*np.sin(az)
        self._NS = np.cos(alt)*np.cos(az)
        self._Time = t
        aIdx = np.where(np.fabs(t-self._TimeObs) < 0.01)[0]
        if len(aIdx) == 0:
            print "No observation at this time ", t
            print self._TimeObs
            return None
        self._idxTime = aIdx[0]
        self._PresRat = pressure*1.0/self._PressureRef        


    def setConstObsPtgRef(self, pCoord):
        """
        alt, az        
        """        
        self.setConstObsPtgRefComp(pCoord[0], pCoord[1])
        
        
    def setConstObsPtgRefComp(self, alt, az):
        self._EWref = np.cos(alt)
        self._NSref = self._EWref*np.cos(az)
        self._EWref *= np.sin(az)
        
# getter

    def computeAtmTrans(self, flagPlot=False):
        if    flagPlot: 
            pl.figure()
            ret = self._transGrayAero()
            pl.plot(ret)
            ret *= self._transMols()
            #pl.plot(ret)
            ret *= self._transMola()
            #pl.plot(ret)
            ret *= self._transO3()
            #pl.plot(ret)
            ret *= self._transH2O()
            #pl.plot(ret)
            #pl.legend(["gray","Mols","Mola","O3","H2O"])
            pl.legend(["gray"])
            pl.title("composant modtran")
            self._CurtTrans = ret

        else:
            ret = self._transGrayAero()
            ret *= self._transMols()
            ret *= self._transMola()
            ret *= self._transO3()
            ret *= self._transH2O()
            self._CurtTrans = ret
        return ret
    
    def computeAtmTransAtConstObs(self, alt, az, t, pressure):
        self.setConstObsComp(alt, az, t, pressure)  
        return self.computeAtmTrans()
    
    def getC_H20(self,IdxObs): 
        return self._Par[self._NbParNoH20+IdxObs]

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
        print "  EW ref   :   ", self._EWref
        print "  SN ref   :   ", self._NSref
        print "  EW       :   ", self._EW
        print "  SN       :   ", self._NS
        print "  pres. rat:   ", self._PresRat
        print "  time     :   ", self._Time
        
    def printBurkeModelParam(self):
        print "Parameters"
        print "  Tgray :   ", self._Par[0]
        print "  Tau0 :   ", self._Par[1]
        print "  Tau1 :   ", self._Par[2]
        print "  Tau2 :   ", self._Par[3]
        print "  alpha :   ", self._Par[4]
        print "  Cmol :   ", self._Par[5]
        print "  C_O3 :   ", self._Par[6]
        print "  C_H2O :   ", self._Par[7:self._NbPar]
                   
    def plotCurrentTrans(self):
        if self._CurtTrans == None:
            print "No current transmission !"
            return None    
        pl.figure()
        pl.plot(self._Tpl._wl, self._CurtTrans)
        pl.xlabel("Angstrom")
        pl.ylabel("%")
        pl.grid()
        pl.title("atm trans at: [alt:%.1f,az:%.1f], pres ratio %.2f, time %.1f"%(np.rad2deg(self._Alt),np.rad2deg(self._Az), self._PresRat, self._Time ))
    
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
        

class BurkeAtmModelv2(BurkeAtmModelv1):
    """
    same thing but for only one instant, add spacial variation C_H20
    """
    
    def __init__(self, ModFileTempl):
        if ModFileTempl == "" : return 
        BurkeAtmModelv1.__init__(self, ModFileTempl, np.array([0]))
        self._NameParam.append('$dC_{H2O}/dEW$')
        self._NameParam.append('$dC_{H2O}/dNS$')
        self._NbPar = 10
                
    def setParamExample1(self):
        BurkeAtmModelv1.setParamExample1(self)
        self._Par[8] = 0.166
        self._Par[9] = -0.2
                        
    def setConstObsComp(self, alt, az, t, pressure):
        #print "setConstObsComp v2"
        self._AirMass = np.fabs(1/np.cos(np.pi/2 - alt))
        self._Alt = alt
        self._Az = az
        self._EW = np.cos(alt)
        self._NS = self._EW*np.cos(az)
        self._EW *= np.sin(az)        
        self._Time = t
        self._idxTime = None
        self._PresRat = pressure*1.0/self._PressureRef
        self._DeltaEW = self._EW - self._EWref
        self._DelatNS = self._NS - self._NSref
           
    def _transH2O(self):
        C_H20 = self._Par[7] + self._Par[8]*self._DeltaEW + self._Par[9]*self._DelatNS
        return (1 - C_H20*self._AbsH2OConst*np.power(self._Tpl._AH2O, self._AirMass))
    
    def getC_H20(self,IdxObs): 
        return None
