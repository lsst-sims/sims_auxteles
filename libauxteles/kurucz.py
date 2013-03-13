'''
Created on 13 nov. 2012

@author: colley
'''

import pickle as pk
import numpy as np
import math
import pyfits as pf
import pylab as pl
import scipy.interpolate as sci
import scipy.optimize as spo
import lmfit as lm



def read_kurucz_all(directory = '/home/colley/projet/lsst/stellar_spectra/k93models/'):
    metal_dir = ['km50', 'km45', 'km40', 'km35', 'km30', 'km25', 'km20', 'km15', 'km10', 'km05', 'km03', 'km02', 'km01', 'kp00', 'kp01', 'kp02', 'kp03', 'kp05', 'kp10']
    step_metal = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, +0.0, +0.1, +0.2, +0.3, +0.5, +1.0]
    step_temperature = [3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 37500, 40000, 42500, 45000, 47500, 50000]
    step_gravity = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

    nb_metal = 19
    nb_temp = 61
    nb_logg = 11
    nb_lambda = 1221

    nb_col = 12684
    nb_lig = 1224

    data = np.zeros((nb_lig, nb_col))
    data[0,0] = -99
    data[1,0] = -99
    data[2,0] = -99

    indice_col = 1
    for m in range(0, nb_metal):
        # cas m = +0.5
        if (m == nb_metal-2):
            nb_temp = 59
        # cas m = +1.0
        if (m == nb_metal-1):
            nb_temp = 57
        for t in range(0, nb_temp):
            nom_file_kurucz = directory + "%s/%s_%s.fits" % (metal_dir[m], metal_dir[m], step_temperature[t])
            for g in range(0, nb_logg):
                loggravity_name = "g%02i" %(int(math.floor(step_gravity[g]*10.)))

                data_k93 = pf.getdata(nom_file_kurucz)
                print "read ",nom_file_kurucz
                wavelength = data_k93.field('wavelength')
                flux = data_k93.field(loggravity_name)

                data[3:1224,0] = wavelength # * 10. # To get wavelength in A
                data[0,indice_col] = step_metal[m]
                data[1,indice_col] = step_temperature[t]
                data[2,indice_col] = step_gravity[g]
                data[3:1224,indice_col] = flux
                
                indice_col += 1

                #print m, t, g, data[0:5,0:5]
    return data


def fits2pickle(fIn, fOut):
    kur = read_kurucz_all(fIn)
    f=open(fOut, "wb")
    pk.dump(kur, f, pk.HIGHEST_PROTOCOL)
    f.close()


def fits2pickle_JMC():
    fIn = '/home/colley/projet/lsst/stellar_spectra/k93models/'
    fOut = '/home/colley/projet/lsst/stellar_spectra/k93.pic'
    fits2pickle(fIn, fOut)
    
    
    
    
class Kurucz(object):
    '''    
    Format Kurucz array :
        column 0 : wavelength angstrom
        line 0 : metallicity log_Z : -5.0 to 1.0
        line 1 : temperature Kelvin: 3000 to 50000
        line 2 : gravity log_g : 0 to 5
    '''

    def __init__(self, filePickle, test=False):
        try:
            f=open(filePickle, "rb")
        except:
            print "can't open file ", filePickle
            return 
        if test:
            print "add _ParamNotDef attribut"
            self._ParamNotDef = []    
        self._Flux = pk.load(f)
        self._oNGP = None
        self._oInterLin = None
        #print self._Flux
        f.close()
        self.setWLuseAll()
        # may be useful for linear interpolated feature  
        self._deleteFluxNotDefined()
        self._CoefUnit = 1
        #self.setCoefUnit(1e-7)
        self._BoundsMet = (-4, 0)
        self._BoundsGra = (0.5, 4.5)
        self._BoundsTemp = (3600.0, 47500.0)
        
        
    def setCoefUnit(self, coef):
        self._CoefUnit *= coef
        self._Flux[3:,1:] *= coef        
        
    def resampleBetween(self, pWLmin, pWLmax, pNb):
        newWL = np.linspace(pWLmin, pWLmax, pNb, True)
        self.resample(newWL)
        
    def resample(self, pWL):
        # interpole only around new WL domain 
        self.setWLInterval(pWL[0]*0.98, pWL[-1]*1.02)
        WLin = self.getWL()
        sizepWL = len(pWL)
        for idx in np.arange(1, len(self._Flux[0,:])): 
            # Bspline fit - interpole                
            tck = sci.splrep(WLin, self.getFluxIdx(idx))
            self._Flux[3:sizepWL+3,idx] = sci.splev(pWL, tck)
        self._Flux[3:sizepWL+3, 0] = pWL 
        # remove obsolete (WL, Flux)     
        idxRemove = np.arange(sizepWL+3, len(self._Flux[0,:]))
        self._Flux = np.delete(self._Flux, idxRemove, 0)
        # use all wl domain and raz interpole object
        self.setWLuseAll()
        self._oInterLin = None
        self._oNGP = None
                
    def _deleteFluxNotDefined(self):
        """
        star flux is not defined everywhere but Kurucz set array flux at zero in this case
        With linear interpolation this zero array may be distort result  ... ?
        so erase all zero array 
        """
        aColNOK = []
        if hasattr(self, '_ParamNotDef'):
            NotDefined = True
            # to visualize where Kurucz model isn't defined, see test
            self._ParamNotDef = []    
        else:
            NotDefined = False
        for idx in np.arange(1, len(self._Flux[0,:])):            
            if self._Flux[3:,idx].sum() == 0.0:
                #print "%d not defined"%idx
                if NotDefined :                   
                    par = self.getParam(idx)
                    self._ParamNotDef.append(par[0])
                    self._ParamNotDef.append(par[1])
                    self._ParamNotDef.append(par[2])
                aColNOK.append(idx)
        if aColNOK != []:
            print "delete %d columns not defined "%len(aColNOK)
            self._Flux = np.delete(self._Flux, aColNOK, 1)
            if NotDefined :                
                self._ParamNotDef = np.array(self._ParamNotDef).reshape(len(aColNOK),3)
                
    def setWLuseAll(self):
        self._IdxMin = 3
        self._IdxMax = len(self._Flux[0])-1
    
    def setWLInterval(self, wlMin , wlMax):   
        idx = np.where(self._Flux[3:,0] >= wlMin)[0]
        if len(idx) == 0:
            print "ERROR [setWLInterval]: wlMin > WL max"
            self.setWLuseAll()
            return 
        self._IdxMin = idx[0]+3
        idx = np.where(self._Flux[3:,0] <= wlMax)[0]
        if len(idx) == 0:
            print "ERROR [setWLInterval]: wlMax < WL min"
            self.setWLuseAll()
            return 
        self._IdxMax = idx[-1]+3
        print wlMin , wlMax
        print "nb wl ", self._IdxMax-self._IdxMin
    
    def getFluxIdx(self, idx):
        return self._Flux[self._IdxMin:self._IdxMax+1, idx]
    
    def getWL(self):
        return self._Flux[self._IdxMin:self._IdxMax+1, 0]
    
    def getParam(self, idx):
        """
        return metallicity, temperature, gravity
        """
        return self._Flux[0:3, idx]
     
    def plotFluxIdx(self,idx):
        pl.figure()
        pl.plot(self.getWL(),self.getFluxIdx(idx))
        pl.xlabel("Angstrom")
        pl.grid()
        pl.title("Kurucz 93 star flux M %.2f T %.2f G %.2f"%(self._Flux[0,idx], self._Flux[1,idx], self._Flux[2,idx]))
                
    def plotMultiFluxesCont(self,idx0, nb):
        pl.figure()
        for idx in range(nb):
            pl.plot(self.getWL(),self.getFluxIdx(idx+idx0))
        strLgd = []
        for idxR in range(nb):
            idx = idxR + idx0
            strLgd.append("M %.2f T %.2f G %.2f"%(self._Flux[0,idx], self._Flux[1,idx], self._Flux[2,idx]))
        pl.legend(strLgd)
        pl.xlabel("Angstrom")
        pl.grid()
        pl.title("Kurucz 93 star fluxes")

    def plotFluxesArrayIdx(self,aIdx):
        pl.figure()
        strLgd = []
        for idx in aIdx:
            pl.plot(self.getWL(),self.getFluxIdx(idx))
            strLgd.append("M %.1f T %.0f G %.2f"%(self._Flux[0,idx], self._Flux[1,idx], self._Flux[2,idx]))
        pl.legend(strLgd)
        pl.xlabel("Angstrom")
        pl.grid()
        pl.title("Kurucz 93 star fluxes")
                
    def getFluxNGP(self, pPar):
        """
        pPar : Met, Temp, Gra
        return nearest grid point flux for given parameters pMet, pTemp, pGra
        """
        idx = self.getIdxNGP(pPar)
        return self.getFluxIdx(idx)
       
    def getIdxNGP(self, pPar):
        """
        pPar : Met, Temp, Gra
        return idx nearest grid point for given parameters pMet, pTemp, pGra
        """        
        if self._oNGP == None:
            a = self._Flux[0:3,1:].transpose()
            self._oNGP = sci.NearestNDInterpolator(a, np.arange(1, len(self._Flux[0]-1)))
        idx = self._oNGP(np.array([pPar]))
        #print idx, self._Flux[0:3,idx]
        return idx
        
    def getFluxInterLin(self, pPar):
        """
        pPar : Met, Temp, Gra
        linear interpolation of flux for given parameters pMet, pTemp, pGra
        """   
        if self._oInterLin == None:
            param = self._Flux[0:3,1:].transpose()
            flux = self._Flux[self._IdxMin: self._IdxMax+1,1:].transpose()      
            self._oInterLin = sci.LinearNDInterpolator(param, flux)
            print "fin LinearNDInterpolator"
        res = self._oInterLin(np.array([pPar]))
        return res.ravel()

    def _fitLeastsqLinAll(self, pFlux, pIter=1, guess=None, pLogT = False):
        """
        fit NOK : probleme avec le domaine de definition metal. et grav.
        2 solutions: 
            ==> solution possible utiliser une methode optimisation sous contrainte.
            ==> fit alternee sur temp puis sur (met, grav) OK  , fitLeastsqLinAll3Step
        """
        return None
        assert (len(pFlux)== len(self.getWL()))
        if guess == None:
            guess = np.array([-2.5, 8000., 2.5])
        if pLogT: 
            guess[1] = np.log10(guess[1])
        def errorModel(par, kur, ydata):
            parI =  par.copy()
            delta = 1
            flag = True
            if par[2] > 4:                
                delta = np.fabs(par[2] - 4)
                flag = False
            elif par[2] <1:                
                delta = np.fabs(par[2] - 1)
                flag = False
            if par[0] < -4:                
                delta *= np.fabs(par[0] +4)
                flag = False
            elif par[0] > 0.5:                
                delta *= np.fabs(par[0] -0.5)
                flag = False
            #print par 
            if pLogT: parI[1]= 10**parI[1]  
            if flag:                     
                res = (ydata - kur.getFluxInterLin(parI))*1e-7
            else:
                print "NGP penalisation"                
                res = (ydata - kur.getFluxNGP(parI))*1e-3
            print par, (res**2).sum()
            return res.ravel()
        # loop iteration
        for idxIter in range(pIter):
            res = spo.leastsq(errorModel, guess, args=(self, pFlux), full_output=True)
            if res[4] >= 0 and res[4]<=4:
                print "FIT %d OK :"%(idxIter+1) + str(res[0])
                guess = res[0]
            else:
                print res
                print "FIT NOK : ",  res[3]
                return None
        ret = res[0]
        if pLogT: ret[1]= 10**ret[1]
        return ret
            
    def fitLeastsqLinTemp(self, pFlux, guess=None):        
        """
        guess : Met, Temp, Gra
        """
        assert (len(pFlux)== len(self.getWL()))
        if guess == None:
            guess = np.array([-2.5, 8000., 2.5])
        def errorModel(par, kur, ydata):
            parI = np.copy(guess)
            parI[1] = par                                  
            res = (ydata - kur.getFluxInterLin(parI))*1e-7
            print parI, (res**2).sum()
            return res.ravel()        
        res = spo.leastsq(errorModel, guess[1], args=(self, pFlux), full_output=True)
        if res[4] >= 0 and res[4]<=4:
            print "FIT OK :"+ str(res[0])
        else:
            print res
            print "FIT NOK : ",  res[3]
            return None
        ret = np.copy(guess)
        ret[1] = res[0]        
        return ret
    
    def fitFminLinTemp(self, pFlux, guess=None):
        """
        guess : Met, Temp, Gra
        """
        assert (len(pFlux)== len(self.getWL()))
        if guess == None:
            guess = np.array([-2.5, 8000., 2.5])
        def errorModel(par, kur, ydata):
            parI = np.copy(guess)
            parI[1] = par                                  
            res = ((ydata - kur.getFluxInterLin(parI))*1e-7)**2
            print parI, (res**2).sum()
            return res.ravel().sum()       
        res = spo.fmin(errorModel, guess[1], args=(self, pFlux))
        ret = np.copy(guess)
        ret[1] = res[0]        
        return ret
    
    def fitWithBounds(self, pFlux, guess=None):
        """
        guess : Met, Temp, Gra
        """
        assert (len(pFlux)== len(self.getWL()))
        if guess == None:
            guess = np.array([-2.5, 8000., 2.5])
        def errorModel(par, kur, ydata):
            fluxTheo = kur.getFluxInterLin(par) 
            if np.isnan(fluxTheo[0]):
                print   "fluxTheo is nan, used NGP "
                fluxTheo = kur.getFluxNGP(par)                          
            res = ((ydata - fluxTheo))**2
            res = res.ravel().sum()
            print par, res
            return res
        myBounds = (self._BoundsMet, self._BoundsTemp, self._BoundsGra)
        # avec  test_fitWithBounds et    np.random.seed(10) 
        #res = spo.minimize(errorModel, guess, args=(self, pFlux), method='SLSQP', bounds=myBounds)
        # donne erreur !
        #res = spo.minimize(errorModel, guess, args=(self, pFlux), method='SLSQP', bounds=myBounds)
        res = spo.minimize(errorModel, guess, args=(self, pFlux), method='SLSQP', bounds=myBounds)
        print res.message
        return res.x
    
    def fitWithBounds2(self, pFlux, guess=None):
        """
        guess : Met, Temp, Gra
        """
        assert (len(pFlux)== len(self.getWL()))
        if guess == None:
            guess = np.array([-2.5, 8000., 2.5])
        def errorModel(par, kur, ydata):
            print "input par", par
            fluxTheo = kur.getFluxInterLin(par) 
            if np.isnan(fluxTheo[0]):
                print   "fluxTheo is nan, used NGP "
                fluxTheo = kur.getFluxNGP(par)                          
            res = ((ydata - fluxTheo))**2
            res = res.ravel().sum()
            print res
            return res
        myBounds = (self._BoundsMet, self._BoundsTemp, self._BoundsGra)
        # avec  test_fitWithBounds et    np.random.seed(10) 
        #res = spo.minimize(errorModel, guess, args=(self, pFlux), method='SLSQP', bounds=myBounds)
        # donne erreur !
        #res = spo.minimize(errorModel, guess, args=(self, pFlux), method='SLSQP', bounds=myBounds)
        res = spo.minimize(errorModel, guess, args=(self, pFlux), method='SLSQP', bounds=myBounds)
        print "sol1:", res.message, res.x
        #res2 = spo.minimize(errorModel,res.x , args=(self, pFlux), method='SLSQP', bounds=myBounds)
        #print "sol2:", res2.message, res2.x
        return res.x
    
    def fitWithBoundsLM(self, pFlux, guess=None):
        """
        guess : Met, Temp, Gra
        """
        assert (len(pFlux)== len(self.getWL()))
        if guess == None:
            guess = np.array([-2.5, 8000., 2.5])
        def errorModel(params, kur, ydata):           
            v=[]
            [v.append(params[a].value) for a in params]
            par = np.array(v)
            print "input par", par
            fluxTheo = kur.getFluxInterLin(par) 
            if np.isnan(fluxTheo[0]):
                print   "fluxTheo is nan, used NGP "
                fluxTheo = kur.getFluxNGP(par)                          
            res = ((ydata - fluxTheo))**2
            res = res.ravel()
            #print res
            return res
        par= lm.Parameters()
        par.add('met', value=guess[0], min=-4, max= 0.0)
        par.add('temp', value=guess[1])
        par.add('gra', value=guess[2], min=0.5, max= 5.0)       
        #myBounds = (self._BoundsMet, self._BoundsTemp, self._BoundsGra)
        res = lm.minimize(errorModel, par, args=(self, pFlux))
        lm.report_errors(par)
        #res2 = spo.minimize(errorModel,res.x , args=(self, pFlux), method='SLSQP', bounds=myBounds)
        #print "sol2:", res2.message, res2.x
        sol = []
        [sol.append(par[a].value) for a in par]
        return np.array(sol)
    
    def fitNoBounds(self, pFlux, guess=None):
        """
        guess : Met, Temp, Gra
        """
        assert (len(pFlux)== len(self.getWL()))
        if guess == None:
            guess = np.array([-2.5, 8000., 2.5])
        def errorModel(par, kur, ydata):                            
            res = ((ydata - kur.getFluxInterLin(par)))**2
            res = res.ravel().sum()
            print par, res
            return res        
        res = spo.minimize(errorModel, guess, args=(self, pFlux), method='SLSQP')
        print res.message
        return res.x
    
    def fitLeastsqLinAll3Step(self, pFlux, guess=None):
        """
        alternated fit: temp => (met, gra) => temp
        """
        # Step 1: fit on temp
        sol = self.fitLeastsqLinTemp(pFlux, guess)
        # Step 2: fit gra, met        
        def errorModel(par, kur, ydata):
            parI = np.copy(sol)
            parI[0] = par[0]
            parI[2] = par[1]                                  
            res = (ydata - kur.getFluxInterLin(parI))*1e-7
            print parI, (res**2).sum()
            return res.ravel()            
        guess =  np.delete(sol, 1)  
        print "guess:", guess
        res = spo.leastsq(errorModel, guess, args=(self, pFlux), full_output=True)
        if res[4] >= 0 and res[4]<=4:
            print "FIT OK :"+ str(res[0])
        else:
            print res
            print "FIT NOK : ",  res[3]
            return None
        ret = np.copy(sol)
        ret[0] = res[0][0]
        ret[2] = res[0][1]
        # Step 3: fit temp again
        sol = self.fitLeastsqLinTemp(pFlux, ret) 
        return sol
        
        
