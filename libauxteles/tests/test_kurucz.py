'''
Created on 14 nov. 2012

@author: colley
'''
import unittest
from kurucz import *
import scipy.optimize as spo


FileKuruczPic = '/home/colley/projet/lsst/stellar_spectra/k93.pic'

def test_plotFlux():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(1000, 12000)
    oKur.plotFluxIdx(1540)
    oKur.plotMultiFluxesCont(1540,5)
    oKur.plotFluxIdx(1544)
    
    
def test_setWLInterval():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(3000, 9000)
    print oKur._IdxMin,  oKur._Flux[oKur._IdxMin-1,0], oKur._Flux[oKur._IdxMin,0], oKur._Flux[oKur._IdxMin+1,0]
    print oKur._IdxMax,  oKur._Flux[oKur._IdxMax-1,0], oKur._Flux[oKur._IdxMax,0], oKur._Flux[oKur._IdxMax+1,0]    
    oKur.setWLInterval(3030, 9010)
    print oKur._IdxMin,  oKur._Flux[oKur._IdxMin-1,0], oKur._Flux[oKur._IdxMin,0], oKur._Flux[oKur._IdxMin+1,0]
    print oKur._IdxMax,  oKur._Flux[oKur._IdxMax-1,0], oKur._Flux[oKur._IdxMax,0], oKur._Flux[oKur._IdxMax+1,0]    
     
def test_resample():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(4500, 10500)
    par = np.array([0.0, 3660., 0.5])
    pl.figure()
    pl.plot(oKur.getWL(), oKur.getFluxInterLin(par))
    oKur.resampleBetween(5000, 10000, 512)
    oKur.setWLInterval(5000, 10000)
    pl.plot(oKur.getWL(), oKur.getFluxInterLin(par))
    pl.legend(["raw",'resample'])
    pl.grid()
    
    
def test_getFluxInterLin():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(2000, 10000)
    pl.figure()
    fl = oKur.getFluxInterLin(np.array([-1.3, 6000, 2.7]))
    pl.plot(oKur.getWL(), fl.T)
    fl = oKur.getFluxInterLin(np.array([-1.3, 6242, 2.7]))
    pl.plot(oKur.getWL(), fl.T)
    fl = oKur.getFluxInterLin(np.array([-1.3, 6300, 2.7]))
    pl.plot(oKur.getWL(), fl.T)
    pl.legend(['6000','6242','6300'])
    pl.grid()
    

def test_fitTempGra():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(2500, 10000)        
    fluxTrue =oKur.getFluxInterLin(np.array([-2, 7125, 3.7]))
    flux = (1+np.random.normal(0, 0.01,len(fluxTrue)))*fluxTrue
    metUsed= -2.
    def errorModel(par, kur, ydata):
        parI = np.array([metUsed, 0. ,0.])
        parI[1:3] = par
        res = (ydata - oKur.getFluxInterLin(parI))*1e-7
        print parI, (res**2).sum()
        return res.ravel()
    print "True parameter"
    print "-2, 6100, 2.7"
    print "\nguess parameter"
    p0 = np.array([ 6100., 2.7])
    p0 += np.random.normal(0, 0.3,2)*p0
    print p0
    res = spo.leastsq(errorModel, p0, args=(oKur, flux), full_output=True)
    if res[4] >= 0 and res[4]<=4:
        print "FIT OK :",  res[3]
        #print res[0], res[2]["nfev"]
        print "\nSolution scipy.leastsq in %d calls "%res[2]["nfev"]
        print res[0]            
        if res[1] != None:                             
            print res[1]
        else: 
            print "no covariance estimated"
    else:
        print res
        print "FIT NOK : ",  res[3]
    pl.figure()
    pl.plot(oKur.getWL(), flux  )
    par = np.array([metUsed, 0., 0.])
    par[1:3] = res[0] 
    pl.plot(oKur.getWL(), oKur.getFluxInterLin(par))
    pl.legend(["True","Fit"])
    
 
class Test(unittest.TestCase):

    def notest_nearestFlux(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(1000, 12000)
        idx= oKur.getIdxNGP(np.array([1.3, 6100, 2.7]))
        oKur.plotFluxIdx(idx)
        
    def notest_getFluxInterLin1(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)
        par = np.array([-2.0, 6250., 2.0])
        parNGP1 = oKur.getParam(oKur.getIdxNGP(par))
        flux1 =oKur.getFluxNGP(par) 
        par = np.array([-2.0, 6000., 2.0])
        parNGP2 = oKur.getParam(oKur.getIdxNGP(par))
        flux2 =oKur.getFluxNGP(par) 
        par = np.array([-2.0, 6020., 2.0])        
        fluxLin = oKur.getFluxInterLin(par)
        pl.figure()
        pl.plot(oKur.getWL(), flux1)
        pl.plot(oKur.getWL(), flux2)
        pl.plot(oKur.getWL(), fluxLin)
        pl.legend(["NGP"+str(parNGP1),"NGP"+str(parNGP2),"Lin" + str(par)])
        pl.grid()
        
        
    def notest_getFluxInterLin2(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)
        par = np.array([-2.0, 7000., 2.0])
        parNGP1 = oKur.getParam(oKur.getIdxNGP(par))
        flux1 =oKur.getFluxNGP(par) 
        par = np.array([-2.5, 7000., 2.0])
        parNGP2 = oKur.getParam(oKur.getIdxNGP(par))
        flux2 =oKur.getFluxNGP(par) 
        par = np.array([-2.1, 7000., 2.0])        
        fluxLin = oKur.getFluxInterLin(par)
        pl.figure()
        pl.plot(oKur.getWL(), flux1)
        pl.plot(oKur.getWL(), flux2)
        pl.plot(oKur.getWL(), fluxLin,'--o')
        pl.legend(["NGP"+str(parNGP1),"NGP"+str(parNGP2),"Lin" + str(par)])
        pl.grid()
        
    def notest_getFluxInterLin3(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)
        par = np.array([-2.0, 5000., 2.0])
        parNGP1 = oKur.getParam(oKur.getIdxNGP(par))
        flux1 =oKur.getFluxNGP(par) 
        par = np.array([-2.0, 5000., 2.5])
        parNGP2 = oKur.getParam(oKur.getIdxNGP(par))
        flux2 =oKur.getFluxNGP(par) 
        par = np.array([-2.0, 5000., 2.1])        
        fluxLin = oKur.getFluxInterLin(par)
        pl.figure()
        pl.plot(oKur.getWL(), flux1)
        pl.plot(oKur.getWL(), flux2)
        pl.plot(oKur.getWL(), fluxLin,'--o')
        pl.legend(["NGP"+str(parNGP1),"NGP"+str(parNGP2),"Lin" + str(par)])
        pl.grid()
       
    def notest_nearestFluxes(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)
        aIdx = []
        temp = 9000
        aIdx.append(oKur.getIdxNGP(np.array([-1, temp, 2.7])))
        aIdx.append(oKur.getIdxNGP(np.array([-2, temp, 2.7])))
        aIdx.append(oKur.getIdxNGP(np.array([-3, temp, 2.7])))
        aIdx.append(oKur.getIdxNGP(np.array([-4, temp, 2.7])))
        oKur.plotFluxesArrayIdx(aIdx)
        aIdx = []
        aIdx.append(oKur.getIdxNGP(np.array([-1, 6000, 2.7])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 7000, 2.7])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 8000, 2.7])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 9000, 2.7])))
        oKur.plotFluxesArrayIdx(aIdx)
        aIdx = []
        aIdx.append(oKur.getIdxNGP(np.array([-1, 6100, 2])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 6100, 3])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 6100, 4])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 6100, 5])))
        oKur.plotFluxesArrayIdx(aIdx)        
        aIdx = []
        aIdx.append(oKur.getIdxNGP(np.array([-3.99999983e+00 ,  7.46518475e+03  , 3.74397055e+00])))
        #aIdx.append(oKur.getIdxNGP(np.array([-1.77619036e+01 ,  6.26773253e+03 ,  2.52832703e+00])))
        aIdx.append(oKur.getIdxNGP(np.array([-4.10844657e+00 ,  7.45306937e+03 ,  3.75453873e+00])))
        aIdx.append(oKur.getIdxNGP(np.array([-2, 6100, 2.7])))
        oKur.plotFluxesArrayIdx(aIdx)  
        
    def notest_fitLeastsqLinAll01(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)
        TruePar = np.array([-2, 6100, 2.7])      
        flux =oKur.getFluxInterLin(TruePar)
        sol = oKur._fitLeastsqLinAll(flux, 1) 
        print sol
        sol = oKur._fitLeastsqLinAll(flux, 5) 
        print sol
         
    def notest_fitLeastsqLinTemp01(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2000, 10000)
        TruePar = np.array([-2.5, 6100, 2.5])      
        flux =oKur.getFluxInterLin(TruePar)
        sol = oKur.fitLeastsqLinTemp(flux)
        print sol
        pl.figure()
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(TruePar))
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(sol))
        pl.legend(["True","Fit"])
        pl.grid()
        
    def notest_fitLeastsqLinTemp02(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2000, 10000)
        TruePar = np.array([-3, 9100, 2])      
        flux =oKur.getFluxInterLin(TruePar)
        sol = oKur.fitLeastsqLinTemp(flux)
        print sol
        pl.figure()
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(TruePar))
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(sol))
        pl.legend(["True","Fit"])
        pl.grid()
        
    def test_fitLeastsqLinAll3Step(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2000, 10000)
        TruePar = np.array([-3.2, 9150, 2.1])      
        fluxTrue =oKur.getFluxInterLin(TruePar)
        flux = fluxTrue*(1+np.random.normal(0, 0.002,len(fluxTrue)))
        sol = oKur.fitLeastsqLinAll3Step(flux)
        print sol
        pl.figure()
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(TruePar))
        pl.plot(oKur.getWL(), flux)
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(sol))
        pl.legend(["True","True noisy","Fit"])
        pl.grid()
        
    def notest_fitFminLinTemp01(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2000, 10000)
        TruePar = np.array([-2.5, 6100, 2.5])      
        flux =oKur.getFluxInterLin(TruePar)
        sol = oKur.fitFminLinTemp(flux)
        print sol
        pl.figure()
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(TruePar))
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(sol))
        sol[1] = 6000
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(sol),'y')
        pl.legend(["True","Fit",'6000'])
        pl.grid()
       
    
    def notest_leastsq(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)        
        flux =oKur.getFluxInterLin(np.array([-2, 6100, 2.7]))        
        def errorModel(par, kur, ydata):
            parI =  par.copy()
            flag = True
            if par[2] > 4:
                parI[2] = 4
                flag = False
            if par[2] <1:
                parI[2] = 1
                flag = False
            if par[0] < -4:
                parI[0] = -4
                flag = False
            if par[0] > 0.5:
                parI[0] = 0.5
                flag = False
            #print par 
            if flag:
                #res = (ydata - kur.getFluxInterLin(parI))*1e-7
                res = (ydata - kur.getFluxInterCubic(parI))*1e-7
            else:
                res = ydata
            print par, (res**2).sum()
            return res.ravel()
        print "True parameter"
        print "-2, 6100, 2.7"
        print "\nguess parameter"
        p0 = np.array([-2.1, 6900., 2.5])
        print p0
        res = spo.leastsq(errorModel, p0, args=(oKur, flux), full_output=True,factor=10, diag=np.array([1., 100,1.]))
        if res[4] >= 0 and res[4]<=4:
            print "FIT OK :",  res[3]
            #print res[0], res[2]["nfev"]
            print "\nSolution scipy.leastsq in %d calls "%res[2]["nfev"]
            print res[0]            
            if res[1] != None:                             
                print res[1]
            else: 
                print "no covariance estimated"
        else:
            print res
            print "FIT NOK : ",  res[3]
        pl.figure()
        pl.plot(oKur.getWL(), oKur.getFluxNGP(np.array([-2, 6100, 2.7]))   )
        par = res[0] 
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(par))
        pl.legend(["True","Fit"])
    
    def notest_fmin(self):
        oKur = Kurucz(FileKuruczPic)
        oKur.setWLInterval(2500, 10000)        
        flux =oKur.getFluxNGP(np.array([-2, 6100, 2.7]))        
        def Model(par, kur, flux):
            parI =  par.copy()
            flag = True
            if par[2] > 4:
                parI[2] = 4
                flag = False
            if par[2] <1:
                parI[2] = 1
                flag = False
            if par[0] < -4:
                parI[0] = -4
                flag = False
            if par[0] > 0.5:
                parI[0] = 0.5
                flag = False
            #print par             
            res = (flux -kur.getFluxInterLin(parI))**2
            return res.ravel().sum()
        print "True parameter"
        print "-2, 6100, 2.7"
        print "\nguess parameter"
        p0 = np.array([-2.1, 6900., 2.5])
        print p0
        res = spo.fmin(Model, p0, args=(oKur, flux))
        print res
        
        pl.figure()
        pl.plot(oKur.getWL(), oKur.getFluxNGP(np.array([-2, 6100, 2.7]))   )
        par = res[0] 
        pl.plot(oKur.getWL(), oKur.getFluxInterLin(par))
        pl.legend(["True","Fit"])
    
    def test_zz(self):
        try:
            pl.show()
        except AttributeError:
            pass
        
        
#test_setWLInterval()
#test_plotFlux()
test_getFluxInterLin()
#test_resample()
test_fitTempGra()
#unittest.main()
try:
    pl.show()
except AttributeError:
    pass
