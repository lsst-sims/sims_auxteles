'''
Created on 14 nov. 2012

@author: colley
'''
import unittest
from kurucz import *
from mpl_toolkits.mplot3d import Axes3D 

FileKuruczPic = '../../data/kurucz93/k93.pic'

def test_plotFlux():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(1000, 12000)
    oKur.plotFluxIdx(1540)
    oKur.plotMultiFluxesCont(1540,5)    
    oKur.plotFluxIdx(1544)
    
def test_plotNotDefined():
    oKur = Kurucz(FileKuruczPic, True)
    fig = pl.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('metalicity')
    ax.set_ylabel('temperature')
    ax.set_zlabel('gravity')
    print oKur._ParamNotDef.shape
    ax.scatter(oKur._ParamNotDef[:,0], oKur._ParamNotDef[:,1], oKur._ParamNotDef[:,2], 'b')
   
    
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
    
    
def test_restrictToWLinterval():
    oKur = Kurucz(FileKuruczPic)    
    par = np.array([0.0, 3660., 0.5])
    pl.plot(oKur.getWL(), oKur.getFluxInterLin(par))
    oKur.restrictToWLinterval(4000, 5000)
    print len(oKur.getWL()),len( oKur.getFluxInterLin(par))
    pl.plot(oKur.getWL(), oKur.getFluxInterLin(par))
    pl.plot(oKur.getWL(), oKur.getFluxInterLin(par),"*")
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
    
    
def test_plotVega():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(2000, 10000)
    pl.figure()
    pl.title('Raw Kurucz Vega')
    fl = oKur.getFluxInterLin(np.array([-0.5, 9550, 3.95]))
    pl.plot(oKur.getWL(), fl.T)
    pl.ylabel(oKur.UnitFluxStr)
    pl.xlabel(oKur.UnitWLStr)
    pl.grid()
    pl.figure()
    pl.title('Vega at Earth')
    flE = fl*5.63e-17
    pl.plot(oKur.getWL(), flE.T)
    pl.ylabel(r"$ergs.cm^{-2}.s^{-1}.A^{-1}$")
    pl.xlabel("Ansgtrom")
    pl.grid()
   

def test_fitTempGra():
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(2500, 10000)        
    fluxTrue =oKur.getFluxInterLin(np.array([-2, 7125, 3.7]))
    flux = (1+np.random.normal(0, 0.0001,len(fluxTrue)))*fluxTrue
    metUsed= -2.
    def errorModel(par, kur, ydata):
        parI = np.array([metUsed, 0. ,0.])
        parI[1:3] = par
        res = (ydata - oKur.getFluxInterLin(parI))
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


def test_fitWithBounds():
    """229878
    Seed: 485374
    
    produced
    
03373.003316
input par [ -1.35454174e+00   1.00209643e+04   1.62860243e+00]
103373.003307
input par [ -1.35454176e+00   1.00209643e+04   1.62860243e+00]
103373.003902
input par [  6.77699007e+00   1.01726894e+04  -6.57282362e+00]
fluxTheo is nan, used NGP 
56349971.6923    
    """
    seed = (int)(np.random.uniform(1,1000000))
    #seed = 229878
    print "Seed:",seed
    np.random.seed(seed)
    oKur = Kurucz(FileKuruczPic)
    oKur.setWLInterval(2500, 10000)    
    oKur.setCoefUnitFlux(1e-7)
    #flux = np.copy(fluxTrue)
    Met = np.random.uniform(oKur._BoundsMet[0], oKur._BoundsMet[1])
    Temp = np.random.uniform(4000, 10000)
    Gra = np.random.uniform(oKur._BoundsGra[0], oKur._BoundsGra[1])
    parTrue = np.array([Met, Temp, Gra])    
    print "True:     ", parTrue
    fluxTrue =oKur.getFluxInterLin(parTrue)
    flux = (1+np.random.normal(0, 0.0001,len(fluxTrue)))*fluxTrue
    print fluxTrue[:10]    
    g0 = np.array([-3, 10000, 3])
    sol = oKur.fitWithBoundsLM(flux, g0)
    print "True:     ", parTrue
    print "estimated:", sol
    print "erreur relativeen %:",(sol - parTrue)*100/parTrue
#    g0 = np.array([-1.0,  5100., 1.7])
#    sol = oKur.fitNoBounds(flux, g0)
#    print "erreur relative:",(sol - parTrue)*100/parTrue
    #idx= oKur.getIdxNGP(np.array([3.65307192e-01  , 1.38417718e+04   ,8.56629378e-01]))
    #oKur.plotFluxIdx(idx)

    
 
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
       
    def test_nearestFluxes(self):
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
        #oKur.plotFluxesArrayIdx(aIdx)
        aIdx = []
        aIdx.append(oKur.getIdxNGP(np.array([-1, 9000, 1])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 9000, 2])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 9000, 3])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 9000, 4])))
        aIdx.append(oKur.getIdxNGP(np.array([-1, 9000, 5])))
        #oKur.plotFluxesArrayIdx(aIdx)        
        aIdx = []
        aIdx.append(oKur.getIdxNGP(np.array([-3.99999983e+00 ,  7.46518475e+03  , 3.74397055e+00])))
        #aIdx.append(oKur.getIdxNGP(np.array([-1.77619036e+01 ,  6.26773253e+03 ,  2.52832703e+00])))
        aIdx.append(oKur.getIdxNGP(np.array([-4.10844657e+00 ,  7.45306937e+03 ,  3.75453873e+00])))
        aIdx.append(oKur.getIdxNGP(np.array([-2, 6100, 2.7])))
        #oKur.plotFluxesArrayIdx(aIdx)  
        
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
        
    def notest_fitLeastsqLinAll3Step(self):
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
        

def test_MKIII2Kurucz():
    # G5III     5150     +2.54       kp00_5250[g25] 
    x = 45
    temp = np.exp(-0.000018*x**3 +  0.002580*x**2  +  ( -0.142389)*x + 11.329507)
    gra = -0.000043*x**3 +   0.003600*x**2   +  (   -0.119226)*x +     4.816304
    met = 0.000003*x**3 +  (   -0.000183 )*x**2  +  (  -0.000820)*x + (  -0.028369)
    print "Kurucz README G5III     T=5150     G=+2.54     "
    print "temp:", temp
    print "gra:",gra
    print "met:",met
    
def test_MKV2Kurucz():
    # A0V       9520     +4.14   
    x = 20
    temp = np.exp(-0.000030*x**3 +  0.003951*x**2  +  ( -0.181333)*x + 11.551453)
    gra = -0.000085*x**3 + 0.006919*x**2   +  (-0.168342)*x +     5.097857
    met = -0.000038*x**3 +  (   0.005342 )*x**2  +  (-0.235097)*x + (2.989520)
    print "Kurucz README A0V       T=9520     G=+4.14   "
    print "temp:", temp
    print "gra:",gra
    print "met:",met


def test_MK2Kurucz():
    oKur = Kurucz(FileKuruczPic)
    print oKur.convertMK("G5", "III")
    print oKur.convertMK("A0", "V")
    print oKur.convertMK("AI", "V")
    print oKur.convertMK("UI", "V")
    print oKur.convertMK("A0", "v")
     
      
#test_setWLInterval()
#test_plotFlux()
#test_getFluxInterLin()
#test_resample()
#test_fitTempGra()
#unittest.main()
#test_fitWithBounds()
#test_restrictToWLinterval()
#test_plotVega()
test_MKIII2Kurucz()
test_MKV2Kurucz()
test_MK2Kurucz()

    

#test_plotFlux()
#test_plotNotDefined()
try:
    pl.show()
except AttributeError:
    pass
