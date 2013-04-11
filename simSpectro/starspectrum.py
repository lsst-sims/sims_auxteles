# -*- coding: utf-8 -*-

#############################################################################
#        starspectrum - a class defining the spectrum of an AT star         #
#                                                                           #
#        author:          Gurvan BAZIN                                      #
#                         Universitäts-Sternwarte-München                   #
#        latest revision: July 2010                                         #
#############################################################################

# include some useful libraries
# and class definitions

#import CalibSys 
from math import *
#from scipy import *
from constants import *
import numpy as np
import scipy as sp
import scipy.special
import pylab as pl


# en m/s
S_SpeedLight = 3e8
S_LevelDbg = 0

# a function to correct for wrong negative values
def notzero(x):
    if x<=0.:
        # if negative, put something very small
        return 1.e-10
    else:
        return x

# gaussian distribution
def gauss(x, sigma, mu=0.):
    return exp(-(x-mu)*(x-mu)/(2.*sigma*sigma))/(sigma*sqrt(2*pi))

def freq2wavelength(a, unit=1e-9):
    """
    a en Hertz
    defaut unit nm
    """
    return (S_SpeedLight/a)/unit

S_LevelPlot = 0

class starspectrum:

    def __init__(self):
        # define some arrays at each step of the calculation
        self.wl = []
        self.nu = []
        self.dEdl = []
        self.dEdn = []

        self.atmdEdl = []
        self.atmdEdn = []

        self.phot = []
        self.photnoise = []

        self.photconv = []
        self.photconvnoise = []

        self.photatm = []
        self.photatmnoise = []

        self.photmirror = []
        self.photmirrornoise = []

        self.photgrism = []
        self.photgrismnoise = []
        
        self.photccd = []
        self.photccdnoise = []
        self.wlccd = []
        self.nuccd = []
        
        self.elecccd = []
        self.elecccdnoise = []
        
        self.caldEdn = []
        self.caldEdnnoise = []
        self.caldEdl = []
        self.caldEdlnoise = []

        self.atmrebin = []

        self.sensdEdn = []
        self.sensdEdl = []

        self.syscaldEdl = []
        self.syscaldEdn = []

        # and characteristic variables for this observation
        self.factorslit = 1.
        self.factorextract = 1.
        
    def setSpectrum(self, wl, spectrum):
        self.wl = wl
        self.nu = clight/wl        
        self.dEdl = spectrum
        self.dEdn = wl*spectrum/self.nu
          
    def readdEdl(self, filename):
        # to read a spectrum in J/m2/s/nm
        facLamba=1
        print "facLamba ", facLamba
        mfile = open(filename, 'r')
        for line in mfile.readlines():
            if line[0] == "#":
                continue
            line = line.split()
            lamb = float(line[0])*facLamba
            self.wl.append(lamb)
            self.nu.append(clight/lamb)
            dedl = float(line[1])
            dedn = lamb*lamb*dedl/clight
            self.dEdn.append(dedn)
            self.dEdl.append(dedl)
        mfile.close()
        self.nu = np.array(self.nu)
        self.wl = np.array(self.wl)
        self.dEdn = np.array(self.dEdn)
        self.dEdl = np.array(self.dEdl)
        self.readverif(filename)
        print "nu read min max", self.nu[0], self.nu[-1]
        print "1/nu read min max", 1/self.nu[-1], 1/self.nu[0]
        print "lbd read min max", self.wl[0], self.wl[-1]
        
        
            
    def plotWL(self,pTitle="star spectrum", wl=None, dEdl=None):
        if dEdl == None: dEdl=self.dEdl
        if wl == None: wl=self.wl
        pl.figure()        
        pl.xlabel("wavelength nm")
        #pl.ylabel("$J.m^{-2}.s^{-1}.nm^{-1}$")
        pl.grid()
        pl.title(pTitle)
        pl.plot(wl, dEdl)
        
        
    def plotNbPhotons(self, pTitle="star spectrum", phot = None):
        if phot == None:
            phot = self.phot
        pl.figure()        
        pl.xlabel("? ")
        pl.ylabel("Number photons")
        pl.grid()
        pl.title(pTitle)
        pl.plot(phot)
        
        
    def plotToWLnm(self, nu, val, pTitle="star spectrum"):
        wl = freq2wavelength(nu)
        pl.figure()        
        pl.xlabel("nm")
        pl.grid()
        pl.title(pTitle)
        pl.plot(wl, val)


    def readverif(self, filename):
        # just for checking... and correctly sort the arrays 
        # in wavelength and frequency 
        if len(self.nu>2):
            if self.nu[0]>self.nu[1]:
                self.nu = np.flipud(self.nu)
                self.wl = np.flipud(self.wl)
                self.dEdn = np.flipud(self.dEdn)
                self.dEdl = np.flipud(self.dEdl)
        else:
            print 'error reading', filename

    def apertureCorrection(self, seeing, slitwidth, factext=0.68):
        # apply an aperture correction with respect to the seeing and slitwidth
        # also apply a correction corresponding to the aperture used to extract the 1D spectrum from the 2D spectrum: 1\sigma by default
        sigma_seeing = seeing/(2.*(2.*log(2))**.5)
        halfwidth_slit = slitwidth/2.
        self.factorslit = .5*( scipy.special.erf(halfwidth_slit/(sigma_seeing*2.**.5)) - scipy.special.erf(-halfwidth_slit/(sigma_seeing*2.**.5)) )
        self.factorextract = factext        
        self.dEdnAper = self.dEdn*self.factorslit*self.factorextract
        self.dEdlAper = self.dEdl*self.factorslit*self.factorextract

    def addPhotonNoise(self):
        # add the calculated photon noise
        self.photccd += np.random.normal(0., self.photccdnoise)

    def computePhoton(self, effarea=1.0, exptime=240.):
        # compute photons from the input spectrum (energy density)
        # assumed effective area of the telescope = 1.0 m (mirror, slit...)
        # todo: take the calypso characteristics
        # exposure of 4 min
        adn =[]
        for i in range(len(self.nu)):
            if i == 0:
                dn = (self.nu[i+1]-self.nu[i])
            elif i == len(self.nu)-1:
                dn = (self.nu[i]-self.nu[i-1])
            else:
                dn = ((self.nu[i+1]-self.nu[i])+(self.nu[i]-self.nu[i-1]))/2.
            adn.append(dn)
            #self.phot.append((self.dEdn[i] * dn ) / (self.nu[i] * hplanck))
            self.phot.append((self.dEdnConv[i] * dn ) / (self.nu[i] * hplanck))
        self.phot = np.array(self.phot) * effarea * exptime
        self.photconv = self.phot
        if S_LevelPlot >= 2 :
            diffNu = np.diff(self.nu)
            pl.figure()
            pl.title("Delta frequency")
            pl.xlabel("Hertz")
            pl.ylabel("Hertz")
            pl.plot(self.nu, np.array(adn))
            pl.plot(self.nu[0:-1], diffNu)
            pl.legend(["diff centred","diff"],loc=0)
            pl.grid()

    def computePhotonNoise(self):
        # calculate the Poisson photon noise
        self.photccdnoise = np.sqrt(self.photccd)
        self.photccdnoise = map(notzero, self.photccdnoise)

    def computePhotonATM(self, atm):
        """
        jmc: 
        atm = [self.nu, self.wl, self.tr]
        """
        # loop on star bins
        # jmc : loop on index nu ccd pixel
        for i in range(len(self.nuccd)):
            # computing distances between the center and the edges of each bin
            # jmc : defined size interval around nuccd[i] left right
            #       size interval is dnccdd + dnccdg
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.

            j = 0
            foundj = False
            # Search for the first bin j of atm spectrum
            # overlapping the bin i of the star
            while j < len(atm[0,:]):
                # jmc : defined size interval around atm[0,j]
                #       size interval is dnatmd + dnatmg
                if j == 0:
                    dnatmd = (atm[0,j+1]-atm[0,j])/2.
                    dnatmg = dnatmd
                elif j == len(atm[0,:])-1:
                    dnatmg = (atm[0,j]-atm[0,j-1])/2.
                    dnatmd = dnatmg
                else:
                    dnatmd = (atm[0,j+1]-atm[0,j])/2.
                    dnatmg = (atm[0,j]-atm[0,j-1])/2.
                    
                if atm[0,j]-dnatmg > self.nuccd[i]-dnccdg:
                    break
                bordGaucheCCD = self.nuccd[i]-dnccdg
                if atm[0,j]-dnatmg <= bordGaucheCCD and atm[0,j]+dnatmd > bordGaucheCCD:
                    foundj = True
                    break
                j+=1
            if foundj:
                k = j
                foundk = False
                # Search for the last bin k of atm spectrum
                # overlapping the bin i of the star
                while k < len(atm[0,:]):
                    # jmc : 
                    if k == 0:
                        dnatmd = (atm[0,k+1]-atm[0,k])/2.
                        dnatmg = dnatmd
                    elif k == len(atm[0,:])-1:
                        dnatmg = (atm[0,k]-atm[0,k-1])/2.
                        dnatmd = dnatmg
                    else:
                        dnatmd = (atm[0,k+1]-atm[0,k])/2.
                        dnatmg = (atm[0,k]-atm[0,k-1])/2.
                        
                    if atm[0,k]-dnatmg > self.nuccd[i]+dnccdd:
                        # jmc atm_k sort du pixel CCD
                        foundk = True
                        break
                    k+=1
            # if everything is ok, compute proportion of absorption for each part the
            # star bin i
            if foundj and foundk:
                nbbin = k-j
                bintmp = 0.
                # jmc: loop on atm touching pixel
                for m in range(nbbin):
                    # jmc: n index atm bin
                    n = j+m
                    if n == 0:
                        # jmc: first bin atm
                        dnatmd = (atm[0,n+1]-atm[0,n])/2.
                        dnatmg = dnatmd
                    elif n == len(atm[0,:])-1:
                        # jmc: last bin atm
                        dnatmg = (atm[0,n]-atm[0,n-1])/2.
                        dnatmd = dnatmg
                    else:
                        # jmc:  
                        dnatmd = (atm[0,n+1]-atm[0,n])/2.
                        dnatmg = (atm[0,n]-atm[0,n-1])/2.
                    
                    if atm[0,n]-dnatmg <= self.nuccd[i]-dnccdg:
                        nu1 = self.nuccd[i]-dnccdg
                    else:
                        nu1 = atm[0,n]-dnatmg
                    if atm[0,n]+dnatmd >= self.nuccd[i]+dnccdd:
                        nu2 = self.nuccd[i]+dnccdd
                    else:
                        nu2 = atm[0,n]+dnatmd
                    # jmc : prop (delta nu ATM )/(delta nu CCD) 
                    prop = (nu2-nu1)/(dnccdg+dnccdd)
                    # bintmp moyenne de la transmission ATM pondéré par la taille relative des bin ATM 
                    bintmp += prop*atm[2,n]
                self.photatm.append(bintmp*self.photccd[i])
                self.photatmnoise.append(bintmp*self.photccdnoise[i])
            # else, put 0
            else:
                self.photatm.append(0.)
                self.photatmnoise.append(0.)
        self.photatm = np.array(self.photatm)
        self.photatmnoise = np.array(self.photatmnoise)
        # if negative bins, put 0
        for i in range(len(self.photatm)):
            # jmc 
            #if self.photatm[i]<1.:
            if self.photatm[i]<0.:
                self.photatm[i] = 0.
                self.photatmnoise[i] = 0.

    def computePhotonMirror(self, mirrorresp):
        # Mirror reflectivity
        ginterp = sp.interpolate.interp1d(mirrorresp.getnu(), mirrorresp.getrenu(), bounds_error=False, fill_value = 0.)
        mirrortr = ginterp(self.nuccd)
        self.mirrortr = mirrortr
        self.photmirror = self.photatm * mirrortr
        self.photmirrornoise = self.photatmnoise * mirrortr
        
    def computePhotonGrism(self, grismresp):
        # Transmission of the grism
        gtr = []
        for i in range(len(self.nuccd)):
            gtr.append(grismresp.getT(clight/self.nuccd[i]))
        self.grismTrans = np.array(gtr)
        self.photgrism = self.photmirror * np.array(gtr)
        self.photgrismnoise = self.photmirrornoise * np.array(gtr)


    def computePhotonCCD(self, grism, lpix1=300.e-9, lpix2=1100.e-9, nbpixel=512):
        # Dispersion of the grism in lambda
        # sin i - sin i' = lambda/a (put i = 0)
        # calculate the position of self.wl on the ccd and interpolate on the linear grid of the 512 pixels
        print  self.wl       
        y = [ tan(grism.outputAngle(e)) for e in self.wl ]
        if S_LevelPlot > 3:
            pl.figure()
            pl.plot(self.wl, y)
            pl.plot(self.wl, y,"*")
            pl.title('output grism')
            pl.grid()
            pl.xlabel(' wave length nm')
        dd = (y[-1]-y[0])/nbpixel
        yccd = [ y[0]+dd*i for i in range(nbpixel) ]
        finterp = sp.interpolate.interp1d(np.flipud(np.array(y)), np.flipud(self.wl), bounds_error=False, fill_value = 0.)
        self.wlccd = np.flipud(finterp(np.flipud(np.array(yccd))))
        self.nuccd = clight/self.wlccd
        print "wlccd size ",np.size(self.wlccd)
        #print "wlccd : ", self.wlccd[0],  self.wlccd[-1]    
        for i in range(len(self.nuccd)):
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.

            # Search for the first bin j of the spectrum
            # overlapping the bin i of the ccd
            j = 0
            foundj = False
            while j < len(self.nu):
                if j == 0:
                    dnd = (self.nu[j+1]-self.nu[j])/2.
                    dng = dnd
                elif j == len(self.nu)-1:
                    dng = (self.nu[j]-self.nu[j-1])/2.
                    dnd = dng
                else:
                    dnd = (self.nu[j+1]-self.nu[j])/2.
                    dng = (self.nu[j]-self.nu[j-1])/2.
                    
                if self.nu[j]-dng <= self.nuccd[i]-dnccdg and self.nu[j]+dnd > self.nuccd[i]-dnccdg:
                    foundj = True
                    break
                elif self.nu[j] > self.nuccd[i]:
                    break
                j+=1
            foundk = False
            if foundj:
                # Search for the last bin k of the spectrum
                # overlapping the bin i of the ccd
                k = j
                while k < len(self.nu):
                    if k == 0:
                        dnd = (self.nu[k+1]-self.nu[k])/2.
                        dng = dnd
                    elif k == len(self.nu)-1:
                        dng = (self.nu[k]-self.nu[k-1])/2.
                        dnd = dng
                    else:
                        dnd = (self.nu[k+1]-self.nu[k])/2.
                        dng = (self.nu[k]-self.nu[k-1])/2.
                    
                    if self.nu[k]-dng > self.nuccd[i]+dnccdd:
                        foundk = True
                        break
                    k+=1
            # if ok
            if foundj and foundk:
                nbbin = k-j
                bintmp = 0
                for m in range(nbbin):
                    n = j+m
                    if n == 0:
                        dnd = (self.nu[n+1]-self.nu[n])/2.
                        dng = dnd
                    elif n == len(self.nu)-1:
                        dng = (self.nu[n]-self.nu[n-1])/2.
                        dnd = dng
                    else:
                        dnd = (self.nu[n+1]-self.nu[n])/2.
                        dng = (self.nu[n]-self.nu[n-1])/2.
                        
                    if self.nu[n]-dng >= self.nuccd[i]-dnccdg:
                        nu1 = self.nu[n]-dng
                    else:
                        nu1 = self.nuccd[i]-dnccdg
                    if self.nu[n]+dnd <= self.nuccd[i]+dnccdd:
                        nu2 = self.nu[n]+dnd
                    else:
                        nu2 = self.nuccd[i]+dnccdd
                    prop = (nu2-nu1)/(dnd+dng)
                    if S_LevelDbg > 2 : print "pixel:", i,prop
                    bintmp += prop*self.photconv[n]
                self.photccd.append(bintmp)
            else:
                self.photccd.append(0)
                self.photccdnoise.append(0)
        self.photccd = np.array(self.photccd)
        #print "photccd size ", size(self.photccd)
                 
    def computeElectronCCD(self, ccdresp, gain=1.):
        # converte in photo-electrons
        # Gain = 1 (assumption)
        # Response of the CCD
        finterp = sp.interpolate.interp1d(ccdresp.getnu(), ccdresp.getrenu(), bounds_error=False, fill_value = 0.)
        ccdr = finterp(self.nuccd)
        self.elecccd = self.photgrism * ccdr
        self.elecccdnoise = self.photgrismnoise * ccdr
        if S_LevelPlot >3:
            pl.figure()
            pl.plot(ccdresp.getnu(), ccdresp.getrenu()*100)
            pl.plot(self.nuccd, ccdr*100)
            pl.plot(self.nuccd, ccdr*100, '*')
            pl.title("Interpolation transmission CCD")
            pl.grid()
            pl.xlabel("Hertz")
            pl.ylabel("%")
            pl.legend(["raw","interpolate","interpolate"])
            
    def addNoiseElectronCCD(self, noise=5.):
        # noise sigma = 5 e-
        for i in range(len(self.elecccd)):
            if self.elecccd[i] > 0.:
                n = np.random.normal(0., noise)
                self.elecccd[i] += n
                self.elecccdnoise[i] = sqrt(self.elecccdnoise[i]*self.elecccdnoise[i] + n*n)

    def rebinatm(self, atm):
        # rebin the calibration atm on the CCD bins
        # loop on star bins
        for i in range(len(self.nuccd)):
            # computing distances between the center and the edges of each bin
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.

            j = 0
            foundj = False
            # Search for the first bin j of atm spectrum
            # overlapping the bin i of the star
            while j < len(atm[0,:]):
                if j == 0:
                    dnatmd = (atm[0,j+1]-atm[0,j])/2.
                    dnatmg = dnatmd
                elif j == len(atm[0,:])-1:
                    dnatmg = (atm[0,j]-atm[0,j-1])/2.
                    dnatmd = dnatmg
                else:
                    dnatmd = (atm[0,j+1]-atm[0,j])/2.
                    dnatmg = (atm[0,j]-atm[0,j-1])/2.
                    
                if atm[0,j]-dnatmg > self.nuccd[i]-dnccdg:
                    break
                if atm[0,j]-dnatmg <= self.nuccd[i]-dnccdg and atm[0,j]+dnatmd > self.nuccd[i]-dnccdg:
                    foundj = True
                    break
                j+=1
            if foundj:
                k = j
                foundk = False
                # Search for the last bin k of atm spectrum
                # overlapping the bin i of the star
                while k < len(atm[0,:]):
                    if k == 0:
                        dnatmd = (atm[0,k+1]-atm[0,k])/2.
                        dnatmg = dnatmd
                    elif k == len(atm[0,:])-1:
                        dnatmg = (atm[0,k]-atm[0,k-1])/2.
                        dnatmd = dnatmg
                    else:
                        dnatmd = (atm[0,k+1]-atm[0,k])/2.
                        dnatmg = (atm[0,k]-atm[0,k-1])/2.
                        
                    if atm[0,k]-dnatmg > self.nuccd[i]+dnccdd:
                        foundk = True
                        break
                    k+=1
            # if everything is ok, compute proportion of absorption for each part the
            # star bin i
            if foundj and foundk:
                nbbin = k-j
                bintmp = 0.
                for m in range(nbbin):
                    n = j+m
                    if n == 0:
                        dnatmd = (atm[0,n+1]-atm[0,n])/2.
                        dnatmg = dnatmd
                    elif n == len(atm[0,:])-1:
                        dnatmg = (atm[0,n]-atm[0,n-1])/2.
                        dnatmd = dnatmg
                    else:
                        dnatmd = (atm[0,n+1]-atm[0,n])/2.
                        dnatmg = (atm[0,n]-atm[0,n-1])/2.
                    
                    if atm[0,n]-dnatmg <= self.nuccd[i]-dnccdg:
                        nu1 = self.nuccd[i]-dnccdg
                    else:
                        nu1 = atm[0,n]-dnatmg
                    if atm[0,n]+dnatmd >= self.nuccd[i]+dnccdd:
                        nu2 = self.nuccd[i]+dnccdd
                    else:
                        nu2 = atm[0,n]+dnatmd
                    prop = (nu2-nu1)/(dnccdg+dnccdd)
                    bintmp += prop*atm[2,n]
                self.atmrebin.append(bintmp)
            # else, put 0
            else:
                self.atmrebin.append(0.)
        self.atmrebin = np.array(self.atmrebin)


    def MakeSensFuncInstruOnly(self, ccdresp, grismresp, mirrorresp, gain=1., effarea=1.0, exptime=240.):
        """
        """
        self._deltaNuCCD = []
        self.sensdEdnJMC = np.zeros(len(self.nuccd), dtype=np.float64)
        self.sensdEdlJMC = np.zeros(len(self.nuccd), dtype=np.float64)
        # CCD
        tck = sp.interpolate.splrep(ccdresp.getnu(), ccdresp.getrenu())
        TransCCDPixel = sp.interpolate.splev(self.nuccd, tck)
  
        # MIROR
        tck = sp.interpolate.splrep(mirrorresp.getnu(), mirrorresp.getrenu())
        TransMIRPixel = sp.interpolate.splev(self.nuccd, tck)
        #print "MakeSensFuncJMC: ", gain, effarea, self.factorslit, self.factorextract
        factorDivers = gain*effarea*self.factorslit*self.factorextract
        # loop on pixel CCD
        for i in range(len(self.nuccd)):
            # -> J.m-2.s-1
            if TransCCDPixel[i] > 0.  and TransMIRPixel[i] > 0.:
                gtr = grismresp.getT(clight/self.nuccd[i])
                TransInstru= gtr*TransCCDPixel[i]*TransMIRPixel[i]                
                calib = (hplanck*self.nuccd[i])/(TransInstru*exptime*factorDivers)
            else:
                calib = 0.
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
            self._deltaNuCCD.append(dnccdg+dnccdd)
            # final estimation
            self.sensdEdnJMC[i] = calib/(dnccdg+dnccdd)
            #self.sensdEdnJMC[i] = calib
        self.sensdEdlJMC = self.sensdEdnJMC*self.nuccd**2/clight
        if S_LevelPlot > 2:
            pl.figure()
            pl.title(" fnc calib JMC")
            pl.grid()
            pl.plot(self.wlccd[50:-50], self.sensdEdlJMC[50:-50])
        self.sensdEdl = self.sensdEdlJMC
        self.sensdEdn = self.sensdEdnJMC 
        self._deltaNuCCD = np.array(self._deltaNuCCD) 

    def MakeSensFuncJMC(self, atm, ccdresp, grismresp, mirrorresp, gain=1., effarea=1.0, exptime=240.):
        self._deltaNuCCD = []
        self.sensdEdnJMC = np.zeros(len(self.nuccd), dtype=np.float64)
        self.sensdEdlJMC = np.zeros(len(self.nuccd), dtype=np.float64)
        # CCD
        tck = sp.interpolate.splrep(ccdresp.getnu(), ccdresp.getrenu())
        TransCCDPixel = sp.interpolate.splev(self.nuccd, tck)
        # ATM
        # smooth atmosphere transmission
        minDeltaWL = np.min(np.diff(atm[1][::-1]))
        xcst = np.arange(min(atm[1]), max(atm[1]), minDeltaWL)
        tck = sp.interpolate.splrep(atm[1][::-1], atm[2][::-1])
        AtmPer = sp.interpolate.splev(xcst, tck, der=0)
        kernel = []        
        minDeltaCCD = np.max(np.fabs(np.diff(self.wlccd)))
        sigma = minDeltaCCD/4
        for e in np.arange(-10.*sigma, 10.*sigma+minDeltaWL, minDeltaWL):
            kernel.append(gauss(e, sigma, 0.))
        kernel = np.array(kernel)
        # convolve
        #convycst = scipy.signal.fftconvolve(ycst, kernel, 'same')
        AtmConv= sp.convolve(AtmPer, kernel, 'same')/kernel.sum()            
       
        tck = sp.interpolate.splrep(xcst,AtmConv )
        TransATMPixel = sp.interpolate.splev(self.wlccd, tck)
        #tck = interpolate.splrep(atm[0],atm[2] )
        #TransATMPixel = interpolate.splev(self.nuccd, tck)
        if S_LevelPlot > 2: 
            pl.figure()
            pl.xlabel("wavelength nm")
            pl.ylabel("?")
            pl.grid()
            pl.title("Atm conv. par %g"%sigma)
            pl.plot(atm[1], atm[2])      
            pl.plot(xcst, AtmConv) 
            pl.plot(xcst, AtmConv) 
            pl.plot(self.wlccd, TransATMPixel) 
            pl.legend(["atm raw","atm conv.","atm@CCD"])           
        # MIROR
        tck = sp.interpolate.splrep(mirrorresp.getnu(), mirrorresp.getrenu())
        TransMIRPixel = sp.interpolate.splev(self.nuccd, tck)
        #print "MakeSensFuncJMC: ", gain, effarea, self.factorslit, self.factorextract
        factorDivers = gain*effarea*self.factorslit*self.factorextract
        # loop on pixel CCD
        for i in range(len(self.nuccd)):
            # -> J.m-2.s-1
            if TransCCDPixel[i] > 0. and TransATMPixel[i] > 0. and TransMIRPixel[i] > 0.:
                gtr = grismresp.getT(clight/self.nuccd[i])
                TransInstru= gtr*TransCCDPixel[i]*TransMIRPixel[i]                
                calib = (hplanck*self.nuccd[i])/(TransInstru*TransATMPixel[i]*exptime*factorDivers)
            else:
                calib = 0.
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
            self._deltaNuCCD.append(dnccdg+dnccdd)
            # final estimation
            self.sensdEdnJMC[i] = calib/(dnccdg+dnccdd)
            #self.sensdEdnJMC[i] = calib
        self.sensdEdlJMC = self.sensdEdnJMC*self.nuccd**2/clight
        if S_LevelPlot > 2:
            pl.figure()
            pl.title(" fnc calib JMC")
            pl.grid()
            pl.plot(self.wlccd[50:-50], self.sensdEdlJMC[50:-50])
        self.sensdEdl = self.sensdEdlJMC
        self.sensdEdn = self.sensdEdnJMC 
        self._deltaNuCCD = np.array(self._deltaNuCCD) 
        
        
    def MakeSensFunc(self, atm, ccdresp, grismresp, mirrorresp, gain=1., effarea=1.0, exptime=240.):
        # calculate the sensitivity function for flux calibration        
        # rebin atm in ccd bin
        self.rebinatm(atm)
        sensdEdntmp = []
        # points outside the main absorption lines of the ATM, from Thomas Matheson IDL procedures, slighly modified for convergence of the interpolation
        # en nm
        WLnm = np.array([ 3080.0, 3160.0, 3240.0, 3320.0, 3400.0, 3480.0, 3560.0, 3640.0, 3680.0, 3720.0, 3760.0, 3780.0, 3820.0, 3860.0, 4020.0, 4200.0, 4400.0, 4560.0, 4760.0, 5000.0, 5120.0, 5240.0, 5400.0, 5560.0, 5760.0, 6020.0, 6420.0, 6780.0, 7100.0, 7460.0, 7780.0, 8100.0, 8380.0, 8780.0, 8900.0, 9900.0, 9940.0,10260.0,10820.0, 11140.0, 12000.0])[::-1]/10.
        nutmp = clight/WLnm
        finterp = sp.interpolate.interp1d(ccdresp.getnu(), ccdresp.getrenu(), bounds_error=False, fill_value = 0.)
        ginterp = sp.interpolate.interp1d(mirrorresp.getnu(), mirrorresp.getrenu(), bounds_error=False, fill_value = 0.)
        hinterp = sp.interpolate.interp1d(self.nuccd, self.atmrebin, bounds_error=False, fill_value = 0.)
        ccdr = finterp(nutmp)
        mirrortr = ginterp(nutmp)
        atmrebintr = hinterp(nutmp)
        for i in range(len(nutmp)):
            # -> J.m-2.s-1
            if ccdr[i] > 0. and mirrortr[i] > 0. and atmrebintr[i] > 0.:
                gtr = grismresp.getT(clight/nutmp[i])
                sensdEdntmp.append(hplanck*nutmp[i]/(ccdr[i]*gain*gtr*mirrortr[i]*effarea*exptime*atmrebintr[i]*self.factorslit*self.factorextract))
            else:
                sensdEdntmp.append(0.)
        nudiff0 = []
        sensdiff0 = []
        for i in range(len(sensdEdntmp)):
            if sensdEdntmp[i] != 0.:
                nudiff0.append(nutmp[i])
                sensdiff0.append(sensdEdntmp[i])
        # cubic splines interpolation
        tck = sp.interpolate.splrep(nudiff0, log10(sensdiff0), s=0)
        self.sensdEdn = 10.**(sp.interpolate.splev(self.nuccd, tck, der=0))
        if True:
            pl.figure()
            pl.title("compare fnc calib")
            pl.plot(WLnm, sensdEdntmp)
            pl.plot(WLnm, sensdEdntmp,"*")
            pl.plot(self.wlccd, self.sensdEdn)
            pl.legend(["raw","raw","interpol CCD"])
            pl.figure()
            pl.title("compare fnc atm")
            pl.plot(WLnm, atmrebintr)
            pl.plot(atm[1], atm[2])
            pl.legend(["atm interpol","raw Atm"])
            
        # next, -> J.m-2.s-1.Hz-1
        for i in range(len(self.nuccd)):
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.

            self.sensdEdn[i] /= (dnccdg+dnccdd)
        self.sensdEdl = np.array(self.sensdEdn)*np.array(self.nuccd)**2/clight
        #self.sensdEdl = self.sensdEdn*self.nuccd**2/clight

    def calibrateSensSpec(self):
        # flux-calibrate the spectrum using the sens function
        for i in range(len(self.nuccd)):
            self.caldEdn.append(self.elecccd[i]*self.sensdEdn[i])
            self.caldEdnnoise.append(self.elecccdnoise[i]*self.sensdEdn[i])
        self.caldEdn = np.array(self.caldEdn)
        self.caldEdnnoise = np.array(self.caldEdnnoise)
        self.caldEdl = self.caldEdn*self.nuccd**2/clight
        self.caldEdlnoise = self.caldEdnnoise*self.nuccd**2/clight


    def calibrateSpec(self, ccdresp, grismresp, mirrorresp, gain=1., effarea=1.0, exptime=240.):
        # Calibrate the spectrum with the perfectly known
        # response of the instrument
        # to use for checking...
        finterp = sp.interpolate.interp1d(ccdresp.getnu(), ccdresp.getrenu(), bounds_error=False, fill_value = 0.)
        ginterp = sp.interpolate.interp1d(mirrorresp.getnu(), mirrorresp.getrenu(), bounds_error=False, fill_value = 0.)
        ccdr = finterp(self.nuccd)
        mirrortr = ginterp(self.nuccd)
        for i in range(len(self.nuccd)):
            # -> J.m-2.s-1
            if ccdr[i] != 0. and mirrortr[i] != 0.:
                gtr = grismresp.getT(clight/self.nuccd[i])
                self.caldEdn.append(self.elecccd[i]*hplanck*self.nuccd[i]/(ccdr[i]*gain*gtr*mirrortr[i]*effarea*exptime))
                self.caldEdnnoise.append(self.elecccdnoise[i]*hplanck*self.nuccd[i]/(ccdr[i]*gain*gtr*mirrortr[i]*effarea*exptime))
            else:
                self.caldEdn.append(0.)
                self.caldEdnnoise.append(0.)
            # next, -> J.m-2.s-1.Hz-1
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.

            self.caldEdn[i] /= (dnccdg+dnccdd)
            self.caldEdnnoise[i] /= (dnccdg+dnccdd)
        self.caldEdl = self.caldEdn*self.nuccd**2/clight
        self.caldEdlnoise = self.caldEdnnoise*self.nuccd**2/clight

    def write(self, spectype, filename):
        print "write file '%s'"%filename
        # what type of spectrum do you want to write ?
        if spectype == 'dEdn':
            spec = [self.nu, self.dEdn, np.zeros(len(self.nu))]
        elif spectype == 'dEdnl':
            spec = [self.wl, self.dEdn, np.zeros(len(self.nu))]
        elif spectype == 'dEdl':
            spec = [self.wl, self.dEdl, np.zeros(len(self.nu))]
        elif spectype == 'phot':
            spec = [self.nu, self.phot, np.zeros(len(self.nu))]
        elif spectype == 'photconv':
            spec = [self.nu, self.photconv, np.zeros(len(self.nu))]
        elif spectype == 'photatm':
            spec = [self.nuccd, self.photatm, self.photatmnoise]
        elif spectype == 'photmirror':
            spec = [self.nuccd, self.photmirror, self.photmirrornoise]
        elif spectype == 'photgrism':
            spec = [self.nuccd, self.photgrism, self.photgrismnoise]
        elif spectype == 'photccd':
            spec = [self.nuccd, self.photccd, self.photccdnoise]
        elif spectype == 'elecccd':
            spec = [self.nuccd, self.elecccd, self.elecccdnoise]
        elif spectype == 'elecccdlambda':
            spec = [self.wlccd[::-1], self.elecccd[::-1], self.elecccdnoise[::-1]]
        elif spectype == 'caldEdn':
            spec = [self.nuccd, self.caldEdn, self.caldEdnnoise]
        elif spectype == 'caldEdl':
            spec = [self.wlccd[::-1], self.caldEdl[::-1], self.caldEdlnoise[::-1]]
        elif spectype == 'syscaldEdn':
            spec = [self.nuccd, self.syscaldEdn, self.caldEdnnoise]
        elif spectype == 'syscaldEdl':
            spec = [self.wlccd[::-1], self.syscaldEdl[::-1], self.caldEdlnoise[::-1]]
        elif spectype == 'atmdEdl':
            spec = [self.wl[::-1], self.atmdEdl[::-1], np.zeros(len(self.wl))]
        elif spectype == 'atmrebin':
            spec = [self.wlccd[::-1], self.atmrebin[::-1], np.zeros(len(self.wlccd))]
        elif spectype == 'sensdEdl':
            spec = [self.wlccd[::-1], self.sensdEdl[::-1], np.zeros(len(self.wlccd))]
        elif spectype == 'sensdEdn':
            spec = [self.nuccd, self.sensdEdn, np.zeros(len(self.nuccd))]
        myfile = open(filename, 'w')
        for i in range(len(spec[0])):
            myfile.write(str(spec[0][i])+'\t'+str(spec[1][i])+'\t'+str(spec[2][i])+'\n')
        myfile.close()


    def convolvePhotondEdl(self, sigma):
        minDeltaWL = np.min(np.diff(self.wl[::-1]))
        xcst = np.arange(min(self.wl), max(self.wl), minDeltaWL)
        tck = sp.interpolate.splrep(self.wl[::-1], self.dEdlAper[::-1])
        ycst = sp.interpolate.splev(xcst, tck, der=0)
        # put the gaussian kernel into an array
        kernel = []
        for e in np.arange(-10.*sigma, 10.*sigma+minDeltaWL, minDeltaWL):
            kernel.append(gauss(e, sigma, 0.))
        kernel = np.array(kernel)
        # convolve
        #convycst = scipy.signal.fftconvolve(ycst, kernel, 'same')
        convycst = sp.convolve(ycst, kernel, 'same')
        # and put it again on the initial grid
        tck = sp.interpolate.splrep(xcst, convycst)
        convy = sp.interpolate.splev(self.wl, tck, der=0)
        self.dEdlConv = convy/kernel.sum()
        self.dEdnConv = (self.wl**2)*self.dEdlConv/clight
        if S_LevelPlot > 0:
            pl.figure()
            pl.xlabel("wavelength nm")
            pl.ylabel("?")
            pl.grid()
            pl.title("Sky spectrum at telescope resolution, conv. par %g"%sigma)
            pl.plot(self.wl, self.dEdl)            
            pl.plot(self.wl, self.dEdlConv) 
            pl.legend(["dEdl","dEdl conv."])           
            pl.figure()
            pl.xlabel("Hertz")
            pl.ylabel("?")
            pl.grid()
            pl.title("Sky spectrum at telescope resolution, conv. par %g"%sigma)
            pl.plot(self.nu, self.dEdn)            
            pl.plot(self.nu, self.dEdnConv) 
            pl.legend(["dEdn","dEdn conv."])           
            

    def convolvePhoton(self, sigma, doPlot=False):
        # Convolve the spectrum to put it at the good resolution
        # determine the minimum step
        if S_LevelPlot > 0:
            pl.figure()
            pl.xlim(200, 1200)
            pl.xlabel("wavelength nm")
            pl.ylabel("nb photon")
            pl.grid()
            pl.title("Sky spectrum at telescope resolution, conv. par %g"%sigma)
            pl.plot(self.wl, self.phot)  
                      
        minstep = 1.e6
        for i in range(len(self.wl)-1):
            tmpstep = self.wl[i]-self.wl[i+1]
            if minstep>tmpstep:
                minstep = tmpstep
        # interpolate the spectrum on a grid with constant step
        xcst = np.arange(min(self.wl), max(self.wl), minstep)
        tck = sp.interpolate.splrep(self.wl[::-1], self.phot[::-1])
        ycst = sp.interpolate.splev(xcst, tck, der=0)
        if S_LevelPlot >1:
            pl.figure()
            pl.plot(self.wl, self.phot)  
            pl.plot(xcst, ycst)
            pl.title("B-Spline approximation photon number by frequency interval")
            pl.legend(["phot","B-spline"])
            pl.xlabel("wavelength nm")
            pl.ylabel("nb photon")
            pl.grid()
            
        # put the gaussian kernel into an array
        kernel = []
        for e in np.arange(-10.*sigma, 10.*sigma+minstep, minstep):
            kernel.append(gauss(e, sigma, 0.))
        kernel = np.array(kernel)
        # convolve
        #convycst = scipy.signal.fftconvolve(ycst, kernel, 'same')
        convycst = sp.convolve(ycst, kernel, 'same')
        # and put it again on the initial grid
        tck =  sp.interpolate.splrep(xcst, convycst)
        #convy = interpolate.splev(self.wl[::-1], tck, der=0)
        #self.photconv = convy[::-1]/kernel.sum()
        self.photconv = self.phot
        #convthreshold = 1.e-4
        #for i in range(len(self.wl)):
        #    print i, '/', len(self.wl)
        #    mysum = 0.
        #    myerr = 0.
        #    mynorm = 0.
        #    for j in range(len(self.wl)-1):
        #        jj = j+1
        #        g =  gauss(self.wl[j], sigma, self.wl[i])
        #        gg =  gauss(self.wl[jj], sigma, self.wl[i])
        #        if g > convthreshold or gg > convthreshold:
        #            mysum += (self.phot[j]*g + self.phot[jj]*gg) * fabs(self.wl[j]-self.wl[jj])
        #            mynorm += (g+gg) * math.fabs(self.wl[j]-self.wl[jj])
        #    self.photconv.append(mysum/mynorm)
        if doPlot:
            pl.plot(self.wl, self.photconv)            
            pl.plot(self.wl, self.photconv/kernel.sum())
            pl.legend(["sky","raw conv. by telescope","conv. normalisation" ])
            print "par là"         
            
            
    def AddCalibSys(self, calibsys):
        # add flux calibration systematics to the calibrated spectrum
        sn = []
        for i in range(len(self.wlccd)):
            if self.caldEdlnoise[i] != 0.:
                sn.append(self.caldEdl[i]/self.caldEdlnoise[i])
            else:
                sn.append(0.)
        sn = np.array(sn)
        dev = calibsys.GetRandomRealization(self.wlccd[::-1], self.caldEdl[::-1], sn[::-1])
        self.syscaldEdl = self.caldEdl + dev[::-1]
        self.syscaldEdn = self.caldEdn + dev
        if S_LevelPlot > 2:
            pl.figure()        
            pl.xlabel("? ")
            pl.ylabel("?")
            pl.grid()
            pl.title("s/n compute by AddCalibSys method")
            pl.subplot(121)
            pl.plot(sn)
            pl.subplot(122)
            pl.title("dev compute by calibsys.GetRandomRealization method")
            pl.plot(dev)
            
            
    def EstimTransmissionATM(self,atm, ccdresp, grismresp, mirrorresp, gain=1., effarea=1.0, exptime=240.):
        EstTransATM = []
        # points outside the main absorption lines of the ATM, from Thomas Matheson IDL procedures, slighly modified for convergence of the interpolation
        finterp = sp.interpolate.interp1d(ccdresp.getnu(), ccdresp.getrenu(), bounds_error=False, fill_value = 0.)
        ginterp = sp.interpolate.interp1d(mirrorresp.getnu(), mirrorresp.getrenu(), bounds_error=False, fill_value = 0.)
        sinterp = sp.interpolate.interp1d(self.nu, self.dEdnConv, bounds_error=False, fill_value = 0.)
        ccdr = finterp(self.nuccd)
        mirrortr = ginterp(self.nuccd)
        StarEC = sinterp(self.nuccd)
        for i in range(len(self.nuccd)): 
        #for i in range(2): 
            if i == 0:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = dnccdd
            elif i == len(self.nuccd)-1:
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.
                dnccdd = dnccdg
            else:
                dnccdd = (self.nuccd[i+1]-self.nuccd[i])/2.
                dnccdg = (self.nuccd[i]-self.nuccd[i-1])/2.           
            num = self.elecccd[i]*hplanck*self.nuccd[i]
            # dem
            gtr = grismresp.getT(clight/self.nuccd[i])
            # pas de prise en compte de self.factorslit*self.factorextract dans dEdnConv
            dem = ccdr[i]*gtr*mirrortr[i]*exptime*StarEC[i]*(dnccdg+dnccdd)
            #print 'res:',num,dem
            EstTransATM.append(num/dem)
        self._EstTransATM = np.array(EstTransATM)
        if S_LevelPlot > 3 : 
            pl.figure()        
            pl.xlabel("nm")
            pl.ylabel("%")
            pl.grid()    
            pl.title("estimate atmosphere transmission with star reference above atmosphere")
            pl.plot(atm[1], atm[2]*100)
            pl.plot(clight/self.nuccd, self._EstTransATM*100)       
            pl.legend(["atm 'true'", "atm estimated"])


        
