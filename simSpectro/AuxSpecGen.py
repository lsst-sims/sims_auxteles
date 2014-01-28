#!/usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
#        AuxSpecGen - Generation of LSST Auxiliary Telescope spectra        #
#                                                                           #
#        author:          Gurvan BAZIN                                      #
#                         Universitäts-Sternwarte-München                   #
#        latest revision: July 2010                                         #
#############################################################################


# include some useful libraries
import getopt
import os
import sys

import ATMCombSpec as atm
import CalibSys as cal
import ccd
from constants import hplanck
import constants as cst
import grism as grm
import mirror as mir
import numpy as np
import matplotlib.pyplot as pl
import scipy as sp
import scipy.interpolate as spi
import starspectrum as star
import tools as tl


# and class definitions
def getDirectory(path):
    """
    path=/path/to/my/file.xx
    return 
    /path/to/my
    """
    ipath = path[::-1]
    idx = ipath.find('/')
    return ipath[idx+1:][::-1]


def getModuleDirectory():
    """
    get directory of the present module
    """
    path = __file__
    return getDirectory(path)



class AuxTeles(object):
    """
    simulation spectrum star acquisition throught atmosphere with telescope spectroscope and CDD
    """
  
    def __init__(self):
        self.ccdqe = ccd.ccd()
        self.gri = grm.grism()
        self.mirrortr = mir.mirror()
        self.calibsys = cal.CalibSys()
        self.atms = atm.ATMCombSpec()    
        self.star = star.starspectrum()
        datadir = getModuleDirectory() + '/data/'        
        self.mirrortr.read(datadir + 'CFHT_Primary_Transmission.txt')
        self.ccdqe.read(datadir    + 'E2V_CCD42-90_QEmodel.txt')
        self.calibsys.read(datadir + 'corr.txt')
        self.doPlot = False
        self.seeing = 1  # [rad]
        self.slitwidth = 1.3 # [rad]
        self.setDiamTelescope(1.2) #[m^2]
        #self.inputres = 1300.0 # [unitless]
        #self.outputres = 600.0 # [unitless]
        self.inputres = 1 # [unitless]
        self.outputres = 2 # [unitless]
        self.nbpixel = 512
        self.WLref = 650.0 # [nm]
        self.instruEff = None
        self.verbose = 0
        
        
    def setResolution(self, inputRes, outputRes):
        self.inputres = inputRes
        self.outputres = outputRes   
        
        
    def ajustNBPixel(self, wl=None):
        if wl==None:
            wl = self.star.wl
        amplWL  = np.fabs(wl[0] - wl[-1])
        deltaWL = self.WLref / self.outputres
        self.nbpixel = 2*int(np.trunc( amplWL / deltaWL ) + 1)
        
        
    def setResolutionAndAjustNBPixel(self, inputRes, outputRes):
        """
        ajust number pixel with star spectrum wavelength definition
        """
        self.setResolution(inputRes, outputRes)
        self.ajustNBPixel(self.star.wl)
        
        
    def setAtmFile(self, atmstar, seeing=None):
        if seeing != None: self.seeing = seeing
        self.atms.read(atmstar)
        
        
    def setAtmTrans(self, wl, trans, seeing=None):
        """
        wl []
        """
        if seeing != None: self.seeing = seeing
        self.atms.setWithTrans(wl, trans)
    
    
    def setSlit(self, slitwidth ):
        """
        rad
        """
        self.slitwidth = slitwidth
                          
    def setTelescope(self, mirorAera):
        """
        m^2
        """
        self.mirAera = mirorAera
        
        
    def setDiamTelescope(self, Diammiror):
        """
        m
        """
        self.mirAera = ((Diammiror/2)**2)*3.14159
                        
                
    def setStarFile(self, starFile, facEC=1.0):
        self.star.readdEdl(starFile, facEC)


    def setStarSpectrum(self, wl, spectrum):
        """
        wl       : [nm] wave length
        spectrum : [J.m^{-2}.s^{-1}.nm^{-1}]
        """        
        self.star.setSpectrum(wl, spectrum)
    
    
    def computePhotonElec(self,  TpsExpo=240):
        """
        with only photon noise
        """
        self.getFluxWithStarLowResolution(TpsExpo, addNoise=False)
        self.sigmaPhotNoise = np.sqrt(self.star.elecccd)  
        self.sigmaPhotNoise = map(star.notzero, self.sigmaPhotNoise)
        noise = np.random.normal(0., self.sigmaPhotNoise)
        self.star.elecccd += noise
        

    def addPhotonNoiseAndSigma(self, aFluxWmnm):
        aNP = self._Wmnm2elec(aFluxWmnm)
        sigmaPhotNoise = np.sqrt(aNP)
        noise = np.random.normal(0., sigmaPhotNoise)
        aSignalWmnm = self._elec2Wmnm(aNP+noise)
        aSigmaWmnm = self._elec2Wmnm(sigmaPhotNoise)
        print 'raw snr ', (aNP+noise).mean()/noise.std()
        return aSignalWmnm, aSigmaWmnm
        
               
    def getCalibFlux(self, TpsExpo=240, addNoise=True):
        # apply the aperture correction with respect to
        # the slit width and the seeing
        self.TpsExpo =  TpsExpo
        if self.doPlot : self.star.plotWL("raw spectrum")
        if self.verbose > 0 : print "apertureCorrection"
        self.star.apertureCorrection(self.seeing, self.slitwidth)
        if self.doPlot : self.star.plotWL("spectrum after apertureCorrection", dEdl=self.star.dEdlAper)        
        # convolution of the spectrum to the resolution of the spectrograph
        sigmaspectrograph = self.WLref/self.outputres
        sigmaspec         = self.WLref/self.inputres
        sigmaconv = (sigmaspectrograph**2 - sigmaspec**2)**.5
        if(sigmaconv<0):
            print 'error: convolution is impossible, the resolution of the input spectrum is less than the resolution of the output spectrum'
            exit(2)
        if self.verbose > 0 : print "convolvePhotondEdl"
        self.star.convolvePhotondEdl(sigmaconv)
        if self.doPlot : self.star.plotWL("dEdl after convolvePhoton", dEdl=self.star.dEdlConv)            
        # conversion to photons
        if self.verbose > 0 : print "computePhoton"
        self.star.computePhoton(self.mirAera, TpsExpo)        
        if self.doPlot : self.star.plotNbPhotons("#photons after computePhoton", self.star.phot, self.star.wl)    
        # rebin on the ccd grid
        if self.verbose > 0 : print "computePhotonCCD ", self.nbpixel        
        self.star.computePhotonCCD(self.gri, nbpixel=self.nbpixel)
        if self.doPlot : self.star.plotNbPhotons("#photon after computePhotonCCD", self.star.photccd)                
        # compute and add photon noise
        if self.verbose > 0 : print "computePhotonNoise"   
        self.star.computePhotonNoise()
        if addNoise :                        
            if self.doPlot : self.star.plotNbPhotons("photon noise computePhotonNoise", self.star.photccdnoise)                        
            print "addPhotonNoise"
            self.star.addPhotonNoise()
            if self.doPlot : self.star.plotNbPhotons("#photon after addPhotonNoise", self.star.photccd)
        # for checking
        #star.write('photccd', starname+'photccd.dat')    
        # through the atmosphere
        if self.verbose > 0 : print "computePhotonATM"
        self.star.computePhotonATM(self.atms.gettrans())
        if self.doPlot : self.star.plotNbPhotons("photatm computePhotonATM", self.star.photatm)        
        # for checking
        #star.write('photatm', starname+'photatm.dat')    
        # reflexion on the mirror
        if self.verbose > 0 : print "computePhotonMirror"
        self.star.computePhotonMirror(self.mirrortr)
        if self.doPlot : 
            self.star.plotNbPhotons("mirrortr by computePhotonMirror", self.star.mirrortr)        
            self.star.plotNbPhotons("photmirorr by computePhotonMirror", self.star.photmirror)        
        # for checking
        #star.write('photmirrsetStarFileor', starname+'photmirror.dat')    
        # through the Grism
        if self.verbose > 0 : print "computePhotonGrism"
        self.star.computePhotonGrism(self.gri)
        if self.doPlot : 
            self.star.plotNbPhotons("grismTrans computePhotonGrism", self.star.grismTrans)
            self.star.plotNbPhotons("photgrism computePhotonGcalibrateSpecrism", self.star.photgrism)    
        # conversion in photo-electrons
        if self.verbose > 0 : print "computeElectronCCD"
        self.star.computeElectronCCD(self.ccdqe)
        if self.doPlot : self.star.plotNbPhotons("elecccd computeElectronCCD", self.star.elecccd)            
        # add a typical read noise       
        if addNoise:
            print "addNoiseElectronCCD"
            self.star.addNoiseElectronCCD()
            if self.doPlot : self.star.plotToWLnm( self.star.nuccd, self.star.elecccd, "Final estimation photon-electon by pixel")                
        # for checking
        #star.write('elecccdlambda', starname+'elecccdlambda.dat')
        #star.write('elecccd', starname+'elecccdnu.dat')    
        # Make the sensitivity function to simulate the flux calibration
        if self.verbose > 0 : print "MakeSensFunc"
        #star.MakeSensFunc(atmc.gettrans(), ccdqe, gri, mirrortr)
        #if True : star.plotWL("sensdEdl MakeSensFunc", star.wlccd[50:-50], star.sensdEdl[50:-50])
        self.star.MakeSensFuncInstruOnly(self.ccdqe, self.gri, self.mirrortr, effarea=self.mirAera, exptime=TpsExpo)        
        # flux-calibrate the spectrum 
        if self.verbose > 0 : print "calibrateSensSpec"
        self.star.calibrateSensSpec()
        if self.doPlot : self.star.plotWL("caldEdl calibrateSensSpec", self.star.wlccd, self.star.caldEdl)        
        #star.write('caldEdn', starname+'caldEdn.dat')
        #star.write('caldEdl', starname+'caldEdl.dat')    
        # add systematic effect due to the flux calibration
        if addNoise:
            print "AddCalibSys"
            self.star.AddCalibSys(self.calibsys)
            if self.doPlot : self.star.plotWL("syscaldEdl AddCalibSys Final", self.star.wlccd, self.star.syscaldEdl)
        else:
            self.star.syscaldEdl = self.star.caldEdl
            self.star.syscaldEdn = self.star.caldEdn                   
        #self.star.plotWL("syscaldEdl AddCalibSys Final", self.star.wlccd, self.star.syscaldEdl)    
        if self.doPlot : self.star.plotWL("raw spectrum")    
        return self.star.wlccd, self.star.syscaldEdl
    
    
    def getFluxWithStarLowResolution(self, TpsExpo=240, addNoise=True):
        """
        like getCalibFlux() but with star low resolution (like Kurucz model done) so the 
        signifying resolution is atmopheric one
        
        """
        # apply the aperture correction with respect to
        # the slit width and the seeing
        self.TpsExpo =  TpsExpo
        if self.doPlot : self.star.plotWL("raw spectrum")
        if self.verbose > 0 : print "apertureCorrection"
        self.star.apertureCorrection(self.seeing, self.slitwidth)
        if self.doPlot : self.star.plotWL("spectrum after apertureCorrection", dEdl=self.star.dEdlAper)        
        # convolution of the atmosphere transmission 
        sigmaspectrograph = self.WLref/self.outputres
        sigmaspec         = self.WLref/self.inputres
        sigmaconv = (sigmaspectrograph**2 - sigmaspec**2)**.5
        if(sigmaconv<0):
            print 'error: convolution is impossible, the resolution of the input spectrum is less than the resolution of the output spectrum'
            exit(2)
        if self.verbose > 0 : print "convolveAtm"
        atmConvol = self.atms.getTransDowngrade( self.inputres, self.outputres )
        #if self.doPlot : self.atms.plotAtm(atmConvol)
        self.star.prodWithConvolArray(self.atms.wl[::-1], atmConvol[::-1], self.doPlot)     
        if self.doPlot : self.star.plotWL("dEdl after convolvePhoton", dEdl=self.star.dEdlConv)    
        # conversion to photons
        if self.verbose > 0 : print "computePhoton"
        self.star.computePhoton(self.mirAera, TpsExpo)        
        if self.doPlot : self.star.plotNbPhotons("#photons after computePhoton", self.star.phot, self.star.wl)    
        # rebin on the ccd grid
        if self.verbose > 0 : print "computePhotonCCD ", self.nbpixel        
        self.star.computePhotonCCD(self.gri, nbpixel=self.nbpixel)
        if self.doPlot : self.star.plotNbPhotons("#photon after computePhotonCCD", self.star.photccd)                
        # compute and add photon noise
        if self.verbose > 0 : print "computePhotonNoise"   
        self.star.computePhotonNoise()
        if addNoise :                        
            if self.doPlot : self.star.plotNbPhotons("photon noise computePhotonNoise", self.star.photccdnoise)                        
            print "addPhotonNoise"
            self.star.addPhotonNoise()
            if self.doPlot : self.star.plotNbPhotons("#photon after addPhotonNoise", self.star.photccd)
        # replace computePhotonATM
        self.star.photatm = self.star.photccd.copy()
        self.star.photatmnoise = np.array(self.star.photccdnoise).copy()      
        if self.verbose > 0 : print "computePhotonMirror"        
        self.star.computePhotonMirror(self.mirrortr)
        if self.doPlot : 
            self.star.plotNbPhotons("mirrortr by computePhotonMirror", self.star.mirrortr)        
            self.star.plotNbPhotons("photmirorr by computePhotonMirror", self.star.photmirror)        
        # for checking
        #star.write('photmirrsetStarFileor', starname+'photmirror.dat')    
        # through the Grism
        if self.verbose > 0 : print "computePhotonGrism"
        self.star.computePhotonGrism(self.gri)
        if self.doPlot : 
            self.star.plotNbPhotons("grismTrans computePhotonGrism", self.star.grismTrans)
            self.star.plotNbPhotons("photgrism computePhotonGcalibrateSpecrism", self.star.photgrism)    
        # conversion in photo-electrons
        if self.verbose > 0 : print "computeElectronCCD"
        self.star.computeElectronCCD(self.ccdqe)
        if self.doPlot : self.star.plotNbPhotons("elecccd computeElectronCCD", self.star.elecccd)            
        # add a typical read noise       
        if addNoise:
            print "addNoiseElectronCCD"
            self.star.addNoiseElectronCCD()
            if self.doPlot : self.star.plotToWLnm( self.star.nuccd, self.star.elecccd, "Final estimation photon-electon by pixel")                
        # for checking
        #star.write('elecccdlambda', starname+'elecccdlambda.dat')
        #star.write('elecccd', starname+'elecccdnu.dat')    
        # Make the sensitivity function to simulate the flux calibration
        if self.verbose > 0 : print "MakeSensFunc"
        #star.MakeSensFunc(atmc.gettrans(), ccdqe, gri, mirrortr)
        #if True : star.plotWL("sensdEdl MakeSensFunc", star.wlccd[50:-50], star.sensdEdl[50:-50])
        self.star.MakeSensFuncInstruOnly(self.ccdqe, self.gri, self.mirrortr, effarea=self.mirAera, exptime=TpsExpo)        
        # flux-calibrate the spectrum 
        if self.verbose > 0 : print "calibrateSensSpec"
        self.star.calibrateSensSpec()
        if self.doPlot : self.star.plotWL("caldEdl calibrateSensSpec", self.star.wlccd, self.star.caldEdl)        
        #star.write('caldEdn', starname+'caldEdn.dat')
        #star.write('caldEdl', starname+'caldEdl.dat')    
        # add systematic effect due to the flux calibration
        if addNoise:
            print "AddCalibSys"
            self.star.AddCalibSys(self.calibsys)
            if self.doPlot : self.star.plotWL("syscaldEdl AddCalibSys Final", self.star.wlccd, self.star.syscaldEdl)
        else:
            self.star.syscaldEdl = self.star.caldEdl
            self.star.syscaldEdn = self.star.caldEdn                   
        #self.star.plotWL("syscaldEdl AddCalibSys Final", self.star.wlccd, self.star.syscaldEdl)    
        #if self.doPlot : self.star.plotWL("raw spectrum")    
    
    def getWLccd(self):
        """
        [nm]
        """
        return self.star.wlccd
    
    
    def getRawFlux_photonElectron(self):
        return self.star.elecccd
    
    
    def plotRawFlux_photonElectron(self):
        self.star.plotNbPhotons("photonElectron", self.star.elecccd) 
    
    
    def getSigmaPhotonNoise(self): 
        return self.sigmaPhotNoise
    
    
    def getSigmaPhotonNoise_Wmnm(self):
        ret = self._elec2Wmnm(self.sigmaPhotNoise) 
        return ret 
    
    
    def getInstruEfficiency(self, doPlot=False):
        ccdresp = self.ccdqe
        grismresp = self.gri
        mirrorresp = self.mirrortr
        if self.instruEff == None:            
            eff = []        
            #
            # Efficiency wavelength depend
            #
            # CCD
            tck = spi.splrep(ccdresp.getnu(), ccdresp.getrenu())
            TransCCDPixel = spi.splev(self.star.nuccd, tck)  
            # MIROR
            tck = spi.splrep(mirrorresp.getnu(), mirrorresp.getrenu())
            TransMIRPixel = spi.splev(self.star.nuccd, tck)
            # glass prism
            for i in range(len(self.star.nuccd)):
                # -> J.m-2.s-1
                if TransCCDPixel[i] > 0.  and TransMIRPixel[i] > 0.:
                    gtr = grismresp.getT(cst.clight/self.star.nuccd[i])
                    TransInstru= gtr*TransCCDPixel[i]*TransMIRPixel[i]                                
                else:
                    TransInstru= 0
                eff.append(TransInstru)
            self.instruEff = np.array( eff )
            #
            # Efficiency with no wavelength depend
            #
            self.instruEff *= self.star.factorextract 
            self.instruEff *= self.star.factorslit 
            self.instruEff *= grismresp.effFirstOrder
        if doPlot:
            eff = []
            for i in range(len(self.star.nuccd)):
                # -> J.m-2.s-1
                if TransCCDPixel[i] > 0.  and TransMIRPixel[i] > 0.:
                    gtr = grismresp.getT(cst.clight/self.star.nuccd[i])
                    TransInstru= gtr                                
                else:
                    TransInstru= 0
                eff.append(TransInstru)
            TransGlassGrism = np.array( eff ) 
            #print len(self.star.nuccd),     len(self.star.wlccd), len(TransGlassGrism)    
            cstSlitSeeing = np.ones(len(self.star.wlccd))*self.star.factorslit
            cstFirstOrderGrism = np.ones(len(self.star.wlccd))*grismresp.effFirstOrder
            pl.figure()
            pl.plot(self.star.wlccd, cstSlitSeeing*100)
            pl.plot(mirrorresp.getwl(), mirrorresp.getre()*100)
            pl.plot(self.star.wlccd, TransGlassGrism*100)
            pl.plot(self.star.wlccd, cstFirstOrderGrism*100)
            pl.plot(ccdresp.getwl(), ccdresp.getre()*100)                        
            pl.plot(self.star.wlccd, self.instruEff*100)
            pl.title("Instrumental efficiency")
            pl.xlabel("nm")
            pl.ylabel("%")
            pl.legend(["slit seeing","mirror",  "glass grism","blaze 1 order", "CCD", "Total"],loc=7)
            pl.grid()
        return self.instruEff
    

    def _elec2Wmnm(self, aElec, doPlot=False):
        Coef = self.star.electron2Jnu()
        ret  = aElec*Coef
        ret *= self.star.nuccd**2/cst.clight         
        #print self.mirAera, self.TpsExpo
        ret /= self.mirAera*self.TpsExpo
        if doPlot:
            pl.figure()
            pl.title("_elec2Wmnm")
            pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')        
            pl.plot(self.getWLccd(), ret)       
            pl.grid()
        return ret 
    
        
    def _Wmnm2elec(self, aWmnm, doPlot=False):
        Coef = self.star.electron2Jnu()
        ret  = (aWmnm*cst.clight)*(self.mirAera*self.TpsExpo)
        ret /= (Coef*(self.star.nuccd**2))         
        return ret 
    
  
    def getRawFlux_Wmnm(self):
        """
        return raw flux in W.m^{-2}.nm^{-1}
        """
        ret = self._elec2Wmnm(self.star.elecccd)
        return ret
        
        
    def getMeanTimeExposForSNR(self, snr, wlmin=None, wlmax =None):
        """
        for photon noise dominant : snr = s/sqrt(s) = sqrt(s)
        snr = sqrt(<elec>*Tsnr/T)
        
        return Tsnr
        """
        if wlmin == None:
            wlmin = self.star.wlccd[::-1][0]
        if wlmax == None:
            wlmax = self.star.wlccd[::-1][-1]
        iMin, iMax = tl.indexInIntervalCheck(self.star.wlccd[::-1], wlmin, wlmax)
#        print self.star.wlccd[::-1][0], self.star.wlccd[::-1][-1]
#        print wlmin, wlmax
#        print iMin, iMax
        return (snr**2*self.TpsExpo)/self.star.elecccd[::-1][iMin: iMax].mean()
    
#
# Plot 
#

        


def finalPlot(pstar, atms, pCCD, pMir, pGrism):
    assert isinstance(pstar, star.starspectrum)
    assert isinstance(atms, atm.ATMCombSpec)
    assert isinstance(pCCD, ccd.ccd )
    assert isinstance(pMir, mir.mirror )    
    assert isinstance(pGrism, grm.grism )
    atms.plotAtm()    
    pl.figure()
    pl.xlabel("wavelength nm")
    pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')
    pl.grid()    
    pl.title("Star flux above atmosphere and observed")
    pl.plot(pstar.wl, pstar.dEdl)    
    pl.plot(pstar.wlccd, pstar.elecccd*cst.hplanck*pstar.nuccd**3/(cst.clight*240*pstar._deltaNuCCD))
    pl.legend(["above atm ", "observed raw"])
    # transmission
    if False:
        pl.figure()
        pl.xlabel("wavelength nm")
        pl.ylabel("%")
        pl.grid()    
        pl.title("Transmission functions")
        pl.plot(atms.wl, 100*atms.tr)
        pl.plot(pCCD.wl, 100*pCCD.re)   
        pl.plot(pMir.wl, 100*pMir.re)
        re = []
        for i in range(len(pstar.nuccd)):
            re.append(pGrism.getT(cst.clight/pstar.nuccd[i]))   
        re = np.array(re)
        pl.plot(cst.clight/pstar.nuccd, re*100)
        pl.legend(["atmosphere", "CCD", "telescope mirror" ,"Grism"], loc=7)
    # index refraction
#    pl.figure()
#    pl.xlabel("wavelength nm")
#    pl.grid()    
#    pl.title("Indices of refraction grism BK7")
#    aLambda = arange(200,1200,10) # nm
#    indice = [pGrism.getN(Lambda) for Lambda in aLambda]
#    pl.plot(aLambda, indice)
    
    # interpolation nuCCD
    tck = sp.interpolate.splrep(pstar.wl[::-1], pstar.dEdl[::-1])
    dEdlCCD = sp.interpolate.splev(pstar.wlccd, tck)
    # Ratio EC vrai EC estimé
    if False:
        pl.figure()
        pl.title("Ratio estimate flux / true flux")
        pl.xlabel("wavelength nm")
        pl.grid()
        # interpolation nuCCD
        ratio =  pstar.syscaldEdl / dEdlCCD 
        print "Ratio:"
        print "min , max : ",pstar.wl[50], pstar.wl[-50]
        print ratio[50:-50].mean(), ratio[50:-50].std()
        pl.plot(pstar.wlccd[10:-10], 1/ratio[10:-10])
        # error relative
        pl.figure()
        pl.title("Relative error true pstar flux and estimated")
        pl.xlabel("wavelength nm")
        pl.ylabel("%")
        pl.grid()
        # interpolation nuCCD
        #errRel = (dEdlCCD-pstar.syscaldEdl)/pstar.syscaldEdl    
        errRel = (dEdlCCD-pstar.syscaldEdl)/dEdlCCD    
        pl.plot(pstar.wlccd[10:-10], errRel[10:-10]*100)
    #
    pstar.plotToWLnm( pstar.nuccd, pstar.elecccd, "Final estimation photon-electon by pixel")    
         
         
# usage message
def usage():
    print 'usage : AuxSpecGen.py --seeing seeing --slitwidth slitwidth --inputres inputres'
    print '                      star starATM calibATM'
    print 'inputres: resolution of the input star spectrum'
    print 'star: spectrum of the star above the atmosphere'
    print 'starATM: atmosphere under which the star is observed'
    print 'calibATM: atmosphere under which the spectrophotometric calibration star is observed'

# main function
def main():

    # process the command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hr:x:w:v", ["help", "inputres=", "slitwidth=", "seeing="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    # default values
    inputres = 0.
    slitwidth = 1.
    seeing = 1.
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            
            usage()
            sys.exit()
        elif o in ("-r", "--inputres"):
            inputres = float(a)
        elif o in ("-x", "--slitwidth"):
            slitwidth = float(a)
        elif o in ("-w", "--seeing"):
            seeing = float(a)
        else:
            assert False, "unhandled option"

    # chech the number of parameters
    if len(args)<3:
        usage()
        exit(2)

    # get the file names
    starname = args[0]
    atmstar = args[1]
    #atmcalib = args[2]

    oTel = AuxTeles()
    oTel.setAtmFile(atmstar, seeing)
    oTel.setStarFile(starname)   
    oTel.setResolutionAndAjustNBPixel(inputres, 500.0)
    oTel.setSlit(slitwidth)
    oTel.getCalibFlux()
    oTel.getInstruEfficiency(True)
    finalPlot(oTel.star,oTel.atms, oTel.ccdqe, oTel.mirrortr, oTel.gri)


if __name__ == "__main__":
    #pl.ioff()
    #pl.matplotlib.use('GTKAgg')
    main()


try:
    pl.show()
except AttributeError:
    pass
