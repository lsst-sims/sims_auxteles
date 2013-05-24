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
import sys
import os

# and class definitions
import starspectrum as star
import ATMCombSpec as atm
import mirror as mir
import grism as grm
import ccd 
import CalibSys as cal 
import scipy as sp 
import constants as cst
import numpy as np
import pylab as pl


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
        self.seeing = 0.1  # [rad]
        self.slitwidth = 1 # [rad]
        self.mirAera = (1.2**2)*3.14159 #[m^2]
        self.inputres = 1300.0 # [unitless]
        self.outputres = 600.0 # [unitless]
        self.nbpixel = 512
        self.WLref = 650.0 # [nm]
        
        
    def setResolution(self, inputRes, outputRes):
        self.inputres = inputRes
        self.outputres = outputRes   
        
        
    def ajustNBPixel(self, wl=None):
        if wl==None:
            wl = self.star.wl
        amplWL  = np.fabs(wl[0] - wl[-1])
        deltaWL = self.WLref / self.outputres
        self.nbpixel = int(np.trunc( amplWL / deltaWL ) + 1)
        
        
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
        if seeing != None: self.seeing = seeing
        self.atms.setWithTrans(wl, trans)
    
    
    def setSpectro(self, slitwidth ):
        """
        rad
        """
        self.slitwidth = slitwidth
                
                
    def setTelescope(self, mirorAera):
        """
        m^2
        """
        self.mirAera = mirorAera
                
                
    def setStarFile(self, starFile):
        self.star.readdEdl(starFile)


    def setStarSpectrum(self, wl, spectrum):
        """
        wl       : [?] wave lenght
        spectrum : [?]
        """        
        self.star.setSpectrum(wl, spectrum)
    
    
    def getCalibFlux(self, TpsExpo=240, addNoise=True):
        # apply the aperture correction with respect to
        # the slit width and the seeing
        if self.doPlot : self.star.plotWL("raw spectrum")
        print "apertureCorrection"
        self.star.apertureCorrection(self.seeing, self.slitwidth)
        if self.doPlot : self.star.plotWL("spectrum after apertureCorrection")        
        # convolution of the spectrum to the resolution of the spectrograph
        sigmaspectrograph = self.WLref/self.outputres
        sigmaspec         = self.WLref/self.inputres
        sigmaconv = (sigmaspectrograph**2 - sigmaspec**2)**.5
        if(sigmaconv<0):
            print 'error: convolution is impossible, the resolution of the input spectrum is less than the resolution of the output spectrum'
            exit(2)
        print "convolvePhotondEdl"
        self.star.convolvePhotondEdl(sigmaconv)
        if self.doPlot : self.star.plotWL("dEdl after convolvePhoton", dEdl=self.star.dEdlConv)            
        # conversion to photons
        print "computePhoton"
        self.star.computePhoton(self.mirAera, TpsExpo)
        if self.doPlot : self.star.plotNbPhotons("#photons after computePhoton", self.star.phot)    
        # rebin on the ccd grid
        print "computePhotonCCD ", self.nbpixel
        self.star.computePhotonCCD(self.gri, nbpixel=self.nbpixel)
        if self.doPlot : self.star.plotNbPhotons("#photon after computePhotonCCD", self.star.photccd)                
        # compute and add photon noise
        print "computePhotonNoise"   
        self.star.computePhotonNoise()
        if self.doPlot : self.star.plotNbPhotons("photon noise computePhotonNoise", self.star.photccdnoise)                        
        if addNoise :
            print "addPhotonNoise"
            self.star.addPhotonNoise()
            if self.doPlot : self.star.plotNbPhotons("#photon after addPhotonNoise", self.star.photccd)        
        # for checking
        #star.write('photccd', starname+'photccd.dat')    
        # through the atmosphere
        print "computePhotonATM"
        self.star.computePhotonATM(self.atms.gettrans())
        if self.doPlot : self.star.plotNbPhotons("photatm computePhotonATM", self.star.photatm)        
        # for checking
        #star.write('photatm', starname+'photatm.dat')    
        # reflexion on the mirror
        print "computePhotonMirror"
        self.star.computePhotonMirror(self.mirrortr)
        if self.doPlot : 
            self.star.plotNbPhotons("mirrortr by computePhotonMirror", self.star.mirrortr)        
            self.star.plotNbPhotons("photmirorr by computePhotonMirror", self.star.photmirror)        
        # for checking
        #star.write('photmirrsetStarFileor', starname+'photmirror.dat')    
        # through the Grism
        print "computePhotonGrism"
        self.star.computePhotonGrism(self.gri)
        if self.doPlot : 
            self.star.plotNbPhotons("grismTrans computePhotonGrism", self.star.grismTrans)
            self.star.plotNbPhotons("photgrism computePhotonGcalibrateSpecrism", self.star.photgrism)    
        # conversion in photo-electrons
        print "computeElectronCCD"
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
        print "MakeSensFunc"
        #star.MakeSensFunc(atmc.gettrans(), ccdqe, gri, mirrortr)
        #if True : star.plotWL("sensdEdl MakeSensFunc", star.wlccd[50:-50], star.sensdEdl[50:-50])
        self.star.MakeSensFuncInstruOnly(self.ccdqe, self.gri, self.mirrortr, effarea=self.mirAera, exptime=TpsExpo)        
        # flux-calibrate the spectrum 
        print "calibrateSensSpec"
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
        
    


def finalPlot(pstar, atms, pCCD, pMir, pGrism):
    assert isinstance(pstar, star.starspectrum)
    assert isinstance(atms, atm.ATMCombSpec)
    assert isinstance(pCCD, ccd.ccd )
    assert isinstance(pMir, mir.mirror )    
    assert isinstance(pGrism, grm.grism )    
    pl.figure()
    pl.xlabel("wavelength nm")
    pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')
    pl.grid()    
    pl.title("Star flux above atmosphere and observed")
    pl.plot(pstar.wl, pstar.dEdl)
    pl.plot(pstar.wlccd, pstar.syscaldEdl)
    pl.plot(pstar.wlccd, pstar.elecccd*cst.hplanck*pstar.nuccd**3/(cst.clight*240*pstar._deltaNuCCD))
    pl.legend(["above atm ", "observed calibrated", "observed raw"])
    # transmission
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
    doPlot = False
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

    # chech if the ATSIM_DATA environment variable is set
    datadir = os.environ.get('ATSIM_DATA')
    if datadir == None:
        print 'you must run ". init.bash" to configure the data directory'
        exit(2)

    # get the file names
    starname = args[0]
    atmstar = args[1]
    atmcalib = args[2]

    # create objects and read the data
    ccdqe = ccd.ccd()
    gri = grm.grism()
    mirrortr = mir.mirror()
    calibsys = cal.CalibSys()
    atms = atm.ATMCombSpec()
    atmc = atm.ATMCombSpec()
    mystar = star.starspectrum()
    diamTeles = 1
    effarea = (np.pi/4)* diamTeles**2
    
    mirrortr.read(datadir+'/'+ os.environ.get('ATSIM_MIRROR'))
    ccdqe.read(datadir+'/'+ os.environ.get('ATSIM_CCDQE'))
    calibsys.read(datadir+'/'+ os.environ.get('ATSIM_SYSCORR'))

    # read the atmspheres
    atms.read(atmstar)
    atmc.read(atmcalib)
    # read the star spectrum given in J/m2/s/nm
    mystar.readdEdl(starname)
    if doPlot : mystar.plotWL("raw spectrum")
    
    # modify the star file name to identify the atmosphere
    starname += '_'+str(atmstar)+'_'
    # for checking
    mystar.write('dEdn', starname+'dEdn.dat')

    # apply the aperture correction with respect to
    # the slit width and the seeing
    print "apertureCorrection"
    mystar.apertureCorrection(seeing, slitwidth)
    if doPlot : mystar.plotWL("spectrum after apertureCorrection")
    

    # convolution of the spectrum to the resolution of the spectrograph
    sigmaspectrograph = 650./400.
    sigmaspec = 650./inputres
    sigmaconv = (sigmaspectrograph**2 - sigmaspec**2)**.5
    if(sigmaconv<0):
        print 'error: convolution is impossible, the resolution of the input spectrum is less than the resolution of the output spectrum'
        exit(2)
    print "convolvePhoton"
    mystar.convolvePhotondEdl(sigmaconv)
    if doPlot : mystar.plotWL("dEdl after convolvePhoton", dEdl=star.dEdlConv)

    
    # conversion to photons
    print "computePhoton"
    mystar.computePhoton(effarea)
    if True : mystar.plotNbPhotons("#photons after computePhoton", mystar.phot/0.68)


    # for checking
    #star.write('phot', starname+'phot.dat')


    # for checking
    #star.write('photconv', starname+'photconv.dat')

    # rebin on the ccd grid
    print "computePhotonCCD"
    mystar.computePhotonCCD(gri)
    if doPlot : mystar.plotNbPhotons("#photon after computePhotonCCD", star.photccd)
    
    
    # compute and add photon noise
    print "computePhotonNoise"   
    mystar.computePhotonNoise()
    if doPlot : mystar.plotNbPhotons("photon noise computePhotonNoise", star.photccdnoise)
    
    
    print "addPhotonNoise"
    mystar.addPhotonNoise()
    if doPlot : mystar.plotNbPhotons("#photon after addPhotonNoise", mystar.photccd)
    
    # for checking
    #star.write('photccd', starname+'photccd.dat')

    # through the atmosphere
    print "computePhotonATM"
    mystar.computePhotonATM(atms.gettrans())
    if doPlot : mystar.plotNbPhotons("photatm computePhotonATM", star.photatm)
    
    # for checking
    #star.write('photatm', starname+'photatm.dat')

    # reflexion on the mirror
    print "computePhotonMirror"
    mystar.computePhotonMirror(mirrortr)
    if doPlot : 
        mystar.plotNbPhotons("mirrortr by computePhotonMirror", mystar.mirrortr)        
        mystar.plotNbPhotons("photmirorr by computePhotonMirror", mystar.photmirror)
    
    # for checking
    #star.write('photmirror', starname+'photmirror.dat')

    # through the Grism
    print "computePhotonGrism"
    mystar.computePhotonGrism(gri)
    if doPlot : 
        mystar.plotNbPhotons("grismTrans computePhotonGrism", mystar.grismTrans)
        mystar.plotNbPhotons("photgrism computePhotonGrism", mystar.photgrism)

    # conversion in photo-electrons
    print "computeElectronCCD"
    mystar.computeElectronCCD(ccdqe)
    if doPlot : mystar.plotNbPhotons("elecccd computeElectronCCD", mystar.elecccd)

    # add a typical read noise
    print "addNoiseElectronCCD"
    mystar.addNoiseElectronCCD()
    if True : mystar.plotToWLnm( mystar.nuccd, mystar.elecccd, "jmc Final estimation photon-electon by pixel")
    
    
    # for checking
    #mystar.write('elecccdlambda', starname+'elecccdlambda.dat')
    #star.write('elecccd', starname+'elecccdnu.dat')

    # Make the sensitivity function to simulate the flux calibration
    print "MakeSensFunc"
    #star.MakeSensFunc(atmc.gettrans(), ccdqe, gri, mirrortr)
    #if True : star.plotWL("sensdEdl MakeSensFunc", star.wlccd[50:-50], star.sensdEdl[50:-50])
    mystar.MakeSensFuncJMC(atmc.gettrans(), ccdqe, gri, mirrortr, effarea= effarea)
    
    
    # modify the star file name to identify the calibration atmosphere
    starname += '_'+str(atmcalib)+'_'
    # for checking
    #star.write('atmrebin', starname+'atmrebin.dat')
    #star.write('sensdEdl', starname+'sensdEdl.dat')
    #star.write('sensdEdn', starname+'sensdEdn.dat')

    # flux-calibrate the spectrum 
    print "calibrateSensSpec"
    mystar.calibrateSensSpec()
    if doPlot : mystar.plotWL("caldEdl calibrateSensSpec", mystar.wlccd, mystar.caldEdl)
    
    #star.write('caldEdn', starname+'caldEdn.dat')
    #star.write('caldEdl', starname+'caldEdl.dat')

    # add systematic effect due to the flux calibration
    print "AddCalibSys"
    mystar.AddCalibSys(calibsys)
    if doPlot : mystar.plotWL("syscaldEdl AddCalibSys", mystar.wlccd, mystar.syscaldEdl)
    
    # output of the simulation
    mystar.write('syscaldEdn', starname+'syscaldEdn.dat')
    mystar.write('syscaldEdl', starname+'syscaldEdl.dat')

    # 
    #star.EstimTransmissionATM(atms.gettrans(), ccdqe, gri, mirrortr, )

    finalPlot(mystar,atms, ccdqe, mirrortr, gri)

if __name__ == "__main__":
    #pl.ioff()
    #pl.matplotlib.use('GTKAgg')
    main()

try:
    pl.show()
except AttributeError:
    pass
