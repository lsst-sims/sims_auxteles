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
        self.seeing = 0.1
        self.slitwidth = 1
        self.mirAera = (1.2**2)*3.14159
        self.inputres = 410
        self.nbpixel  = 512
        
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
    
    def getCalibFlux(self, TpsExpo=240):
        # apply the aperture correction with respect to
        # the slit width and the seeing
        print "apertureCorrection"
        self.star.apertureCorrection(self.seeing, self.slitwidth)
        if self.doPlot : self.star.plotWL("spectrum after apertureCorrection")        
        # convolution of the spectrum to the resolution of the spectrograph
        sigmaspectrograph = 650./400.
        sigmaspec = 650./self.inputres
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
        print "computePhotonCCD"
        self.star.computePhotonCCD(self.gri, self.nbpixel)
        if self.doPlot : self.star.plotNbPhotons("#photon after computePhotonCCD", self.star.photccd)                
        # compute and add photon noise
        print "computePhotonNoise"   
        self.star.computePhotonNoise()
        if self.doPlot : self.star.plotNbPhotons("photon noise computePhotonNoise", self.star.photccdnoise)                
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
        print "AddCalibSys"
        self.star.AddCalibSys(self.calibsys)
        if self.doPlot : self.star.plotWL("syscaldEdl AddCalibSys Final", self.star.wlccd, self.star.syscaldEdl)        
        #self.star.plotWL("syscaldEdl AddCalibSys Final", self.star.wlccd, self.star.syscaldEdl)        
        return self.star.wlccd, self.star.syscaldEdl
        
        

def simuSpectroAuxTeles(inputres, slitwidth, seeing, starname, atmstar):
    # create objects and read the data
    ccdqe = ccd.ccd()
    gri = grm.grism()
    mirrortr = mir.mirror()
    calibsys = cal.CalibSys()
    atms = atm.ATMCombSpec()    
    ostar = star.starspectrum()
    doPlot = True
    
    datadir = getDirectory(__file__) + '/data/'
    
    mirrortr.read(datadir + 'CFHT_Primary_Transmission.txt')
    ccdqe.read(datadir    + 'E2V_CCD42-90_QEmodel.txt')
    calibsys.read(datadir + 'corr.txt')
    ccdqe.plotResponse()
    calibsys.plotCorrel()
    
    # read the atmspheres
    atms.read(atmstar)   
    # read the ostar spectrum given in J/m2/s/nm
    ostar.readdEdl(starname)
    if doPlot : ostar.plotWL("raw spectrum")
    
    # modify the ostar file name to identify the atmosphere
    starname += '_'+str(atmstar)+'_'
    # for checking
    #ostar.write('dEdn', starname+'dEdn.dat')

    # apply the aperture correction with respect to
    # the slit width and the seeing
    print "apertureCorrection"
    ostar.apertureCorrection(seeing, slitwidth)
    if doPlot : ostar.plotWL("spectrum after apertureCorrection")
    

    # convolution of the spectrum to the resolution of the spectrograph
    sigmaspectrograph = 650./400.
    sigmaspec = 650./inputres
    sigmaconv = (sigmaspectrograph**2 - sigmaspec**2)**.5
    if(sigmaconv<0):
        print 'error: convolution is impossible, the resolution of the input spectrum is less than the resolution of the output spectrum'
        exit(2)
    print "convolvePhoton"
    ostar.convolvePhotondEdl(sigmaconv)
    if doPlot : ostar.plotWL("dEdl after convolvePhoton", dEdl=ostar.dEdlConv)

    
    # conversion to photons
    print "computePhoton"
    ostar.computePhoton()
    if doPlot : ostar.plotNbPhotons("#photons after computePhoton", ostar.phot)


    # for checking
    #ostar.write('phot', starname+'phot.dat')


    # for checking
    #ostar.write('photconv', starname+'photconv.dat')

    # rebin on the ccd grid
    print "computePhotonCCD"
    ostar.computePhotonCCD(gri)
    if doPlot : ostar.plotNbPhotons("#photon after computePhotonCCD", ostar.photccd)
    
    
    # compute and add photon noise
    print "computePhotonNoise"   
    ostar.computePhotonNoise()
    if doPlot : ostar.plotNbPhotons("photon noise computePhotonNoise", ostar.photccdnoise)
    
    
    print "addPhotonNoise"
    ostar.addPhotonNoise()
    if doPlot : ostar.plotNbPhotons("#photon after addPhotonNoise", ostar.photccd)
    
    # for checking
    #ostar.write('photccd', starname+'photccd.dat')

    # through the atmosphere
    print "computePhotonATM"
    ostar.computePhotonATM(atms.gettrans())
    if doPlot : ostar.plotNbPhotons("photatm computePhotonATM", ostar.photatm)
    
    # for checking
    #ostar.write('photatm', starname+'photatm.dat')

    # reflexion on the mirror
    print "computePhotonMirror"
    ostar.computePhotonMirror(mirrortr)
    if doPlot : 
        ostar.plotNbPhotons("mirrortr by computePhotonMirror", ostar.mirrortr)        
        ostar.plotNbPhotons("photmirorr by computePhotonMirror", ostar.photmirror)
    
    # for checking
    #ostar.write('photmirror', starname+'photmirror.dat')

    # through the Grism
    print "computePhotonGrism"
    ostar.computePhotonGrism(gri)
    if doPlot : 
        ostar.plotNbPhotons("grismTrans computePhotonGrism", ostar.grismTrans)
        ostar.plotNbPhotons("photgrism computePhotonGcalibrateSpecrism", ostar.photgrism)

    # conversion in photo-electrons
    print "computeElectronCCD"
    ostar.computeElectronCCD(ccdqe)
    if doPlot : ostar.plotNbPhotons("elecccd computeElectronCCD", ostar.elecccd)

    # add a typical read noise
    print "addNoiseElectronCCD"
    ostar.addNoiseElectronCCD()
    if doPlot : ostar.plotToWLnm( ostar.nuccd, ostar.elecccd, "Final estimation photon-electon by pixel")
    
    
    # for checking
    #ostar.write('elecccdlambda', starname+'elecccdlambda.dat')
    #ostar.write('elecccd', starname+'elecccdnu.dat')

    # Make the sensitivity function to simulate the flux calibration
    print "MakeSensFunc"
    #ostar.MakeSensFunc(atmc.gettrans(), ccdqe, gri, mirrortr)
    #if True : ostar.plotWL("sensdEdl MakeSensFunc", ostar.wlccd[50:-50], ostar.sensdEdl[50:-50])
    ostar.MakeSensFuncInstruOnly(ccdqe, gri, mirrortr)


    # flux-calibrate the spectrum 
    print "calibrateSensSpec"
    ostar.calibrateSensSpec()
    if doPlot : ostar.plotWL("caldEdl calibrateSensSpec", ostar.wlccd, ostar.caldEdl)
    
    #ostar.write('caldEdn', starname+'caldEdn.dat')
    #ostar.write('caldEdl', starname+'caldEdl.dat')

    # add systematic effect due to the flux calibration
    print "AddCalibSys"
    ostar.AddCalibSys(calibsys)
    if doPlot : ostar.plotWL("syscaldEdl AddCalibSys Final", ostar.wlccd, ostar.syscaldEdl)
    




def finalPlot(star, atms, pCCD, pMir, pGrism):
    assert isinstance(star, star.starspectrum)
    assert isinstance(atms, atm.ATMCombSpec)
    assert isinstance(pCCD, ccd )
    assert isinstance(pMir, mir.mirror )    
    assert isinstance(pGrism, grm.grism )    
    pl.figure()
    pl.xlabel("wavelength nm")
    #pl.ylabel(r'$J.m^{-2}.s^{-1}.nm^{-1}$')
    pl.grid()    
    pl.title("Star flux above atmosphere and observed")
    pl.plot(star.wl, star.dEdl)
    pl.plot(star.wlccd, star.syscaldEdl)
    pl.plot(star.wlccd, star.elecccd*cst.hplanck*star.nuccd**3/(cst.clight*240*star._deltaNuCCD))
    pl.legend(["above atm ","observed calibrated","observed raw"])
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
    for i in range(len(star.nuccd)):
        re.append(pGrism.getT(cst.clight/star.nuccd[i]))   
    re = np.array(re)
    pl.plot(cst.clight/star.nuccd, re*100)
    pl.legend(["atmosphere","CCD","telescope mirror" ,"Grism"], loc=7)
    # index refraction
#    pl.figure()
#    pl.xlabel("wavelength nm")
#    pl.grid()    
#    pl.title("Indices of refraction grism BK7")
#    aLambda = arange(200,1200,10) # nm
#    indice = [pGrism.getN(Lambda) for Lambda in aLambda]
#    pl.plot(aLambda, indice)
    
    # interpolation nuCCD
    tck = sp.interpolate.splrep(star.wl[::-1], star.dEdl[::-1])
    dEdlCCD = sp.interpolate.splev(star.wlccd, tck)
    # Ratio EC vrai EC estimé
    if False:
        pl.figure()
        pl.title("Ratio estimate flux / true flux")
        pl.xlabel("wavelength nm")
        pl.grid()
        # interpolation nuCCD
        ratio =  star.syscaldEdl / dEdlCCD 
        print "Ratio:"
        print "min , max : ",star.wl[50], star.wl[-50]
        print ratio[50:-50].mean(), ratio[50:-50].std()
        pl.plot(star.wlccd[10:-10], 1/ratio[10:-10])
    # error relative
    pl.figure()
    pl.title("Relative error true star flux and estimated")
    pl.xlabel("wavelength nm")
    pl.ylabel("%")
    pl.grid()
    # interpolation nuCCD
    #errRel = (dEdlCCD-star.syscaldEdl)/star.syscaldEdl    
    errRel = (dEdlCCD-star.syscaldEdl)/dEdlCCD    
    pl.plot(star.wlccd[10:-10], errRel[10:-10]*100)
    #
    star.plotToWLnm( star.nuccd, star.elecccd, "Final estimation photon-electon by pixel")    
         
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
    ccdqe = ccd()
    gri = grm.grism()
    mirrortr = mir.mirror()
    calibsys = cal.CalibSys()
    atms = atm.ATMCombSpec()
    atmc = atm.ATMCombSpec()
    star = star.starspectrum()

    mirrortr.read(datadir+'/'+ os.environ.get('ATSIM_MIRROR'))
    ccdqe.read(datadir+'/'+ os.environ.get('ATSIM_CCDQE'))
    calibsys.read(datadir+'/'+ os.environ.get('ATSIM_SYSCORR'))

    # read the atmspheres
    atms.read(atmstar)
    atmc.read(atmcalib)
    # read the star spectrum given in J/m2/s/nm
    star.readdEdl(starname)
    if doPlot : star.plotWL("raw spectrum")
    
    # modify the star file name to identify the atmosphere
    starname += '_'+str(atmstar)+'_'
    # for checking
    star.write('dEdn', starname+'dEdn.dat')

    # apply the aperture correction with respect to
    # the slit width and the seeing
    print "apertureCorrection"
    star.apertureCorrection(seeing, slitwidth)
    if doPlot : star.plotWL("spectrum after apertureCorrection")
    

    # convolution of the spectrum to the resolution of the spectrograph
    sigmaspectrograph = 650./400.
    sigmaspec = 650./inputres
    sigmaconv = (sigmaspectrograph**2 - sigmaspec**2)**.5
    if(sigmaconv<0):
        print 'error: convolution is impossible, the resolution of the input spectrum is less than the resolution of the output spectrum'
        exit(2)
    print "convolvePhoton"
    star.convolvePhotondEdl(sigmaconv)
    if doPlot : star.plotWL("dEdl after convolvePhoton", dEdl=star.dEdlConv)

    
    # conversion to photons
    print "computePhoton"
    star.computePhoton()
    if doPlot : star.plotNbPhotons("#photons after computePhoton", star.phot)


    # for checking
    #star.write('phot', starname+'phot.dat')


    # for checking
    #star.write('photconv', starname+'photconv.dat')

    # rebin on the ccd grid
    print "computePhotonCCD"
    star.computePhotonCCD(gri)
    if doPlot : star.plotNbPhotons("#photon after computePhotonCCD", star.photccd)
    
    
    # compute and add photon noise
    print "computePhotonNoise"   
    star.computePhotonNoise()
    if doPlot : star.plotNbPhotons("photon noise computePhotonNoise", star.photccdnoise)
    
    
    print "addPhotonNoise"
    star.addPhotonNoise()
    if doPlot : star.plotNbPhotons("#photon after addPhotonNoise", star.photccd)
    
    # for checking
    #star.write('photccd', starname+'photccd.dat')

    # through the atmosphere
    print "computePhotonATM"
    star.computePhotonATM(atms.gettrans())
    if doPlot : star.plotNbPhotons("photatm computePhotonATM", star.photatm)
    
    # for checking
    #star.write('photatm', starname+'photatm.dat')

    # reflexion on the mirror
    print "computePhotonMirror"
    star.computePhotonMirror(mirrortr)
    if doPlot : 
        star.plotNbPhotons("mirrortr by computePhotonMirror", star.mirrortr)        
        star.plotNbPhotons("photmirorr by computePhotonMirror", star.photmirror)
    
    # for checking
    #star.write('photmirror', starname+'photmirror.dat')

    # through the Grism
    print "computePhotonGrism"
    star.computePhotonGrism(gri)
    if doPlot : 
        star.plotNbPhotons("grismTrans computePhotonGrism", star.grismTrans)
        star.plotNbPhotons("photgrism computePhotonGrism", star.photgrism)

    # conversion in photo-electrons
    print "computeElectronCCD"
    star.computeElectronCCD(ccdqe)
    if doPlot : star.plotNbPhotons("elecccd computeElectronCCD", star.elecccd)

    # add a typical read noise
    print "addNoiseElectronCCD"
    star.addNoiseElectronCCD()
    if doPlot : star.plotToWLnm( star.nuccd, star.elecccd, "Final estimation photon-electon by pixel")
    
    
    # for checking
    #star.write('elecccdlambda', starname+'elecccdlambda.dat')
    #star.write('elecccd', starname+'elecccdnu.dat')

    # Make the sensitivity function to simulate the flux calibration
    print "MakeSensFunc"
    #star.MakeSensFunc(atmc.gettrans(), ccdqe, gri, mirrortr)
    #if True : star.plotWL("sensdEdl MakeSensFunc", star.wlccd[50:-50], star.sensdEdl[50:-50])
    star.MakeSensFuncJMC(atmc.gettrans(), ccdqe, gri, mirrortr)
    
    
    # modify the star file name to identify the calibration atmosphere
    starname += '_'+str(atmcalib)+'_'
    # for checking
    #star.write('atmrebin', starname+'atmrebin.dat')
    #star.write('sensdEdl', starname+'sensdEdl.dat')
    #star.write('sensdEdn', starname+'sensdEdn.dat')

    # flux-calibrate the spectrum 
    print "calibrateSensSpec"
    star.calibrateSensSpec()
    if doPlot : star.plotWL("caldEdl calibrateSensSpec", star.wlccd, star.caldEdl)
    
    #star.write('caldEdn', starname+'caldEdn.dat')
    #star.write('caldEdl', starname+'caldEdl.dat')

    # add systematic effect due to the flux calibration
    print "AddCalibSys"
    star.AddCalibSys(calibsys)
    if doPlot : star.plotWL("syscaldEdl AddCalibSys", star.wlccd, star.syscaldEdl)
    
    # output of the simulation
    star.write('syscaldEdn', starname+'syscaldEdn.dat')
    star.write('syscaldEdl', starname+'syscaldEdl.dat')

    # 
    #star.EstimTransmissionATM(atms.gettrans(), ccdqe, gri, mirrortr, )

    finalPlot(star,atms, ccdqe, mirrortr, gri)

if __name__ == "__main__":
    #pl.ioff()
    #pl.matplotlib.use('GTKAgg')
    main()

try:
    pl.show()
except AttributeError:
    pass
