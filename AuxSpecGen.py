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
from starspectrum import *
from ATMCombSpec import *
from mirror import *
from grism import *
from ccd import *
from CalibSys import *
import pylab as pl
import scipy as sp 

def finalPlot(star, atms, pCCD, pMir, pGrism):
    assert isinstance(star, starspectrum)
    assert isinstance(atms, ATMCombSpec)
    assert isinstance(pCCD, ccd )
    assert isinstance(pMir, mirror )    
    assert isinstance(pGrism, grism )    
    pl.figure()
    pl.xlabel("wavelength nm")
    pl.ylabel("$J.m^{-2}.s^{-1}.nm^{-1}$")
    pl.grid()    
    pl.title("Star flux above atmosphere and estimated")
    pl.plot(star.wl, star.dEdl)
    pl.plot(star.wlccd, star.syscaldEdl)
    pl.legend(["above atm ","estimated"])
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
        re.append(pGrism.getT(clight/star.nuccd[i]))   
    re = array(re)
    pl.plot(clight/star.nuccd, re*100)
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
    errRel = (dEdlCCD-star.syscaldEdl)/star.syscaldEdl    
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
    gri = grism()
    mirrortr = mirror()
    calibsys = CalibSys()
    atms = ATMCombSpec()
    atmc = ATMCombSpec()
    star = starspectrum()

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
    star.EstimTransmissionATM(atms.gettrans(), ccdqe, gri, mirrortr, )

    finalPlot(star,atms, ccdqe, mirrortr, gri)

if __name__ == "__main__":
    #pl.ioff()
    #pl.matplotlib.use('GTKAgg')
    main()

try:
    pl.show()
except AttributeError:
    pass
