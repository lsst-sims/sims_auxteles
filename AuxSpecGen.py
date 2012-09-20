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

    # modify the star file name to identify the atmosphere
    starname += '_'+str(atmstar)+'_'
    # for checking
    star.write('dEdn', starname+'dEdn.dat')

    # apply the aperture correction with respect to
    # the slit width and the seeing
    star.apertureCorrection(seeing, slitwidth)

    # conversion to photons
    star.computePhoton()
    # for checking
    #star.write('phot', starname+'phot.dat')
    
    # convolution of the spectrum to the resolution of the spectrograph
    sigmaspectrograph = 650./400.
    sigmaspec = 650./inputres
    sigmaconv = (sigmaspectrograph**2 - sigmaspec**2)**.5
    if(sigmaconv<0):
        print 'error: convolution is impossible, the resolution of the input spectrum is less than the resolution of the output spectrum'
        exit(2)
    star.convolvePhoton(sigmaconv)
    # for checking
    #star.write('photconv', starname+'photconv.dat')

    # rebin on the ccd grid
    star.computePhotonCCD(gri)

    # compute and add photon noise
    star.computePhotonNoise()
    star.addPhotonNoise()
    # for checking
    #star.write('photccd', starname+'photccd.dat')

    # through the atmosphere
    star.computePhotonATM(atms.gettrans())
    # for checking
    #star.write('photatm', starname+'photatm.dat')

    # reflexion on the mirror
    star.computePhotonMirror(mirrortr)
    # for checking
    #star.write('photmirror', starname+'photmirror.dat')

    # through the Grism
    star.computePhotonGrism(gri)

    # conversion in photo-electrons
    star.computeElectronCCD(ccdqe)

    # add a typical read noise
    star.addNoiseElectronCCD()
    # for checking
    #star.write('elecccdlambda', starname+'elecccdlambda.dat')
    #star.write('elecccd', starname+'elecccdnu.dat')

    # Make the sensitivity function to simulate the flux calibration
    star.MakeSensFunc(atmc.gettrans(), ccdqe, gri, mirrortr)
    # modify the star file name to identify the calibration atmosphere
    starname += '_'+str(atmcalib)+'_'
    # for checking
    #star.write('atmrebin', starname+'atmrebin.dat')
    #star.write('sensdEdl', starname+'sensdEdl.dat')
    #star.write('sensdEdn', starname+'sensdEdn.dat')

    # flux-calibrate the spectrum 
    star.calibrateSensSpec()
    #star.write('caldEdn', starname+'caldEdn.dat')
    #star.write('caldEdl', starname+'caldEdl.dat')

    # add systematic effect due to the flux calibration
    star.AddCalibSys(calibsys)
    # output of the simulation
    star.write('syscaldEdn', starname+'syscaldEdn.dat')
    star.write('syscaldEdl', starname+'syscaldEdl.dat')


if __name__ == "__main__":
    main()

