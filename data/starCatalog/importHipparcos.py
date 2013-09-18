#!/usr/bin/env python

#
# JM Colley , University Denis Diderot, laboratory APC
# Sept 2013
#

import os
import sys
import getopt


sys.path.append('../../libauxteles')


import hipparcos as hip

def usage():
    print 'NAME'    
    print '  importHipparcos - import selected hipparcos catalog in python format pickle '
    print ''    
    print 'SYNTAX '
    print '  importHipparcos <path to hip_main.dat>'
    print ' '
    print 'DESCRIPTION '
    print '  create file pickle hipparcos.pic contents selected hipparcos catalog'
    print ''
    
MyArgs = sys.argv
NbArgs = len(MyArgs)
#print MyArgs

# aucun argument
if NbArgs == 1:
    usage()
    sys.exit(2)

obj = hip.Hipparcos(MyArgs[1])
oCat= StarCatalog(obj.nStar)
obj.extractFieldForAuxteles(oCat)
obj.close()
oCat.convertSpectralToKurucz()
oCat.removeStarNOK()    
mfile='hipparcos.pic'
oCat.save(mfile)

