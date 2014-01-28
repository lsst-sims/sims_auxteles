#!/usr/bin/env python

#
# JM Colley , University Denis Diderot, laboratory APC
# Sept 2013
#

import getopt
import os
import sys

import hipparcos as hip


sys.path.append('../../libauxteles')



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
oCat=  hip.StarCatalog(obj.nStar)
obj.extractFieldForAuxteles(oCat)
obj.close()
#oCat.convertSpectralToKurucz()
oCat.removeNotUsed()
oCat.printCat(0, 20) 
mfile='hipparcos.pic'
oCat.save(mfile)

