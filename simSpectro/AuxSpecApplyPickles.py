#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ATMCombSpec import *
from ccd import *
from filter import *
from grism import *
from mirror import *
from starspectrum import *


liststars = []
liststarsfile = open('/groups/LSST/bazin/Catalogs/Pickles/list.Pickles')
for line in liststarsfile.readlines():
    line = line.split()
    liststars.append(line[0])
liststarsfile.close()


#nbatm = 16
nbatm = 1 # for test

for atmid in range(nbatm):
    atm = ATMCombSpec()
    atmname = 'Modtran_'+str(atmid)+'.dat'
    print 'Lecture de ', atmname
    atm.read('/groups/LSST/bazin/T05.00/'+atmname)

    for starname in liststars:
        print 'Lecture de ', starname
        star = starspectrum()
        star.readdEdlPickles('/groups/LSST/bazin/Catalogs/Pickles/'+starname+'.dat', 0.1)
        star.extrapol()
        starname += '_'+str(atmid)+'_'

        print 'Passage de l\'atmosphère'
        star.computeEnergyATM(atm.gettrans())
        star.write('atmdEdl', starname+'atmdEdl.dat')
