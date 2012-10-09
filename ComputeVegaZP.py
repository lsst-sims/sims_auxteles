#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

from filter import *
from scipy.integrate import *

v = filter()
v.readl('/Users/gurvan/Models/Filters/Landolt/sv_-27A.dat', 1.)

# Vega spectrum in erg/cm2/s/A
vega_l = []
vega_f = []
file = open('/Users/gurvan/Models/alpha_lyr_stis_004.ascii', 'r')
for line in file.readlines():
    line = line.split()
    vega_l.append(float(line[0]))
    vega_f.append(float(line[1]))
file.close()

vega_l = array(vega_l)
vega_f = array(vega_f)

print 2.5*log10(trapz(vega_f*v.getTl(vega_l), vega_l))
