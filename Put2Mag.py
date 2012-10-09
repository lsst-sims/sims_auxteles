#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

from filter import *
from scipy.integrate import *
import sys

if len(sys.argv) != 4:
    print 'usage: Put2Mag.py file.in mag file.out'
    exit(2)

v = filter()
v.readl('/Users/gurvan/Models/Filters/Landolt/sv_-27A.dat', 1.)

# Leborgne spectrum in erg/cm2/s/A
hd_l = []
hd_f = []
file = open(sys.argv[1], 'r')
for line in file.readlines():
    line = line.split()
    hd_l.append(float(line[0]))
    hd_f.append(float(line[1]))
file.close()

hd_l = array(hd_l)
hd_f = array(hd_f)

#-2.5log(f*alpha)+zp = mag
#alpha = 10**((zp-mag)/(2.5))/f

flux = trapz(hd_f*v.getTl(hd_l), hd_l)
alpha = 10**((-13.7458349352-float(sys.argv[2]))/2.5)/flux

fileout = open(sys.argv[3], 'w')
for i in range(len(hd_l)):
    fileout.write(str(hd_l[i])+'\t'+str(hd_f[i]*alpha)+'\n')
fileout.close()

