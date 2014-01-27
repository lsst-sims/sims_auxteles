# -*- coding: utf-8 -*-

#############################################################################
#       CalibSys - a class defining the systematic effects of the flux      #
#                  calibration of the star spectrum                         #
#        author:          Gurvan BAZIN                                      #
#                         Universitäts-Sternwarte-München                   #
#        latest revision: July 2010                                         #
#############################################################################

# include some useful libraries
# and class definitions
import math
import numpy
import numpy.fft
import numpy.random
import scipy.interpolate

import pylab as pl


S_LevelPlot = 0


class CalibSys:

    def __init__(self):
        # distance and correlation function
        self.dist = []
        self.corr = []
        self.step = 0.

    def read(self, filename):
        # read the file
        file = open(filename, 'r')
        for line in file.readlines():
            line = line.split()
            if line[0]=="#":
                continue
            self.dist.append(float(line[0]))
            self.corr.append(float(line[1]))
        file.close()
        self.dist = numpy.array(self.dist)
        self.corr = numpy.array(self.corr)
        self.step = self.dist[1] - self.dist[0]
        #self.plotCorrel()
        
        
    def plotCorrel(self):
        pl.figure()
        pl.xlabel("nm")
        pl.ylabel("?")
        pl.grid()
        pl.title("correlation function")
        pl.plot(self.dist, self.corr)
        
    def plotNoise(self,x,y,title=''):
        pl.figure()
        pl.xlabel("nm")
        pl.ylabel("?")
        pl.grid()
        pl.title(title)
        pl.plot(x,y)
      
        
    def GetRandomRealization(self, wl, s, sn):
        distmax = max(wl) - min(wl)
        c = []
        d = []
        dtmp = 0.
        i = 0
        # set a correlation function matching the spectrum
        while(dtmp <= distmax):
            if i < len(self.dist):
                c.append(self.corr[i])
                d.append(self.dist[i])
                i += 1
            else:
                c.append(0.)
                d.append(d[-1] + self.step)
            dtmp = d[-1]
        c = numpy.array(c)
        # for symetry reason...
        c = numpy.concatenate([c, c[::-1][0:len(c)-1]])
        # Generate the random realization of the deviation
        PowerSpec = abs(numpy.fft.fft(c))
        noise = numpy.random.normal(numpy.zeros(len(c)))
        fourierdeviation = numpy.fft.fft(noise)*numpy.sqrt(PowerSpec)
        deviation = numpy.real(numpy.fft.ifft(fourierdeviation))
        
        # interpolate the deviation on the wl of the spectrum
        tck = scipy.interpolate.splrep(d+min(wl), deviation[0:len(d)])
        intdev = scipy.interpolate.splev(wl, tck, der=0)
        if S_LevelPlot > 3: self.plotNoise(wl, intdev, "raw noise")
        # return the deviation times the rms
        # given that the rms is at a good approximation propto the inverse of the s/n
        rms = []
        for e in sn:
            if e != 0.:
                rms.append(1./(2.*e))
            else:
                rms.append(0.)
        rms = numpy.array(rms)
        ret = intdev*rms*s
        if S_LevelPlot > 3: self.plotNoise(wl, ret, "noise")
        return ret


