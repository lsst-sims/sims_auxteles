# -*- coding: utf-8 -*-

#############################################################################
#           grism - a class defining the properties of the grism            #
#                                                                           #
#        author:          Gurvan BAZIN                                      #
#                         Universitäts-Sternwarte-München                   #
#        latest revision: June 2010                                         #
#############################################################################

# include some useful libraries
# and class definitions
import math
from constants import *
import numpy as np



class grism:
    
    def __init__(self):
        # define some useful quantities
        # Sellmeier coefficients
        self.B = np.array([1.03961212, 2.31792344e-1, 1.01046945])
        self.C = np.array([6.00069867e3, 2.00179144e4, 1.03560653e8])
        # grating characteristic
        self.a = 1./100.e3 # m
        # prism angle
        self.prismAngle = math.pi/180. # 1 degree

    def getN(self, l):
        # return the indice of the prism as a function of the wavelength
        return np.sqrt(1. + sum( [ self.B[i]*l*l/(l*l-self.C[i]) for i in range(3) ] ))

    def getT(self, l):
        # return the transmission of the prism
        # normal indicence for the input beam
        kair = 2.*np.pi/l
        kverre = 2.*np.pi*self.getN(l)/l
        t = 2.*kair/(kair+kverre)
        # output beam
        theta1 = self.gratingAngle(l) + self.prismAngle
        theta2 = math.asin(math.sin(theta1)*self.getN(l))
        t *= 2.*self.getN(l)*math.cos(theta1)/(self.getN(l)*math.cos(theta1) + math.cos(theta2))
        return t

    def gratingAngle(self, l):
        """
        l [nm] wavelength
        
        Compute angle diffraction with grating relation for order 0
        """
        # return the angle at the ouput of the grating
        #print l, l*1.e-9, l*1.e-9/self.a, self.a
        return math.asin(l*1.e-9/self.a)

    def outputAngle(self, l):
        """
        l [nm] wavelength
        
        Compute deviation angle of grism
        """
        # return the total angle at the ouput of the grism
        return math.asin(self.getN(l) * math.sin(self.gratingAngle(l) + self.prismAngle)) - self.prismAngle
