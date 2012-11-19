'''
Created on 13 nov. 2012

@author: colley
'''

import pickle as pk
import numpy as np
import math
import pyfits as pf
import pylab as pl
import scipy.interpolate as sci



def read_kurucz_all(directory = '/home/colley/projet/lsst/stellar_spectra/k93models/'):
    metal_dir = ['km50', 'km45', 'km40', 'km35', 'km30', 'km25', 'km20', 'km15', 'km10', 'km05', 'km03', 'km02', 'km01', 'kp00', 'kp01', 'kp02', 'kp03', 'kp05', 'kp10']
    step_metal = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, +0.0, +0.1, +0.2, +0.3, +0.5, +1.0]
    step_temperature = [3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 37500, 40000, 42500, 45000, 47500, 50000]
    step_gravity = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

    nb_metal = 19
    nb_temp = 61
    nb_logg = 11
    nb_lambda = 1221

    nb_col = 12684
    nb_lig = 1224

    data = np.zeros((nb_lig, nb_col))
    data[0,0] = -99
    data[1,0] = -99
    data[2,0] = -99

    indice_col = 1
    for m in range(0, nb_metal):
        # cas m = +0.5
        if (m == nb_metal-2):
            nb_temp = 59
        # cas m = +1.0
        if (m == nb_metal-1):
            nb_temp = 57
        for t in range(0, nb_temp):
            nom_file_kurucz = directory + "%s/%s_%s.fits" % (metal_dir[m], metal_dir[m], step_temperature[t])
            for g in range(0, nb_logg):
                loggravity_name = "g%02i" %(int(math.floor(step_gravity[g]*10.)))

                data_k93 = pf.getdata(nom_file_kurucz)
                print "read ",nom_file_kurucz
                wavelength = data_k93.field('wavelength')
                flux = data_k93.field(loggravity_name)

                data[3:1224,0] = wavelength # * 10. # To get wavelength in A
                data[0,indice_col] = step_metal[m]
                data[1,indice_col] = step_temperature[t]
                data[2,indice_col] = step_gravity[g]
                data[3:1224,indice_col] = flux
                
                indice_col += 1

                #print m, t, g, data[0:5,0:5]
    return data



def fits2pickle(fIn, fOut):
    kur = read_kurucz_all(fIn)
    f=open(fOut, "wb")
    pk.dump(kur, f, pk.HIGHEST_PROTOCOL)
    f.close()



def fits2pickle_JMC():
    fIn = '/home/colley/projet/lsst/stellar_spectra/k93models/'
    fOut = '/home/colley/projet/lsst/stellar_spectra/k93.pic'
    fits2pickle(fIn, fOut)
    
    
    
class Kurucz(object):
    '''    
    Format Kurucz array :
        column 0 : wavelength
        line 0 : metallicity
        line 1 : temperature
        line 2 : gravity
    '''

    def __init__(self, filePickle):
        try:
            f=open(filePickle, "rb")
        except:
            print "can't open file ", filePickle
            return 
        self._Flux = pk.load(f)
        #print self._Flux
        f.close()
    
    
    def plotFlux(self,idx):
        pl.figure()
        pl.plot(self._Flux[3:700,0],self._Flux[3:700,idx])
        pl.xlabel("Angstrom")
        pl.grid()
        pl.title("Kurucz 93 star flux M %.2f T %.2f G %.2f"%(self._Flux[0,idx], self._Flux[1,idx], self._Flux[2,idx]))
        
        
    def plotMultiFlux(self,idx0, nb):
        pl.figure()
        for idx in range(nb):
            pl.plot(self._Flux[3:700,0],self._Flux[3:700,idx+idx0])
        strLgd = []
        for idxR in range(nb):
            idx = idxR + idx0
            strLgd.append("M %.2f T %.2f G %.2f"%(self._Flux[0,idx], self._Flux[1,idx], self._Flux[2,idx]))
        pl.legend(strLgd)
        pl.xlabel("Angstrom")
        pl.grid()
        pl.title("Kurucz 93 star flux")
    
    
    def getIdxNearestFlux(self, pMet, pTemp, pGra):
        a = self._Flux[0:3,1:].transpose()
        oNGP = sci.NearestNDInterpolator(a, np.arange(1, len(self._Flux[0]-1)))
        b = np.array([[pMet, pTemp, pGra]])
        idx = oNGP(b)
        print idx, self._Flux[0:3,idx]
        return idx