'''
Created on 6 nov. 2012

@author: colley
'''

import numpy as np
import math
import pyfits as pf


S_verbose = 1


def removeAllEltlike(myList,elt):
    while(True):
        try:
            myList.remove(elt)
        except:
            return 
        
        
def getFirstCharNotBlanck(myStr):
    for elt in  myStr:
        if elt != ' ' : return elt
        
        
def removeBlanckInString(myStr):
    for idx, elt in  enumerate(myStr):
        if   elt != ' ':
            return myStr[idx:]
    return ''
            
            
def readTextFileColumn(fileName, sep = None):
    """
    read text file with N columns and P lines, so return a real array with dimension PxN
    sep = None : process blank and tabulation
    """
    # test existance du fichier
    listReal = []
    cptLine = 0
    cptLineTot = 0
    fd = open(fileName, 'r')
    for line in fd.readlines():
        cptLineTot +=1
        #print "'%s'"%line 
        if line[0] == "#":
            continue
        if sep ==None:
            listCol = line.split()
        else:
            listCol = line.split(sep)
        if S_verbose >0 : print listCol
        removeAllEltlike(listCol,'')  
        removeAllEltlike(listCol,'\n')  
        print listCol
        if len(listCol) == 0:
            if S_verbose >0 : print "ligne vide ... continue"
            continue        
        if getFirstCharNotBlanck(listCol[0]) == "#":
            if S_verbose >0 : print "ligne commente ... continue"
            continue
        if cptLine == 0:
            removeAllEltlike(listCol,'\n')
            nbCol = len(listCol)
        elif nbCol != len(listCol):
            if len(listCol) >= 1:
                clearStr = removeBlanckInString(listCol[0])
                if clearStr[0] == '\n':
                    if S_verbose >0 : print "ligne vide ... continue"
                    continue                    
            print "ERROR line %d, read %d elements instead %d"%(cptLineTot, len(listCol),nbCol )
            return False
        for Elt in listCol:
            try:
                _RealVal = float(Elt)                    
            except:
                print "ERROR line %d, can't convert '%s' in real number"%(cptLineTot, Elt)
                return False
            listReal.append(_RealVal)
        cptLine += 1
        # to numpy array
    outArray = np.array(listReal)
    outArray = np.reshape(outArray,(cptLine, nbCol))
    return outArray



def read_kurucz_all():
    directory = '/home/colley/projet/lsst/stellar_spectra/k93models/'

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


if __name__ == "__main__":
    a=read_kurucz_all()
    print a