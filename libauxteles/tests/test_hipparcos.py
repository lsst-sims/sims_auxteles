'''
Created on 13 sept. 2013

@author: colley
'''
from hipparcos import *


S_PathHip='/home/colley/temp/lsst/hipparcos'


def test_Line():
    obj = Hipparcos(S_PathHip)
    obj.read()
    obj.close()
    
    
def test_read():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractFieldForAuxteles(oCat)
    oCat.printCat(10000, 10200)
    obj.close()
    
def test_readKurucz():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractFieldForAuxteles(oCat)
    oCat.printCat(10000, 10200)
    oCat.convertSpectralToKurucz()
    obj.close()
    
    
def test_readKuruczRemoveNOK():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractFieldForAuxteles(oCat)
    oCat.printCat(0, 20)
    oCat.convertSpectralToKurucz()
    oCat.removeStarNOK()
    oCat.printCat(0, 50)
    obj.close()
    
    
def test_saveRead():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractFieldForAuxteles(oCat)
    obj.close()
    oCat.convertSpectralToKurucz()
    oCat.removeStarNOK()
    oCat.printCat(50, 60)
    mfile='/home/colley/temp/toto.pic'
    oCat.save(mfile)
    obj = oCat.read(mfile)
    print "restore object"
    obj.printCat(50, 60)
    
    
   
#test_Line()

#test_read()
#test_readKurucz()
#test_readKuruczRemoveNOK()
test_saveRead()