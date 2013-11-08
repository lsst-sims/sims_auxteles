'''
Created on 13 sept. 2013

@author: colley
'''
from hipparcos import *


S_PathHip='/home/colley/temp/lsst/hipparcos'
S_pathCat = "../../data/starCatalog/hipparcos.pic"


#
# CLASS Hipparcos
#

def test_Line():
    obj = Hipparcos(S_PathHip)
    obj.read()
    obj.close()
    
    
def test_read():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractSingleStarForAuxteles(oCat)
    oCat.printCat(10000, 12000)
    obj.close()
    
def test_readKurucz():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractSingleStarForAuxteles(oCat)
    oCat.printCat(10000, 10200)
    oCat.convertSpectralToKurucz()
    oCat.removeStarKuruczNOK()
    obj.close()
    
    
def test_readKuruczRemoveNOK():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractSingleStarForAuxteles(oCat)
    oCat.printCat(0, 20)
    oCat.convertSpectralToKurucz()
    oCat.removeStarNOK()
    oCat.printCat(0, 50)
    obj.close()
    
    
def test_saveRead():
    obj = Hipparcos(S_PathHip)
    oCat= StarCatalog(obj.nStar)
    obj.extractSingleStarForAuxteles(oCat)
    obj.close()
    oCat.convertSpectralToKurucz()
    oCat.removeStarNOK()
    oCat.printCat(50, 60)
    mfile='/home/colley/temp/toto.pic'
    oCat.save(mfile)
    obj = oCat.read(mfile)
    print "restore object"
    obj.printCat(50, 60)
    
#
# TEST method 
#

#test_Line()
#test_read()
test_readKurucz()
#test_readKuruczRemoveNOK()
#test_saveRead()


#
# CLASS StarCatalog
#

def testStarCatalog_findTypeStar01():
    oTemp = StarCatalog(1)
    oCat = StarCatalog(1)
    oCat = oTemp.read(S_pathCat)
    del oTemp
    #oCat.printCat(1000, 1500)
    idx =oCat.findTypeStar("G",pLum='V')
    print "naine jaune ", len(idx)
    oCat.printCat(0,0, idx[:20])
    
    
def testStarCatalog_RedDwarf():
    oTemp = StarCatalog(1)
    oCat = StarCatalog(1)
    oCat = oTemp.read(S_pathCat)
    del oTemp
    #oCat.printCat(1000, 1500)
    idx =oCat.findTypeStar("M",pLum='V')
    print "naine rouge ", len(idx)
    oCat.printCat(0,0, idx[:20])
    
    
def testStarCatalog_OrangeDwarf():
    oTemp = StarCatalog(1)
    oCat = StarCatalog(1)
    oCat = oTemp.read(S_pathCat)
    del oTemp
    #oCat.printCat(1000, 1500)
    idx =oCat.findTypeStar("K",pLum='V')
    print "naine orange ", len(idx)
    oCat.printCat(0,0, idx[:20])
       
#
# TEST method 
#
#testStarCatalog_findTypeStar01()
#testStarCatalog_RedDwarf()
#testStarCatalog_OrangeDwarf()