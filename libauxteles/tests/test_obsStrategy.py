'''
Created on 20 sept. 2013

@author: colley
'''

from obsStrategy import *


S_pathCat = "../../data/starCatalog/hipparcos.pic"


def test_init():
    oObs = ObsStrategyReal01(S_pathCat)
    oObs.oCat.printCat(0, 50)
    
    
def test_selectMag():
    oObs = ObsStrategyReal01(S_pathCat)
    print "nb star ", oObs.oCat.nbStar
    oObs.oCat.printCat(0, 20)
    oObs.selectMag(7, 9)
    print "nb star ", oObs.oCat.nbStar
    oObs.oCat.printCat(0, 20)
    
    
def test_computCoordLoc():
    oObs = ObsStrategyReal01(S_pathCat)
    oObs.selectMag(7, 8)
    oObs.computCoordRefObs([2013, 9,20, 12, 30])
    oObs.histoDistZen()
    

def test_plotVisible(): 
    oObs = ObsStrategyReal01(S_pathCat)
    oObs.selectMag(7, 7.015)
    oObs.computCoordRefObs([2013, 9,20, 13, 50])
    oObs.selectVisibleStar(60)
    oObs.plotSelectStarsz()
    
    
def test_multiplot():
    oObs = ObsStrategyReal01(S_pathCat)
    oObs.selectMag(7, 7.2)
    oObs.doMultiplot([2013, 9,20, 13, 00], 3, 10)
    
    
def test_selectAlgo01():
    oObs = ObsStrategyReal01(S_pathCat)
    oObs.selectMag(7, 7.5)
    oObs.selectAlgo01([2013, 9,20, 13, 00], [2013, 9,20, 18, 00])
    
    
def test_doMultiplotSelect():
    oObs = ObsStrategyReal01(S_pathCat)
    oObs.selectMag(7, 7.5)
    oObs.selectAlgo01([2013, 9,20, 13, 00], [2013, 9,20, 18, 00])
    oObs.doMultiplotSelect([2013, 9,20, 13, 00], 3, 3)  
    
    
def test_getRisingSettingSun():
    oObs = ObsStrategyReal01(S_pathCat)
    print oObs.getRisingSettingSun([2013, 10, 8])


def test_selectAlgo02():
    oObs = ObsStrategyReal01(S_pathCat)
    oObs.selectMag(7, 7.2)
    oObs.selectAlgo02([2013, 6, 21])
    oObs.doMultiplotAlgo02()
      
#
# 
#

#test_init()
#test_selectMag()
#test_computCoordLoc()
#test_plotVisible()
#test_multiplot()
#test_selectAlgo01()
#test_doMultiplotSelect()
#test_getRisingSettingSun()
test_selectAlgo02()


#try:
#    pl.show()
#except AttributeError:
#    pass
