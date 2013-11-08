'''
Created on 3 juin 2013

@author: colley
'''
import unittest
from starTargetSimu import *
import kurucz as kz
import obsStrategy as obS


S_FileKuruczPic = '../../data/kurucz93/k93.pic'
S_StarPos = "../../data/simuRef/starPos2013-12-21.pkl"



def test_plot():
    oKur = kz.Kurucz(S_FileKuruczPic)    
    oKur.restrictToWLinterval(3000, 9000)
    cat = StarTargetSimuAll(oKur)
    cat.plotAll()
    
    
    
def test_SetAllStars():
    oKur = kz.Kurucz(S_FileKuruczPic)    
    oKur.restrictToWLinterval(3000, 9000)
    cat = StarTargetSimuAll(oKur)
    tmp = obS.DataPositionStarTarget()
    d = tmp.read(S_StarPos)
    print d
    cat.setAllStars(d.Kur[:,0], d.Kur[:,2], d.Kur[:,1], d.Mag, d.Par, d.IdCst)
    cat.plotAll()
    print cat.IDstar2Idx(32778)
    print cat.IDstar2Idx(60078)
    print cat.IDstar2Idx(3277)
  
   

#test_plot()
test_SetAllStars()

try:
    pl.show()
except AttributeError:
    pass

    
