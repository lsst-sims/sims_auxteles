'''
Created on 3 juin 2013

@author: colley
'''
import unittest
from starTargetSimu import *
import kurucz as kz

FileKuruczPic = '../../data/kurucz93/k93.pic'


def test_plot():
    oKur = kz.Kurucz(FileKuruczPic)    
    oKur.restrictToWLinterval(3000, 9000)
    cat = StarTargetSimuAll(oKur)
    cat.plotAll()
    pl.show()
    
    

test_plot()


