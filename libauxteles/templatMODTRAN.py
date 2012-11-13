'''
Created on 6 nov. 2012

@author: colley
'''

import tools as tl 
import pylab as pl


class TemplateMODTRAN(object):
    '''
    classdocs
    '''


    def __init__(self, fileName):
        '''
        Constructor
        '''
        out = tl.readTextFileColumn(fileName)
        self._wl = out[:,0]
        self._Tr03 = out[:,3]
        self._Trmols  = out[:,5]
        self._Trmola  = out[:,7]
        self._TrH2O  = out[:,2]*out[:,4]
    
    
    def plotTemplate(self):
        pl.figure()
        pl.plot(self._wl, self._Tr03)
        pl.plot(self._wl, self._Trmols)
        pl.plot(self._wl, self._Trmola)
        pl.plot(self._wl, self._TrH2O)
        TrAll = self._Tr03*self._Trmols*self._Trmola* self._TrH2O
        pl.plot(self._wl, TrAll,'y')
        pl.title("Template MODTRAN transmission")
        pl.legend(["03","mols","mola/02","H2O","all"],loc=4)
        pl.grid()
        pl.xlabel("nm")
        pl.ylabel("%")
        pl.ylim([0,1.02])
        
if __name__ == "__main__":
    Tpl = TemplateMODTRAN('/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt')
    Tpl.plotTemplate()
    pl.show()

        