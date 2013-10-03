'''
Created on 20 sept. 2013

@author: colley
'''


import hipparcos as hip
import MJDtools as mjd
import Astrotools as astro
import numpy as np
import pylab as pl
import copy
from numpy.core.numeric import concatenate

S_LevelPlot = 5


class ObsStrategyReal01(object):
    def __init__(self, fileCat):
        oTemp = hip.StarCatalog(0)
        self.oCat =  hip.StarCatalog(0)      
        self.oCat =  oTemp.read(fileCat)
        self.MaxDZ = 60
        
        
    def selectMag(self, pMin , pMax):
        self.oCat.selectMag(pMin, pMax)
        
        
    def computCoordLoc(self, pDate, pSite=None):
        """
        pDate: format [YYYY, MM, DD, HH, MM]
        pSite: format [longitude, latitude]
        """       
        dateJD = mjd.tomjd([pDate[0], pDate[1], pDate[2], (pDate[3] + pDate[4]/60.)/24.])
        self.coodLoc = np.empty_like(self.oCat.coord)
        if pSite == None:
            for idx in range (self.oCat.nbStar):
                self.coodLoc[idx] = astro.eq2loc([self.oCat.coord[idx,0], self.oCat.coord[idx,1],dateJD] )                                   
        else:
            for idx in range (self.oCat.nbStar):
                self.coodLoc[idx] = astro.eq2loc([self.oCat.coord[idx,0], self.oCat.coord[idx,1],dateJD, pSite[0], pSite[1] ] )


    def selectVisibleStar(self, pElv):
        self.lIdxVis = []
        for idx in range (self.oCat.nbStar):
            if self.coodLoc[idx, 1] < pElv: 
                self.lIdxVis.append(idx)


    def selectAlgo01(self, pTbeg, pTend):
        """
        selection 4 stars
        """
        # compute position a pTbeg
        self.computCoordLoc(pTbeg)
        self.selectVisibleStar(self.MaxDZ)        
        # possible stars near horizon    
        print  self.coodLoc[self.lIdxVis]  
        y = np.sin(np.deg2rad(self.coodLoc[self.lIdxVis, 1]))*np.sin(np.deg2rad(self.coodLoc[self.lIdxVis, 0]))        
        lIdxHor1 = np.where(np.logical_and(np.logical_and(self.coodLoc[self.lIdxVis, 0]>-135,self.coodLoc[self.lIdxVis, 0]<-90), self.coodLoc[self.lIdxVis, 1] > 55 ))[0]
        lIdxHor2 = np.where(np.logical_and(np.logical_and(self.coodLoc[self.lIdxVis, 0]>-80,self.coodLoc[self.lIdxVis, 0]<-20), self.coodLoc[self.lIdxVis, 1] > 55 ))[0]
        if  S_LevelPlot > 2:
            print lIdxHor1, lIdxHor2
            idx = np.concatenate((lIdxHor1, lIdxHor2))
            self.plotStarWithIdx(idx, "Star possible")            
        # possible stars near zenith
        lIdxZen3 = np.where(np.logical_and(np.logical_and(np.logical_and(self.coodLoc[self.lIdxVis, 0]>90,self.coodLoc[self.lIdxVis, 0]<180), y < 0.08 ), self.coodLoc[self.lIdxVis, 1] <30))[0]
        lIdxZen4 = np.where(np.logical_and(np.logical_and(np.logical_and(self.coodLoc[self.lIdxVis, 0]>0,self.coodLoc[self.lIdxVis, 0]<90),  y < 0.08), self.coodLoc[self.lIdxVis, 1] <30))[0]
        if  S_LevelPlot > 2:
            print lIdxHor1, lIdxHor2
            idx = np.concatenate(( lIdxZen3, lIdxZen4))
            self.plotStarWithIdx(idx, "Star possible")
        # selection star
        selectStarVis=  [lIdxHor1[0], lIdxHor2[0], lIdxZen3[0], lIdxZen4[0]]   
        print  selectStarVis
        self.selectStar = [self.lIdxVis[idx] for idx in selectStarVis]
        if  S_LevelPlot > 2:                        
            self.plotStarWithIdx(selectStarVis, "Star select")
        
#
# Plot
#         
    def plotStarWithIdx(self, pPidx, pTit, pFile=None):
        pl.figure()
        print pFile       
        ax = pl.subplot(111, polar=True)
        ax.set_title(pTit)
        for idx in range(len(pPidx)):
            idxVis = self.lIdxVis[pPidx[idx]]            
            ax.plot(np.deg2rad(self.coodLoc[idxVis, 0]), np.sin(np.deg2rad(self.coodLoc[idxVis, 1])),'*b')
        ax.set_rmax(1.0)
        ax.grid(True)
        if pFile != None:
            pl.savefig(pFile)
            
           
    def plotVisible(self,pTit, pFile=None):
        pl.figure()
        print pFile       
        ax = pl.subplot(111, polar=True)
        ax.set_title(pTit)
        for idx in range(len(self.lIdxVis)):
            idxVis = self.lIdxVis[idx]
            if self.coodLoc[idxVis, 1] < self.MaxDZ: 
                ax.plot(np.deg2rad(self.coodLoc[idxVis, 0]), np.sin(np.deg2rad(self.coodLoc[idxVis, 1])),'.b')
            else:
                ax.plot(np.deg2rad(self.coodLoc[idxVis, 0]), np.sin(np.deg2rad(self.coodLoc[idxVis, 1])),'.r')
        ax.set_rmax(1.0)
        ax.grid(True)
        if pFile != None:
            pl.savefig(pFile)
        return ax
         
         
    def doMultiplot(self, pDate, pPasMn, pNpas):
        mDate = copy.copy(pDate)
        for idx in range(pNpas):
            mDate[4] = pDate[4]+idx*pPasMn
            print mDate
            self.computCoordLoc(mDate)
            self.selectVisibleStar(90)
            mfile="/home/colley/temp/montage/starView%03d.png"%idx
            self.plotVisible( "Cerro Pachon %d mn"%(idx*pPasMn), mfile)


    def doMultiplotSelect(self, pDate, pPasMn, pNpas):
        mDate = copy.copy(pDate)
        for idx in range(pNpas):
            mDate[4] = pDate[4]+idx*pPasMn
            print mDate
            self.computCoordLoc(mDate)
            self.selectVisibleStar(90)
            mfile="/home/colley/temp/montage/starViewSelect%03d.png"%idx
            ax = self.plotVisible( "Cerro Pachon %d mn"%(idx*pPasMn))
            for idx in self.selectStar:                
                ax.plot(np.deg2rad(self.coodLoc[idx, 0]), np.sin(np.deg2rad(self.coodLoc[idx, 1])),'*y',  markersize=10)
            pl.savefig(mfile) 
        
#
# STAT
#
    def histoDistZen(self):
        pl.figure()
        pl.hist(self.coodLoc[:,1], 20)
        
        
        