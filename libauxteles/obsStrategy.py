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
import ephem as eph 

S_LevelPlot = 5


class ObsStrategyReal01(object):
    def __init__(self, fileCat):
        oTemp = hip.StarCatalog(0)
        self.oCat =  hip.StarCatalog(0)      
        self.oCat =  oTemp.read(fileCat)
        self.MaxDZ = 60 # distance zenital max
        self.deltaTimeObs = 20 # in minute
        
    
    def getRisingSettingSun(self, pDate):
        """
        pDate = [2012, 12, 30]
        
        return: pyepehem date object
          * UT time setting for [2012, 12, 29]or [2012, 12, 30]
          * UT time rising for [2012, 12, 30]
          
        coherent with 
        http://www.imcce.fr/fr/ephemerides/phenomenes/rts/rts.php
        """
        sun = eph.Sun() 
        lsst = eph.Observer() 
        lsst.lon="-70:44:58" 
        lsst.lat="-30:14:40" 
        # test with Paris
        #lsst.lon="2:21:07" 
        #lsst.lat="48:51:24" 
        lsst.horizon = '-12'
        lsst.date='%d/%d/%d'%(pDate[0], pDate[1], pDate[2])        
        try:
            rising  = lsst.next_rising(sun, use_center=True)            
            setting = lsst.next_setting(sun, use_center=True)
            if setting > rising:
                # get setting day before
                setting = lsst.previous_setting(sun, use_center=True)    
        except eph.AlwaysUpError:
            print "ephem error"
            return None
#        print lsst.previous_rising(sun, use_center=True)
#        print lsst.previous_setting(sun, use_center=True)                
#        print lsst.next_rising(sun, use_center=True)
#        print lsst.next_setting(sun, use_center=True)
        print  setting,  rising
        return [setting, rising]
          
          
    def selectMag(self, pMin , pMax):
        self.oCat.selectMag(pMin, pMax)
         
    
    def _ymdhm(self, pDate):
        """
        pDate: format [YYYY, MM, DD, HH, MM]
        convert hm in decimal day
        return  [YYYY, MM, DD, .DD]
        """
        return  [pDate[0], pDate[1], pDate[2], (pDate[3] + pDate[4]/60.)/24.]
    
     

    def selectAlgo01(self, pTbeg, pTend):
        """
        selection 4 stars
        """
        # compute position a pTbeg
        self.oCat.computCoordRefObs(pTbeg)
        self.oCat.selectVisibleStar(self.MaxDZ)        
        # possible stars near horizon    
        print  self.coodLoc[self.lIdxSelect]
        coodLoc = self.oCat.coodLoc 
        lIdxSelect = self.oCat.lIdxSelect
        y = np.sin(np.deg2rad(coodLoc[lIdxSelect, 1]))*np.sin(np.deg2rad(coodLoc[lIdxSelect, 0]))        
        lIdxHor1 = np.where(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>-135,coodLoc[lIdxSelect, 0]<-90), coodLoc[lIdxSelect, 1] > 55 ))[0]
        lIdxHor2 = np.where(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>-80,coodLoc[lIdxSelect, 0]<-20), coodLoc[lIdxSelect, 1] > 55 ))[0]
        if  S_LevelPlot > 2:
            print lIdxHor1, lIdxHor2
            idx = np.concatenate((lIdxHor1, lIdxHor2))
            self.plotStarWithIdx(idx, "Star possible")            
        # possible stars near zenith
        lIdxZen3 = np.where(np.logical_and(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>90,coodLoc[lIdxSelect, 0]<180), y < 0.08 ), coodLoc[lIdxSelect, 1] <30))[0]
        lIdxZen4 = np.where(np.logical_and(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>0,coodLoc[lIdxSelect, 0]<90),  y < 0.08), coodLoc[lIdxSelect, 1] <30))[0]
        if  S_LevelPlot > 2:
            print lIdxHor1, lIdxHor2
            idx = np.concatenate(( lIdxZen3, lIdxZen4))
            self.plotStarWithIdx(idx, "Star possible")
        # selection star
        selectStarVis=  [lIdxHor1[0], lIdxHor2[0], lIdxZen3[0], lIdxZen4[0]]   
        print  selectStarVis
        self.selectStar = [lIdxSelect[idx] for idx in selectStarVis]
        if  S_LevelPlot > 2:                        
            self.plotStarWithIdx(selectStarVis, "Star select")


    def _selectNearZenith(self):
        coodLoc = self.oCat.coodLoc
        lIdxSelect = self.oCat.lIdxSelect
        y = np.sin(np.deg2rad(coodLoc[lIdxSelect, 1]))*np.sin(np.deg2rad(coodLoc[lIdxSelect, 0]))  
        self.lIdxZen3 = np.where(np.logical_and(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>95,coodLoc[lIdxSelect, 0]<145), y < 0.08 ), coodLoc[lIdxSelect, 1] <30))[0]
        self.lIdxZen4 = np.where(np.logical_and(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>0,coodLoc[lIdxSelect, 0]<80),  y < 0.08), coodLoc[lIdxSelect, 1] <30))[0]
        if  S_LevelPlot > 2:            
            idx = np.concatenate(( self.lIdxZen3, self.lIdxZen4))
            self.oCat.plotIdxSelectStar(idx, "possible star")
    
    
    def _selectNearHorizon(self):
        coodLoc = self.oCat.coodLoc 
        lIdxSelect = self.oCat.lIdxSelect            
        self.lIdxHor1 = np.where(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>-115,coodLoc[lIdxSelect, 0]<-90), coodLoc[lIdxSelect, 1] > 55 ))[0]
        self.lIdxHor2 = np.where(np.logical_and(np.logical_and(coodLoc[lIdxSelect, 0]>-70,coodLoc[lIdxSelect, 0]<-30), coodLoc[lIdxSelect, 1] > 55 ))[0]
        if  S_LevelPlot > 2:
            #print self.lIdxHor1, self.lIdxHor2
            idx = np.concatenate((self.lIdxHor1, self.lIdxHor2))
            self.oCat.plotIdxSelectStar(idx, "possible star")            
    
    
    def _choiceStars(self):
        pass
    
    
    def _computeRefObsStarTarget(self):
        for idx, tMJD in enumerate(self.timeObs):
            self.computCoordRefObsMJD(tMJD, self.idxStarTarget[idx])
    
    
    def selectAlgo02(self, pDate):
        # define number of observation
        # =============================
        Tephem= self.getRisingSettingSun(pDate)
        if Tephem == None:
            print "Can't observe with condition !"
            return 
        print Tephem
        #print list(Tephem[0].triple())
        self.tDeb = Tephem[0]
        self.tEnd = Tephem[1]
        tDeb = mjd.tomjd_day(list(Tephem[0].triple()))
        tEnd = mjd.tomjd_day(list(Tephem[1].triple()))
        # coherence test
        duree = (tEnd-tDeb)*24.0
        if duree<0 or duree > 14: 
            print "duration NOK ? %fh"%duree
            return False
        print "duree de la nuit:", (tEnd-tDeb)*24.0
        # alloc data structure
        self.nbPeriodObs = int(duree*60/self.deltaTimeObs)
        self.coordStarTarget = np.empty((self.nbPeriodObs,4,2), dtype=np.float64)
        self.timeObs = tDeb + (np.arange(self.nbPeriodObs)*20.)/(60*24)
        print self.timeObs, tEnd
        self.idxStarTarget = np.ones((self.nbPeriodObs,4), dtype=np.int64)        
        print "nbPeriodObs:", self.nbPeriodObs
        # First selection at tDeb
        self.oCat.computCoordRefObsMJD(tDeb)
        self.oCat.selectVisibleStar(self.MaxDZ)        
        self._selectNearZenith()
        self._selectNearHorizon()
        # choice ref star random
        selectStarVis = []
        idx = np.random.randint(len(self.lIdxHor1))
        selectStarVis.append(self.lIdxHor1[idx])
        idx = np.random.randint(len(self.lIdxHor2))
        selectStarVis.append(self.lIdxHor2[idx])
        idx = np.random.randint(len(self.lIdxZen3))
        selectStarVis.append(self.lIdxZen3[idx])
        idx = np.random.randint(len(self.lIdxZen4))
        selectStarVis.append(self.lIdxZen4[idx])        
        selectStar = [self.oCat.lIdxSelect[idx] for idx in selectStarVis]        
        # compute position of ref. star for each observation period
        self.idxStarTarget *= np.array(selectStar)
        print self.idxStarTarget
        for idx, tMJD in enumerate(self.timeObs):
            self.coordStarTarget[idx] = self.oCat.computCoordRefObsMJD(tMJD, self.idxStarTarget[idx])
        # test risind selectin star
        for idxObs, tMJD in enumerate(self.timeObs):
            for idxStar in range(4):
                if self.coordStarTarget[idxObs, idxStar, 1] > self.MaxDZ:
                    # must replace this star 
                    print "Must replace star Id %d at %f"%(self.idxStarTarget[idxObs, idxStar],tMJD)
                    # select possible star
                    self.oCat.computCoordRefObsMJD(tMJD)
                    self.oCat.selectVisibleStar(self.MaxDZ)        
                    self._selectNearHorizon()
                    # select sky zone
                    Cpt1 = 0
                    Cpt2 = 0
                    for idxS in range(4):
                        if idxS != idxStar:
                            # test only other star
                            lgt = self.coordStarTarget[idxObs, idxS, 0]
                            print lgt
                            if lgt >= 0:                             
                                if lgt > 90 : Cpt1 +=1
                                else : Cpt2 +=1
                            else:
                                if lgt < -90 : Cpt1 +=1
                                else : Cpt2 +=1
                            print  Cpt1 , Cpt2                    
                    if Cpt1 > Cpt2:                    
                        lIdxHor = self.lIdxHor2
                        print "select zone 2"
                    else: 
                        lIdxHor = self.lIdxHor1
                        print "select zone 1"
                    idxStarSel = np.random.randint(len(lIdxHor))                    
                    idxNewStar = self.oCat.lIdxSelect[lIdxHor[idxStarSel]]
                    # replace idx star
                    self.idxStarTarget[idxObs:, idxStar] = idxNewStar
                    # update coordloc all , not optimal but easy ..
                    for idx2, tMJD2 in enumerate(self.timeObs):
                        self.coordStarTarget[idx2] = self.oCat.computCoordRefObsMJD(tMJD2, self.idxStarTarget[idx2])
                    
            
        
        
        
##############################
# Plot
##############################
        
         
         
    def doMultiplot(self, pDate, pPasMn, pNpas):
        mDate = copy.copy(pDate)
        for idx in range(pNpas):
            mDate[4] = pDate[4]+idx*pPasMn
            print mDate
            self.oCat.computCoordRefObs(mDate)
            self.selectVisibleStar(90)
            mfile="/home/colley/temp/montage/starView%03d.png"%idx
            self.plotSelectStars( "Cerro Pachon %d mn"%(idx*pPasMn), mfile)


    def doMultiplotSelect(self, pDate, pPasMn, pNpas):
        mDate = copy.copy(pDate)
        for idx in range(pNpas):
            mDate[4] = pDate[4]+idx*pPasMn
            print mDate
            self.oCat.computCoordRefObs(mDate)
            self.selectVisibleStar(90)
            mfile="/home/colley/temp/montage/starViewSelect%03d.png"%idx
            ax = self.plotSelectStarsz( "Cerro Pachon %d mn"%(idx*pPasMn))
            for idx in self.selectStar:                
                ax.plot(np.deg2rad(self.coodLoc[idx, 0]), np.sin(np.deg2rad(self.coodLoc[idx, 1])),'*y',  markersize=10)
            pl.savefig(mfile)
            
            
    def doMultiplotAlgo02(self):
        t0 = self.timeObs[0]
        for idxTime, tMJD in enumerate(self.timeObs):
            self.oCat.computCoordRefObsMJD(tMJD)
            self.oCat.selectVisibleStar(90)
            mfile="/home/colley/temp/montage/starViewAlgo02-%03d.png"%idxTime
            d0 = (tMJD-t0)*24
            hour = int(d0)
            mnt = (d0-hour)*60
            sTime = "%dh %dmn"%(hour, mnt)
            print sTime
            ax = self.oCat.plotSelectStars( "Cerro Pachon %s"%sTime, self.MaxDZ)
            for idx in self.idxStarTarget[idxTime]:                
                ax.plot(np.deg2rad(self.oCat.coodLoc[idx, 0]), np.sin(np.deg2rad(self.oCat.coodLoc[idx, 1])),'*y',  markersize=10)
            pl.savefig(mfile)
            print mfile
        
        
#
# STAT
#
    def histoDistZen(self):
        pl.figure()
        pl.hist(self.coodLoc[:,1], 20)
        
        
        