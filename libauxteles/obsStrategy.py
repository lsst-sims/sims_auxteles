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
import pickle as pk
import tools as tl


S_LevelPlot = 5



class DataPositionStarTarget(object):
    """
    class to store result of simulation of position of stars target
    """
    
    def setPeriodObs(self, pDate, pBeg, pEnd):
        """
        pDate: like [2013, 6, 21]
        pBeg:[mjd]
        pEnd:[mjd]
        """
        self.Date = pDate
        self.Beg = pBeg
        self.End = pEnd
        
        
    def setStarPosition(self, pId, pPos, pTime):
        """        
        pId:   (n,4)   Identificator hipparcos catalog
        pPos:  (n,4,2) position in horizontal coordinate [deg] (azimuth, dist_zenital)
        pTime: (n,1)   [mjd] time of observation
        
        with: n number of period observation 
        """
        self.IdPos   = pId
        self.Pos  = pPos
        self.Time = pTime
        
        
    def setStarConstant(self, pId, pKur, pPar, pMag, pMK):
        """
        pId: (s,1) Identificator hipparcos catalog
        pKur: (s,3) Kurucz parameters (temp, grav, met)
        pPar: (s,1) parallaxe gives by Hipparcos
        pMag: (s,1) magnitude gives by Hipparcos
        pMK : (s,8) spectral type Morgan-keeman
        
        with: s number of star observed
        """
        self.IdCst = pId
        self.Kur = pKur
        self.Par = pPar
        self.Mag = pMag
        self.MK = pMK
        
        
    def save(self, filename):
        f=open(filename, "wb")
        pk.dump(self, f, pk.HIGHEST_PROTOCOL)
        f.close()
        
        
    def read(self, filename):
        f=open(filename, "rb")
        obj = pk.load(f)
        f.close()
        return obj        


    def __str__(self):
        ret ='ID\t(temp,\t\tgra,\tmet)\tparal\tmag\tTypeSpect\n'
        for idx, ID in enumerate(self.IdCst):
            mkStr = hip.arrayInt2str(self.MK[idx])
            ret += '%d\t(%.1f,\t%.2f,\t%.2f)\t%.1f\t%.1f\t%s\n'%(ID, self.Kur[idx,0], self.Kur[idx,1],self.Kur[idx,2], self.Par[idx], self.Mag[idx], mkStr)
        ret += '\nHipparcos ID star observed by period:\n'
        for idx, star in enumerate(self.IdPos):
            ret += '%03d %s\n'%(idx+1, str(star))
        return ret 


###############################################################################
class ObsStrategyReal(object):
###############################################################################    
    def __init__(self, fileCat):
        oTemp = hip.StarCatalog(0)
        self.oCat =  hip.StarCatalog(0)      
        self.oCat =  oTemp.read(fileCat)
        self.oCat.convertSpectralToKurucz()
        self.oCat.removeStarKuruczNOK()
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
        """
        pDate like [2013, 6, 21]
        
        defined:
         * self.tDeb, self.tEnd : time begin end observation
         
        return:
           Position of 4 star at each instant, manage rising and setting star
        """
        # define number of observation
        # =============================
        Tephem= self.getRisingSettingSun(pDate)
        if Tephem == None:
            print "Can't observe with condition !"
            return 
        print Tephem
        #print list(Tephem[0].triple())
        self.Date = pDate
        self.tDeb = Tephem[0]
        self.tEnd = Tephem[1]
        tDeb = mjd.tomjd_day(list(Tephem[0].triple()))
        tEnd = mjd.tomjd_day(list(Tephem[1].triple()))
        print "mjd ", tDeb, tEnd
        # coherence test
        duree = (tEnd-tDeb)*24.0
        if duree<0 or duree > 14: 
            print "duration NOK ? %fh"%duree
            return False
        print "duree de la nuit:", (tEnd-tDeb)*24.0
        # alloc data structure
        #
        # number of period observation
        self.nbPeriodObs = int(duree*60/self.deltaTimeObs)
        # star target local coordinates
        self.coordStarTarget = np.empty((self.nbPeriodObs,4,2), dtype=np.float64)
        # time period observation
        self.timeObs = tDeb + (np.arange(self.nbPeriodObs)*20.)/(60*24)
        print self.timeObs, tEnd
        # index star target in catalog self.oCat
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
        print self.oCat.kurucz[selectStar]
        for idx, tMJD in enumerate(self.timeObs):
            self.coordStarTarget[idx] = self.oCat.computCoordRefObsMJD(tMJD, self.idxStarTarget[idx])
        # test rising selection stars
        for idxObs, tMJD in enumerate(self.timeObs):
            for idxStar in range(4):
                if self.coordStarTarget[idxObs, idxStar, 1] > self.MaxDZ:                                    
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
                    print self.oCat.kurucz[idxNewStar]
                    # replace idx star
                    self.idxStarTarget[idxObs:, idxStar] = idxNewStar
                    # update coordloc all , not optimal but easy ..
                    for idx2, tMJD2 in enumerate(self.timeObs):
                        self.coordStarTarget[idx2] = self.oCat.computCoordRefObsMJD(tMJD2, self.idxStarTarget[idx2])
                        
                        
    def selectAlgo03(self, pDate, nbPeriodObs=24):
        """
        pDate like [2013, 6, 21]
        
        defined:
         * self.tDeb, self.tEnd : time begin end observation
         
        return:
          only 4 stars and fix nbPeriodObs
        """
        # define number of observation
        # =============================
        Tephem= self.getRisingSettingSun(pDate)
        if Tephem == None:
            print "Can't observe with condition !"
            return 
        print Tephem
        #print list(Tephem[0].triple())
        self.Date = pDate
        self.tDeb = Tephem[0]
        self.tEnd = Tephem[1]
        tDeb = mjd.tomjd_day(list(Tephem[0].triple()))
        tEnd = mjd.tomjd_day(list(Tephem[1].triple()))
        print "mjd ", tDeb, tEnd
        # coherence test
        dureeHour = (tEnd-tDeb)*24.0
        if dureeHour<0 or dureeHour > 14: 
            print "duration NOK ? %fh"%dureeHour
            return False
        print "dureeHour de la nuit:", (tEnd-tDeb)*24.0
        # alloc data structure
        #
        # number of period observation
        self.deltaTimeObs = 5.0
        self.nbPeriodObs = int(dureeHour*60/self.deltaTimeObs)
        # star target coordinates
        self.coordStarTarget = np.empty((self.nbPeriodObs,4,2), dtype=np.float64)
        # time period observation
        self.timeObs = tDeb + np.arange(self.nbPeriodObs)*(self.deltaTimeObs/(60*24))
        # index star target in catalog self.oCat
        self.idxStarTarget = np.ones((self.nbPeriodObs,4), dtype=np.int64)        
        print "nbPeriodObs:", self.nbPeriodObs
        print self.timeObs
        # First selection at tDeb        
        self.oCat.computCoordRefObsMJD(tDeb)
        self.oCat.selectVisibleStar(self.MaxDZ)        
        self._selectNearZenith()
        self._selectNearHorizon()
        # choice ref star random
        selectStarVis = []        
        selectStarVis.append(self.lIdxHor1[0])        
        selectStarVis.append(self.lIdxHor2[0])
        selectStarVis.append(self.lIdxZen3[0])        
        selectStarVis.append(self.lIdxZen4[0])        
        selectStar = [self.oCat.lIdxSelect[idx] for idx in selectStarVis]            
        # compute position of ref. star for each observation period
        self.idxStarTarget *= np.array(selectStar)
        print self.idxStarTarget
        print self.oCat.kurucz[selectStar]
        for idx, tMJD in enumerate(self.timeObs):
            self.coordStarTarget[idx] = self.oCat.computCoordRefObsMJD(tMJD, self.idxStarTarget[idx])
        # test rising selection stars
        rising = False
        for idxObs, tMJD in enumerate(self.timeObs):            
            for idxStar in range(4):
                print tMJD, idxStar
                if self.coordStarTarget[idxObs, idxStar, 1] > self.MaxDZ:                                    
                    print "Must replace star Id %d at %f"%(self.idxStarTarget[idxObs, idxStar],tMJD)
                    rising = True
                if rising: break
            if rising: break
        # fix end time when first star rising
        tEnd = tMJD
        print "mjd ", tDeb, tEnd
        # coherence test
        dureeHour = (tEnd-tDeb)*24.0
        if dureeHour<0 or dureeHour > 14: 
            print "duration NOK ? %fh"%dureeHour
            return False
        print "dureeHour de la nuit:", (tEnd-tDeb)*24.0
        # alloc data structure
        #
        # number of period observation
        self.nbPeriodObs = nbPeriodObs*1.0
        deltaTpsObsHour =  dureeHour/self.nbPeriodObs
        self.deltaTimeObs =  deltaTpsObsHour*60
        # star target coordinates
        self.coordStarTarget = np.empty((self.nbPeriodObs,4,2), dtype=np.float64)
        # time period observation
        self.timeObs = tDeb + np.arange(self.nbPeriodObs)*deltaTpsObsHour/24        
        print "new time observation: ", self.timeObs, tEnd
        # index star target in catalog self.oCat
        self.idxStarTarget = np.ones((self.nbPeriodObs,4), dtype=np.int64)        
        print "nbPeriodObs:", self.nbPeriodObs
        # First selection at tDeb        
        self.idxStarTarget *= np.array(selectStar)
        print "idxStarTarget: ", self.idxStarTarget
        print self.oCat.kurucz[selectStar]
        for idx, tMJD in enumerate(self.timeObs):
            self.coordStarTarget[idx] = self.oCat.computCoordRefObsMJD(tMJD, self.idxStarTarget[idx])
        print "Result:"
        print '  nb obs: ', self.nbPeriodObs
        print '  Time observation:\n', self.timeObs
        print '  position: \n', self.coordStarTarget
                        
                        
    def saveResultAlgo02(self, pFile):
        oData = DataPositionStarTarget()
        oData.setPeriodObs(self.Date, self.tDeb, self.tEnd)
        IdStar = self.oCat.Id[self.idxStarTarget]
        oData.setStarPosition(IdStar, self.coordStarTarget, self.timeObs)
        # alloc
        nbStar=0
        lId = [] # identificator Hipparcos
        lIdx = [] # index 
        print "idxStarTarget:",  self.idxStarTarget  
        print "IdStar: ",IdStar
        rIdxStarTarget = self.idxStarTarget.ravel()
        for idx,iStar in enumerate(IdStar.ravel()):
            print idx,iStar, lId
            if not (iStar in lId):
                lId.append(iStar)
                lIdx.append(rIdxStarTarget[idx])
                nbStar += 1
        aId = np.array(lId, dtype=np.int64)
        aKur = np.empty((nbStar,3), dtype=np.float32)
        aPar = np.empty(nbStar, dtype=np.float32)
        aMag = np.empty(nbStar, dtype=np.float32)
        aSpect = np.zeros((nbStar,8), dtype=np.uint8)
        for idx, iStar in enumerate(aId):
            i = lIdx[idx]
            aKur[idx]   = self.oCat.kurucz[i]
            aPar[idx]   = self.oCat.paral[i]
            aMag[idx]   = self.oCat.mag[i]
            aSpect[idx] = self.oCat.spectral[i]    
        oData.setStarConstant(aId, aKur, aPar, aMag, aSpect)
        oData.save(pFile)
        print oData
        
        
        
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
            
            
    def doMultiplotAlgo(self, nameFile="starView"):
        t0 = self.timeObs[0]
        for idxTime, tMJD in enumerate(self.timeObs):
            self.oCat.computCoordRefObsMJD(tMJD)
            self.oCat.selectVisibleStar(90)
            mfile=tl.getRootPackage()+"/output/%s-%03d.png"%(nameFile,idxTime)
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
        
        
    def doPlotStarTraj(self):
        """
        trajectory of selected stars on one plot
        """
        pl.figure()
        ax = pl.subplot(111, polar=True)
        ax.set_title("Trajectory of stars selected")
        for aPos in self.coordStarTarget:
            for pos in aPos:                
                ax.plot(np.deg2rad(pos[0]), np.sin(np.deg2rad(pos[1])),'*y', markersize=10)

        
#
# STAT
#
    def histoDistZen(self):
        pl.figure()
        pl.hist(self.coodLoc[:,1], 20)
        
        
        