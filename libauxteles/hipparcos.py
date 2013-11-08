'''
Created on 13 sept. 2013

@author: colley
'''
import numpy as np 
import kurucz as kur
import pickle as pk
import pylab as pl
import Astrotools as astro
import MJDtools as mjd


def arrayInt2str(aInt):
    mystr = ''
    for s in aInt:
        mystr += chr(s) 
    return  mystr  


def MorganKeeman2Kurucz(sMK):
    """
    ok only if sMK in  [BAFGKM][0..9]['I','II','III','IV','V']
    """
    spec = sMK[0:2]
    lum = sMK[2:]
    a=lum.find(" ")
    if a > 0 : lum = lum[:a]
    ret = kur.convertMK(spec, lum)
    if ret == None:
        print "Can't convert: ", sMK
    return ret        

        
def MorganKeeman2Kurucz_robust(sMK):
    """
    manage 2 problems in sMK string
     * presence of / character, => try substring 
     * if sMK is only [BAFGKM][0..9], => add random luninosity   
     
    return:
     * None
     * temp, grav, met    
    """    
    aLum = ["I", "II","III", "III", "IV", "IV", "V", "V", "V"]
    sLum = 9
    a = sMK.find("/") 
    if a > 0:
        if a <= 2: 
            print " %s try with %s"%(sMK, sMK[a+1:])
            return MorganKeeman2Kurucz_robust(sMK[a+1:])         
        else:
            print " %s try with %s"%(sMK, sMK[:a])
            return MorganKeeman2Kurucz_robust(sMK[:a])            
    spec = sMK[0:2]
    lum = sMK[2:]
    if lum[0] == " ":
        #print " no lim ",sMK
        #return None
        lum = aLum[np.random.randint(sLum)]
    else:
        a=lum.find(" ")
        if a > 0 : lum = lum[:a]
    ret = kur.convertMK(spec, lum)
    if ret == None:
        print "Can't convert: ", sMK
#    else:
#        print "ok for ", sMK
    return ret        
        
   
    
    
class Hipparcos(object):
    
    def __init__(self, pathHip):
        """
        """
        self.fd = open(pathHip+'/hip_main.dat', 'r')
        self.nStar = 118218
    
    
    def close(self):
        self.fd.close()
    
        
    def getnbLine(self):
        cptLine = 0
        for line in self.fd:
            cptLine +=1
        print cptLine


    def extractSingleStarForAuxteles(self, oCat):
        print 'Select single star'
        cptLine = 0
        cptNOK = 0
        cptSingle = 0
        for line in self.fd:
            #print line
            if line[346] == ' ' and line[350] != 'S': 
                cptSingle += 1
            else:
                # multiple systems case or ambigus
                continue
            try: 
                Id = int(line[8:14])
                dec = float(line[64:76])
                ra = float(line[51:63])
                mag = float(line[41:46])
                para = float(line[79:86])                
            except:
                cptNOK +=1
                continue
            #print line[346],  line[350]  
            oCat.Id[cptLine]=Id
            oCat.coord[cptLine] = np.array([ra, dec])
            oCat.mag[cptLine] = mag
            oCat.paral[cptLine] = para
            spec = line[435:447]                     
            #oCat.spectral[cptLine,0] = ord(line[346])         
            for idx, sym in enumerate(spec):
                if idx < 8:
                    oCat.spectral[cptLine,idx] = ord(sym)
            cptLine += 1
            #if spec[0] == 'O' : print spec
            #if cptLine > 10: break
        #print 'NOK ', cptNOK
        print 'number single star OK ', cptLine
        #print 'single ', cptSingle
                
                    

class StarCatalog(object):
    
    def __init__(self, nbStar):
        self.nbStar = nbStar
        self.Id = np.zeros(nbStar, dtype=np.int32) # ID hipparcos
        self.coord = np.zeros((nbStar,2), dtype=np.float64)
        self.mag = np.zeros(nbStar, dtype=np.float32)
        self.paral = np.zeros(nbStar, dtype=np.float32)
        self.spectral = np.zeros((nbStar,8), dtype=np.uint8)
        
        
    def printCat(self, beg, end, pIdx=None):
        if pIdx == None:
            myIdx = range(beg, end)
        else:
            myIdx = pIdx
        print "Id\tcoord\t\t\t\tmag\tparax\tspectrale type"
        for idx in myIdx:
            spec = ''
            for s in self.spectral[idx]:
                spec += chr(s) 
            print "%d\t(%.8f, %.8f)\t%.1f\t%.1f\t%s"%(self.Id[idx],self.coord[idx,0], self.coord[idx,1], self.mag[idx], self.paral[idx], spec)
        
    
    def convertSpectralToKuruczOld(self):
        """
        ok even if []
        """
        self.kurucz = np.zeros((self.nbStar,3), dtype=np.float32)
        nOK = 0
        nNOK = 0
        for idx in range(self.nbStar):
            if self.Id[idx] != 0:
                mystr = ''
                for s in self.spectral[idx]:
                    mystr += chr(s)
                spec = mystr[0:2]
                lum = mystr[2:]
                a=lum.find(" ")
                if a > 0 : lum = lum[:a]
                ret = kur.convertMK(spec, lum)
                if ret != None:
                    self.kurucz[idx,:] = ret
                    nOK += 1
                else:
                    nNOK += 1
        print "OK %d on %d"%(nOK, nOK+nNOK)
        
        
    def convertSpectralToKurucz(self):
        """
        manage / character (split in sub-string) and add a random luminosity is not available 
        """
        self.kurucz = np.zeros((self.nbStar,3), dtype=np.float32)
        nOK = 0
        nNOK = 0
        for idx in range(self.nbStar):
            if self.Id[idx] != 0:
                mystr = ''
                for s in self.spectral[idx]:
                    mystr += chr(s)
                #ret = MorganKeeman2Kurucz(mystr)
                ret = MorganKeeman2Kurucz_robust(mystr)
                if ret != None:
                    self.kurucz[idx,:] = ret
                    nOK += 1
                else:
                    nNOK += 1        
        print "OK %d on %d"%(nOK, nOK+nNOK)

            
    def _removeIdx(self, lIdx):
        self.Id=np.delete(self.Id, lIdx, 0)
        self.coord=np.delete(self.coord, lIdx, 0)
        if hasattr(self, 'kurucz'):
            self.kurucz=np.delete(self.kurucz, lIdx, 0)
        self.mag=np.delete(self.mag, lIdx, 0)
        self.paral=np.delete(self.paral, lIdx, 0)
        self.spectral=np.delete(self.spectral, lIdx, 0)
        self.nbStar = len(self.Id)
        

        
    def removeStarKuruczNOK(self):
        lIdx=[]
        for idx in range(self.nbStar):
            if np.sum(self.kurucz[idx]) == 0 :
                lIdx.append(idx)
                #print "%d is nok"%idx
        # delete line
        self._removeIdx(lIdx)
        #print "nb star ok ", self.nbStar
            
            
    def removeNotUsed(self):
        lIdx=[]
        for idx in range(self.nbStar):
            if np.sum(self.coord[idx]) == 0 :
                lIdx.append(idx)
                #print "%d is nok"%idx
        # delete line
        self._removeIdx(lIdx)
        #print "nb star ok ", self.nbStar
            
    
    def findTypeStar(self, pType, pSsType=None, pLum=None, pIdx=None):    
        """
        problem with luminosity selection , but ok for V
        """    
        idx = []
        if pIdx == None:
            myIdx = np.arange(self.nbStar)
        else:
            myIdx = pIdx
        #print myIdx
        typeS = ord(pType)
        if pLum != None:
            lum = [ord(l) for l in pLum]
        if pSsType != None:    
            ssT = ord(pSsType)
        if  pSsType == None:
            if pLum == None:             
                idx = np.where(self.spectral[myIdx,0] == typeS)[0]
            else:
                sLum = len(pLum)
                idx = np.where(np.logical_and(self.spectral[myIdx,0] == typeS, np.all(self.spectral[myIdx,2:2+sLum] == lum,1)))[0]               
        else:
            if pLum == None:                
                idx = np.where(np.logical_and(self.spectral[myIdx,0] == typeS, self.spectral[myIdx,1] == ssT))[0]
            else:
                sLum = len(pLum)
                idx = np.where(np.logical_and(np.logical_and(self.spectral[myIdx,0] == typeS, np.all(self.spectral[myIdx,2:2+sLum] == lum,1)), self.spectral[myIdx,1] == ssT))[0]            
        return idx
                   
                   
    def selectMag(self, pMin , pMax):
        lIdx=[]
        for idx in range(self.nbStar):
            mag = self.mag[idx]
            if mag > pMax or  mag < pMin:
                lIdx.append(idx)
        self._removeIdx(lIdx)


    def selectVisibleStar(self, pElv):
        """
        fill array self.lIdxSelect index of star with zenith distance compatible with pElv
        """
        self.lIdxSelect = []
        for idx in range (self.nbStar):
            if self.coodLoc[idx, 1] < pElv: 
                self.lIdxSelect.append(idx)
          
#
# PLOT
#

    def plotIdxSelectStar(self, pPidx, pTit, pFile=None):
        pl.figure()
        print pFile       
        ax = pl.subplot(111, polar=True)
        ax.set_title(pTit)
        for idxS in pPidx:
            idx = self.lIdxSelect[idxS]
            ax.plot(np.deg2rad(self.coodLoc[idx, 0]), np.sin(np.deg2rad(self.coodLoc[idx, 1])),'*b')
        ax.set_rmax(1.0)
        ax.grid(True)
        if pFile != None:
            pl.savefig(pFile)
            
           
    def plotSelectStars(self, pTit, MaxDZ= 90, pFile=None):
        pl.figure()
        ax = pl.subplot(111, polar=True)
        ax.set_title(pTit)
        for idxSel in self.lIdxSelect:            
            if self.coodLoc[idxSel, 1] < MaxDZ: 
                ax.plot(np.deg2rad(self.coodLoc[idxSel, 0]), np.sin(np.deg2rad(self.coodLoc[idxSel, 1])),'.b')
            else:
                ax.plot(np.deg2rad(self.coodLoc[idxSel, 0]), np.sin(np.deg2rad(self.coodLoc[idxSel, 1])),'.r')
        ax.set_rmax(1.0)
        ax.grid(True)
        if pFile != None:
            pl.savefig(pFile)
        return ax

#
# observer referential
#
    def computCoordRefObsMJD(self, dateMJD, pLidx=None, pSite=None):
        """
        pSite: format [longitude, latitude]
        """        
        if pLidx == None:
            Lidx = range (self.nbStar)
        else:
            Lidx = pLidx
        coodLoc = np.empty((len(Lidx),2), dtype=np.float64)
        if pSite == None:
            for idx, idxStar in enumerate(Lidx):
                #print idxStar, self.nbStar, self.coord.shape
                coodLoc[idx] = astro.eq2loc([self.coord[idxStar,0], self.coord[idxStar,1],dateMJD] )                                   
        else:
            for idx, idxStar in enumerate(Lidx):
                coodLoc[idx] = astro.eq2loc([self.coord[idxStar,0], self.coord[idxStar,1], dateMJD, pSite[0], pSite[1] ] )
        if pLidx == None:
            self.coodLoc = coodLoc
            self.dateMJD = dateMJD
        else:
            return coodLoc
        
                
    def computCoordRefObs(self, pDate, pLidx=None, pSite=None):
        """
        pDate: format [YYYY, MM, DD, HH, MM]
        pSite: format [longitude, latitude]
        """       
        dateJD = mjd.tomjd(self._ymdhm(pDate))
        return self.computCoordRefObsMJD(dateJD, pLidx, pSite)
        
        
#
# I/O
#           
    def save(self, filename):
        f=open(filename, "wb")
        pk.dump(self, f, pk.HIGHEST_PROTOCOL)
        f.close()
        
        
    def read(self, filename):
        f=open(filename, "rb")
        obj = pk.load(f)
        f.close()
        return obj

