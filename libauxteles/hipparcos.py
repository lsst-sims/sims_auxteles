'''
Created on 13 sept. 2013

@author: colley
'''
import numpy as np 
import kurucz as kur
import pickle as pk


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
        self.Id = np.zeros(nbStar, dtype=np.int32)
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
        
    
    def convertSpectralToKurucz(self):
        self.kurucz = np.zeros((self.nbStar,3), dtype=np.float32)
        nOK = 0
        nNOK = 0
        for idx in range(self.nbStar):
            if self.Id[idx] != 0:
                mystr = ''
                for s in self.spectral[idx]:
                    mystr += chr(s)
                spec = mystr[1:3]
                lum = mystr[3:]
                a=lum.find(" ")
                if a > 0 : lum = lum[:a]
                ret = kur.convertMK(spec, lum)
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
            
            
    def save(self, filename):
        f=open(filename, "wb")
        pk.dump(self, f, pk.HIGHEST_PROTOCOL)
        f.close()
        
        
    def read(self, filename):
        f=open(filename, "rb")
        obj = pk.load(f)
        f.close()
        return obj

