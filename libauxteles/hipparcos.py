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


    def extractFieldForAuxteles(self, oCat):
        cptLine = 0
        cptNOK = 0
        for line in self.fd:
            #print line
            try: 
                Id = int(line[8:14])
                dec = float(line[64:76])
                ra = float(line[51:63])
                mag = float(line[41:46])
                para = float(line[79:86])
                
            except:
                cptNOK +=1
                continue
            oCat.Id[cptLine]=Id
            oCat.coord[cptLine] = np.array([ra, dec])
            oCat.mag[cptLine] = mag
            oCat.paral[cptLine] = para
            spec = line[435:447] 
            #print line[346]  
            oCat.spectral[cptLine,0] = ord(line[346])         
            for idx, sym in enumerate(spec):
                if idx < 7:
                    oCat.spectral[cptLine,idx+1] = ord(sym)
            cptLine += 1
            if spec[0] == 'O' : print spec
            #if cptLine > 10: break
        print 'NOK ', cptNOK
        print 'OK ', cptLine
                
                    

class StarCatalog(object):
    
    def __init__(self, nbStar):
        self.nbStar = nbStar
        self.Id = np.zeros(nbStar, dtype=np.int32)
        self.coord = np.zeros((nbStar,2), dtype=np.float64)
        self.mag = np.zeros(nbStar, dtype=np.float32)
        self.paral = np.zeros(nbStar, dtype=np.float32)
        self.spectral = np.zeros((nbStar,8), dtype=np.uint8)
        
        
    def printCat(self, beg, end):
        for idx in range(beg, end):
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
                str = ''
                for s in self.spectral[idx]:
                    str += chr(s)
                spec = str[1:3]
                lum = str[3:]
                a=lum.find(" ")
                if a > 0 : lum = lum[:a]
                ret = kur.convertMK(spec, lum)
                if ret != None:
                    self.kurucz[idx,:] = ret
                    nOK += 1
                else:
                    nNOK += 1
        print "OK %d on %d"%(nOK, nOK+nNOK)
            
    
    def removeStarNOK(self):
        lIdx=[]
        for idx in range(self.nbStar):
            if np.sum(self.kurucz[idx]) == 0 :
                lIdx.append(idx)
                #print "%d is nok"%idx
        # delete line
        self.Id=np.delete(self.Id, lIdx, 0)
        self.coord=np.delete(self.coord, lIdx, 0)
        self.kurucz=np.delete(self.kurucz, lIdx, 0)
        self.mag=np.delete(self.mag, lIdx, 0)
        self.paral=np.delete(self.paral, lIdx, 0)
        self.spectral=np.delete(self.spectral, lIdx, 0)
        self.nbStar = len(self.Id)
        print "nb star ok ", self.nbStar
            
            
    def save(self, filename):
        f=open(filename, "wb")
        pk.dump(self, f, pk.HIGHEST_PROTOCOL)
        f.close()
        
        
    def read(self, filename):
        f=open(filename, "rb")
        obj = pk.load(f)
        f.close()
        return obj
