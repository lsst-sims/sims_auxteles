'''
Created on 6 nov. 2012

@author: colley
'''

import numpy as np
import math
import pyfits as pf
import scipy.interpolate as sci
import scipy as sp
import pylab as pl



S_verbose = 0


def removeAllEltlike(myList,elt):
    while(True):
        try:
            myList.remove(elt)
        except:
            return 
        
        
def getFirstCharNotBlanck(myStr):
    for elt in  myStr:
        if elt != ' ' : return elt
        
        
def removeBlanckInString(myStr):
    for idx, elt in  enumerate(myStr):
        if   elt != ' ':
            return myStr[idx:]
    return ''
            
            
def readTextFileColumn(fileName, sep = None):
    """
    read text file with N columns and P lines, so return a real array with dimension PxN
    sep = None : process blank and tabulation
    """
    # test existance du fichier
    listReal = []
    cptLine = 0
    cptLineTot = 0
    fd = open(fileName, 'r')
    for line in fd.readlines():
        cptLineTot +=1
        #print "'%s'"%line 
        if line[0] == "#":
            continue
        if sep ==None:
            listCol = line.split()
        else:
            listCol = line.split(sep)
        if S_verbose >0 : print listCol
        removeAllEltlike(listCol,'')  
        removeAllEltlike(listCol,'\n')  
        #print listCol
        if len(listCol) == 0:
            if S_verbose >0 : print "ligne vide ... continue"
            continue        
        if getFirstCharNotBlanck(listCol[0]) == "#":
            if S_verbose >0 : print "ligne commente ... continue"
            continue
        if cptLine == 0:
            removeAllEltlike(listCol,'\n')
            nbCol = len(listCol)
        elif nbCol != len(listCol):
            if len(listCol) >= 1:
                clearStr = removeBlanckInString(listCol[0])
                if clearStr[0] == '\n':
                    if S_verbose >0 : print "ligne vide ... continue"
                    continue                    
            print "ERROR line %d, read %d elements instead %d"%(cptLineTot, len(listCol),nbCol )
            return False
        for Elt in listCol:
            try:
                _RealVal = float(Elt)                    
            except:
                print "ERROR line %d, can't convert '%s' in real number"%(cptLineTot, Elt)
                return False
            listReal.append(_RealVal)
        cptLine += 1
        # to numpy array
    outArray = np.array(listReal)
    outArray = np.reshape(outArray,(cptLine, nbCol))
    return outArray


def interpolBSpline(pXin, pYin, pXout, pFlagPlot = False):
    """
    pXin : increasing ordered array 
    pYin : same size as pXin
    pXout :increasing ordered array, completely in pXin
    """
    assert len(pXin) == len(pYin)
    # test coherence pXin pXout
    if  pXin[0] >  pXout[0]: raise
    if  pXin[-1] <  pXout[-1]: raise       
    tck = sci.splrep(pXin, pYin)
    Yout = sci.splev(pXout, tck)
    # test Yout ???    
    if pFlagPlot:
        pl.figure()
        pl.plot(pXin, pYin)
        pl.plot(pXout, Yout)
        pl.grid()
        pl.title("interpolBSpline function ")
        pl.legend(["Raw","Interpol"])
    return Yout


def interpolLinear(pXin, pYin, pXout, pFlagPlot = False):
    """
    pXin : increasing ordered array 
    pYin : same size as pXin
    pXout :increasing ordered array, completely in pXin
    """
    assert len(pXin) == len(pYin)
    # test coherence pXin pXout
    if  pXin[0] >  pXout[0]: raise
    if  pXin[-1] <  pXout[-1]: raise       
    oInter = sci.interp1d(pXin, pYin)
    Yout = oInter(pXout)
    # test Yout ???    
    if pFlagPlot:
        pl.figure()
        pl.plot(pXin, pYin)
        pl.plot(pXout, Yout)
        pl.grid()
        pl.title("interpolBSpline function ")
        pl.legend(["Raw","Interpol"])
    return Yout


def indexInIntervalCheck(pX, xMin, xMax):
    """
    pX numpy array sorted small to tall
    xMin, xMax in same unit than pX and min(pX) <= xMin < xMax  <=  max(pX) 
    
    return :
      * None if min(pX) <= xMin < xMax  <=  max(pX) is False
      * idxMin, idxMax as  xMin <= pX[idxMin:idxMax] <= xMax      
    """
    idx = np.where(pX >= xMin)[0]
    if len(idx) == 0:
        print "ERROR [indexInInterval]: xMin <  min(pX)"
        return None
    IdxMin = idx[0]
    idx = np.where(pX <= xMax)[0]
    if len(idx) == 0:
        print "ERROR [indexInInterval]: xMax > max(pX)"
        return None
    IdxMax = idx[-1]
    return IdxMin, IdxMax


def indexInInterval(pX, xMin, xMax):
    """
    pX numpy array sorted small to tall
    xMin, xMax in same unit than pX 
    
    return :
      * idxMin, idxMax as  xMin <= pX[idxMin:idxMax] <= xMax      
    """
    idx = np.where(pX >= xMin)[0]
    if len(idx) == 0:
        IdxMin = 0
    else:
        IdxMin = idx[0]
    idx = np.where(pX <= xMax)[0]
    if len(idx) == 0:
        IdxMax = len(pX)
    else:
        IdxMax = idx[-1]
    return IdxMin, IdxMax
    
    
def gauss(x, sigma, mu=0.):
    return np.exp(-(x-mu)*(x-mu)/(2.*sigma*sigma))/(sigma*np.sqrt(2*np.pi))
    
    
def downgradeResol(xIn, yIn, resIn, resOut, ref= 650, xOut=None):
    """
    xIn: can be not regular
    yIn :  signal value
    resIn :  resolution In 
    resOut : resolution out
    ref : resolution = ref / delta
    xOut : optional 
    OUTPUT: 
      convolve signal by gaussian function with sigma ((ref/resOut)**2 - (ref/resIn)**2)**.5
    """
    sigma = ( (ref/resOut)**2 - (ref/resIn)**2)
    if sigma <=0 :
        # no change
        print "[convolGauss] no change"
        return yIn        
    sigma = sigma**.5
    minDelta = np.min(np.diff(xIn))
    xcst = np.arange(xIn[0], xIn[-1], minDelta)
    tck = sci.splrep(xIn, yIn)
    ycst = sci.splev(xcst, tck, der=0)
    # put the gaussian kernel into an array
    kernel = []
    for e in np.arange(-10.*sigma, 10.*sigma+minDelta, minDelta):
        kernel.append(gauss(e, sigma, 0.))
    kernel = np.array(kernel)
    # convolve    
    convycst = sp.convolve(ycst, kernel, 'same')
    # and put it again on the initial grid
    tck = sci.splrep(xcst, convycst)
    if xOut != None:        
        yOut = sci.splev(xOut, tck, der=0)/kernel.sum()
    else:
        yOut = sci.splev(xIn, tck, der=0)/kernel.sum()
    return yOut
    
   
    
