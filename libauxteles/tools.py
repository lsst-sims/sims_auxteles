'''
Created on 6 nov. 2012

@author: colley
'''

import math

import numpy as np
import pyfits as pf
import pylab as pl
import scipy as sp
import scipy.interpolate as sci
import os.path as osp


S_verbose = 0
S_doPlot = False

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
    """
    remove first blanck
    """
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


def productOf2array(x1, y1, x2, y2, interpolMeth = interpolLinear):
    """
    product of 2 array not defined for the same value
    return product is defined at intersection x1, x2 with x1 value
        array x1:        |****|       
        array x2:   |******|
        array xInter:    |*|
    """   
    x1Inter = None
    if x2[0] < x1[0]:
        """
        3 cases: A, B, C
        array 1:        |****|
        array 2: A |***|
        array 2: B |******|
        array 2: C |*************|
        """        
        if x2[-1] < x1[0]:
            """
            cas A
            """
            print "no intersection"
            return None
        xMin = x1[0]
        if x2[-1] < x1[-1]:
            """
            cas B
            """            
            xMax = x2[-1]
        else:
            """
            cas C
            """            
            xMax = x1[-1]
            x1Inter = x1
            y1Inter = y1
    else:
        """
        3 cases: D, E, F
        array 1:        |*****|
        array 2: D        |**|
        array 2: E        |******|
        array 2: F              |***|
        """        
        if x1[-1] < x2[0]:
            """
            cas F
            """
            print "no intersection"
            return None
        xMin = x2[0]
        if x2[-1] < x1[-1]:
            """
            cas D
            """            
            xMax = x2[-1]
        else:
            """
            cas E
            """
            xMax = x1[-1]
    #
    #print "xMin, xMax ",xMin, xMax
    if x1Inter == None:
        IdxMin, IdxMax = indexInInterval(x1, xMin, xMax)
        x1Inter = x1[IdxMin: IdxMax]
        y1Inter = y1[IdxMin: IdxMax]
    y2Inter = interpolMeth(x2, y2, x1Inter)
    return x1Inter, y2Inter*y1Inter  
        


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
    sigma = ( (ref/resOut)**2 - (ref/resIn)**2 )
    if sigma <=0 :
        # no change
        print "[convolGauss] resOut > resIn , no downgrad", resIn,resOut
        if xOut != None:
            tck  = sci.interp1d(xIn, yIn)      
            yOut = tck(xOut)
        else:
            yOut = yIn
        return yOut      
    sigma = sigma**.5
    minDelta = np.min(np.diff(xIn))
    nbPoint = 2*int((xIn[-1] - xIn[0])/minDelta)
    xcst = np.linspace(xIn[0], xIn[-1], nbPoint, True)
    #xcst = np.arange(xIn[0], xIn[-1], minDelta)
    oInter = sci.interp1d(xIn, yIn)
    ycst = oInter(xcst)
    # put the gaussian kernel into an array
    kernel = []
    for e in np.arange(-10.*sigma, 10.*sigma+minDelta, minDelta):
        kernel.append(gauss(e, sigma, 0.))
    kernel = np.array(kernel)
    # convolve    
    convycst = sp.convolve(ycst, kernel, 'same')
    # and put it again on the initial grid
    oInter = sci.interp1d(xcst, convycst)
    if S_doPlot:
        pl.figure()
        pl.plot(xIn, yIn)
    if xOut != None:        
        yOut = oInter(xOut)/kernel.sum()
        if S_doPlot:pl.plot(xOut, yOut)
    else:
        yOut = oInter(xIn)/kernel.sum()
        if S_doPlot:pl.plot(xIn, yOut)    
    return yOut
    
    
def downgradeResolSpline(xIn, yIn, resIn, resOut, ref= 650, xOut=None):
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
        print "[convolGauss] resOut > resIn , no downgrad", resIn,resOut
        if xOut != None:
            tck = sci.splrep(xIn, yIn)      
            yOut = sci.splev(xOut, tck, der=0)
        else:
            yOut = yIn
        return yOut      
    sigma = sigma**.5
    minDelta = np.min(np.diff(xIn))
    xcst = np.arange(xIn[0], xIn[-1], minDelta/2)
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
    if S_doPlot:
        pl.figure()
        pl.plot(xIn, yIn)
    if xOut != None:        
        yOut = sci.splev(xOut, tck, der=0)/kernel.sum()
        if S_doPlot:pl.plot(xOut, yOut)
    else:
        yOut = sci.splev(xIn, tck, der=0)/kernel.sum()
        if S_doPlot:pl.plot(xIn, yOut)    
    return yOut

def parallaxeToKm(parallaxe):
    """
    parallaxe [mas]
    """
    return 3.085e16/parallaxe


def stellarRadius(temp, magApp, parallaxe):   
    """
    IN: 
        temp [K]
        magApp [noUnit]
        parallaxe [mas]    
    OUT: 
        radius [km]
    
    source : http://cas.sdss.org/dr6/en/proj/advanced/hr/radius1.asp
    """ 
    absMag = magApp - 5*np.log10(1000.0/parallaxe) + 5
    # 4.83 abs magnitude of sun
    rap = ((5800.0/temp)**2)*(2.512**(4.83-absMag))**0.5
    # 696342 km sun radius
    return 696342 * rap


def coefKuruczEarth(temp, magApp, parallaxe):
    stelRad = stellarRadius(temp, magApp, parallaxe)
    distEarth = parallaxeToKm(parallaxe)
    #print stelRad, distEarth
    return (stelRad/distEarth)**2
    

def plotFile(nameFile):
    pl.figure()
    pl.title(nameFile)
    ret = readTextFileColumn(nameFile)
    pl.plot(ret[:,0], ret[:,1])
    
    
def plotListFile(lNameFile):
    pl.figure()
    lgd = []  
    for mfile in lNameFile:
        ret = readTextFileColumn(mfile)
        pl.plot(ret[:,0], ret[:,1])
        lgd.append(mfile.split('/')[-1])
    pl.legend(lgd)


# linear interpolation between (x1,y1) and (x2,y2) at point x = x_val
# by G. Blanc
def x_linear_interpolate_1D(x_val, x1, x2, y1, y2):
    if (x_val > x2 or x_val < x1): 
        print "Pb in x_linear_interpolate_1D : ", x_val, x1, x2, y1, y2
        return -99
    a = (y1-y2)/(x1-x2)
    b = (y2*x1-y1*x2)/(x1-x2)
    return a * x_val + b


# by G. Blanc
def integre_trapeze(x, fn):
    if (len(x) != len(fn)):
        print "integre_trapeze : x et fn n'ont pas la meme dimension !"
        return -99.
    som = 0.
    for i in range(0,len(x)-1):
        som = som + (x[i+1]-x[i]) * (fn[i] + fn[i+1]) / 2.

    return som

# Return index of tab just below value
# by G. Blanc
def get_index_tab_below(value, tab):
    for i in range(0, len(tab)):
        if (i == len(tab)-1):
            index = i
            break
        if (value > tab[i]):
            continue
        else:
            if (i == 0):
                index = 0
                break
            else:
                index = i-1
                break
    return index


# Interpolate on spectrum flux...
# by G. Blanc
def sumSFDforABmag(flux_spectrum, wavelength_spectrum, filter_trans, wavelength_filter):

    if ((len(flux_spectrum) != len(wavelength_spectrum)) or (len(filter_trans) != len(wavelength_filter))):
        print "Houston!, We have a problem; size of tab non OK in convolution"
        return -99

    convol = np.zeros(len(wavelength_filter))

    #ws = index of wave spectrum just below the value of the wave of filter
    ws = get_index_tab_below(wavelength_filter[0], wavelength_spectrum)
    for wf in range(0,len(wavelength_filter)-1):
        # To faster the code
        # ONLY for filters...
        if (filter_trans[wf] < 0.001):
            continue

        # increase ws to get lambda spectrum close to lambda of filter 
        while (wavelength_filter[wf] > wavelength_spectrum[ws+1]):
            ws += 1

        # Interpolate spectrum flux in x = wavelength_filter[wf]
        flux_interp = x_linear_interpolate_1D(wavelength_filter[wf], wavelength_spectrum[ws], wavelength_spectrum[ws+1], \
                                                  flux_spectrum[ws], flux_spectrum[ws+1])

        #print "DEBUG convol : wf = ", wf, wavelength_filter[wf], wavelength_spectrum[ws], wavelength_filter[wf+1], wavelength_spectrum[ws+1],  flux_spectrum[ws], flux_interp, flux_spectrum[ws+1]

        #convol[wf] = flux_interp * filter_trans[wf] #* wavelength_filter[wf]
        # for photon counting detectors
        convol[wf] = flux_interp * filter_trans[wf] * wavelength_filter[wf]

    flux_tot = integre_trapeze(wavelength_filter, convol)
    return flux_tot



def getFileNameLSSTfilter(letter):
    lLetter=["u","g","r","i","z","y"]
    lfilter = letter.lower()
    if lfilter not in lLetter:
        return None
    return '../../data/filter/short_total_%s.dat'%lfilter



def deltaMagnitudeAB(SFDstar, wl_SFDstar, TransTrue, TransEst, wl_Trans, letterFilter):
    """
    SFDstar : spectral flux density star,  W.x{-2}.nm{-1}
     wl_SFDstar: wavelength in nm
     TransTrue : true transmission 
     TransEst : estimated transmission
     wl_Trans :  wavelength in nm
     letterFilter : lettter in {ugrizy}
     
     return -2.5log_{10}(int_SFD.TransTrue.Filter/(int_SFD.TransEst.Filter))         
    """
    fileFilter = getFileNameLSSTfilter(letterFilter)
    if fileFilter == None:
        print '[deltaMagnitudeAB] ERROR letter filter "%s" ', letterFilter
        return None
    ret = readTextFileColumn(fileFilter)
    wl = ret[:,0] # in nm
    trFilter = ret[:,1]
    print fileFilter
    #print trFilter[:10]
    # True transmission 
    wl_prod , trTot = productOf2array(wl, trFilter, wl_Trans, TransTrue)
    sumTranstrue = sumSFDforABmag(SFDstar, wl_SFDstar, trTot, wl_prod)
    if sumTranstrue < 0:
        print '[deltaMagnitudeAB] ERROR sumSFDforABmag(TransTrue) return ',sumTranstrue
        return None
    # estimated transmission 
    wl_prod , trTot = productOf2array(wl, trFilter, wl_Trans, TransEst)
    sumTransEst  = sumSFDforABmag(SFDstar, wl_SFDstar, trTot, wl_prod)
    if sumTranstrue < 0:
        print '[deltaMagnitudeAB] ERROR sumSFDforABmag(TransEst) return ',sumTransEst
        return None
    return -2.5*np.log10(sumTranstrue/sumTransEst)
    

def getRootPackage():
    return osp.normpath(getDirectory(__file__)+'/..')

  
def getDirectory(path):
    """
    path=/path/to/my/file.xx
    return 
    /path/to/my
    """
    ipath = path[::-1]
    idx = ipath.find('/')
    return ipath[idx+1:][::-1]

def getOutputDir():
    return osp.normpath(getDirectory(__file__)+'/../output')


if __name__ == "__main__":
    #plotFile('/home/colley/temp/lsst/auxteles/short_total_g.dat')
    lFile = ['/home/colley/temp/lsst/auxteles/short_total_u.dat']
    lFile.append('/home/colley/temp/lsst/auxteles/short_total_g.dat')
    lFile.append('/home/colley/temp/lsst/auxteles/short_total_r.dat')
    lFile.append('/home/colley/temp/lsst/auxteles/short_total_i.dat')
    lFile.append('/home/colley/temp/lsst/auxteles/short_total_z.dat')
    lFile.append('/home/colley/temp/lsst/auxteles/short_total_y.dat')
    plotListFile(lFile)
    pl.xlabel('nm')  
    pl.ylabel('transmission')    
    pl.title("LSST Filter")
    pl.gri
    pl.show()