'''
Created on 17 mars 2014

@author: colley
'''
import matplotlib.pyplot as pl
import numpy as np
import tools as tl
import copy 


#
# Static Constant
#

S_DoPlot = 0




def getKeywordTP7():
    """
    return list of concatenated keyword with underscore, like "O3_TRANS" to 
    identified a colmun in TP7 file
    """
    _key1TP7 = "FREQ COMBIN    H2O   UMIX     O3  TRACE     N2    H2O MOLEC AER+CLD  HNO3 AER+CLD    -LOG    CO2     CO    CH4    N2O     O2    NH3     NO    NO2    SO2  CFC11  CFC12  CFC13  CFC14  CFC22 CFC113 CFC114 CFC115 CLONO2   HNO4 CHCL2F   CCL4   N2O5"
    _key2TP7 = "CM-1  TRANS  TRANS  TRANS  TRANS  TRANS   CONT   CONT   SCAT  TRANS  TRANS abTRNS  COMBIN  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS  TRANS"
    keyTP7 = [a+'_'+b for a, b in zip(_key1TP7.split(), _key2TP7.split())]
    return keyTP7


def printCardLine(lCard):
    sOut = ""
    for line in lCard:
        sOut += line
    print sOut


def readTP7_multiCard(pFileTP7, pListCol):
    """
    read TP7 multicard file, extract card and a column selection given by  
    keyword in pListCol list
    return:
     * numpy array dim 3 (number Card x number line wavelength x len(pListCol))
     * list on line card 
    """
    KeyTP7 = getKeywordTP7()
    try:
        lSelectIdx = [KeyTP7.index(key) for key in pListCol]
    except:
        print "Unknown keyword", key
        raise
    try:
        fTP7 = open(pFileTP7, 'r')
    except:
        print "can't open file ", pFileTP7
        raise
    # status line : "DATA"  or "CARD"
    status = "CARD"
    lArray = []
    lCard = []
    glArray = []
    glCard = []
    for line in fTP7:
        if status == "CARD":            
            if line.find('CM-1') >=0:
                status = "DATA"
                # remove last line FREQ ...
                lCard.pop()
                glCard.append(lCard)
                lCard = []
            else:
                lCard.append(line)
        else:
            if line.find('-9999.') >=0:
                glArray.append(lArray)
                lArray = []
                status = "CARD"
            else:
                lCol = line.split()
                #print line             
                sa = [float(lCol[idx]) for idx in lSelectIdx]
                lArray.append(sa)
    # to numpy array
    # same size 
    sizeArray = [len(lA) for lA in glArray]
    if len(set(sizeArray)) !=1:
        print "number of element if diffrent", sizeArray
        raise
    return np.array(glArray), glCard
    

def readTP7(pFile, pListCol, pWLMin=None, pWLmax=None):
    """
    read tp7 output modtran and select column
    WARNING : must add # to catch header !
    parameters: 
     * pFile :  string with path and name file
     * pListCol :  list of selected column
     * pWLMin, pWLmax : min max wavelength in Angstrom
     
    return:
     * numpy array where column 0 is  wavelength in Angstrom and other columns 
       are selected by pListCol
    """  
    _keyTP7 = "FREQ COMBIN    H2O   UMIX     O3  TRACE     N2    H2Ob MOLEC AER+CLD  HNO3 AER+CLD    -LOG    CO2     CO    CH4    N2O     O2    NH3     NO    NO2    SO2  CFC11  CFC12  CFC13  CFC14  CFC22 CFC113 CFC114 CFC115 CLONO2   HNO4 CHCL2F   CCL4   N2O5"
    KeyTP7 = _keyTP7.split()
    # reorder keyword like tp7 file sense, import for pyplot.legend()
    lSelectIdx = [KeyTP7.index(key) for key in pListCol]
    # read tp7 file
    aRet = tl.readTextFileColumn(pFile)
    if aRet is None: 
        print "ERROR: ", aRet
        return
    # selected columns
    lSelectIdx.insert(0,0) # add wavele,gth colmun
    aRet = np.flipud(aRet.T[lSelectIdx].T)
    print aRet
    # cm^-1 to angstrom
    wlAngs = 1.0e8/aRet[:,0]
    if pWLmax != None:        
        IdxRet = tl.indexInIntervalCheck(wlAngs, pWLMin, pWLmax)
        if IdxRet[0] is None:
            raise
            return None
        # suppress line
        aRet = aRet[IdxRet[0]: IdxRet[1]]
    # replace wl in angstrom
    aRet[:,0] = 1.0e8/aRet[:,0]
    # growing sense
    # plot
    if S_DoPlot:
        pl.figure()
        pl.title('MODTRAN TP7 transmission selected')
        pl.xlabel("Angstrom")
        for col in aRet[:,1:].T:
            pl.plot(aRet[:,0], col)
        pl.grid()
        pl.legend(pListCol, loc="best")
    return aRet


def extractTemplateFromTP7(pFileTP7, pDirOut="../output/"):
    """
    create from tp7 numpy file template
    """
    lKey = ["FREQ_CM-1" , "H2O_TRANS", "H2O_CONT", "O3_TRANS", "O2_TRANS",\
            "MOLEC_SCAT", "NO2_TRANS"]
    aTP7, aCard = readTP7_multiCard(pFileTP7, lKey)
    lAM = []
    lTplt = []
    for idx, aTemp in enumerate(aTP7):
        pl.figure()
        angle = float(aCard[idx][5].split()[2])
        am = 1/np.cos(np.deg2rad(angle))
        lAM.append(am)
        pl.title("Air mass : %3.2f"%am)
        pl.ylim(0.1, 1.0)
        pl.plot(aTemp[:,0], aTemp[:,1])
        pl.grid()
        aRet = np.flipud(aTemp)
        wlAngs = 1.0e7/aRet[:,0]
        tO3 = aRet[:, lKey.index("O3_TRANS")]
        tMols = aRet[:, lKey.index("MOLEC_SCAT")]
        tMola = aRet[:, lKey.index("O2_TRANS")]*\
                aRet[:, lKey.index("NO2_TRANS")]
        tH2O = aRet[:, lKey.index("H2O_TRANS")]*\
               aRet[:, lKey.index("H2O_CONT")]
        tplMOD = np.array([tO3, tMols, tMola, tH2O], dtype=np.float32)
        #amStr = "%2.1f"%am
        lTplt.append(tplMOD)
    print wlAngs
    np.save(pDirOut+'tpltMODTRAN.npy', np.array(lTplt))
    np.save(pDirOut+'tpltMODTRAN_AM.npy', np.array(lAM))    
    np.save(pDirOut+'tpltMODTRAN_wlnm.npy', wlAngs)
        
    