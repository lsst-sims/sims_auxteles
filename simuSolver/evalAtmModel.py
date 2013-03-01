'''
Created on 21 fevr. 2013

@author: colley
'''
import basicSimuV01 as bsimu
import pylab as pl
import numpy as np


def plotStdERtransVersusSNR():
    #lSNR = [200.0, 400.0, 800.0, 1600.0, 3200.0]
    lSNR = [400.0, 800.0, 1600.0]
    lEC = []
    laeroEC =[] 
    lnoAeroEC =[] 
    for snr in lSNR:
        np.random.seed(103)    
        oSim = bsimu.SimuAtmStarSolve2(1,16)
        oSim.oObs.addNoisebySNR (snr )
        guess = oSim.oObs.getGuessDefault()
        oSim.oSol.solveOnlyAtm(guess[4:])
        tempStar = oSim.oStarCat.getAllTemperature()
        oSim.oSol._parEst = np.concatenate((tempStar, oSim.oSol._parEst))
        errRelTot = oSim.oSol.transmisErrRelAtmAll()
        print errRelTot.mean(), errRelTot.std()
        lEC.append(errRelTot.std())
        errRel = oSim.oSol.aeroErrRel(oSim.oSol._parEst)
        laeroEC.append(errRel.std())
        errRel = oSim.oSol.noAeroErrRel(oSim.oSol._parEst)
        lnoAeroEC.append(errRel.std())
    print laeroEC
    print lnoAeroEC
    pl.figure()
    pl.plot(np.array(lSNR), np.array(lEC))
    pl.grid()
    pl.figure()
    pl.title("Standard deviation true relative error atmopsheric transmission estimated")
    pl.loglog(np.array(lSNR), np.array(lEC))
    pl.loglog(np.array(lSNR), np.array(lEC),'*')
    pl.grid()
    pl.ylabel("standard deviation")
    pl.xlabel("SNR noise")
    
    pl.figure()   
    pl.loglog(np.array(lSNR), np.array(laeroEC))
    pl.loglog(np.array(lSNR), np.array(lnoAeroEC))
    pl.legend(["aero",'noAero'])
    pl.grid()
    pl.ylabel("standard deviation")
    pl.xlabel("SNR noise")
#    pl.figure()   
#    #pl.loglog(np.array(lSNR), np.array(laeroEC))
#    pl.loglog(np.array(lSNR), np.array(lnoAeroEC))
#    #pl.legend(["aero",'noAero'])
#    pl.grid()
#    pl.ylabel("standard deviation")
#    pl.xlabel("SNR noise")
    pl.show()
   
if __name__ == '__main__':
    plotStdERtransVersusSNR()