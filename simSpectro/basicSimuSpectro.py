#
# basic compute photon-electron for star 
#

import numpy as np


C_speedLigth_m_s = 3.0e8
C_Planck_J_s = 6.63e-34
C_ergs2J = 1.0e-7


def wl2freq(wl):
    """
    wl [m]
    
    return [Hz]
    """
    return C_speedLigth_m_s/wl

def freq2wl(fq):
    """
    fq [Hz]
    
    return [m]
    """
    return C_speedLigth_m_s/fq


def stellarRadius(temp, magApp, parallaxe):   
    """
    IN: 
        temp [K]
        magApp : visual magnitude [noUnit]
        parallaxe [mas]    
    OUT: 
        radius [km]
    
    source : http://cas.sdss.org/dr6/en/proj/advanced/hr/radius1.asp
    """ 
    absMag = magApp - 5*np.log10(1000.0/parallaxe) + 5
    # 4.83 abs magnitude of sun
    rap = ((5800.0/temp)**2)*(2.51**(4.83-absMag))**0.5
    # 696342 km sun radius
    return 696342 * rap

def parallaxeToKm(parallaxe):
    """
    parallaxe [mas]
    """
    return 3.085e16/parallaxe


useVega = False


if useVega:
    #  flux density Vega ( M -0.5, T 9600, G 3.95) from Kurucz 93 at 5550 A
    #  ergs.cm^{-2}.s^{-1}.A^{-1}
    fdVega=5.51e7
    print "Flux density Vega Kurucz93 : %g ergs.cm^{-2}.s^{-1}.A^{-1}"%fdVega 
    
    #To convert to observed
    #flux at Earth, multiply by a factor of (R/D)^2 where R is the stellar
    #radius, and D is the distance to Earth.
    # 
    # vega constants
    #  temperatue = 9600K
    #  visual magnitude =  0.03
    #  parallaxe 129 mas
    #
    vegaRadius = stellarRadius(9550, 0.03, 129)
    vegaDist = parallaxeToKm(129)
    print "vegaRadius(R): %g km"%vegaRadius
    print "vegaDist(D)  : %g km"% vegaDist
    
    coefEarth = (vegaRadius/vegaDist)**2
    print "(R/D)^2: ", coefEarth
    
    
    fdStarEarth = fdVega*coefEarth
else:
    fdStarEarth =  1.70E-14    
    # 100 coeff to ergs.cm^{-2}.s^{-1}.A^{-1} from J.m{-2}.s^{-1}.nm^{-1}
    #fdStarEarth =  (2.75152e-11)*100
       
       
print "Flux density star at Earth :  %g  ergs.cm^{-2}.s^{-1}.A^{-1}"%fdStarEarth
#
# Convert en J.m{-2}.s^{-1}.m^{-1}
#
#ergs = 1e-7 J
#cm^{-2} = 1e4 m^{-2}
#A^{-1} = 1e10 m^{-1}
fdStarEarthSI = fdStarEarth*C_ergs2J*1.0e4*1.0e10
print "Flux density star at Earth :  %g  J.m{-2}.s^{-1}.m^{-1}"%fdStarEarthSI

    
wl = 600e-9
print "estimation at %g m, or %g Hz"%(wl,wl2freq(wl) )
demiLargeurWL = (3.86e-10)/2
#demiLargeurWL = 0.7e-9
#diamTeles = 4
diamTeles = 4
SurfaceMir = (np.pi/4)*(diamTeles**2)
tpsExposure = 240
DeltaFreq = wl2freq(wl - demiLargeurWL) - wl2freq(wl + demiLargeurWL)
print "DeltaFreq:", DeltaFreq
energyPhot = C_Planck_J_s * wl2freq(wl)
coef2Freq =  (wl**2) / C_speedLigth_m_s
print coef2Freq

NphotonRaw = fdStarEarthSI*SurfaceMir*tpsExposure*coef2Freq*DeltaFreq/energyPhot
print "nb photon tot %g"%NphotonRaw

# transmission instru and atm
nAtm = 0.783
nMir = 0.911
nSlit = 0.78
#nSlit = 0.68
nGrism = 0.96
nGrismOrder = 0.13
#nGrismOrder = 1
nCCD = 0.79

nTot = nAtm*nMir*nSlit*nGrism*nGrismOrder*nCCD
print 'Trans Total %g'%nTot

# 
NphotonTrans = NphotonRaw*nTot
print "nb photon trans %g"%NphotonTrans



ref = 43244
#ref = 8.042e9

print NphotonTrans/ref, ref/NphotonTrans
