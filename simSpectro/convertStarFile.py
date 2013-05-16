# -*- coding: utf-8 -*-


def convertStarFile(fileIn, fileOut, coefWL = 1.0 , coefFlux= 1.0):
    # to read a spectrum in J/m2/s/nm   
    fileIn = open(fileIn, 'r')
    fileOut = open(fileOut, 'w')
    for line in fileIn.readlines():
        if line[0] == "#":
            fileOut.write(line)
        line = line.split()
        lamb = float(line[0])*coefWL
        dedl = float(line[1])*coefFlux
        fileOut.write("%f %g\n"%(lamb, dedl))
    fileIn.close()
    fileOut.close()
    print 'Done.'
        
            
#convertStarFile('example/alpha_lyr_stis_004.txt','example/alpha_lyr_stis_004nm.txt',0.1)
convertStarFile('example/alpha_lyr_stis_004nm.txt', 'example/alpha_lyr_stis_004nmBis.txt',1, 1e-2)
