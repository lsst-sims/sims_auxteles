import pyfits

class Spec2D:

    def __init__(self):
        hdulist = 0
        trimsec = 'TRIMSEC'
        xmin = xmax = ymin = ymax = 0
        
    def open(self, filename):
        hdulist = pyfits.open(filename)
        tmp = hdulist[0].header[trimsec]
        xmin = int(tmp[tmp.find('[')+1:tmp.find(':')])
        xmax = int(tmp[tmp.find(':')+1:tmp.find(',')])
        ymin = int(tmp[tmp.find(',')+1:tmp.find(':', 3)])
        ymax = int(tmp[tmp.find(':', 3)+1:tmp.find(']')])

    def AddSpectrum(self, starspec, seeing):
        sigma = seeing/(2.*sqrt(2.*log(2.)))
        for i in range(xmin, xmax):
            for j in range(ymin, ymax):
                hdulist.data[i, j] += starspec.get2Delecccd(xmax-xmin, ymax-ymin, y, sigma)

    def write(self):
        hdulist.flush()

    def write2fits(self, filename):
        hdulist.writeto(filename)


