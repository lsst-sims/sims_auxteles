from scipy import *

from constants import *


class filter:

    def __init__(self):
        self.wl = []
        self.nu = []
        self.T = []

    def readl(self, filename, wlfactor):
        file = open(filename, 'r')
        for line in file.readlines():
            line = line.split()
            lamb = float(line[0])*wlfactor
            trans = float(line[1])
            self.wl.append(lamb)
            self.nu.append(clight/lamb)
            self.T.append(trans)
        file.close()
        self.wl = array(self.wl)
        self.nu = array(self.nu)
        self.T = array(self.T)

    def getTl(self, lamb):
        f = interpolate.interp1d(self.wl, self.T, bounds_error=False, fill_value = 0.)
        return f(lamb)

    def getTn(self, freq):
        f = interpolate.interp1d(flipud(self.nu), flipud(self.T), bounds_error=False, fill_value = 0.)
        return f(freq)
