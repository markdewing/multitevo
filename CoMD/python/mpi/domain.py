
import parallel

import numpy as np

class Domain:
    def __init__(self, xproc, yproc, zproc, globalExtent):
        assert(xproc * yproc * zproc == parallel.NRanks)
        self.procGrid = [xproc, yproc, zproc]

        self.procCoord = np.zeros(3, dtype=np.int)
        myrank = parallel.MyRank
        self.procCoord[0] = myrank % self.procGrid[0]
        myrank /= self.procGrid[0]
        self.procCoord[1] = myrank % self.procGrid[1]
        self.procCoord[2] = myrank / self.procGrid[1]

        self.globalMin = np.zeros(3, dtype=np.float)
        self.globalMax = np.zeros(3, dtype=np.float)
        self.globalMax[:] = globalExtent[:]
        self.globalExtent = self.globalMax - self.globalMin


        self.localExtent = self.globalExtent / self.procGrid
        self.localMin = self.globalMin + self.procCoord     * self.localExtent
        self.localMax = self.globalMin + (self.procCoord+1) * self.localExtent


    def processorNum(self, dix, diy, diz):
        ix = (self.procCoord[0] + dix + self.procGrid[0]) % self.procGrid[0]
        iy = (self.procCoord[1] + diy + self.procGrid[1]) % self.procGrid[1]
        iz = (self.procCoord[2] + diz + self.procGrid[2]) % self.procGrid[2]

        return ix + self.procGrid[0]*(iy + self.procGrid[1]*iz)
