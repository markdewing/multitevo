
import parallel

import numpy as np

class Halo:
    HALO_X_MINUS = 0
    HALO_X_PLUS = 1
    HALO_Y_MINUS = 2
    HALO_Y_PLUS = 3
    HALO_Z_MINUS = 4
    HALO_Z_PLUS = 5

    HALO_X_AXIS = 0
    HALO_Y_AXIS = 1
    HALO_Z_AXIS = 2

    def __init__(self, boxes, domain):
        self.bufCapacity = 0
        self.nCells = np.zeros(6, dtype=np.int)
        self.pbcFactor = np.zeros((6,3))
        self.nbrRank = np.zeros(6, dtype=np.int)

        self.nCells[self.HALO_X_MINUS] = 2*(boxes.gridSize[1] + 2) * (boxes.gridSize[2] + 2)
        self.nCells[self.HALO_Y_MINUS] = 2*(boxes.gridSize[0] + 2) * (boxes.gridSize[2] + 2)
        self.nCells[self.HALO_Z_MINUS] = 2*(boxes.gridSize[0] + 2) * (boxes.gridSize[1] + 2)
        self.nCells[self.HALO_X_PLUS] = self.nCells[self.HALO_X_MINUS]
        self.nCells[self.HALO_Y_PLUS] = self.nCells[self.HALO_Y_MINUS]
        self.nCells[self.HALO_Z_PLUS] = self.nCells[self.HALO_Z_MINUS]

        self.cellList = []
        for ii in range(6):
            self.cellList.append(self.mkAtomCellList(boxes, ii, self.nCells[ii]))

        if domain.procCoord[self.HALO_X_AXIS] == 0:
            self.pbcFactor[self.HALO_X_MINUS, self.HALO_X_AXIS] = 1.0

        if domain.procCoord[self.HALO_X_AXIS] == domain.procGrid[self.HALO_X_AXIS]-1:
            self.pbcFactor[self.HALO_X_PLUS, self.HALO_X_AXIS] = -1.0

        if domain.procCoord[self.HALO_Y_AXIS] == 0:
            self.pbcFactor[self.HALO_Y_MINUS, self.HALO_Y_AXIS] = 1.0

        if domain.procCoord[self.HALO_Y_AXIS] == domain.procGrid[self.HALO_Y_AXIS]-1:
            self.pbcFactor[self.HALO_Y_PLUS, self.HALO_Y_AXIS] = -1.0

        if domain.procCoord[self.HALO_Z_AXIS] == 0:
            self.pbcFactor[self.HALO_Z_MINUS, self.HALO_Z_AXIS] = 1.0

        if domain.procCoord[self.HALO_Z_AXIS] == domain.procGrid[self.HALO_Z_AXIS]-1:
            self.pbcFactor[self.HALO_Z_PLUS, self.HALO_Z_AXIS] = -1.0

        self.nbrRank[self.HALO_X_MINUS] = domain.processorNum(-1, 0, 0)
        self.nbrRank[self.HALO_X_PLUS]  = domain.processorNum(+1, 0, 0)
        self.nbrRank[self.HALO_Y_MINUS] = domain.processorNum(0, -1, 0)
        self.nbrRank[self.HALO_Y_PLUS]  = domain.processorNum(0, +1, 0)
        self.nbrRank[self.HALO_Z_MINUS] = domain.processorNum(0, 0, -1)
        self.nbrRank[self.HALO_Z_PLUS]  = domain.processorNum(0, 0, +1)


    def haloExchange(self, sim):
        for iAxis in range(3):
            self.exchangeData(iAxis, sim)

    def exchangeData(self, iAxis, sim):
        faceM = 2*iAxis
        faceP = faceM + 1

        bufM = self.loadAtomsBuffer(sim, faceM)
        bufP = self.loadAtomsBuffer(sim, faceP)

        nbrRankM = self.nbrRank[faceM]
        nbrRankP = self.nbrRank[faceP]

        recvBufP = parallel.sendReceiveParallel(bufM, nbrRankM, nbrRankP)
        recvBufM = parallel.sendReceiveParallel(bufP, nbrRankP, nbrRankM)

        self.unloadAtomsBuffer(sim, faceM, recvBufP)
        self.unloadAtomsBuffer(sim, faceP, recvBufM)

    def loadAtomsBuffer(self, sim, face):
        shift = self.pbcFactor[face,:] * sim.atoms.bounds
        nCells = self.nCells[face]
        cellList = self.cellList[face]

        buf = []
        for iCell in range(nCells):
            iBox = cellList[iCell]
            for ii in range(sim.boxes.nAtoms[iBox]):
                iOff = iBox*sim.boxes.MAXATOMS + ii
                entry = {}
                entry['gid'] = sim.atoms.gid[iOff]
                entry['type'] = sim.atoms.species[iOff]
                entry['rx'] = sim.atoms.r[iOff][0] + shift[0]
                entry['ry'] = sim.atoms.r[iOff][1] + shift[1]
                entry['rz'] = sim.atoms.r[iOff][2] + shift[2]
                entry['px'] = sim.atoms.p[iOff][0]
                entry['py'] = sim.atoms.p[iOff][1]
                entry['pz'] = sim.atoms.p[iOff][2]
                buf.append(entry)

        return buf

    def unloadAtomsBuffer(self, sim, face, buf):
        for entry in buf:
            gid = entry['gid']
            species = entry['type']
            rx = entry['rx']
            ry = entry['ry']
            rz = entry['rz']
            px = entry['px']
            py = entry['py']
            pz = entry['pz']
            sim.boxes.putAtomInBox(sim.atoms, gid, species, rx, ry, rz, px, py, pz)

    def mkAtomCellList(self, boxes, iFace, nCells):
        xBegin = -1
        xEnd = boxes.gridSize[0] + 1
        yBegin = -1
        yEnd = boxes.gridSize[1] + 1
        zBegin = -1
        zEnd = boxes.gridSize[2] + 1

        if iFace == self.HALO_X_MINUS: xEnd = xBegin + 2
        if iFace == self.HALO_X_PLUS:  xBegin = xEnd - 2
        if iFace == self.HALO_Y_MINUS: yEnd = yBegin + 2
        if iFace == self.HALO_Y_PLUS:  yBegin = yEnd - 2
        if iFace == self.HALO_Z_MINUS: zEnd = zBegin + 2
        if iFace == self.HALO_Z_PLUS:  zBegin = zEnd - 2

        cells = []
        for ix in range(xBegin,xEnd):
            for iy in range(yBegin,yEnd):
                for iz in range(zBegin,zEnd):
                    cells.append(boxes.getBoxFromTuple(ix,iy,iz))

        assert(nCells == len(cells))
        return cells


    def sortAtomsInCell(self, atoms, boxes, iBox):
        nAtoms = boxes.nAtoms[iBox]

        sortType = [('gid', np.int),
                    ('type', np.int),
                    ('rx', np.float), ('ry', np.float), ('rz', np.float),
                    ('px', np.float), ('py', np.float), ('pz', np.float)
                   ]

        tmp = np.zeros(nAtoms, dtype=sortType)

        begin = iBox * boxes.MAXATOMS
        iTmp = 0
        for ii in range(begin, begin + nAtoms):
            tmp[iTmp]['gid'] = atoms.gid[ii]
            tmp[iTmp]['type'] = atoms.species[ii]
            tmp[iTmp]['rx'] = atoms.r[ii, 0]
            tmp[iTmp]['ry'] = atoms.r[ii, 1]
            tmp[iTmp]['rz'] = atoms.r[ii, 2]
            tmp[iTmp]['px'] = atoms.p[ii, 0]
            tmp[iTmp]['py'] = atoms.p[ii, 1]
            tmp[iTmp]['pz'] = atoms.p[ii, 2]
            iTmp += 1

        tmp.sort(order='gid')

        begin = iBox * boxes.MAXATOMS
        iTmp = 0
        for ii in range(begin, begin + nAtoms):
            atoms.gid[ii] = tmp[iTmp]['gid']
            atoms.species[ii] = tmp[iTmp]['type'] = atoms.species[ii]
            atoms.r[ii, 0] = tmp[iTmp]['rx']
            atoms.r[ii, 1] = tmp[iTmp]['ry']
            atoms.r[ii, 2] = tmp[iTmp]['rz']
            atoms.p[ii, 0] = tmp[iTmp]['px']
            atoms.p[ii, 1] = tmp[iTmp]['py']
            atoms.p[ii, 2] = tmp[iTmp]['pz']
            iTmp += 1


