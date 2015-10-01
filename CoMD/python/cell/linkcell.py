from __future__ import print_function, division

import numpy as np
import math

def initLinkCells(sim, box_size):
    gridSize = np.array([int(b/sim.pot.cutoff) for b in box_size], dtype=np.int)
    boxSize = box_size*1.0/gridSize
    nLocalBoxes = gridSize[0]*gridSize[1]*gridSize[2]

    nHaloBoxes = 2*((gridSize[0] + 2) * (gridSize[1] + gridSize[2] + 2) + (gridSize[1] * gridSize[2]))

    nTotalBoxes = nLocalBoxes + nHaloBoxes

    return linkcell(nLocalBoxes, nTotalBoxes, gridSize, 1.0/boxSize) 

class linkcell:
    def __init__(self, nLocalBoxes, nTotalBoxes, gridSize, invBoxSize):
        self.MAXATOMS = 64
        self.nLocalBoxes = nLocalBoxes
        self.nTotalBoxes = nTotalBoxes
        self.nAtoms = np.zeros(nTotalBoxes, dtype=np.int)
        self.gridSize = gridSize
        self.invBoxSize = invBoxSize

    def getBoxFromTuple(self, ix, iy, iz):
        iBox = 0
        # Halo in Z+
        if iz == self.gridSize[2]:
            iBox = self.nLocalBoxes + 2*self.gridSize[2]*self.gridSize[1] + 2*self.gridSize[2]*(self.gridSize[0]+2)  \
                   + (self.gridSize[0]+2)*(self.gridSize[1]+2) + (self.gridSize[0] + 2)*(iy+1) + (ix+1)

        # Halo in Z-
        elif iz == -1:
            iBox = self.nLocalBoxes + 2*self.gridSize[2]*self.gridSize[1] + 2*self.gridSize[2]*(self.gridSize[0]+2) \
                    + (self.gridSize[0]+2)*(iy+1) + (ix+1)

        # Halo in Y+
        elif iy == self.gridSize[1]:
            iBox = self.nLocalBoxes + 2*self.gridSize[2]*self.gridSize[1] + self.gridSize[2]*(self.gridSize[0]+2) \
                    + (self.gridSize[0]+2)*iz + (ix+1)

        # Halo in Y-
        elif iy == -1:
            iBox = self.nLocalBoxes + 2*self.gridSize[2]*self.gridSize[1] + iz*(self.gridSize[0] + 2) + (ix+1)

        # Halo in X+
        elif ix == self.gridSize[0]:
            iBox = self.nLocalBoxes + self.gridSize[1]*self.gridSize[2] + iz*self.gridSize[1] + iy

        # Halo in X-
        elif ix == -1:
            iBox = self.nLocalBoxes + iz*self.gridSize[1] + iy

        # Local link cell
        else:
            iBox = ix + self.gridSize[0]*iy + self.gridSize[0]*self.gridSize[1]*iz

        if iBox >= self.nTotalBoxes:
            print('iBox too large',iBox,self.nTotalBoxes)
            print('ix,iy,iz',ix,iy,iz)
        assert iBox >=0
        assert iBox < self.nTotalBoxes
        return iBox

    def getBoxFromCoord(self, r):
        ix = int(math.floor(r[0]*self.invBoxSize[0]))
        iy = int(math.floor(r[1]*self.invBoxSize[1]))
        iz = int(math.floor(r[2]*self.invBoxSize[2]))
        return self.getBoxFromTuple(ix, iy, iz)

    def putAtomInBox(self, atoms, gid, iType, x, y, z, px=0.0, py=0.0, pz=0.0):
        iBox = self.getBoxFromCoord([x, y, z])
        iOff = iBox*self.MAXATOMS
        iOff += self.nAtoms[iBox]

        if iBox < self.nLocalBoxes:
            atoms.nLocal += 1

        self.nAtoms[iBox] += 1
        atoms.gid[iOff] = gid
        atoms.species[iOff] = iType
        atoms.r[iOff][0] = x
        atoms.r[iOff][1] = y
        atoms.r[iOff][2] = z
        atoms.p[iOff][0] = px
        atoms.p[iOff][1] = py
        atoms.p[iOff][2] = pz

    def getTuple(self, iBox):
        ix,iy,iz = 0,0,0
        # local box
        if iBox < self.nLocalBoxes:
            ix = iBox % self.gridSize[0]
            iBox //= self.gridSize[0]
            iy = iBox % self.gridSize[1]
            iz = iBox // self.gridSize[1]
        else:   # halo box
            ink = iBox - self.nLocalBoxes
            if ink < 2*self.gridSize[1]*self.gridSize[2]:
                if ink < self.gridSize[1]*self.gridSize[2]:
                    ix = 0
                else:
                    ink -= self.gridSize[1]*self.gridSize[2]
                    ix = self.gridSize[0] + 1
                iy = 1 + ink % self.gridSize[1]
                iz = 1 + ink // self.gridSize[1]

            elif ink < 2*self.gridSize[2]*(self.gridSize[1] + self.gridSize[0] + 2):
                ink -= 2 * self.gridSize[2] * self.gridSize[1]
                if ink < (self.gridSize[0] + 2)*self.gridSize[2]:
                    iy = 0
                else:
                    ink -= (self.gridSize[0] + 2)*self.gridSize[2]
                    iy = self.gridSize[1] + 1
                ix = ink % (self.gridSize[0] + 2)
                iz = 1 + ink // (self.gridSize[0] + 2)
            else:
                ink -= 2*self.gridSize[2]*(self.gridSize[1] + self.gridSize[0] + 2)
                if ink < (self.gridSize[0] + 2)*(self.gridSize[1] + 2):
                    iz = 0
                else:
                    ink -= (self.gridSize[0] + 2)*(self.gridSize[1] + 2)
                    iz = self.gridSize[2] + 1
                ix = ink % (self.gridSize[0] + 2)
                iy = ink // (self.gridSize[0] + 2)

            ix -= 1
            iy -= 1
            iz -= 1

        return ix,iy,iz

    def getNeighborBoxes(self, iBox, nbrBoxes):
        ix,iy,iz = self.getTuple(iBox)
        count = 0
        for i in range(ix-1, ix+2):
            for j in range(iy-1, iy+2):
                for k in range(iz-1, iz+2):
                    nbrBoxes[count] = self.getBoxFromTuple(i, j, k)
                    count += 1
        return count

    def copyAtom(self, atoms, iAtom, iBox, jAtom, jBox):
        iOff = self.MAXATOMS*iBox + iAtom
        jOff = self.MAXATOMS*jBox + jAtom
        atoms.gid[jOff] = atoms.gid[iOff]
        atoms.species[jOff] = atoms.species[iOff]
        atoms.r[jOff] = atoms.r[iOff]
        atoms.p[jOff] = atoms.p[iOff]
        atoms.f[jOff] = atoms.f[iOff]

    def moveAtom(self, atoms, iId, iBox, jBox):
        nj = self.nAtoms[jBox]
        self.copyAtom(atoms, iId, iBox, nj, jBox)
        self.nAtoms[jBox] += 1

        self.nAtoms[iBox] -= 1
        ni = self.nAtoms[iBox]
        if ni:
            self.copyAtom(atoms, ni, iBox, iId, iBox)
        if jBox > self.nLocalBoxes:
            atoms.nLocal -= 1

    def emptyHaloCells(self):
        for ii in range(self.nLocalBoxes, self.nTotalBoxes):
            self.nAtoms[ii] = 0

    def updateLinkCells(self, atoms):
        # empty halo cells
        self.emptyHaloCells()
        for iBox in range(self.nLocalBoxes):
            iOff = iBox * self.MAXATOMS
            ii = 0
            while ii < self.nAtoms[iBox]:
                jBox = self.getBoxFromCoord(atoms.r[iOff + ii])
                if jBox != iBox:
                    self.moveAtom(atoms, ii, iBox, jBox)
                else:
                    ii += 1




