from __future__ import print_function, division

import numpy as np
import constants

class LJ_Pot(object):
    def __init__(self):
        self.sigma = 2.315
        self.epsilon = 0.167
        self.mass = 63.55 * constants.amuToInternalMass
        self.lat = 3.615
        self.lattice_type = 'FCC'
        self.cutoff = 2.5*self.sigma
        self.name = "Cu"
        self.atomic_no = 29

    def output(self):
        pass

    def computeForce(self, atoms, sim):
        POT_SHIFT = 1.0
        #POT_SHIFT = 0.0
        rCut2 = self.cutoff * self.cutoff

        #for i in range(atoms.nAtoms):
        for iBox in range(sim.boxes.nLocalBoxes):
            iOff = sim.boxes.MAXATOMS * iBox
            for ii in range(sim.boxes.nAtoms[iBox]):
                atoms.f[iOff,:] = 0.0
                iOff += 1

        ePot = 0.0

        s6 = self.sigma * self.sigma * self.sigma * self.sigma * self.sigma * self.sigma

        rCut6 = s6/(rCut2*rCut2*rCut2)
        eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0)

        nbrBoxes = np.zeros(27, dtype=np.int)
        # loop over atoms

        iPot = 0
        #for i in range(atoms.nAtoms):
        #    for j in range(atoms.nAtoms):
        for iBox in range(sim.boxes.nLocalBoxes):
            nIBox = sim.boxes.nAtoms[iBox]
            nNbrBoxes = sim.boxes.getNeighborBoxes(iBox, nbrBoxes)
            for jTmp in range(nNbrBoxes):
                jBox = nbrBoxes[jTmp]
                assert jBox >= 0

                nJBox = sim.boxes.nAtoms[jBox]
                if nJBox == 0:
                    continue

                for ii in range(nIBox):
                    iOff = sim.boxes.MAXATOMS * iBox + ii
                    iId = atoms.gid[iOff]
                    for jj in range(nJBox):
                        jOff = sim.boxes.MAXATOMS * jBox + jj
                        jId = atoms.gid[jOff]
                        if jBox < sim.boxes.nLocalBoxes and jId <= iId:
                            continue


                        r2 = 0.0
                        dx = atoms.r[iOff,0] - atoms.r[jOff,0]
                        r2 += dx*dx

                        dy = atoms.r[iOff,1] - atoms.r[jOff,1]
                        r2 += dy*dy

                        dz = atoms.r[iOff,2] - atoms.r[jOff,2]
                        r2 += dz*dz

                        if r2 > rCut2:
                            continue

                        ir2 = 1.0/r2

                        r6 = s6 * (ir2*ir2*ir2)
                        eLocal = r6*(r6-1.0) - eShift

                        if jBox < sim.boxes.nLocalBoxes:
                            ePot += eLocal
                        else:
                            ePot += 0.5 * eLocal

                        iPot += 1


                        fr = -4.0*self.epsilon * r6 * ir2*(12.0*r6 - 6.0)

                        atoms.f[iOff,0] -= dx*fr
                        atoms.f[jOff,0] += dx*fr
                        atoms.f[iOff,1] -= dy*fr
                        atoms.f[jOff,1] += dy*fr
                        atoms.f[iOff,2] -= dz*fr
                        atoms.f[jOff,2] += dz*fr

        ePot = ePot*4.0*self.epsilon
        return ePot
