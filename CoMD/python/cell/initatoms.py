from __future__ import print_function, division

import constants
import numpy as np
import math
import random

class atoms:
    def __init__(self, n):
        self.r = np.zeros( (n,3) )
        self.f = np.zeros( (n,3) )
        self.p = np.zeros( (n,3) )
        self.species = np.zeros(n, dtype=np.int)
        self.gid = np.zeros(n, dtype=np.int)
        self.bounds = np.zeros(3)
        self.maxTotalAtoms = n
        self.nLocal = 0
        self.nGlobal = 0

    def setBounds(self, nx, ny, nz, lat):
        self.bounds[0] = nx*lat
        self.bounds[1] = ny*lat
        self.bounds[2] = nz*lat

def initAtoms(maxatoms):
    return atoms(maxatoms)


def createFccLattice(nx, ny, nz, lat, atoms, boxes):
    """Create atom positions on a face centered cubic (FCC) lattice
       with nx * ny * nz unit cells and lattice constance 'lat'
    """
    nb = 4 # number of atoms in this basis

    basis  = [ (0.25, 0.25, 0.25),
               (0.25, 0.75, 0.75),
               (0.75, 0.25, 0.75),
               (0.75, 0.75, 0.25)
             ]

    idx = 0
    # loop over ix,iy,iz
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for ib in range(nb):
                    rx = (ix+basis[ib][0])*lat
                    ry = (iy+basis[ib][1])*lat
                    rz = (iz+basis[ib][2])*lat

                    #atoms.r[idx,0] = rx
                    #atoms.r[idx,1] = ry
                    #atoms.r[idx,2] = rz
                    gid = ib + nb*(iz + nz*(iy + ny*ix))
                    boxes.putAtomInBox(atoms, gid, 0, rx, ry, rz)
                    idx += 1

    atoms.nGlobal = idx
    assert idx == nb*nx*ny*nz
    #idx should equal nx*ny*nz*nb

def kineticEnergy(sim):
    ke = 0.0
    #for i in range(sim.atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        iOff = sim.boxes.MAXATOMS * iBox
        for ii in range(sim.boxes.nAtoms[iBox]):
            iSpecies = sim.atoms.species[iOff]
            invMass = 0.5/sim.species[iSpecies].mass
            ke += invMass * (sim.atoms.p[iOff,0]*sim.atoms.p[iOff,0] +  \
                    sim.atoms.p[iOff,1]*sim.atoms.p[iOff,1] + sim.atoms.p[iOff,2]*sim.atoms.p[iOff,2])
            iOff += 1
    return ke

def computeVcm(sim):
    vcm = np.zeros(3)
    totalMass = 0.0
    #for i in range(sim.atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        iOff = sim.boxes.MAXATOMS * iBox
        for ii in range(sim.boxes.nAtoms[iBox]):
            vcm += sim.atoms.p[iOff,:]
            iSpecies = sim.atoms.species[iOff]
            totalMass += sim.species[iSpecies].mass
            iOff += 1
    vcm = vcm/totalMass
    return vcm

def setVcm(sim, newVcm):
    oldVcm = computeVcm(sim)
    shift = newVcm - oldVcm
    #for i in range(sim.atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        iOff = sim.boxes.MAXATOMS * iBox
        for ii in range(sim.boxes.nAtoms[iBox]):
            iSpecies = sim.atoms.species[iOff]
            mass = sim.species[iSpecies].mass
            sim.atoms.p[iOff,:] += mass*shift
            iOff += 1


def setTemperature(sim, temperature):
    #for i in range(sim.atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        iOff = sim.boxes.MAXATOMS * iBox
        for ii in range(sim.boxes.nAtoms[iBox]):
            iSpecies = sim.atoms.species[iOff]
            mass = sim.species[iSpecies].mass
            sigma = math.sqrt(constants.kB_eV * temperature/mass)
            sim.atoms.p[iOff,0] = mass * sigma * random.gauss(1.0, 1.0)
            sim.atoms.p[iOff,1] = mass * sigma * random.gauss(1.0, 1.0)
            sim.atoms.p[iOff,2] = mass * sigma * random.gauss(1.0, 1.0)
            iOff += 1

    if temperature == 0.0:
        return

    vZero = np.zeros(3)
    setVcm(sim, vZero)
    ke = kineticEnergy(sim)
    temp = (ke/sim.atoms.nGlobal)/constants.kB_eV/1.5
    scaleFactor = math.sqrt(temperature/temp)

    #for i in range(sim.atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        iOff = sim.boxes.MAXATOMS * iBox
        for ii in range(sim.boxes.nAtoms[iBox]):
            sim.atoms.p[iOff,:] *= scaleFactor
            iOff += 1
