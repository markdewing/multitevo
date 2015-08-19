
import constants
import parallel

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


def createFccLattice(nx, ny, nz, lat, atoms, boxes, domain):
    """Create atom positions on a face centered cubic (FCC) lattice
       with nx * ny * nz unit cells and lattice constance 'lat'
    """
    nb = 4 # number of atoms in this basis

    basis  = [ (0.25, 0.25, 0.25),
               (0.25, 0.75, 0.75),
               (0.75, 0.25, 0.75),
               (0.75, 0.75, 0.25)
             ]

    begin = [int(x) for x in np.floor(domain.localMin/lat)]
    end = [int(x) for x in np.ceil(domain.localMax/lat)]

    idx = 0
    for ix in range(begin[0], end[0]):
        for iy in range(begin[1], end[1]):
            for iz in range(begin[2], end[2]):
                for ib in range(nb):
                    rx = (ix+basis[ib][0])*lat
                    ry = (iy+basis[ib][1])*lat
                    rz = (iz+basis[ib][2])*lat
                    if rx < domain.localMin[0] or rx >= domain.localMax[0]:
                        continue
                    if ry < domain.localMin[1] or ry >= domain.localMax[1]:
                        continue
                    if rz < domain.localMin[2] or rz >= domain.localMax[2]:
                        continue

                    gid = ib + nb*(iz + nz*(iy + ny*ix))
                    boxes.putAtomInBox(atoms, gid, 0, rx, ry, rz)
                    idx += 1

    nlocal = np.zeros(1, dtype=np.int)
    nglobal = np.zeros(1, dtype=np.int)
    nlocal[0] = idx

    parallel.addParallel(nlocal, nglobal)

    atoms.nGlobal = nglobal[0]
    if atoms.nGlobal != nb*nx*ny*nz:
        print 'nGlobal = ',atoms.nGlobal
        print 'nb,nx,ny,nz,product',nb,nx,ny,nz,nb*nx*ny*nz
    assert atoms.nGlobal == nb*nx*ny*nz


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

    sim.setEKinetic(ke)
    sim.collectGlobalEnergy()

    return sim.getEKinetic()


def computeVcm(sim):
    vcm = np.zeros(3)
    vcmGlobal = np.zeros(3)
    totalMass = 0.0
    #for i in range(sim.atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        iOff = sim.boxes.MAXATOMS * iBox
        for ii in range(sim.boxes.nAtoms[iBox]):
            vcm += sim.atoms.p[iOff,:]
            iSpecies = sim.atoms.species[iOff]
            totalMass += sim.species[iSpecies].mass
            iOff += 1

    parallel.addParallel(vcm, vcmGlobal)

    vcm = vcmGlobal/totalMass
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
