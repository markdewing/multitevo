
import constants
import numpy as np
import math
import random

ftype = np.float

class atoms:
    def __init__(self, n):
        self.r = np.zeros( (n,3), dtype=ftype )
        self.f = np.zeros( (n,3), dtype=ftype )
        self.p = np.zeros( (n,3), dtype=ftype )
        self.species = np.zeros(n, dtype=np.int)
        self.nAtoms = n
        self.bounds = np.zeros(3, dtype=ftype)

    def setBounds(self, nx, ny, nz, lat):
        self.bounds[0] = nx*lat
        self.bounds[1] = ny*lat
        self.bounds[2] = nz*lat

def initAtoms(nx, ny, nz):
    natom = 4*nx*ny*nz
    return atoms(natom)


def createFccLattice(nx, ny, nz, lat, atoms):
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

                    atoms.r[idx,0] = rx
                    atoms.r[idx,1] = ry
                    atoms.r[idx,2] = rz
                    idx += 1

    #idx should equal nx*ny*nz*nb

def kineticEnergy(sim):
    ke = 0.0
    for i in range(sim.atoms.nAtoms):
        iSpecies = sim.atoms.species[i]
        invMass = 0.5/sim.species[iSpecies].mass
        ke += invMass * (sim.atoms.p[i,0]*sim.atoms.p[i,0] +  \
                    sim.atoms.p[i,1]*sim.atoms.p[i,1] + sim.atoms.p[i,2]*sim.atoms.p[i,2])
    return ke

def computeVcm(sim):
    vcm = np.zeros(3, dtype=ftype)
    totalMass = 0.0
    for i in range(sim.atoms.nAtoms):
        vcm += sim.atoms.p[i,:]
        iSpecies = sim.atoms.species[i]
        totalMass += sim.species[iSpecies].mass
    vcm = vcm/totalMass
    return vcm

def setVcm(sim, newVcm):
    oldVcm = computeVcm(sim)
    shift = newVcm - oldVcm
    for i in range(sim.atoms.nAtoms):
        iSpecies = sim.atoms.species[i]
        mass = sim.species[iSpecies].mass
        sim.atoms.p[i,:] += mass*shift


def setTemperature(sim, temperature):
    for i in range(sim.atoms.nAtoms):
        iSpecies = sim.atoms.species[i]
        mass = sim.species[iSpecies].mass
        sigma = math.sqrt(constants.kB_eV * temperature/mass)
        sim.atoms.p[i,0] = mass * sigma * random.gauss(1.0, 1.0)
        sim.atoms.p[i,1] = mass * sigma * random.gauss(1.0, 1.0)
        sim.atoms.p[i,2] = mass * sigma * random.gauss(1.0, 1.0)

    if temperature == 0.0:
        return

    vZero = np.zeros(3, dtype=ftype)
    setVcm(sim, vZero)
    ke = kineticEnergy(sim)
    temp = (ke/sim.atoms.nAtoms)/constants.kB_eV/1.5
    scaleFactor = math.sqrt(temperature/temp)
    for i in range(sim.atoms.nAtoms):
        sim.atoms.p[i,:] *= scaleFactor
