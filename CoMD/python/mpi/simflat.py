
import domain
import ljforce
import initatoms
import linkcell
import halo
import parallel

import numpy as np

class SimFlat:
    def __init__(self, nSteps, printRate, dt):
        self.nSteps = nSteps
        self.printRate = printRate
        self.dt = dt
        self.atoms = None
        self.pot = ljforce.lj_pot()
        self.species = initSpecies(self.pot)
        self.ePotLocal = 0.0
        self.eKineticLocal = 0.0
        self._ePot = 0.0
        self._eKinetic = 0.0
        self.globalEnergyGood = False
        self.boxes = None
        self.atomHalo = None
        self.domain = None

    def collectGlobalEnergy(self):
        eLocal = np.zeros(2)
        eSum = np.zeros(2)

        eLocal[0] = self.eKineticLocal
        eLocal[1] = self.ePotLocal

        parallel.addParallel(eLocal, eSum)

        self._eKinetic = eSum[0]
        self._ePot = eSum[1]

        self.globalEnergyGood = True

    def getEPot(self):
        if not self.globalEnergyGood:
            self.collectGlobalEnergy()
        return self._ePot

    def getEKinetic(self):
        if not self.globalEnergyGood:
            self.collectGlobalEnergy()
        return self._eKinetic

    def setEPot(self, epot):
        self.ePotLocal = epot
        self.globalEnergyGood = False

    def setEKinetic(self, ekin):
        self.eKineticLocal = ekin
        self.globalEnergyGood = False


class species:
    def __init__(self, mass):
        self.mass = mass

def initSpecies(pot):
    return [species(pot.mass)]


def redistributeAtoms(sim, atoms):
    sim.boxes.updateLinkCells(atoms)
    sim.atomHalo.haloExchange(sim)
    for ii in range(sim.boxes.nTotalBoxes):
        sim.atomHalo.sortAtomsInCell(sim.atoms, sim.boxes, ii)

def initSimulation(cmd):
    sim = SimFlat(cmd.nSteps, cmd.printRate, cmd.dt)

    latticeConstant = cmd.lat

    if cmd.lat < 0.0:
        latticeConstant = sim.pot.lat


    box_size = np.zeros(3)
    box_size[0] = cmd.nx*latticeConstant
    box_size[1] = cmd.ny*latticeConstant
    box_size[2] = cmd.nz*latticeConstant

    sim.domain = domain.Domain(cmd.xproc, cmd.yproc, cmd.zproc, box_size)

    sim.boxes = linkcell.initLinkCells(sim, sim.domain)

    sim.atoms = initatoms.initAtoms(sim.boxes.nTotalBoxes * sim.boxes.MAXATOMS)

    sim.atoms.setBounds(cmd.nx, cmd.ny, cmd.nz, latticeConstant)

    initatoms.createFccLattice(cmd.nx, cmd.ny, cmd.nz, latticeConstant, sim.atoms, sim.boxes, sim.domain)

    initatoms.setTemperature(sim, cmd.temp)

    sim.atomHalo = halo.Halo(sim.boxes, sim.domain)

    redistributeAtoms(sim, sim.atoms)



    return sim
