from __future__ import print_function, division

import ljforce
import initatoms
import linkcell
import halo
import numpy as np

class SimFlat:
    def __init__(self, nSteps, printRate, dt):
        self.nSteps = nSteps
        self.printRate = printRate
        self.dt = dt
        self.atoms = None
        self.pot = ljforce.lj_pot()
        self.species = initSpecies(self.pot)
        self.ePot = 0.0
        self.eKinetic = 0.0
        self.boxes = None
        self.atomHalo = None


class species:
    def __init__(self, mass):
        self.mass = mass

def initSpecies(pot):
    return [species(pot.mass)]


def redistributeAtoms(sim, atoms):
    sim.boxes.updateLinkCells(atoms)
    sim.atomHalo.haloExchange(sim)

def initSimulation(cmd):
    sim = SimFlat(cmd.nSteps, cmd.printRate, cmd.dt)

    latticeConstant = cmd.lat

    if cmd.lat < 0.0:
        latticeConstant = sim.pot.lat




    box_size = np.zeros(3)
    box_size[0] = cmd.nx*latticeConstant
    box_size[1] = cmd.ny*latticeConstant
    box_size[2] = cmd.nz*latticeConstant
    sim.boxes = linkcell.initLinkCells(sim, box_size)

    sim.atoms = initatoms.initAtoms(sim.boxes.nTotalBoxes * sim.boxes.MAXATOMS)

    sim.atoms.setBounds(cmd.nx, cmd.ny, cmd.nz, latticeConstant)

    initatoms.createFccLattice(cmd.nx, cmd.ny, cmd.nz, latticeConstant, sim.atoms, sim.boxes)

    initatoms.setTemperature(sim, cmd.temp)

    sim.atomHalo = halo.Halo(sim.boxes)

    redistributeAtoms(sim, sim.atoms)



    return sim
