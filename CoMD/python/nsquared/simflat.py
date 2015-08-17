
import ljforce
import initatoms

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

class species:
    def __init__(self, mass):
        self.mass = mass

def initSpecies(pot):
    return [species(pot.mass)]

def initSimulation(cmd):
    sim = SimFlat(cmd.nSteps, cmd.printRate, cmd.dt)

    latticeConstant = cmd.lat

    if cmd.lat < 0.0:
        latticeConstant = sim.pot.lat

    sim.atoms = initatoms.initAtoms(cmd.nx, cmd.ny, cmd.nz)

    sim.atoms.setBounds(cmd.nx, cmd.ny, cmd.nz, latticeConstant)

    initatoms.createFccLattice(cmd.nx, cmd.ny, cmd.nz, latticeConstant, sim.atoms)

    initatoms.setTemperature(sim, cmd.temp)

    return sim
