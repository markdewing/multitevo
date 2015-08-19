
# Translation of CoMD to Python

import simflat
import initatoms
import time
import constants
import parallel

import argparse
import sys


def advanceVelocity(sim, atoms, dt):
    #for i in range(atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        for ii in range(sim.boxes.nAtoms[iBox]):
            iOff = sim.boxes.MAXATOMS * iBox + ii
            atoms.p[iOff,0] += dt*atoms.f[iOff,0]
            atoms.p[iOff,1] += dt*atoms.f[iOff,1]
            atoms.p[iOff,2] += dt*atoms.f[iOff,2]

def advancePosition(sim, atoms, dt):
    #for i in range(atoms.nAtoms):
    for iBox in range(sim.boxes.nLocalBoxes):
        for ii in range(sim.boxes.nAtoms[iBox]):
            iOff = sim.boxes.MAXATOMS * iBox + ii
            iSpecies = atoms.species[iOff]
            invMass = 1.0/sim.species[iSpecies].mass
            atoms.r[iOff,0] += dt*atoms.p[iOff,0]*invMass
            atoms.r[iOff,1] += dt*atoms.p[iOff,1]*invMass
            atoms.r[iOff,2] += dt*atoms.p[iOff,2]*invMass


#def redistributeAtoms(sim, atoms): #    sim.boxes.updateLinkCells(atoms) #    sim.atomHalo.haloExchange(sim)

def timestep(sim, nSteps, dt):
    for ii in range(nSteps):
        advanceVelocity(sim, sim.atoms, 0.5*dt)
        advancePosition(sim, sim.atoms, dt)
        simflat.redistributeAtoms(sim, sim.atoms)
        sim.pot.computeForce(sim.atoms, sim)
        advanceVelocity(sim, sim.atoms, 0.5*dt)
    #print 'ke,pe = ',initatoms.kineticEnergy(sim)/sim.atoms.nAtoms,sim.ePot/sim.atoms.nAtoms
    initatoms.kineticEnergy(sim)


def initValidate(sim):
    ke = initatoms.kineticEnergy(sim)
    pe = sim.getEPot()

    if parallel.printRank():
        print "Initial PE (cohesive energy) : ",pe/sim.atoms.nGlobal
        print "Initial energy : ",(ke+pe)/sim.atoms.nGlobal," atom count : ",sim.atoms.nGlobal

firstCall = True
iStepPrev = -1
def printInfo(sim, iStep, elapsedTime):
    global firstCall, iStepPrev
    if firstCall:
        firstCall = False
        #print "# Loop  Time(fs)   Total Energy  Potential Energy   Kinetic Energy  Temperature  Performance (us/atom)  # Atoms"
        if parallel.printRank():
            print "#{:>100s}".format("Performance")
            print "#{:>6s} {:>10s} {:>18s} {:>18s} {:>18s} {:>12s} {:^12s} {:^12s}".\
                format("Loop", "Time(fs)", "Total Energy", "Potential Energy", "Kinetic Energy", "Temperature",  "(us/atom)", "# Atoms")

    nEval = iStep - iStepPrev
    iStepPrev = iStep

    simtime = iStep * sim.dt
    n = sim.atoms.nGlobal
    eTotal = (sim.getEPot() + sim.getEKinetic())/n
    eK = sim.getEKinetic() / n
    eU = sim.getEPot() / n
    Temp = (sim.getEKinetic()/n)/(constants.kB_eV * 1.5)
    timePerAtom = 1.0e6 * elapsedTime/(n*nEval)

    if parallel.printRank():
        print " {:6d} {:10.2f} {:18.12f} {:18.12f} {:18.12f} {:12.4f} {:10.4f} {:12d}".format(iStep, simtime, eTotal, eU, eK, Temp, timePerAtom, n)
    #print iStep, simtime, eTotal, eU, eK, Temp, timePerAtom, n


def parseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x","--nx", type=int, default=20, help="number of unit cells in x")
    parser.add_argument("-y","--ny", type=int, default=20, help="number of unit cells in y")
    parser.add_argument("-z","--nz", type=int, default=20, help="number of unit cells in z")
    parser.add_argument("-i","--xproc", type=int, default=1, help="number of ranks in x direction")
    parser.add_argument("-j","--yproc", type=int, default=1, help="number of ranks in y direction")
    parser.add_argument("-k","--zproc", type=int, default=1, help="number of ranks in z direction")
    parser.add_argument("-N","--nSteps", type=int, default=100, help="total number of time steps")
    parser.add_argument("-n","--printRate", type=int, default=10, help="number of steps between output")
    parser.add_argument("-D","--dt", type=float, default=1, help="time step (in fs)")
    parser.add_argument("-l","--lat", type=float, default=-1, help="lattice parameter (Angstroms)")
    parser.add_argument("-T","--temp", type=float, default=600, help="initial temperature (K)")

    args = parser.parse_args()
    return args

def sanityCheck(cmd):
    nProcs = cmd.xproc * cmd.yproc * cmd.zproc
    if nProcs != parallel.NRanks:
        if parallel.printRank():
            print "Number of MPI ranks (%d) must match xproc * yproc * zproc (%d)"%(parallel.NRanks, nProcs)
            return False

    return True


def run_comd():
    parallel.initParallel()
    cmd = parseCommandLine()

    if not sanityCheck(cmd):
        sys.exit(1)
    # init subsystems
    sim = simflat.initSimulation(cmd)
    #print 'nAtoms = ',sim.atoms.nGlobal
    #for i in range(min(sim.atoms.nAtoms, 10)):
    #for i in range(sim.atoms.nAtoms):
    #    print i,sim.atoms.r[i,:]

    sim.pot.computeForce(sim.atoms, sim)

    initatoms.kineticEnergy(sim)

    initValidate(sim)

    for iStep in range(0,sim.nSteps+1,sim.printRate):
        start = time.clock()
        timestep(sim, sim.printRate, sim.dt)
        end = time.clock()

        printInfo(sim, iStep, end-start)

    # validate result

if __name__ == '__main__':
    run_comd()

