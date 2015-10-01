from __future__ import print_function, division

import simflat
import initatoms
import time
import constants

import argparse

# Translation of CoMD to Python


def advanceVelocity(atoms, dt):
    for i in range(atoms.nAtoms):
        atoms.p[i,0] += dt*atoms.f[i,0]
        atoms.p[i,1] += dt*atoms.f[i,1]
        atoms.p[i,2] += dt*atoms.f[i,2]

def advancePosition(sim, atoms, dt):
    for i in range(atoms.nAtoms):
        iSpecies = atoms.species[i]
        invMass = 1.0/sim.species[iSpecies].mass
        atoms.r[i,0] += dt*atoms.p[i,0]*invMass
        atoms.r[i,1] += dt*atoms.p[i,1]*invMass
        atoms.r[i,2] += dt*atoms.p[i,2]*invMass


def timestep(sim, nSteps, dt):
    for ii in range(nSteps):
        advanceVelocity(sim.atoms, 0.5*dt)
        advancePosition(sim, sim.atoms, dt)
        # redistributeAtoms
        sim.ePot = sim.pot.computeForce(sim.atoms, sim)
        advanceVelocity(sim.atoms, 0.5*dt)
    #print('ke,pe = ',initatoms.kineticEnergy(sim)/sim.atoms.nAtoms,sim.ePot/sim.atoms.nAtoms)
    sim.eKinetic = initatoms.kineticEnergy(sim)


def initValidate(sim):
    ke = initatoms.kineticEnergy(sim)
    pe = sim.ePot

    print("Initial PE (cohesive energy) : ",pe/sim.atoms.nAtoms)
    print("Initial energy : ",(ke+pe)/sim.atoms.nAtoms," atom count : ",sim.atoms.nAtoms)

firstCall = True
iStepPrev = -1
def printInfo(sim, iStep, elapsedTime):
    global firstCall, iStepPrev
    if firstCall:
        firstCall = False
        #print("# Loop  Time(fs)   Total Energy  Potential Energy   Kinetic Energy  Temperature  Performance (us/atom)  # Atoms")
        print("#{:>100s}".format("Performance"))
        print("#{:>6s} {:>10s} {:>18s} {:>18s} {:>18s} {:>12s} {:^12s} {:^12s}".\
            format("Loop", "Time(fs)", "Total Energy", "Potential Energy", "Kinetic Energy", "Temperature",  "(us/atom)", "# Atoms"))

    nEval = iStep - iStepPrev
    iStepPrev = iStep

    simtime = iStep * sim.dt
    n = sim.atoms.nAtoms
    eTotal = (sim.ePot + sim.eKinetic)/n
    eK = sim.eKinetic / n
    eU = sim.ePot / n
    Temp = (sim.eKinetic/n)/(constants.kB_eV * 1.5)
    timePerAtom = 1.0e6 * elapsedTime/(n*nEval)


    print(" {:6d} {:10.2f} {:18.12f} {:18.12f} {:18.12f} {:12.4f} {:10.4f} {:12d}".format(iStep, simtime, eTotal, eU, eK, Temp, timePerAtom, n))
    #print(iStep, simtime, eTotal, eU, eK, Temp, timePerAtom, n)

def printPerformanceResults(sim, elapsedTime):
    n = sim.atoms.nAtoms
    nEval = sim.nSteps
    timePerAtom = 1.0e6 * elapsedTime/(n*nEval)
    print()
    print("Average all atom update rate: %10.2f us/atom" % timePerAtom)
    print()


def parseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x","--nx", type=int, default=4, help="number of unit cells in x")
    parser.add_argument("-y","--ny", type=int, default=4, help="number of unit cells in y")
    parser.add_argument("-z","--nz", type=int, default=4, help="number of unit cells in z")
    parser.add_argument("-N","--nSteps", type=int, default=100, help="total number of time steps")
    parser.add_argument("-n","--printRate", type=int, default=10, help="number of steps between output")
    parser.add_argument("-D","--dt", type=float, default=1, help="time step (in fs)")
    parser.add_argument("-l","--lat", type=float, default=-1, help="lattice parameter (Angstroms)")
    parser.add_argument("-T","--temp", type=float, default=600, help="initial temperature (K)")

    args = parser.parse_args()
    return args


def run_comd():
    cmd = parseCommandLine()
    # init subsystems
    sim = simflat.initSimulation(cmd)
    print('nAtoms = ',sim.atoms.nAtoms)
    #for i in range(min(sim.atoms.nAtoms, 10)):
    #for i in range(sim.atoms.nAtoms):
    #    print i,sim.atoms.r[i,:]

    sim.ePot = sim.pot.computeForce(sim.atoms, sim)
    sim.eKinetic = initatoms.kineticEnergy(sim)

    initValidate(sim)


    timestepTime = 0.0
    timestepTimeOneIteration = 0.0

    iStep = 0
    for jStep in range(0, sim.nSteps, sim.printRate):
        printInfo(sim, iStep, timestepTimeOneIteration)

        start = time.clock()
        timestep(sim, sim.printRate, sim.dt)
        end = time.clock()

        timestepTimeOneIteration = end - start
        timestepTime += timestepTimeOneIteration

        iStep += sim.printRate


    printInfo(sim, iStep, timestepTimeOneIteration)

    printPerformanceResults(sim, timestepTime)

    # validate result

if __name__ == '__main__':
    run_comd()

