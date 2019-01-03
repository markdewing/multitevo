

include("initAtoms.jl")
include("simflat.jl")
include("ljforce.jl")

using Printf
using ArgParse

# Translation of CoMD to Julia

function initValidate(sim)
    ke = sim.eKinetic
    pe = sim.ePot

    natoms = sim.atoms.nAtoms
    @printf "Initial PE (cohesive energy) : %f\n" pe/natoms
    @printf "Initial energy : %f   atom count : %d\n" (ke+pe)/natoms natoms
end

firstCall = true
iStepPrev = -1
function printInfo(sim, iStep, elapsedTime)
    global firstCall, iStepPrev

    if firstCall
        firstCall = false
        @printf "%100s\n" "Perfomance"
        @printf "#%6s %10s %18s %18s %18s %12s %12s %10s\n" "Loop"  "Time(fs)"  "Total Energy"  "Potential Energy"  "Kinetic Energy" "Temperature" "(us/atom)"  "# Atoms"
    end

    nEval = iStep - iStepPrev
    iStepPrev = iStep

    simtime = iStep * sim.dt
    n = sim.atoms.nAtoms
    eTotal = (sim.ePot + sim.eKinetic)/n
    eK = sim.eKinetic/n
    eU = sim.ePot /n
    Temp = (sim.eKinetic/n)/(kB_eV * 1.5)
    timePerAtom = 1.0e6 * elapsedTime/(n*nEval)

    @printf "%6d  %10.2f %18.12f %18.12f %18.12f %12.4f %10.4f %10d\n" iStep simtime eTotal eU eK Temp timePerAtom n

end

function printPerformanceResults(sim, elapsedTime)
    n = sim.atoms.nAtoms
    nEval = sim.nSteps
    timePerAtom = 1.0e6 * elapsedTime/(n*nEval)
    @printf "\nAverage all atom update rate: %10.f us/atom\n" timePerAtom
end

function run_comd()
    cmd = parseCommandLine()
    sim = initSimulation(cmd)

    sim.pot.forceFunc(sim.pot, sim)
    sim.eKinetic = kineticEnergy(sim.atoms, sim.species)

    initValidate(sim)

    timestepTimeOneIteration = 0.0
    timestepTime= 0.0
    iStep = 0
    for jStep in 1:sim.printRate:sim.nSteps
        printInfo(sim, iStep, timestepTimeOneIteration)
        start = time_ns()
        timestep(sim, sim.printRate, sim.dt)
        stop = time_ns()
        timestepTimeOneIteration = (stop - start)*1e-9
        timestepTime += timestepTimeOneIteration
        iStep += sim.printRate
    end

    printInfo(sim, iStep, timestepTimeOneIteration)

    printPerformanceResults(sim, timestepTime)
end

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nx", "-x"
          help = "number of unit cells in x"
          arg_type = Int
          default = 4
        "--ny", "-y"
          help = "number of unit cells in y"
          arg_type = Int
          default = 4
        "--nz", "-z"
          help = "number of unit cells in z"
          arg_type = Int
          default = 4
        "--nSteps", "-N"
          help = "total number of time steps"
          arg_type = Int
          default = 100
        "--printRate", "-n"
          help = "number of steps between output"
          arg_type = Int
          default = 100
        "--dt", "-D"
          help = "time step (in fs)"
          arg_type = Float
          default = 1.0
        "--lat", "-l"
          help = "lattice parameter (Angstroms)"
          arg_type = Float
          default = -1.0
        "--temp", "-T"
          help = "initial temperature (K)"
          arg_type = Float
          default = 600.0
    end
    return parse_args(s)
end

#function parseCommandLine()
#    cmd = Dict()
#    cmd["nx"] = 10
#    cmd["ny"] = 10
#    cmd["nz"] = 10
#    cmd["nSteps"] = 100
#    cmd["printRate"] = 10
#    cmd["dt"] = 1.0
#    cmd["lat"] = -1.0
#    cmd["temp"] = 600.0
#
#    return cmd
#end

run_comd()
