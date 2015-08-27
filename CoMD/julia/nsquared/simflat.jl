

include("initAtoms.jl")
include("constants.jl")
include("ljforce.jl")

type Species
    mass
end

Species() = Species(1.0)

function initSpecies(pot)
    [Species(pot.mass)]
end


type SimFlat
    nSteps::Int
    printRate::Int
    dt::Float
    atoms::Atoms
    species::Array{Species}
    pot::LJ_Pot
    ePot::Float
    eKinetic::Float
end


function dumpPos(atoms)
    for i in 1:atoms.nAtoms
        @printf "%d:  %f %f %f\n" i atoms.r[1, i] atoms.r[2,i] atoms.r[3,i]
    end
end


function initSimulation(cmd)
    nSteps = cmd["nSteps"]
    dt = cmd["dt"]


    nx = cmd["nx"]
    ny = cmd["ny"]
    nz = cmd["nz"]

    temp = cmd["temp"]

    atoms = initAtoms(4*nx*ny*nz)
    pot = init_lj_pot()
    species = initSpecies(pot)

    latticeConstant = cmd["lat"]
    if latticeConstant < 0.0
        latticeConstant = pot.lat
    end

    setBounds(atoms, nx, ny, nz, latticeConstant)

    createFccLattice(nx, ny, nz, latticeConstant, atoms)

    setTemperature(atoms, species, temp)

    return  SimFlat(nSteps, 10, dt, atoms, species, pot, 0.0, 0.0)
end

function kineticEnergy(atoms, species)
    ke = 0.0
    for i in 1:atoms.nAtoms
        iSpecies = atoms.species[i]
        invMass = 0.5/species[iSpecies].mass
        ke += invMass * dot(atoms.p[:,i], atoms.p[:,i])
    end
    return ke
end

function computeVcm(atoms, species)
    vcm = zeros(Float, 3)
    totalMass = 0.0
    for i in 1:atoms.nAtoms
        vcm += atoms.p[:,i]
        iSpecies = atoms.species[i]
        totalMass += species[iSpecies].mass
    end
    vcm = vcm/totalMass
end

function setVcm(atoms, species, newVcm)
    oldVcm = computeVcm(atoms, species)
    shift = newVcm - oldVcm
    for i in 1:atoms.nAtoms
        iSpecies = atoms.species[i]
        mass = species[iSpecies].mass
        atoms.p[:,i] += mass*shift
    end
end

function setTemperature(atoms, species, temperature)
    for i in 1:atoms.nAtoms
        iSpecies = atoms.species[i]
        mass = species[iSpecies].mass
        sigma = sqrt(kB_eV * temperature/mass)
        atoms.p[:,i] = mass * sigma * randn(3)
    end

    if temperature == 0.0
        return
    end

    vZero = zeros(Float, 3)
    setVcm(atoms, species, vZero)
    ke = kineticEnergy(atoms, species)
    temp = (ke/atoms.nAtoms)/kB_eV/1.5
    scaleFactor = sqrt(temperature/temp)
    for i in 1:atoms.nAtoms
        atoms.p[:,i] *= scaleFactor
    end

end


function advanceVelocity(atoms, dt)
    for i in 1:atoms.nAtoms
        atoms.p[:,i] += dt*atoms.f[:,i]
    end
end

function advancePosition(sim, dt)
    for i in 1:sim.atoms.nAtoms
        iSpecies = sim.atoms.species[i]
        invMass = 1.0/sim.species[iSpecies].mass
        sim.atoms.r[:,i] += dt*sim.atoms.p[:,i] * invMass
    end
end

function timestep(sim, nSteps, dt)
    for ii in 1:nSteps
        advanceVelocity(sim.atoms, 0.5*dt)
        advancePosition(sim, dt)
        sim.pot.forceFunc(sim.pot, sim)
        advanceVelocity(sim.atoms, 0.5*dt)
    end
    sim.eKinetic = kineticEnergy(sim.atoms, sim.species)
end
