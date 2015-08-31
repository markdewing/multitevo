

include("initAtoms.jl")
include("constants.jl")
include("ljforce.jl")
include("linkcell.jl")
include("halo.jl")

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
    boxes::LinkCell
    halo::HaloExchange
end


function dumpPos(atoms, boxes)
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            @printf "%d:  %20.10f %20.10f %20.10f\n" iOff atoms.r[1, iOff] atoms.r[2,iOff] atoms.r[3,iOff]
            iOff += 1
        end
    end
end

function dumpVel(atoms, boxes)
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            @printf "%d:  %20.10f %20.10f %20.10f\n" iOff atoms.p[1, iOff] atoms.p[2,iOff] atoms.p[3,iOff]
            iOff += 1
        end
    end
end

function dumpForce(atoms, boxes)
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            @printf "%d:  %20.10f %20.10f %20.10f\n" iOff atoms.f[1, iOff] atoms.f[2,iOff] atoms.f[3,iOff]
            iOff += 1
        end
    end
end


function initSimulation(cmd)
    nSteps = cmd["nSteps"]
    dt = cmd["dt"]
    printRate = cmd["printRate"]


    nx = cmd["nx"]
    ny = cmd["ny"]
    nz = cmd["nz"]

    temp = cmd["temp"]

    pot = init_lj_pot()
    species = initSpecies(pot)

    latticeConstant = cmd["lat"]
    if latticeConstant < 0.0
        latticeConstant = pot.lat
    end

    box_size = [nx, ny, nz]*latticeConstant
    boxes = initLinkCells(pot, box_size)

    halo = initHalo(boxes)

    atoms = initAtoms(boxes.nTotalBoxes  * boxes.MAXATOMS)
    setBounds(atoms, nx, ny, nz, latticeConstant)

    createFccLattice(nx, ny, nz, latticeConstant, atoms, boxes)

    setTemperature(atoms, boxes, species, temp)

    redistributeAtoms(atoms, boxes, halo)

    return  SimFlat(nSteps, printRate, dt, atoms, species, pot, 0.0, 0.0, boxes, halo)
end

function kineticEnergy(atoms, boxes, species)
    ke = 0.0
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            iSpecies = atoms.species[iOff]
            invMass = 0.5/species[iSpecies].mass
            ke += invMass * dot(atoms.p[:,iOff], atoms.p[:,iOff])
            iOff += 1
        end
    end
    return ke
end

function computeVcm(atoms, boxes, species)
    vcm = zeros(Float, 3)
    totalMass = 0.0
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            vcm += atoms.p[:,iOff]
            iSpecies = atoms.species[iOff]
            totalMass += species[iSpecies].mass
            iOff += 1
        end
    end
    vcm = vcm/totalMass
end

function setVcm(atoms, boxes, species, newVcm)
    oldVcm = computeVcm(atoms, boxes, species)
    shift = newVcm - oldVcm
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            iSpecies = atoms.species[iOff]
            mass = species[iSpecies].mass
            atoms.p[:,iOff] += mass*shift
            iOff += 1
        end
    end
end

function setTemperature(atoms, boxes, species, temperature)
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            iSpecies = atoms.species[iOff]
            mass = species[iSpecies].mass
            sigma = sqrt(kB_eV * temperature/mass)
            atoms.p[:,iOff] = mass * sigma * randn(3)
            iOff += 1
        end
    end

    if temperature == 0.0
        return
    end

    vZero = zeros(Float, 3)
    setVcm(atoms, boxes, species, vZero)
    ke = kineticEnergy(atoms, boxes, species)
    temp = (ke/atoms.nGlobal)/kB_eV/1.5
    scaleFactor = sqrt(temperature/temp)
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            atoms.p[:,iOff] *= scaleFactor
            iOff += 1
        end
    end

end


function advanceVelocity(atoms, boxes, dt)
    for iBox in 1:boxes.nLocalBoxes
        iOff = boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:boxes.nAtoms[iBox]
            atoms.p[:,iOff] += dt*atoms.f[:,iOff]
            iOff += 1
        end
    end
end

function advancePosition(sim, dt)
    for iBox in 1:sim.boxes.nLocalBoxes
        iOff = sim.boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:sim.boxes.nAtoms[iBox]
            iSpecies = sim.atoms.species[iOff]
            invMass = 1.0/sim.species[iSpecies].mass
            sim.atoms.r[:,iOff] += dt*sim.atoms.p[:,iOff] * invMass
            iOff += 1
        end
    end
end

function timestep(sim, nSteps, dt)
    for ii in 1:nSteps
        advanceVelocity(sim.atoms, sim.boxes, 0.5*dt)
        advancePosition(sim, dt)
        redistributeAtoms(sim.atoms, sim.boxes, sim.halo)
        sim.pot.forceFunc(sim.pot, sim)
        advanceVelocity(sim.atoms, sim.boxes, 0.5*dt)
    end
    sim.eKinetic = kineticEnergy(sim.atoms, sim.boxes, sim.species)
end

function redistributeAtoms(atoms, boxes, halo)
    updateLinkCells(atoms, boxes)
    haloExchange(halo, atoms, boxes)

    natom = 0
    for iBox in 1:boxes.nLocalBoxes
        natom += boxes.nAtoms[iBox]
    end
    if natom != atoms.nGlobal
        println("Number of atoms incorrect = $natom, should be $(atoms.nGlobal)")
    end

end
