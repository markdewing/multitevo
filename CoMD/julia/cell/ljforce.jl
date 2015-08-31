
include("constants.jl")

type LJ_Pot
    sigma::Float
    epsilon::Float
    mass::Float
    lat::Float
    lattice_type::String
    cutoff::Float
    name::String
    atomic_no::Int
    forceFunc
end

function init_lj_pot()
    sigma = 2.315
    epsilon = 0.167
    mass = 63.55 * amuToInternalMass
    lat = 3.615
    lattice_type = "FCC"
    cutoff = 2.5*sigma
    name = "Cu"
    atomic_no = 29
    return LJ_Pot(sigma, epsilon, mass, lat, lattice_type, cutoff, name, atomic_no, computeForce)
end


function computeForce(pot, sim)
    POT_SHIFT = 1.0

    atoms = sim.atoms

    rCut2 = pot.cutoff * pot.cutoff

    for iBox in 1:sim.boxes.nLocalBoxes
        iOff = sim.boxes.MAXATOMS * (iBox-1) + 1
        for ii in 1:sim.boxes.nAtoms[iBox]
            atoms.f[:, iOff] = 0.0
            iOff += 1
        end
    end

    ePot = 0.0

    s6 = pot.sigma *  pot.sigma *  pot.sigma *   pot.sigma *  pot.sigma  * pot.sigma

    rCut6 = s6/(rCut2 * rCut2 * rCut2)

    eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0)

    nbrBoxes = zeros(Int, 27)

    for iBox in 1:sim.boxes.nLocalBoxes
        nIBox = sim.boxes.nAtoms[iBox]
        nNbrBoxes = getNeighborBoxes(sim.boxes, iBox, nbrBoxes)
        for jTmp in 1:nNbrBoxes
            jBox = nbrBoxes[jTmp]
            nJBox = sim.boxes.nAtoms[jBox]
            if nJBox == 0
                continue
            end

            for ii in 1:nIBox
                iOff = sim.boxes.MAXATOMS * (iBox-1) + ii
                iId = atoms.gid[iOff]
                for jj in 1:nJBox
                    jOff = sim.boxes.MAXATOMS * (jBox-1) + jj
                    jId = atoms.gid[jOff]
                    if (jBox <= sim.boxes.nLocalBoxes) && (jId <= iId)
                        continue
                    end


                    dd = atoms.r[:, iOff] - atoms.r[:,jOff]
                    r2 = dot(dd, dd)


                    if r2 > rCut2
                        continue
                    end

                    ir2 = 1.0/r2

                    r6 = s6 * (ir2*ir2*ir2)

                    eLocal = r6*(r6 - 1.0) - eShift

                    if jBox <= sim.boxes.nLocalBoxes
                        ePot += eLocal
                    else
                        ePot += 0.5*eLocal
                    end

                    fr = -4.0 * pot.epsilon * r6 * ir2*(12.0*r6 - 6.0)

                    atoms.f[:,iOff] -= dd*fr
                    atoms.f[:,jOff] += dd*fr
                end
            end
        end
    end

    ePot = ePot*4.0*pot.epsilon
    sim.ePot = ePot
end
