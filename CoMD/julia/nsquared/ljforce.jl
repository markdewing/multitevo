
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

    for i in 1:atoms.nAtoms
        atoms.f[:,i] = 0.0
    end

    ePot = 0.0

    s6 = pot.sigma *  pot.sigma *  pot.sigma *   pot.sigma *  pot.sigma  * pot.sigma

    rCut6 = s6/(rCut2 * rCut2 * rCut2)

    eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0)

    for i in 1:atoms.nAtoms
        for j in 1:atoms.nAtoms
            if i < j
                dd = atoms.r[:,i] - atoms.r[:,j]
                for k in 1:3
                    if dd[k] < -0.5*atoms.bounds[k]
                        dd[k] += atoms.bounds[k]
                    end
                    if dd[k] > 0.5*atoms.bounds[k]
                        dd[k] -= atoms.bounds[k]
                    end
                end
                r2 = dot(dd, dd)

                if r2 > rCut2
                    continue
                end

                ir2 = 1.0/r2

                r6 = s6 * (ir2*ir2*ir2)

                eLocal = r6*(r6 - 1.0) - eShift

                ePot += eLocal

                fr = -4.0 * pot.epsilon * r6 * ir2*(12.0*r6 - 6.0)

                atoms.f[:,i] -= dd*fr
                atoms.f[:,j] += dd*fr
            end
        end
    end

    ePot = ePot*4.0*pot.epsilon
    sim.ePot = ePot
end
