
include("constants.jl")

mutable struct Atoms
    r::Array{Float,2}
    f::Array{Float,2}
    p::Array{Float,2}
    species
    nAtoms::Int
    bounds::Array{Float}
end

function initAtoms(n)
    r = zeros(Float, 3, n)
    f = zeros(Float, 3, n)
    p = zeros(Float, 3, n)
    species = ones(Int32, n)
    nAtoms = n
    bounds = zeros(Float, 3)
    return Atoms(r, f, p, species, nAtoms, bounds)
end

function setBounds(atoms, nx, ny, nz, lat)
    atoms.bounds[1] = nx*lat
    atoms.bounds[2] = ny*lat
    atoms.bounds[3] = ny*lat
end

function createFccLattice(nx, ny, nz, lat, atoms)
    nb = 4  # number of atoms in this basis

    basis = Array[ [0.25, 0.25, 0.25],
                   [0.25, 0.75, 0.75],
                   [0.75, 0.25, 0.75],
                   [0.75, 0.75, 0.25] ]

    idx = 1
    for ix in 0:nx-1
        for iy in 0:ny-1
            for iz in 0:nz-1
                for ib in 1:nb
                    rx = (ix + basis[ib][1])*lat
                    ry = (iy + basis[ib][2])*lat
                    rz = (iz + basis[ib][3])*lat

                    atoms.r[1, idx] = rx
                    atoms.r[2, idx] = ry
                    atoms.r[3, idx] = rz
                    idx += 1
                end
            end
        end
    end
end


