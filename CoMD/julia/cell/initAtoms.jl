
include("constants.jl")

type Atoms
    r::Array{Float,2}
    f::Array{Float,2}
    p::Array{Float,2}
    species::Array{Int}
    gid::Array{Int}
    maxTotalAtoms::Int
    bounds::Array{Float}
    nLocal::Int
    nGlobal::Int
end

function initAtoms(n)
    r = zeros(Float, 3, n)
    f = zeros(Float, 3, n)
    p = zeros(Float, 3, n)
    species = ones(Int, n)
    gid = zeros(Int, n)
    maxTotalAtoms = n
    bounds = zeros(Float, 3)
    return Atoms(r, f, p, species, gid, maxTotalAtoms, bounds, 0, 0)
end

function setBounds(atoms, nx, ny, nz, lat)
    atoms.bounds[1] = nx*lat
    atoms.bounds[2] = ny*lat
    atoms.bounds[3] = ny*lat
end

function createFccLattice(nx, ny, nz, lat, atoms, boxes)
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

                    gid = ib + nb*(iz + nz*(iy + ny*ix))
                    putAtomInBox(boxes, atoms, gid, 1, [rx, ry, rz])
                    idx += 1
                end
            end
        end
    end
    atoms.nGlobal = idx - 1
end


