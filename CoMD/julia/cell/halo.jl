
include("constants.jl")

HALO_X_MINUS = 1
HALO_X_PLUS = 2
HALO_Y_MINUS = 3
HALO_Y_PLUS = 4
HALO_Z_MINUS = 5
HALO_Z_PLUS = 6

HALO_X_AXIS = 1
HALO_Y_AXIS = 2
HALO_Z_AXIS = 3

mutable struct HaloExchange
    nCells::Array{Int, 1}
    pbcFactor::Array{Float, 2}
    cellList::Array{Array{Int, 1}, 1}
end

function initHalo(boxes)
    nCells = zeros(Int, 6)
    pbcFactor = zeros(Float, 6, 3)

    gs = boxes.gridSize

    nCells[HALO_X_MINUS] = 2*(gs[2]+2) * (gs[3] + 2)
    nCells[HALO_Y_MINUS] = 2*(gs[1]+2) * (gs[3] + 2)
    nCells[HALO_Z_MINUS] = 2*(gs[1]+2) * (gs[2] + 2)
    nCells[HALO_X_PLUS] = nCells[HALO_X_MINUS]
    nCells[HALO_Y_PLUS] = nCells[HALO_Y_MINUS]
    nCells[HALO_Z_PLUS] = nCells[HALO_Z_MINUS]

    cellList = Array{Any}(undef, 6)
    for ii in 1:6
        cellList[ii] = mkAtomCellList(boxes, ii, nCells[ii])
    end

    pbcFactor[HALO_X_MINUS, HALO_X_AXIS] = 1.0
    pbcFactor[HALO_X_PLUS, HALO_X_AXIS] = -1.0
    pbcFactor[HALO_Y_MINUS, HALO_Y_AXIS] = 1.0
    pbcFactor[HALO_Y_PLUS, HALO_Y_AXIS] = -1.0
    pbcFactor[HALO_Z_MINUS, HALO_Z_AXIS] = 1.0
    pbcFactor[HALO_Z_PLUS, HALO_Z_AXIS] = -1.0

    HaloExchange(nCells, pbcFactor, cellList)
end

function mkAtomCellList(boxes, iFace, nCell)
    xBegin = -1
    xEnd = boxes.gridSize[1] + 1
    yBegin = -1
    yEnd = boxes.gridSize[2] + 1
    zBegin = -1
    zEnd = boxes.gridSize[3] + 1

    if iFace == HALO_X_MINUS
        xEnd = xBegin + 2
    end
    if iFace == HALO_X_PLUS
        xBegin = xEnd - 2
    end
    if iFace == HALO_Y_MINUS
        yEnd = yBegin + 2
    end
    if iFace == HALO_Y_PLUS
        yBegin = yEnd - 2
    end
    if iFace == HALO_Z_MINUS
        zEnd = zBegin + 2
    end
    if iFace == HALO_Z_PLUS
        zBegin = zEnd - 2
    end

    list = Array{Int}(undef,nCell)

    ii = 1
    for ix in xBegin:xEnd-1
        for iy in yBegin:yEnd-1
            for iz in zBegin:zEnd-1
                iBox = getBoxFromTuple(boxes, ix, iy, iz)
                list[ii] = iBox
                ii += 1
            end
        end
    end

    return list
end

function haloExchange(halo, atoms, boxes)
    for iAxis in 1:3
        exchangeData(halo, iAxis, atoms, boxes)
    end
end

function exchangeData(halo, iAxis, atoms, boxes)
    faceM = 2*(iAxis-1) + 1
    faceP = faceM + 1
    bufM = loadAtomsBuffer(halo, atoms, boxes, faceM)
    bufP = loadAtomsBuffer(halo, atoms, boxes, faceP)
    unloadAtomsBuffer(halo, atoms, boxes, faceM, bufP)
    unloadAtomsBuffer(halo, atoms, boxes, faceP, bufM)
end

mutable struct TransferBuffer
    gid::Int
    species::Int
    r::Array{Float}
    p::Array{Float}
end

function loadAtomsBuffer(halo, atoms, boxes, face)
    shift = halo.pbcFactor[face,:] .* atoms.bounds

    nCells = halo.nCells[face]
    cellList = halo.cellList[face]
    size = 0
    for iCell in 1:nCells
        iBox = cellList[iCell]
        for ii in 1:boxes.nAtoms[iBox]
            size += 1
        end
    end

    buf = Array{TransferBuffer}(undef, size)
    idx = 1
    for iCell in 1:nCells
        iBox = cellList[iCell]
        for ii in 1:boxes.nAtoms[iBox]
            iOff = (iBox-1)*boxes.MAXATOMS + ii
            pos = atoms.r[:, iOff] + shift
            entry = TransferBuffer(atoms.gid[iOff], atoms.species[iOff], pos, atoms.p[:,iOff])
            buf[idx] = entry
            idx += 1
        end
    end
    return buf
end

function unloadAtomsBuffer(halo, atoms, boxes, face, buf)
    for entry in buf
        putAtomInBox(boxes, atoms, entry.gid, entry.species, entry.r, entry.p)
    end
end


