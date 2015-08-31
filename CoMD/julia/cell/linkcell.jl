
include("constants.jl")

type LinkCell
    MAXATOMS::Int
    nLocalBoxes::Int
    nTotalBoxes::Int
    nAtoms::Array{Int}
    gridSize::Array{Int}
    invBoxSize::Array{Float}
end

function initLinkCells(pot, box_size)
    gridSize = zeros(Int, 3)
    for i in 1:3
        gridSize[i] = floor(Int, box_size[i]/pot.cutoff)
    end
    boxSize = [box_size[i]*1.0/gridSize[i] for i in 1:3]
    nLocalBoxes = gridSize[1] * gridSize[2] * gridSize[3]
    nHaloBoxes = 2*((gridSize[1] + 2) * (gridSize[2] + gridSize[3] + 2) + (gridSize[2] * gridSize[3]))
    nTotalBoxes = nLocalBoxes + nHaloBoxes


    invBoxSize = zeros(Float, 3)
    for i in 1:3
        invBoxSize[i] = 1.0/boxSize[i]
    end

    nAtoms = zeros(Int, nTotalBoxes)
    LinkCell(64, nLocalBoxes, nTotalBoxes, nAtoms, gridSize, invBoxSize)
end

function getBoxFromTuple(boxes, ix, iy, iz)
    gs = boxes.gridSize
    n = boxes.nLocalBoxes
    iBox = 0
    # Halo in Z+
    if iz == gs[3]
        iBox = n + 2*gs[3]*gs[2] + 2*gs[3]*(gs[1]+2) + (gs[1]+2)*(gs[2]+2) + (gs[1]+2)*(iy+1) + (ix+1)

    # Halo in Z-
    elseif iz == -1
        iBox = n + 2*gs[3]*gs[2] + 2*gs[3]*(gs[1]+2) + (gs[1]+2)*(iy+1) + (ix+1)

    # Halo in Y+
    elseif iy == gs[2]
        iBox = n + 2*gs[3]*gs[2] + gs[3]*(gs[1]+2) + (gs[1]+2)*iz + (ix+1)

    # Halo in Y-
    elseif iy == -1
        iBox = n + 2*gs[3]*gs[2] + iz*(gs[1]+2) + (ix+1)

    # Halo in X+
    elseif ix == gs[1]
        iBox = n +  gs[2]*gs[3] + iz*gs[2] + iy

    # Halo in X-
    elseif ix == -1
        iBox = n + iz*gs[2] + iy

    # Local link cell
    else
        iBox = ix + gs[1]*iy  + gs[1]*gs[2]*iz
    end

    # One-based indexing
    iBox += 1

    if iBox > boxes.nTotalBoxes
        println("iBox too large: $iBox, Total Boxes = $(boxes.nTotalBoxes)")
        println(" ix,iy,iz $ix, $iy, $iz")
        println(" gridSize $(boxes.gridSize)")
    end

    return iBox
end


function getBoxFromCoord(boxes, r)
    ix = floor(Int, r[1]*boxes.invBoxSize[1])
    iy = floor(Int, r[2]*boxes.invBoxSize[2])
    iz = floor(Int, r[3]*boxes.invBoxSize[3])


    getBoxFromTuple(boxes, ix, iy, iz)
end

function putAtomInBox(boxes, atoms, gid, iType, r, p=[0.0, 0.0, 0.0])
    #println("x,y,z = $x,$y,$z")
    iBox = getBoxFromCoord(boxes, r)
    iOff = (iBox-1)*boxes.MAXATOMS
    iOff += boxes.nAtoms[iBox] + 1

    if iBox < boxes.nLocalBoxes
        atoms.nLocal += 1
    end

    boxes.nAtoms[iBox] += 1

    atoms.gid[iOff] = gid
    atoms.species[iOff] = iType
    atoms.r[:,iOff] = r
    atoms.p[:,iOff] = p
end

function getTuple(boxes, iBox)
    ix = 0
    iy = 0
    iz = 0
    gs = boxes.gridSize

    iBox = iBox-1  # One-based indexing

    # local box
    if iBox < boxes.nLocalBoxes
        ix = iBox % gs[1]
        iBox = div(iBox, gs[1])
        iy = iBox % gs[2]
        iz = div(iBox, gs[2])
    else # halo box
        ink = iBox - boxes.nLocalBoxes
        if ink < 2*gs[2]*gs[3]
            if ink < gs[2]*gs[3]
                ix = 0
            else
                ink -= gs[2]*gs[3]
                ix = gs[1] + 1
            end
            iy = 1 + ink % gs[2]
            iz = 1 + div(ink, gs[2])
        elseif ink < 2*gs[3]*(gs[2] + gs[1] + 2)
            ink -= 2*gs[3]*gs[2]
            if ink < (gs[1]+2)*gs[3]
                iy = 0
            else
                ink -= (gs[1]+2)*gs[3]
                iy = gs[2]+1
            end
            ix = ink % (gs[2] + 2)
            iz = 1 + div(ink, gs[1] + 2)
        else
            ink -= 2*gs[3]*(gs[2] + gs[1] + 2)
            if ink < (gs[1] + 2)*(gs[2] + 2)
                iz = 0
            else
                ink -= (gs[1] + 2)*(gs[2] + 2)
                iz = gs[3] + 1
            end
            ix = ink %(gs[1] + 2)
            iy = div(ink, gs[1] + 2)
        end
        ix -= 1
        iy -= 1
        iz -= 1
    end
    return ix,iy,iz

end

function getNeighborBoxes(boxes, iBox, nbrBoxes)
    ix, iy, iz = getTuple(boxes, iBox)
    count = 0
    for i in ix-1:ix+1
        for j in iy-1:iy+1
            for k in iz-1:iz+1
                count += 1
                nbrBoxes[count] = getBoxFromTuple(boxes, i, j, k)
            end
        end
    end
    return count
end

function copyAtom(boxes, atoms, iAtom, iBox, jAtom, jBox)
    iOff = boxes.MAXATOMS * (iBox-1) + iAtom
    jOff = boxes.MAXATOMS * (jBox-1) + jAtom
    atoms.gid[jOff] = atoms.gid[iOff]
    atoms.species[jOff] = atoms.species[iOff]
    atoms.r[:,jOff] = atoms.r[:,iOff]
    atoms.p[:,jOff] = atoms.p[:,iOff]
    atoms.f[:,jOff] = atoms.f[:,iOff]
end


function moveAtom(boxes, atoms, iId, iBox, jBox)
    nj = boxes.nAtoms[jBox] + 1
    copyAtom(boxes, atoms, iId, iBox, nj, jBox)
    boxes.nAtoms[jBox] += 1

    boxes.nAtoms[iBox] -= 1
    ni = boxes.nAtoms[iBox] + 1
    if ni != 0
        copyAtom(boxes, atoms, ni, iBox, iId, iBox)
    end
    if jBox >= boxes.nLocalBoxes
        atoms.nLocal -= 1
    end
end

function emptyHaloCells(boxes)
    for ii in boxes.nLocalBoxes+1 : boxes.nTotalBoxes
        boxes.nAtoms[ii] = 0
    end
end

function updateLinkCells(atoms, boxes)
    emptyHaloCells(boxes)
    for iBox in 1:boxes.nLocalBoxes
        iOff = (iBox-1) * boxes.MAXATOMS
        ii = 1
        while ii <= boxes.nAtoms[iBox]
            jBox = getBoxFromCoord(boxes, atoms.r[:,iOff + ii])
            if jBox != iBox
                moveAtom(boxes, atoms, ii, iBox, jBox)
            else
               ii += 1
            end
        end
    end
end

