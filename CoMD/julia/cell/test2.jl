
include("ljforce.jl")
include("simflat.jl")
include("initAtoms.jl")
include("linkcell.jl")
include("halo.jl")

#sim = SimFlat(10, 0.1)

#atoms.r[1,2] = 6.0
pot = init_lj_pot()
species = initSpecies(pot)

lat = 3.615
nx = 4
ny = 4
nz = 4

box_size = [nx,ny,nz]*lat

boxes = initLinkCells(pot, box_size)

atoms = initAtoms(boxes.nTotalBoxes * boxes.MAXATOMS)
setBounds(atoms,nx,ny,nz,lat)

halo = initHalo(boxes)

nSteps = 10
printRate = 10
dt = 0.1


createFccLattice(nx, ny, nz, lat, atoms, boxes)

sim = SimFlat(nSteps, printRate, dt, atoms, species, pot, 0.0, 0.0, boxes, halo)


println("box_size $box_size")
println("gridSize $(boxes.gridSize)")
println("nLocalBoxes $(boxes.nLocalBoxes)")
println("nTotalBoxes $(boxes.nTotalBoxes)")
println("invBoxSize $(boxes.invBoxSize)")

#for ix in -1:6
#    iBox = getBoxFromTuple(boxes, ix, 0, 0)
#    println("ix = $ix, iBox = $iBox")
#end

function test_tuple(boxes)
    ix = 0
    iy = 0
    iz = 0
    #println("ix,iy,iz = $ix,$iy,$iz")
    #iBox = getBoxFromTuple(boxes, ix, iy, iz)
    #println("iBox = $iBox")
    for ix in -1:boxes.gridSize[1]
        for iy in -1:boxes.gridSize[2]
            for iz in -1:boxes.gridSize[2]
                iBox = getBoxFromTuple(boxes, ix, iy, iz)
                #iBox -= 1
                println("$ix $iy $iz $iBox")
                nix, niy, niz = getTuple(boxes, iBox)
                if nix != ix || niy != iy || niz != iz
                    println("diff: $nix $niy $niz")
                end
            end
        end
    end
end

#x = 0.90374
#y = 9.94125
#z = 6.32525
#println("x,y,z = $x,$y,$z")
#iBox = getBoxFromCoord(boxes, [x,y,z])
#println("iBox = $iBox")

#iBox = 7
#nbr=zeros(Int, 27)
#nBox = getNeighborBoxes(boxes, iBox, nbr)
#println("nBox = $nBox")
#println("nbr = $nbr")


#computeForce(sim.pot, sim)
#println("ke = ",sim.eKinetic)
#println("pe = ",sim.ePot)
#for i in 1:sim.atoms.nAtoms
#    println("i, ",i," force = ",sim.atoms.f[:,i])
#end

#updateLinkCells(boxes, sim.atoms)
#haloExchange(halo, sim)
redistributeAtoms(sim.atoms, sim.boxes, sim.halo)

pot.forceFunc(pot, sim)
e_pot = sim.ePot/sim.atoms.nGlobal
println("pot e = $e_pot")

