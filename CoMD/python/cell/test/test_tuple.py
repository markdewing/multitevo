
import sys
sys.path.append('..')
import simflat
import initatoms
import linkcell

import numpy as np

sim = simflat.SimFlat(10, 10, 0.1)

lat = 3.615
nx = 4
ny = 4
nz = 4

box_size = np.zeros(3)
box_size[0] = nx*lat
box_size[1] = ny*lat
box_size[2] = nz*lat

boxes = linkcell.initLinkCells(sim, box_size)
sim.atoms = initatoms.atoms(boxes.nTotalBoxes * boxes.MAXATOMS)
sim.atoms.setBounds(nx,ny,nz,lat)
sim.boxes = boxes

initatoms.createFccLattice(nx, ny, nz, lat, sim.atoms, sim.boxes)


def test_tuple(boxes):
    """
        Test the conversion from 3-D indices to a linear index and back
    """
    ix = 0
    iy = 0
    iz = 0
    for ix in range(-1,boxes.gridSize[0]+1):
        for iy in range(-1,boxes.gridSize[1]+1):
            for iz in range(-1,boxes.gridSize[2]+1):
                iBox = boxes.getBoxFromTuple(ix, iy, iz)
                nix, niy, niz = boxes.getTuple(iBox)
                if ix != nix or iy != niy or iz != niz:
                    print "Difference found, iBox = ",iBox
                    print "  Original indices: ",ix,iy,iz
                    print "      New indices: ",nix, niy, niz


test_tuple(boxes)
