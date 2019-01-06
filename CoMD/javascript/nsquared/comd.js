
function createFccLattice(nx, ny, nz, lat, atoms)
{
    const nb = 4;
    const basis = [ [0.25, 0.25, 0.25],
                    [0.25, 0.75, 0.75],
                    [0.75, 0.25, 0.75],
                    [0.75, 0.75, 0.25]
                  ];

    let idx = 0;
    for (let ix = 0; ix < nx; ix++) {
        for (let iy = 0; iy < ny; iy++) {
            for (let iz = 0; iz < nz; iz++) {
                for (let ib = 0; ib < nb; ib++) {
                    rx = (ix + basis[ib][0]) * lat;
                    ry = (iy + basis[ib][1]) * lat;
                    rz = (iz + basis[ib][2]) * lat;

                    atoms[idx][0] = rx;
                    atoms[idx][1] = ry;
                    atoms[idx][2] = rz;
                    idx += 1;
                }
            }
        }
    }

    // idx should equal nx*ny*nz*nb
}

function run_comd()
{
    console.log("Running CoMD");
}

run_comd()

