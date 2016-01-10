

Import["constants.m"] 

(* Create atoms positions on a face centered cubic (FCC) lattice
   with nx*ny*nz unit cells and lattice constant 'lat' *)

createFccLattice[nx_, ny_, nz_, lat_] := Module[
  {basis, atoms},
  (basis = { {0.25, 0.25, 0.25},
            {0.25, 0.75, 0.75},
            {0.75, 0.25, 0.75},
            {0.75, 0.75, 0.25} };
  atoms =  {};
  Do[ 
    AppendTo[atoms, ({ix, iy, iz} + basis[[ib]])*lat],
       {ib, Length[basis]}, {ix, nx}, {iy, ny}, {iz, nz}
    ];
  atoms)
]


kineticEnergy[p_, mass_] := 0.5 * Total[Map[Dot[#, #] &, p]]/mass 

setTemperature[nAtoms_, mass_, temperature_] := Module[
    {sigma, p, Vcm, ke, scale},
   (
    sigma = Sqrt[kBeV * temperature / mass];
    p = mass * sigma * RandomVariate[NormalDistribution[], {nAtoms, 3}];

    (* compute center of mass velocity *)

    Vcm = Total[p]/Length[p];

    (* adjust center of mass velocity to be zero *)
    p -= Table[Vcm, {nAtoms}];

    ke = kineticEnergy[p, mass];
    temp = ke/nAtoms/kBeV/1.5;
    scale = Sqrt[temperature/temp];
    p *= scale
   )
]

