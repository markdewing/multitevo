
Import["constants.m"]

Begin["LJPot`"]

sigma = 2.315
epsilon = 0.167
mass = 63.55 * amuToInternalMass
lat = 3.615

cutoff = 2.5*sigma


(* Compute force and potential energy in same loop *)

computeForce[atoms_, bounds_] := Module[
 {F, rCut2, rCut6, eShift, ePot, dx, r6, eLocal, fr},
 (
  F = Table[0.0, {Length[atoms]},{3}];
  rCut2 = cutoff*cutoff;
  rCut6 = sigma^6/cutoff^6;
  eShift = rCut6*(rCut6-1.0);
  ePot = 0.0;

  For[i = 1, i <= Length[atoms], i++,
    For[j = 1, j < i, j++, (
        dx = atoms[[i]] - atoms[[j]];
        Do[
            (dx[[k]] = If[dx[[k]] < -0.5*bounds[[k]],
                          dx[[k]] + bounds[[k]], dx[[k]]];
             dx[[k]] = If[dx[[k]] > 0.5*bounds[[k]],
                          dx[[k]] - bounds[[k]], dx[[k]]];),
            {k, 3}];

        r2 = Dot[dx, dx];
        If[r2 <= rCut2, (
            r6 = (sigma*sigma/r2)^3;
            eLocal = r6*(r6-1.0) - eShift;
            ePot += eLocal;

            fr = -4.0*epsilon*r6*(12.0*r6 - 6.0)/r2;

            F[[i]] -= dx*fr;
            F[[j]] += dx*fr;
            )
        ];
      )
    ]
  ];
  {F, 4*epsilon*ePot}
)]

End[]

