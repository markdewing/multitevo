
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



(* Use symbolic form for the potential *)

V[r_] := 4*epsilon*((sigma/r)^12 - (sigma/r)^6)

(* Optimization for this potential - create a function that takes r^2 as an argument - avoid taking sqrt *)
V2[s_] := V[Sqrt[s]]

(* Compute the force from the potential *)

Force = Simplify[-D[V[r], r]/r]
Force3[r_] = Force

(* Again, optimization to avoid sqrt *)
Force2[s_] := Force3[Sqrt[s]]

computeForceSym[atoms_, bounds_] := Module[
 {F, rCut2, eShift, ePot, dx, eLocal, fr},
 (
  F = Table[0.0, {Length[atoms]},{3}];
  rCut2 = cutoff*cutoff;
  eShift = V2[rCut2];
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
            eLocal = V2[r2] - eShift;
            ePot += eLocal;

            fr = Force2[r2];

            F[[i]] += dx*fr; (* be aware of change in sign from non-symbolic version *)
            F[[j]] -= dx*fr;
            )
        ];
      )
    ]
  ];
  {F, ePot}
)]


End[]

