#!/usr/local/bin/MathematicaScript -script

(* Translation of CoMD to Wolfram Language (Mathematica) *)

Import["ljforce.m"]
Import["initatoms.m"]

nx = 4
ny = 4
nz = 4
temperature = 600 (* in Kelvin *)
dt = 1.0 (* in fs *)
printRate = 10 (* steps between printing output *)
nBlock = 10


(* Setup *)

bounds = {nx, ny, nz} * LJPot`lat
atoms = createFccLattice[nx, ny, nz, LJPot`lat]
p = setTemperature[Length[atoms], LJPot`mass, temperature]

ke = kineticEnergy[p, LJPot`mass];
{F, pe} = LJPot`computeForceSym[atoms, bounds]


(* Formatted output *)

Print[OutputForm[StringForm["Initial PE (cohesive energy) : ``", pe/Length[atoms]]]]
Print[OutputForm[StringForm["Initial energy : `` atom count: ``", (ke+pe)/Length[atoms], Length[atoms]]]]

PadString[s_,nchar_] := StringJoin[PadLeft[Characters[s], nchar, " "]]

PrintHeader := 
(
 Print[OutputForm[PadString["Performance",103]]];
 Print[OutputForm[
    StringJoin[
        "#",
        PadString["Loop",6],
        PadString["Time(fs)",11],
        PadString["Total Energy",19],
        PadString["Potential Energy",19],
        PadString["Kinetic Energy",19],
        PadString["Temperature",13],
        PadString["(us/atom)",13],
        PadString["# Atoms",12]
    ]]]
)

PrintHeader[]

PrintInfo[ke_,pe_,time_,nAtoms_,nStep_,iStep_] := 
(
  eTotal = (ke + pe)/nAtoms;
  timePerAtom = 10^6 * time/(nAtoms*nStep);
  temp = (ke/nAtoms)/(kBeV*1.5);
  (*Print["ke",ke/nAtoms,"pe",pe/nAtoms,"total e",eTotal,"time",OutputForm[time],"us per atom",OutputForm[timePerAtom]]*)
  Print[OutputForm[
   StringJoin[
    ToString[PaddedForm[iStep, 6]],
    ToString[PaddedForm[iStep*dt,{9,2}]],
    ToString[PaddedForm[eTotal,{17,12}]],
    ToString[PaddedForm[pe/nAtoms,{17,12}]],
    ToString[PaddedForm[ke/nAtoms,{17,12}]],
    ToString[PaddedForm[temp,{11,4}]],
    ToString[PaddedForm[timePerAtom,{10,4}]],
    ToString[PaddedForm[nAtoms,11]]
   ]
  ]]
)


(* Advance position and velocity for several timesteps *)

timestep[nSteps_, dt_, atoms_, p_, bounds_, F_] := Module[
{p2=p,atoms2=atoms,F2=F},
(
 Do[
 (
  p2 += 0.5*dt*F2; (* advanceVelocity *)
  atoms2 += dt*p2/LJPot`mass; (* advancePosition *)
  {F2, ePot} = LJPot`computeForceSym[atoms2, bounds];
  p2 += 0.5*dt*F2; (* advanceVelocity *)
 ), {nSteps}];
 {atoms2, p2, F2, ePot})
]

timestepTimeOneIteration = 0.0

(* Main code *)

Do[(
     PrintInfo[ke,pe,timestepTimeOneIteration,Length[atoms],printRate,jStep*printRate];
     {timestepTimeOneIteration, res} =
     Timing[
        ({atoms, p, F, pe} = timestep[printRate, dt, atoms, p, bounds, F];
         ke = kineticEnergy[p, LJPot`mass];
        )
     ]
    ),
    {jStep, nBlock}
  ]

