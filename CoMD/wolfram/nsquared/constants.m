
(* could (and should, eventually) get all these from built-in unit conversions:  UnitConvert[Quantity[1, "eV"], "Joules"] *)

(* 1 amu in kilograms *)
amuInKilograms = 1.660538921 10^-27

(* 1 femtosecond in seconds *)
fsInSeconds = 1.0 10^-15

(* 1 Angstrom in meters *)
angsInMeters = 1.0 10^-10

(* 1 eV (electron Volt) in Joules *)
eVInJoules = 1.602176565 10^-19

amuToInternalMass = amuInKilograms * angsInMeters * angsInMeters / (fsInSeconds * fsInSeconds * eVInJoules)

(* Boltzmann constant in eV's *)
kBeV = 8.6173324 10^-5

(* Hartrees to eVs *)
hartreeToEv = 27.21138505

(* Bohrs to Angstroms *)
bohrToAngs = 0.52917721092

