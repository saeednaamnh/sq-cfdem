# Known Issues

## RESOLVED: DEM script must not run its own time loop under CFDEM
Symptom: coupled runs produced results byte-identical to dry runs; no CFD
fields written; OpenFOAM "Time =" never printed; no "Starting time loop".

Root cause: the DEM input script contained `run ${N_STEPS}`. Under CFDEM the
CFD solver owns the time loop and issues `run <couplingInterval>` via the
runLiggghts command model. A `run N` in the DEM script executes the entire
simulation during cfdemCloud construction, so the solver never reaches
`while (runTime.loop())`. Result: no PISO solve, no drag on particles, no
field output -- a silently dry simulation.

Fix: DEM script ends with `run 0` only (matches CFDEM tutorial pattern).
Baked into sqfoam new.

Validation: compare granular_temp.txt against the dry reference AT THE SAME
STEP. Differ = coupled. Identical = broken. (md5 of dumps is not reliable
mid-run; partial writes give false divergence.)
