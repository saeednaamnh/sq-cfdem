# sq-cfdem — Superquadric (prolate spheroid) extension for CFDEMcoupling-PUBLIC

Companion code project to:
Naamneh (2026), "Orientation-dependent transport coefficients of sheared
prolate granular gases from DEM", and the planned CFD-DEM follow-up.

## What this is
CFDEMcoupling-PUBLIC (OpenFOAM-6 branch) ships the *interfaces* for
superquadric coupling — virtual accessors `quat()/shapeArray()/blockiness()`
in `cfdemCloud`, a `<model>Superquadric` dispatch in `newVoidFractionModel`,
`requiresQuaternion()/requiresSuperquadric()` hooks in `forceModel`, and the
`registerFieldsToDEM()` exchange registry — but the concrete implementations
are premium-only. This project implements the missing classes under GPLv3.

## Layout
```
install/        five numbered scripts: prereqs -> OpenFOAM-6 -> sources -> env -> build
src/sqCfdem/    the library (wmake libso)
  common/spheroidGeometry.H            fully implemented math core (drag,
                                       projected areas, sphericities, quats,
                                       Jeffery/Bretherton helpers)
  cfdemCloudRotationSuperquadric/      cloud subclass: field registration +
                                       accessors + hdtorque push  [skeleton]
  voidFractionModel/dividedVoidFractionSuperquadric/
                                       satellite-point voidage    [core done;
                                       parallel scatter-gather TODO]
  forceModel/HolzerSommerfeldDrag/     orientation-dependent drag [near done]
  forceModel/spheroidLiftTorque/       lift + pitching torque     [core done;
                                       DNS coefficient calibration via dict]
patches/        LIGGGHTS 3.8.0 verification notes (NO patch needed: hdtorque
                exchange + torque application already in PUBLIC source)
verification/   V0..V7 ladder, one README each with pass tolerances
docs/           full design document
```

## Build order (after install/ scripts)
1. LIGGGHTS as library, `USE_SUPERQUADRICS = "ON"` in MAKE/Makefile.user
2. `cfdemCompCFDEMsrc && cfdemCompCFDEMsol`
3. `cd src/sqCfdem && wmake libso`

## couplingProperties (user-facing selection)
```
particleShapeType   superquadric;      // activates *Superquadric dispatch
voidFractionModel   divided;           // -> dividedSuperquadric (this lib)
forceModels ( HolzerSommerfeldDrag spheroidLiftTorque gradPForce viscForce Archimedes );
```

## Honest status
Code is written against interfaces read directly from the upstream OF-6
branch (signatures of `setvoidFraction`, `setForce`, `registerFieldsToDEM`,
registration macros mirrored from `DiFeliceDrag`/`dividedVoidFraction`).
It has NOT been compiled — expect one focused iteration of wmake fixes,
concentrated in: forceSubModel accessor names (rho/nu fields, partToArray
signature), the locateModel single-point call, and the exact LIGGGHTS
property names in the exchange registry. Grep markers: `TODO(build-1)`
(first-compile items), `TODO(build-2)` (parallel scatter-gather,
periodic images), `PLACEHOLDER` (lift/torque DNS coefficients).

## License
GPLv3, inherited from CFDEMcoupling/LIGGGHTS. Cite via CITATION.cff (add
on first release; mint a Zenodo DOI).
