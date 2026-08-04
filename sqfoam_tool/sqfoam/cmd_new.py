import re
import os, re, shutil, math
from pathlib import Path
from config import load, SHAPES, contact_dt, fluid_props

def _coupling(): return COUPLING
def _blockmesh(): return BLOCKMESH

def cmd_new(a):
    cfg=load(); s=SHAPES[a.shape]
    name=a.name or f"{a.shape}_St{int(a.St) if a.St==int(a.St) else a.St}"
    parent=Path(a.out) if a.out else Path.cwd()
    cdir=parent/name
    if cdir.exists():
        print(f"!! {cdir} exists — pick another --name or remove it"); return
    rho_f,nu_f,mu_f,tau_p=fluid_props(s,a.St)

    # --- DEM script: transform the user's in_shear.lmp template ---
    tmpl=None
    for cand in [Path(cfg["DRY_RESTARTS"])/a.shape/"in_shear.lmp",
                 Path(cfg["DRY_RESTARTS"])/"prolate_3"/"in_shear.lmp"]:
        if cand.is_file(): tmpl=cand; break
    if not tmpl:
        print("!! no in_shear.lmp template found under DRY_RESTARTS"); return
    txt=tmpl.read_text()
    txt=txt.replace("variable    a           equal",
        f"variable    N_STEPS     equal   {a.steps}\nvariable    a           equal",1)
    txt=txt.replace("run  5000000","run  ${N_STEPS}")
    # coupling needs an atom map (fix couple/cfd requirement)
    if "atom_modify" not in txt:
        txt=txt.replace("atom_style      superquadric",
                        "atom_style      superquadric\natom_modify     map array")
    txt=re.sub(r'print "  Steps: \d+"','print "  Steps: ${N_STEPS}"',txt)
    for pat in ["restarts/","results/"]:
        txt=txt.replace(pat,"../DEM/"+pat)
    txt=txt.replace(
"fix  integrate layerGroup nve/superquadric integration_scheme 1",
"""fix  cfd  all       couple/cfd couple_every 100 mpi
fix  cfd2 layerGroup couple/cfd/force/implicit transfer_torque yes

fix  integrate layerGroup nve/superquadric integration_scheme 1""",1)

    # CRITICAL: under CFDEM the CFD solver drives time stepping via
    # runLiggghts. A "run N" in the DEM script blocks cloud construction
    # -> CFD loop never iterates -> no coupling, no fields, dry results.
    txt = re.sub(r"(?m)^run\s+\$\{N_STEPS\}.*$",
                 "# run ${N_STEPS}  <-- disabled: CFD drives stepping", txt)
    txt = re.sub(r"(?m)^run\s+5000000.*$",
                 "# run 5000000  <-- disabled: CFD drives stepping", txt)

    (cdir/"DEM"/"restarts").mkdir(parents=True)
    (cdir/"DEM"/"post"/"restart").mkdir(parents=True, exist_ok=True)
    (cdir/"DEM"/"results"/"shear"/"profiles").mkdir(parents=True)
    (cdir/"DEM"/"in.liggghts_run").write_text(txt)

    # symlink placement restart (particle-for-particle with Table I)
    src_restart=Path(cfg["DRY_RESTARTS"])/a.shape/"restarts"/"restart_placed.liggghts"
    link=cdir/"DEM"/"restarts"/"restart_placed.liggghts"
    if src_restart.is_file():
        os.symlink(src_restart, link); rlink=f"linked -> {src_restart}"
    else:
        rlink=f"MISSING: put restart_placed.liggghts at {link}"

    # --- CFD side ---
    cfd=cdir/"CFD"; (cfd/"0").mkdir(parents=True)
    (cfd/"system").mkdir(); (cfd/"constant").mkdir()
    tut=Path(cfg["TUTORIAL_CFD"])
    for f in ["fvSchemes","fvSolution","decomposeParDict"]:
        shutil.copy(tut/"system"/f, cfd/"system"/f)
    for f in ["g","liggghtsCommands","turbulenceProperties","RASProperties"]:
        if (tut/"constant"/f).is_file(): shutil.copy(tut/"constant"/f, cfd/"constant"/f)
    dp=(cfd/"system"/"decomposeParDict"); t=dp.read_text()
    t=re.sub(r'n\s*\(\s*\d+\s+\d+\s+\d+\s*\)','n               (2 2 1)',t); dp.write_text(t)

    # _disable_writeLiggghts: needs a post/restart path; not needed here
    _lc = cfd/"constant"/"liggghtsCommands"
    if _lc.is_file():
        _t = _lc.read_text()
        _t = _t.replace("\n    writeLiggghts", "\n    //writeLiggghts")
        _lc.write_text(_t)

    (cfd/"constant"/"transportProperties").write_text(
f"""FoamFile {{ version 2.0; format ascii; class dictionary; object transportProperties; }}
transportModel  Newtonian;
nu              nu [ 0 2 -1 0 0 0 0 ] {nu_f:.6e};
""")
    # use V6's proven couplingProperties template (has modelType, full config)
    _tmpl=Path(__file__).parent/"couplingProperties.template"
    if _tmpl.is_file():
        (cfd/"constant"/"couplingProperties").write_text(_tmpl.read_text())
    else:
        (cfd/"constant"/"couplingProperties").write_text(COUPLING)
    coupl=100*contact_dt(s); dtC=coupl/2.0; end=a.steps*contact_dt(s)
    (cfd/"system"/"controlDict").write_text(
f"""FoamFile {{ version 2.0; format ascii; class dictionary; object controlDict; }}
application     cfdemSolverPisoSQ;
startFrom startTime; startTime 0; stopAt endTime;
endTime         {end:.4f};
deltaT          {dtC:.16e};
writeControl adjustableRunTime; writeInterval {100*dtC:.16f}; purgeWrite 0;
writeFormat ascii; writePrecision 7; timeFormat general; timePrecision 6;
runTimeModifiable yes;
""")
    (cfd/"system"/"blockMeshDict").write_text(BLOCKMESH)
    _mkfields(cfd/"0", rho_f)

    # resumable runner
    (cdir/"run.sh").write_text(RUNNER); os.chmod(cdir/"run.sh",0o755)

    print(f"created {cdir}")
    print(f"  shape={a.shape} (c/a={s['ca']})  St={a.St}  nu_f={nu_f:.3e}  tau_p={tau_p:.1f}")
    print(f"  restart: {rlink}")
    print(f"  -> sqfoam run {cdir}")

def _mkfields(z, rho_f):
    def w(n,cls,dim,intern,wb,wt):
        (z/n).write_text(
f"""FoamFile {{ version 2.0; format ascii; class {cls}; object {n}; }}
dimensions {dim}; internalField uniform {intern};
boundaryField {{ cycXa {{ type cyclic; }} cycXb {{ type cyclic; }}
 cycYa {{ type cyclic; }} cycYb {{ type cyclic; }}
 wallBottom {{ {wb} }} wallTop {{ {wt} }} }}
""")
    w("U","volVectorField","[0 1 -1 0 0 0 0]","(0 0 0)",
      "type fixedValue; value uniform (-1 0 0);","type fixedValue; value uniform (1 0 0);")
    w("p","volScalarField","[0 2 -2 0 0 0 0]","0","type zeroGradient;","type zeroGradient;")
    w("voidfraction","volScalarField","[0 0 0 0 0 0 0]","1","type zeroGradient;","type zeroGradient;")
    w("Ksl","volScalarField","[1 -3 -1 0 0 0 0]","0","type zeroGradient;","type zeroGradient;")
    w("Us","volVectorField","[0 1 -1 0 0 0 0]","(0 0 0)","type zeroGradient;","type zeroGradient;")
    w("rho","volScalarField","[1 -3 0 0 0 0 0]",f"{rho_f:.4e}","type zeroGradient;","type zeroGradient;")
    w("k","volScalarField","[0 2 -2 0 0 0 0]","0","type zeroGradient;","type zeroGradient;")
    w("epsilon","volScalarField","[0 2 -3 0 0 0 0]","0","type zeroGradient;","type zeroGradient;")
    w("nut","volScalarField","[0 2 -1 0 0 0 0]","0","type zeroGradient;","type zeroGradient;")

COUPLING='''FoamFile { version 2.0; format ascii; class dictionary; object couplingProperties; }
couplingInterval 100;
voidFractionModel divided;
particleShapeType superquadric;
locateModel engine;
meshMotionModel noMeshMotion;
IOModel off; probeModel off;
dataExchangeModel twoWayMPI;
averagingModel dilute;
clockModel off; smoothingModel off;
forceModels ( HolzerSommerfeldDrag spheroidRotationTorque gradPForce viscForce );
momCoupleModels ( implicitCouple );
turbulenceModelType turbulenceProperties;
implicitCoupleProps { velFieldName "U"; granVelFieldName "Us"; voidfractionFieldName "voidfraction"; }
dividedProps { alphaMin 0.1; }
dividedSuperquadricProps { alphaMin 0.10; nSatellites 64; maxCellsPerParticle 30; }
HolzerSommerfeldDragProps { velFieldName "U"; voidfractionFieldName "voidfraction"; voidageCorrection off; }
spheroidRotationTorqueProps { velFieldName "U"; }
gradPForceProps { pFieldName "p"; velocityFieldName "U"; voidfractionFieldName "voidfraction"; }
viscForceProps { divTauFieldName "divTau"; velocityFieldName "U"; voidfractionFieldName "voidfraction"; }
engineProps { treeSearch true; }
'''

BLOCKMESH='''FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
convertToMeters 1;
vertices ( (0 0 0)(15 0 0)(15 10 0)(0 10 0)(0 0 10)(15 0 10)(15 10 10)(0 10 10) );
blocks ( hex (0 1 2 3 4 5 6 7) (8 5 5) simpleGrading (1 1 1) );
edges ();
boundary
( cycXa { type cyclic; neighbourPatch cycXb; faces ((0 4 7 3)); }
  cycXb { type cyclic; neighbourPatch cycXa; faces ((1 2 6 5)); }
  cycYa { type cyclic; neighbourPatch cycYb; faces ((0 1 5 4)); }
  cycYb { type cyclic; neighbourPatch cycYa; faces ((3 7 6 2)); }
  wallBottom { type wall; faces ((0 3 2 1)); }
  wallTop    { type wall; faces ((4 5 6 7)); } );
mergePatchPairs ();
'''

RUNNER='''#!/bin/bash
# resumable coupled runner (generated by sqfoam)
set -e
cd "$(dirname "$0")"
RS="DEM/restarts/restart_shear.liggghts"
DEMIN="DEM/in.liggghts_run"
if [ -f "$RS" ]; then
  echo "=== RESUMING from shear checkpoint ==="
  sed -i 's#^read_restart.*#read_restart    ../DEM/restarts/restart_shear.liggghts#' "$DEMIN"
else
  echo "=== FRESH START from placement ==="
  sed -i 's#^read_restart.*#read_restart    ../DEM/restarts/restart_placed.liggghts#' "$DEMIN"
fi
cd CFD
[ -d processor0 ] || { blockMesh > ../log_mesh.txt 2>&1; decomposePar -force > ../log_decomp.txt 2>&1; }
mpirun -np ${NCORES:-4} cfdemSolverPisoSQ -parallel 2>&1 | tee -a ../log_solve.txt \\
  | grep -E "satellite points|Shear rate|FATAL|Courant"
'''
