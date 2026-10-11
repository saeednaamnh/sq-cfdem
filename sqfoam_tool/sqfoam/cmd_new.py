import re
import os, re, shutil, math
from pathlib import Path
from config import load, SHAPES, contact_dt, fluid_props

def _coupling(): return COUPLING
def _blockmesh(): return BLOCKMESH

import case_checks as C


def _start_restart(tdir, start):
    """Restart file the case starts from: the placement, or the dry baseline's final state."""
    if start == "placed":
        return tdir/"restarts"/"restart_placed.liggghts"
    rdir = tdir/"restarts"
    final = sorted(rdir.glob("restart_shear*final*.liggghts"))
    if len(final) == 1:
        return final[0]                      # the explicitly final state of the dry run
    if len(final) > 1:
        raise C.CaseError("--start dry: several final restarts in %s: %s"
                          % (rdir, ", ".join(p.name for p in final)))
    cands = sorted(rdir.glob("restart_shear*.liggghts"), key=lambda p: p.stat().st_mtime)
    if not cands:
        raise C.CaseError("--start dry: no restart_shear*.liggghts in %s" % rdir)
    return cands[-1]                         # otherwise the most recent shear restart


def _template(tdir):
    """The DEM shear script of a baseline folder: in_shear.lmp, or the only in_shear*.lmp."""
    exact = tdir/"in_shear.lmp"
    if exact.is_file():
        return exact
    cands = sorted(tdir.glob("in_shear*.lmp"))
    if len(cands) == 1:
        return cands[0]
    if not cands:
        raise C.CaseError("no in_shear*.lmp in %s (set TEMPLATE_DIRS in ~/.sqfoam.conf)" % tdir)
    raise C.CaseError("several DEM scripts in %s (%s): keep one, or name the right one in_shear.lmp"
                      % (tdir, ", ".join(p.name for p in cands)))


def cmd_new(a):
    cfg=load(); s=SHAPES[a.shape]
    name=a.name or f"{a.shape}_St{int(a.St) if a.St==int(a.St) else a.St}"
    parent=Path(a.out) if a.out else Path.cwd()
    cdir=parent/name
    if cdir.exists():
        print(f"!! {cdir} exists — pick another --name or remove it"); return
    rho_f,nu_f,mu_f,tau_p=fluid_props(s,a.St)

    # --- DEM script: the shape's own published dry baseline, no fallback ---
    tdir=Path(cfg["DRY_RESTARTS"])/cfg["TEMPLATE_DIRS"][a.shape]
    try:
        tmpl=_template(tdir)
    except C.CaseError as e:
        print(f"!! {e}"); return 1
    txt=tmpl.read_text()
    ax=C.semi_axes(txt)
    if ax is None or max(abs(x-y) for x,y in zip(ax,(s["a"],s["b"],s["c"])))>1e-6:
        print(f"!! {tmpl} has semi-axes {ax}, not those of {a.shape}"); return 1
    uw=C.wall_velocity(txt); up,lo=C.wall_motion_signs(txt)
    if uw is None or up is None or lo is None:
        print(f"!! cannot read wall_vel / wall 'move linear' fixes in {tmpl}"); return 1
    if abs(uw-s["U"])>1e-9:
        print(f"!! {tmpl} moves the walls at {uw}, but the {a.shape} baseline is at U={s['U']}"); return 1
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
    try:
        txt, group = C.inject_coupling(txt, couple_every=100)
    except C.CaseError as e:
        print(f"!! {e}"); return 1

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

    # starting state: placement (cold start) or the dry baseline's final state.
    # run.sh resumes from restart_shear.liggghts if present, else restart_placed.
    try:
        src_restart=_start_restart(tdir, a.start)
    except C.CaseError as e:
        print(f"!! {e}"); return 1
    lname="restart_placed.liggghts" if a.start=="placed" else "restart_shear.liggghts"
    link=cdir/"DEM"/"restarts"/lname
    if src_restart.is_file():
        os.symlink(src_restart, link); rlink=f"linked -> {src_restart}"
    else:
        print(f"!! missing restart {src_restart}"); return 1

    # --- CFD side ---
    cfd=cdir/"CFD"; (cfd/"0").mkdir(parents=True)
    (cfd/"system").mkdir(); (cfd/"constant").mkdir()
    tut=Path(cfg["TUTORIAL_CFD"])
    for f in ["fvSchemes","fvSolution","decomposeParDict"]:
        shutil.copy(tut/"system"/f, cfd/"system"/f)
    for f in ["g","liggghtsCommands","turbulenceProperties","RASProperties"]:
        if (tut/"constant"/f).is_file(): shutil.copy(tut/"constant"/f, cfd/"constant"/f)
    # fluid gravity = DEM gravity (zero if the DEM script has none): no one-sided buoyancy
    try:
        gvec=C.dem_gravity(txt) or (0.0,0.0,0.0)
    except C.CaseError as e:
        print(f"!! {e}"); shutil.rmtree(cdir); return 1
    gf=cfd/"constant"/"g"
    gf.write_text(re.sub(r"value\s*\([^)]*\)\s*;", "value           ( %r %r %r );" % gvec, gf.read_text()))
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
    cp=_tmpl.read_text() if _tmpl.is_file() else COUPLING
    cp, nv = re.subn(r"voidageCorrection\s+\S+?\s*;", f"voidageCorrection       {a.voidage};", cp)
    if nv != 1:
        print("!! couplingProperties template: expected one voidageCorrection entry"); return 1
    # void fraction: periodic directions from the mesh written below (cyclic patches)
    per=C.mesh_periodic(BLOCKMESH)
    cp, nd = re.subn(r"(dividedSuperquadricProps\s*\{)", "\\1\n    periodicDirections  (%d %d %d);" % per, cp)
    if nd != 1:
        print("!! couplingProperties template: no dividedSuperquadricProps block"); return 1
    (cfd/"constant"/"couplingProperties").write_text(cp)
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
    _mkfields(cfd/"0", rho_f, top=up*uw, bottom=lo*uw)

    # resumable runner
    (cdir/"run.sh").write_text(RUNNER); os.chmod(cdir/"run.sh",0o755)

    # static audit of what was just written; an invalid case is renamed, not left runnable
    R=C.audit((cdir/"DEM"/"in.liggghts_run").read_text(), gf.read_text(),
              (cfd/"0"/"U").read_text(), (cfd/"constant"/"couplingProperties").read_text(),
              expected_axes=(s["a"],s["b"],s["c"]), blockmesh_txt=BLOCKMESH)
    for i,st,m in R: print(f"  [{st:4}] {i} {m}")
    if any(st=="FAIL" for _,st,_ in R):
        bad=cdir.with_name(cdir.name+".INVALID"); cdir.rename(bad)
        print(f"!! audit failed: case moved to {bad}"); return 1

    print(f"created {cdir}")
    print(f"  shape={a.shape} (c/a={s['ca']})  St={a.St} (isolated particle)  nu_f={nu_f:.3e}  tau_p={tau_p:.1f}")
    if a.voidage=="on" and "nu_core" in s:
        eps=1.0-s["nu_core"]                      # Stokes regime: K ~ eps^(1-chi), chi -> 3.7
        print(f"  with the voidage correction the drag in the core is ~{eps**(1-3.7):.1f}x the isolated value:"
              f" effective St ~ {a.St*eps**(3.7-1):.3g}")
    print(f"  template: {tmpl}")
    print(f"  start: {a.start}")
    print(f"  restart: {rlink}")
    print(f"  -> sqfoam run {cdir}")

def _mkfields(z, rho_f, top, bottom):
    def w(n,cls,dim,intern,wb,wt):
        (z/n).write_text(
f"""FoamFile {{ version 2.0; format ascii; class {cls}; object {n}; }}
dimensions {dim}; internalField uniform {intern};
boundaryField {{ cycXa {{ type cyclic; }} cycXb {{ type cyclic; }}
 cycYa {{ type cyclic; }} cycYb {{ type cyclic; }}
 wallBottom {{ {wb} }} wallTop {{ {wt} }} }}
""")
    w("U","volVectorField","[0 1 -1 0 0 0 0]","(0 0 0)",
      f"type fixedValue; value uniform ({bottom!r} 0 0);",f"type fixedValue; value uniform ({top!r} 0 0);")
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
