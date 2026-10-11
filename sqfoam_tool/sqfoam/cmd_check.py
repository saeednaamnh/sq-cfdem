import os
from pathlib import Path
from config import load, SHAPES
import case_checks as C

def _read_T(path, step=None):
    if not os.path.isfile(path): return None
    last = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"): continue
            p = line.split()
            if len(p) < 3: continue
            if step is None:
                last = (p[0], p[-1])
            elif p[0] == str(step):
                return (p[0], p[-1])
    return last

def cmd_check(a):
    cfg = load(); case = Path(a.case).resolve()
    shape = getattr(a, "shape", None)
    if not shape:                       # from the particle shape in the case's own DEM script
        dem = case/"DEM"/"in.liggghts_run"
        shape = C.shape_of(dem.read_text(), SHAPES) if dem.is_file() else None
    if not shape:
        print("!! cannot infer the shape from DEM/in.liggghts_run; pass --shape prolate_3"); return

    run_T = case/"DEM"/"results"/"shear"/"granular_temp.txt"
    dry_T = Path(cfg["DRY_RESTARTS"])/cfg["TEMPLATE_DIRS"][shape]/"results"/"shear"/"granular_temp.txt"

    print("="*64)
    print("  CONFIGURATION AUDIT — %s  (shape %s)" % (case.name, shape))
    print("="*64)
    files = [case/"DEM"/"in.liggghts_run", case/"CFD"/"constant"/"g",
             case/"CFD"/"0"/"U", case/"CFD"/"constant"/"couplingProperties"]
    missing = [str(f) for f in files if not f.is_file()]
    if missing:
        print("  cannot audit, missing: %s" % ", ".join(missing))
    else:
        s = SHAPES[shape]
        bm = case/"CFD"/"system"/"blockMeshDict"
        R = C.audit(*(f.read_text() for f in files), expected_axes=(s["a"], s["b"], s["c"]),
                    blockmesh_txt=bm.read_text() if bm.is_file() else None)
        for i, st, m in R:
            print("  [%-4s] %s %s" % (st, i, m))
        # A11: no restart may be a link out of the case (a checkpoint written through it
        # would overwrite another case's restart, e.g. the dry baseline's final state)
        links = [p for p in (case/"DEM"/"restarts").glob("*") if p.is_symlink()]
        outside = [p for p in links if case.resolve() not in p.resolve().parents]
        if outside:
            R.append(("A11", "FAIL", "restart(s) linked out of the case: %s; replace each link by a copy"
                      " (cp --remove-destination \"$(readlink -f X)\" X)" % ", ".join("%s -> %s" % (p.name, p.resolve()) for p in outside)))
        else:
            R.append(("A11", "OK", "no restart linked out of the case"))
        i, st, msg = R[-1]
        print("  [%-4s] %s %s" % (st, i, msg))
        if any(st == "FAIL" for _, st, _ in R):
            print("\n  *** do not run this case until every FAIL is fixed ***")
    print()
    print("  COUPLING LIVENESS (wet vs dry granular temperature)")

    if not run_T.is_file():
        print("  no granular_temp.txt yet"); return
    if not dry_T.is_file():
        print("  !! no dry reference at %s" % dry_T); return

    last = _read_T(str(run_T))
    if last is None:
        print("  no data rows yet"); return
    step, T_wet = last
    ref = _read_T(str(dry_T), step=step)

    print("  step %s" % step)
    print("    this run (wet) : T = %s" % T_wet)
    if ref is None:
        print("    dry reference  : (no matching step)")
        print("\n  INCONCLUSIVE at this step."); return
    T_dry = ref[1]
    print("    dry reference  : T = %s" % T_dry)
    print()
    if T_wet == T_dry:
        print("  *** BROKEN: identical to dry — fluid doing NO work. ***")
        print("  Check the DEM script's 'run ${N_STEPS}' is disabled.")
    else:
        try: d = abs(float(T_wet)-float(T_dry))/float(T_dry)*100
        except Exception: d = float("nan")
        print("  OK: diverged from dry by %.3f%% — real fluid-solid coupling." % d)
    print("="*64)
