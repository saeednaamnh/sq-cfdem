import os
from pathlib import Path
from config import load, SHAPES

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
    if not shape:
        for s in SHAPES:
            if case.name.startswith(s): shape = s; break
    if not shape:
        print("!! cannot infer shape; pass --shape prolate_3"); return

    run_T = case/"DEM"/"results"/"shear"/"granular_temp.txt"
    dry_T = Path(cfg["DRY_RESTARTS"])/shape/"results"/"shear"/"granular_temp.txt"

    print("="*64)
    print("  COUPLING CHECK — %s  (shape %s)" % (case.name, shape))
    print("="*64)

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
