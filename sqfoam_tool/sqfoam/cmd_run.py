import os, subprocess, glob
from pathlib import Path
from config import load

def cmd_run(a):
    cfg=load(); case=Path(a.case).resolve()
    runner=case/"run.sh"
    if not runner.is_file():
        print(f"!! no run.sh in {case} (generate with sqfoam new)"); return
    env=dict(os.environ, NCORES=cfg["NCORES"])
    log=case/"run.log"
    if a.fg:
        subprocess.run(["bash",str(runner)], cwd=case, env=env)
    else:
        with open(log,"a") as lf:
            p=subprocess.Popen(["nohup","bash",str(runner)],
                cwd=case, env=env, stdout=lf, stderr=subprocess.STDOUT,
                preexec_fn=os.setpgrp)
        print(f"launched in background (PID {p.pid}); log -> {log}")
        print(f"  monitor: sqfoam status {case}")

def cmd_status(a):
    case=Path(a.case).resolve()
    # alive?
    out=subprocess.run(["pgrep","-af","cfdemSolverPisoSQ"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout
    alive=str(case) in out or "cfdemSolverPisoSQ" in out
    gt=case/"DEM"/"results"/"shear"/"granular_temp.txt"
    step="(no data yet)"; T="?"
    if gt.is_file():
        lines=[l for l in gt.read_text().splitlines() if l and not l.startswith("#")]
        if lines:
            parts=lines[-1].split()
            step=parts[0]; T=parts[-1] if len(parts)>=3 else "?"
    orient=case/"DEM"/"results"/"shear"/"dump.orient"
    frames=0
    if orient.is_file():
        frames=sum(1 for l in open(orient) if l.strip()=="ITEM: TIMESTEP")
    print(f"case: {case.name}")
    print(f"  process: {'RUNNING' if alive else 'not running'}")
    print(f"  last step: {step}   T={T}")
    print(f"  orient frames: {frames}" +
          ("  (enough for extraction)" if frames>=100 else "  (need ~100+ for converged fit)"))
    if not alive:
        print(f"  resume: sqfoam run {case}")
