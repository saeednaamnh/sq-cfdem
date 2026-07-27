import sys, os, subprocess
from pathlib import Path
from config import load, SHAPES

def cmd_extract(a):
    cfg=load(); case=Path(a.case).resolve()
    shear=case/"DEM"/"results"/"shear"
    if not (shear/"dump.orient").is_file():
        print(f"!! no dump.orient in {shear}"); return
    # infer shape from dir name
    shape=None
    for s in SHAPES:
        if case.name.startswith(s): shape=s; break
    if not shape:
        print("!! cannot infer shape from case name; rename to <shape>_..."); return

    pipe=Path(cfg["EXTRACT_PIPELINE"])
    v6cmp=pipe/"v6_compare.py"
    if v6cmp.is_file():
        truth=a.truth or str(pipe/"figures"/"coefficients_v6.json")
        cmd=["python3",str(v6cmp),"--v6-shear",str(shear),
             "--case",shape,"--truth",truth]
        print("running:"," ".join(cmd)); print()
        subprocess.run(cmd, cwd=str(pipe))
    else:
        print(f"!! v6_compare.py not found in {pipe}")
        print("   copy v6_compare.py there, or set EXTRACT_PIPELINE in ~/.sqfoam.conf")
