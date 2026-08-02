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

    # ---- DRY-IDENTITY GUARD (catch silently-failed coupling) ----
    import hashlib
    def _md5(fp):
        h=hashlib.md5()
        with open(fp,'rb') as f:
            for chunk in iter(lambda: f.read(1<<20), b''): h.update(chunk)
        return h.hexdigest()
    run_dump=os.path.join(shear,"dump.orient")
    dry_dump=os.path.join(cfg["DRY_RESTARTS"], shape, "results","shear","dump.orient")
    if os.path.isfile(run_dump) and os.path.isfile(dry_dump):
        if _md5(run_dump)==_md5(dry_dump):
            print("="*60)
            print("  !!! DRY-IDENTITY GUARD TRIPPED !!!")
            print(f"  This run's dump.orient is BYTE-IDENTICAL to the dry")
            print(f"  reference for {shape}. The fluid coupling had NO effect")
            print(f"  -- this is NOT a valid wet result. Do not trust it.")
            print(f"  Likely cause: coupling failed to transfer force.")
            print("="*60)
            return
        else:
            print(f"  [dry-guard OK] dump differs from dry reference -- real coupling")
    # ------------------------------------------------------------

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
