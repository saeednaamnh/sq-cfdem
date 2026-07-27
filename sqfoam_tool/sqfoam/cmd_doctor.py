import os, shutil, subprocess
from pathlib import Path
from config import load, save, CONF, DEFAULTS

def _has(cmd): return shutil.which(cmd) is not None

def cmd_doctor(a):
    cfg = load()
    if not CONF.is_file():
        save(cfg); print(f"created {CONF} (edit paths there if needed)\n")

    print("="*60); print("  sqfoam doctor — environment check"); print("="*60)

    ok=True
    def chk(label, cond, hint=""):
        nonlocal ok
        mark = "OK " if cond else "XX "
        if not cond: ok=False
        print(f"  [{mark}] {label}" + (f"   -> {hint}" if (not cond and hint) else ""))
        return cond

    # OpenFOAM
    chk("OpenFOAM sourced ($WM_PROJECT_VERSION)",
        bool(os.environ.get("WM_PROJECT_VERSION")),
        "source /opt/openfoam6/etc/bashrc")
    # CFDEM env
    chk("CFDEM env ($CFDEM_PROJECT_DIR)",
        bool(os.environ.get("CFDEM_PROJECT_DIR")),
        "source the CFDEM bashrc block")
    # LIGGGHTS binary
    lmp = Path(cfg["CFDEM_LIGGGHTS_SRC_DIR"])/"lmp_auto"
    chk(f"LIGGGHTS lmp_auto ({lmp})", lmp.is_file(),
        "build LIGGGHTS: make auto")
    # our solver
    chk(f"solver {cfg['SOLVER']} on PATH", _has(cfg["SOLVER"]),
        "wmake the cfdemSolverPisoSQ application")
    # our library
    libhits=[]
    for base in [os.environ.get("FOAM_USER_LIBBIN",""),
                 os.environ.get("CFDEM_LIB_DIR","")]:
        if base and (Path(base)/"libsqCfdem.so").is_file(): libhits.append(base)
    chk("libsqCfdem.so built", bool(libhits),
        "cd src/sqCfdem && wmake libso")
    # dry restarts (for grid / dry-consistent runs)
    dr=Path(cfg["DRY_RESTARTS"])
    chk(f"dry restarts dir ({dr})", dr.is_dir(),
        "point DRY_RESTARTS in ~/.sqfoam.conf at your placement restarts")
    # extraction pipeline
    ep=Path(cfg["EXTRACT_PIPELINE"])/"extract_coefficients_v5.py"
    chk(f"extraction pipeline ({ep})", ep.is_file(),
        "set EXTRACT_PIPELINE in ~/.sqfoam.conf")
    # postprocess script
    pp=Path(cfg["POSTPROCESS_PY"])
    chk(f"postprocess_shear.py ({pp})", pp.is_file(),
        "set POSTPROCESS_PY in ~/.sqfoam.conf")

    print("="*60)
    print("  ALL GREEN — ready to run" if ok else
          "  some checks failed — fix the hints above, edit ~/.sqfoam.conf")
    print("="*60)
