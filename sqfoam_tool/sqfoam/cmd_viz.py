import subprocess
from pathlib import Path
from config import load

def cmd_viz(a):
    cfg=load(); case=Path(a.case).resolve()
    pp=Path(cfg["POSTPROCESS_PY"])
    if not pp.is_file():
        print(f"!! postprocess_shear.py not found at {pp}")
        print("   set POSTPROCESS_PY in ~/.sqfoam.conf"); return
    dump_dir=case/"DEM"/"results"/a.which
    pattern="dump.full" if a.which=="shear" else "dump.place"
    if not dump_dir.is_dir():
        print(f"!! {dump_dir} not found"); return
    vtk_dir=case/"DEM"/"results"/f"vtk_{a.which}"
    cmd=["python3",str(pp),"--pvd","--dump_dir",str(dump_dir),
         "--vtk_dir",str(vtk_dir),"--pattern",pattern,
         "--vtk_res",str(a.res)]
    print("running:"," ".join(cmd)); print()
    subprocess.run(cmd)
    pvd=vtk_dir/"collection.pvd"
    if pvd.is_file():
        print(f"\nParaView: open {pvd}")
        print("  (glyphs are real ellipsoid meshes; colour by "
              "orientation_angle / speed / particle_type)")
