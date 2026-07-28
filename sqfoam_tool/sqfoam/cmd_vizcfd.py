import os, subprocess, glob
from pathlib import Path
from config import load

def cmd_vizcfd(a):
    cfg=load(); case=Path(a.case).resolve()
    cfd=case/"CFD"
    if not cfd.is_dir():
        print(f"!! no CFD/ directory in {case}"); return

    procdirs=sorted(glob.glob(str(cfd/"processor*")))
    foam=cfd/"case.foam"

    if a.decomposed:
        # open decomposed directly (no reconstruct) - fastest
        foam.touch()
        print(f"prepared DECOMPOSED case for ParaView (no reconstruct).")
        print(f"\nParaView:")
        print(f"  paraview {foam} &")
        print(f"  -> in Properties panel set 'Case Type' = Decomposed Case, click Apply")
    else:
        # reconstruct parallel -> single, then stub
        if not procdirs:
            print("!! no processor*/ dirs found - was the run parallel?")
            print("   (if it was serial, just: touch CFD/case.foam and open it)")
            return
        print(f"reconstructing {len(procdirs)} processor dirs -> single case ...")
        r=subprocess.run(["reconstructPar","-case",str(cfd)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         universal_newlines=True)
        tail="\n".join(r.stdout.splitlines()[-6:])
        print(tail)
        if r.returncode!=0:
            print("\n!! reconstructPar failed - try:  sqfoam viz-cfd <case> --decomposed")
            return
        foam.touch()
        print(f"\nParaView:")
        print(f"  paraview {foam} &")
        print(f"  -> click Apply, then Color by:  voidfraction  (or U, p, Ksl)")

    # list what fields/times are available
    times=sorted([d.name for d in cfd.iterdir()
                  if d.is_dir() and d.name.replace('.','').isdigit()],
                 key=lambda x: float(x))
    if times:
        print(f"\n  time steps available: {times[0]} ... {times[-1]}  ({len(times)} total)")
    print(f"\n  TIP: to overlay the PARTICLES too, in the same ParaView window also")
    print(f"       open  {case}/DEM/results/vtk_shear/collection.pvd")
    print(f"       (generate it first with:  sqfoam viz {case.name})")
    print(f"\n  NOTE: unresolved CFD-DEM uses cells bigger than particles, so the")
    print(f"        continuum fields look blocky by design (volume-averaged), not smooth.")
