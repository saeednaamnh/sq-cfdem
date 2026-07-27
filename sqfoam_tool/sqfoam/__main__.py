#!/usr/bin/env python3
"""
sqfoam - unified command-line driver for the sq-cfdem superquadric CFD-DEM solver.

Turns the multi-step install/build/run/postprocess workflow into simple commands:

  sqfoam doctor                     check the environment / build status
  sqfoam new  --shape prolate_3 --St 100 [--name myrun]   generate a case
  sqfoam run  <case_dir>            launch (resumable) coupled run in background
  sqfoam status <case_dir>          is it running? how far?
  sqfoam extract <case_dir>         run coefficient extraction (phi, psi, i1...)
  sqfoam viz <case_dir>             generate ParaView ellipsoid VTK/PVD
  sqfoam grid --full                generate the whole St x aspect-ratio grid

Config lives in ~/.sqfoam.conf (auto-created by `sqfoam doctor`).
"""
import sys, os, argparse
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from cmd_doctor  import cmd_doctor
from cmd_new     import cmd_new
from cmd_run     import cmd_run, cmd_status
from cmd_extract import cmd_extract
from cmd_viz     import cmd_viz
from cmd_grid    import cmd_grid

def main():
    p = argparse.ArgumentParser(
        prog="sqfoam",
        description="Superquadric CFD-DEM workflow driver (sq-cfdem)")
    sub = p.add_subparsers(dest="cmd")
    sub.required = True

    sub.add_parser("doctor", help="check environment & build status")

    n = sub.add_parser("new", help="generate a single case")
    n.add_argument("--shape", required=True,
                   choices=["sphere","prolate_2","prolate_3","prolate_5"])
    n.add_argument("--St", type=float, required=True, help="Stokes number")
    n.add_argument("--name", default=None, help="case dir name (default: shape_StN)")
    n.add_argument("--steps", type=int, default=5000000, help="DEM shear steps")
    n.add_argument("--out", default=None, help="parent dir (default: cwd)")

    r = sub.add_parser("run", help="launch resumable coupled run (background)")
    r.add_argument("case"); r.add_argument("--fg", action="store_true",
                   help="run in foreground instead of background")

    s = sub.add_parser("status", help="is the run alive & how far?")
    s.add_argument("case")

    e = sub.add_parser("extract", help="extract transport coefficients")
    e.add_argument("case")
    e.add_argument("--truth", default=None, help="Table I json for comparison")

    v = sub.add_parser("viz", help="generate ParaView ellipsoid VTK/PVD")
    v.add_argument("case")
    v.add_argument("--which", default="shear", choices=["shear","placement"])
    v.add_argument("--res", type=int, default=8, help="ellipsoid mesh resolution")

    g = sub.add_parser("grid", help="generate the full parameter grid")
    g.add_argument("--full", action="store_true")
    g.add_argument("--out", default="physics_grid")

    a = p.parse_args()
    {"doctor":cmd_doctor, "new":cmd_new, "run":cmd_run, "status":cmd_status,
     "extract":cmd_extract, "viz":cmd_viz, "grid":cmd_grid}[a.cmd](a)

if __name__ == "__main__":
    main()
