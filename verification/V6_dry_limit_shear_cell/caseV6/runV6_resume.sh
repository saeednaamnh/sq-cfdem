#!/bin/bash
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
grep -n "read_restart" "$DEMIN"
cd CFD
if [ ! -d processor0 ]; then
    blockMesh > ../log_mesh.txt 2>&1
    decomposePar -force > ../log_decomp.txt 2>&1
fi
mpirun -np 4 cfdemSolverPisoSQ -parallel 2>&1 | tee -a ../log_solve.txt \
  | grep -E "satellite points|Shear rate|FATAL|Courant" | head -20
cd ..
