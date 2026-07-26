#!/bin/bash
set -e
cd "$(dirname "$0")"
if [ ! -f DEM/restarts/restart_placed.liggghts ]; then
  echo "=== PLACEMENT (dry, your in_place.lmp; ~their original runtime) ==="
  cd DEM
  mpirun -np 4 $CFDEM_LIGGGHTS_SRC_DIR/lmp_auto < in.liggghts_init_place.lmp 2>&1 | tail -5
  cd ..
fi
echo "=== COUPLED SHEAR (dry limit) ==="
cd CFD
blockMesh > ../log_mesh.txt 2>&1
decomposePar -force > ../log_decomp.txt 2>&1
mpirun -np 4 cfdemSolverPisoSQ -parallel 2>&1 | tee ../log_solve.txt | grep -E "satellite points|SHEAR|Vol frac|Shear rate|FATAL" | head -20
cd ..
echo "=== outputs in DEM/results/shear/ (your dry pipeline format) ==="
ls DEM/results/shear 2>/dev/null
