#!/bin/bash
set -e
cd "$(dirname "$0")"
cd DEM
mpirun -np 4 $CFDEM_LIGGGHTS_SRC_DIR/lmp_auto < in.liggghts_init > ../log_dem_init.txt 2>&1
cd ..
grep -i "inserted" log_dem_init.txt | head -2
tail -1 DEM/quat_check.txt
cd CFD
blockMesh > ../log_mesh.txt 2>&1
decomposePar -force > ../log_decomp.txt 2>&1
mpirun -np 4 cfdemSolverPisoSQ -parallel 2>&1 | tee ../log_solve.txt | grep -E "satellite points|^Time = 1[02]$|^Time = [2468]$|FATAL" 
cd ..
echo "=== quat series frames: $(grep -c TIMESTEP DEM/post/quat_v3.txt) ==="
