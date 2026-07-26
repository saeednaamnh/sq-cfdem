#!/bin/bash
# minimal V1 runner: no cleanup traps
set -e
cd "$(dirname "$0")"

# 1) DEM init (writes restart + quat_check)
cd DEM
mpirun -np 4 $CFDEM_LIGGGHTS_SRC_DIR/lmp_auto < in.liggghts_init > ../log_dem_init.txt 2>&1
cd ..
grep -i "inserted" log_dem_init.txt | head -2
tail -1 DEM/quat_check.txt

# 2) mesh + decompose (fresh every time)
cd CFD
blockMesh > ../log_mesh.txt 2>&1
decomposePar -force > ../log_decomp.txt 2>&1
# 3) coupled solve
mpirun -np 4 cfdemSolverPisoSQ -parallel 2>&1 | tee ../log_solve.txt | grep -E "satellite points|registering field 'quaternion'|^Time =|FATAL" | tail -8
cd ..
echo "=== forces tail ==="
tail -3 DEM/post/forces_v1.txt
