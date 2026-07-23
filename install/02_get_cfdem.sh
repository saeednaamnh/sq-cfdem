#!/bin/bash
# Fetch CFDEMcoupling-PUBLIC and LIGGGHTS-PUBLIC source side by side.
set -e
mkdir -p $HOME/CFDEM && cd $HOME/CFDEM
[ -d CFDEMcoupling-PUBLIC ] || git clone https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC.git
mkdir -p $HOME/LIGGGHTS && cd $HOME/LIGGGHTS
# If your existing LIGGGHTS install still has its source tree, symlink it here
# instead of recloning, so DEM binary == coupled DEM library:
[ -d LIGGGHTS-PUBLIC ] || git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git
echo "IMPORTANT: LIGGGHTS must be REBUILT as a library with superquadrics ON:"
echo "  edit \$CFDEM_LIGGGHTS_SRC_DIR/MAKE/Makefile.user:  USE_SUPERQUADRICS = \"ON\""
echo "  (sets -DSUPERQUADRIC_ACTIVE_FLAG -DNONSPHERICAL_ACTIVE_FLAG, which the"
echo "   CFDEM superquadric dispatch and sqCfdem both depend on)"
