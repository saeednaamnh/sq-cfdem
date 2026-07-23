#!/bin/bash
# Build order matters: LIGGGHTS lib (superquadrics ON) -> CFDEM -> sqCfdem
set -e
# 1) LIGGGHTS as coupled library with superquadrics
cd $CFDEM_LIGGGHTS_SRC_DIR
grep -q 'USE_SUPERQUADRICS *= *"ON"' MAKE/Makefile.user 2>/dev/null || \
  { echo '>>> EDIT MAKE/Makefile.user: USE_SUPERQUADRICS = "ON" — then rerun'; exit 1; }
make auto -j $(nproc)
# 2) CFDEMcoupling libs + solvers (also compiles the LIGGGHTS coupling fix)
cfdemCompCFDEMsrc
cfdemCompCFDEMsol
# 3) our library
cd $SQCFDEM_DIR/src/sqCfdem && wmake libso
# 4) smoke test
cfdemSysTest
