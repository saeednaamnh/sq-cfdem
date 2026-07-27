# Building the sq-cfdem stack (one-time)

Ubuntu 18.04/20.04/22.04. Captures every fix from the original bring-up.

## 1. System packages
    sudo apt-get install -y build-essential flex bison zlib1g-dev \
      libboost-system-dev libboost-thread-dev libopenmpi-dev openmpi-bin \
      gnuplot libreadline-dev libncurses-dev libxt-dev cmake git

## 2. OpenFOAM-6
Ubuntu 18.04: packaged (fast):
    sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
    sudo add-apt-repository http://dl.openfoam.org/ubuntu
    sudo apt-get update && sudo apt-get install -y openfoam6
Otherwise build from source (version-6 branch).

## 3. LIGGGHTS-PUBLIC + CFDEMcoupling-PUBLIC
    mkdir -p ~/CFDEM ~/LIGGGHTS
    cd ~/CFDEM   && git clone https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC.git
    cd ~/LIGGGHTS && git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git
Add the CFDEM env block to ~/.bashrc (see repo docs/ubuntu1804_walkthrough.md),
open a new shell, then:
    cfdemSysTest

## 4. Build order (critical)
    # LIGGGHTS as library + executable, superquadrics ON:
    cd $CFDEM_LIGGGHTS_SRC_DIR
    grep USE_SUPERQUADRICS MAKE/Makefile.user   # must be "ON"
    make auto -j$(nproc)
    make makeshlib && make -f Makefile.shlib auto -j$(nproc)   # <-- shared lib!
    # CFDEM libs + solvers:
    cfdemCompCFDEMsrc && cfdemCompCFDEMsol
    # our library + solver:
    cd $SQCFDEM_DIR/src/sqCfdem && wmake libso
    cd $SQCFDEM_DIR/applications/cfdemSolverPisoSQ && wmake

## 5. Install the sqfoam driver
    cd sqfoam_tool && ./install_sqfoam.sh
    sqfoam doctor      # verifies all of the above

Gotchas captured here:
- LIGGGHTS needs BOTH `make auto` (exe) AND the shlib build (CFDEM links -llmp_auto)
- CFDEM targets OpenFOAM-6 (not 5.x, not ESI)
- solver must be compiled with -DMS -Dsuperquadrics_flag (the SQ Make/options does this)
- couplingProperties must use `particleShapeType superquadric` + `voidFractionModel divided`
