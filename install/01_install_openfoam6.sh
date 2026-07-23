#!/bin/bash
# OpenFOAM-6 (openfoam.org) from source — required by current CFDEMcoupling-PUBLIC.
# Upstream README: repo is compatible with OpenFOAM-6; older versions
# "may prove difficult". Do NOT use ESI v-releases (v1806 etc.).
set -e
mkdir -p $HOME/OpenFOAM && cd $HOME/OpenFOAM
[ -d OpenFOAM-6 ]   || git clone --depth 1 --branch version-6 https://github.com/OpenFOAM/OpenFOAM-6.git
[ -d ThirdParty-6 ] || git clone --depth 1 --branch version-6 https://github.com/OpenFOAM/ThirdParty-6.git
source $HOME/OpenFOAM/OpenFOAM-6/etc/bashrc || true
cd ThirdParty-6 && ./Allwmake -j $(nproc)
cd ../OpenFOAM-6 && ./Allwmake -j $(nproc) 2>&1 | tee build.log
echo "Sanity: run 'simpleFoam -help' and 'foamInstallationTest'"
# NOTE (Ubuntu 22.04+/gcc>=10): OF-6 is a 2018 release; if the build trips on
# newer gcc, install g++-9 (apt) and export WM_CC=gcc-9 WM_CXX=g++-9 before
# Allwmake, or build inside an Ubuntu 20.04 container.
