#!/bin/bash
# Prerequisite check for the sq-cfdem stack (Ubuntu/Debian assumed)
echo "== compilers/tools =="
for t in g++ gcc mpirun mpicxx cmake git flex make; do
  command -v $t >/dev/null && echo "OK  $t" || echo "MISSING  $t"
done
echo "== LIGGGHTS =="
echo "CFDEMcoupling compiles against the LIGGGHTS *source tree*, not the binary."
echo "Locate your LIGGGHTS-PUBLIC 3.8.0 source; it will be rebuilt as a library."
echo "== suggested apt packages =="
echo "sudo apt-get install build-essential flex bison zlib1g-dev libboost-system-dev \
 libboost-thread-dev libopenmpi-dev openmpi-bin gnuplot libreadline-dev \
 libncurses-dev libxt-dev cmake git"
