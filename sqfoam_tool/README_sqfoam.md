# sqfoam — one-command driver for sq-cfdem

After the stack is built (docs/INSTALL_STACK.md) and `./install_sqfoam.sh` run:

    sqfoam doctor                          # check everything is in place
    sqfoam new --shape prolate_3 --St 100  # generate a case
    sqfoam run  prolate_3_St100            # launch (resumable, background)
    sqfoam status prolate_3_St100          # alive? how far? how many frames?
    sqfoam extract prolate_3_St100         # phi, psi, i1... vs Table I
    sqfoam viz  prolate_3_St100            # ParaView ellipsoid VTK/PVD
    sqfoam grid --full                     # generate all 24 St x aspect cases

Config: ~/.sqfoam.conf (paths to LIGGGHTS, dry restarts, extraction pipeline,
postprocess_shear.py). `sqfoam doctor` creates it; edit if your paths differ.

What each case contains: a coupling-transformed copy of your in_shear.lmp,
a matched OpenFOAM case (mesh, fields, dictionaries), a symlinked placement
restart (particle-for-particle with Table I), and a resumable run.sh.

Runs survive interrudavep: `sqfoam run <case>` resumes from the last shear
checkpoint automatically. One coupled job per 4 cores — don't run two at once.
