# Known Issues

## CFD field output not written during coupled runs
The DEM-driven coupling advances particles inside particleCloud.evolve(),
but OpenFOAM's runTime outer loop does not iterate, so runTime.write() is
never reached and no CFD time directories are written (processor*/[time]/).
DEM output (dump.orient, dump.full, stats) is unaffected and correct.

Impact: cannot visualize CFD continuum fields (U, voidfraction, p) in
ParaView. Particle visualization (sqfoam viz) works fully.

Does NOT affect the science: all transport-coefficient extraction uses the
DEM orientation data, which is written correctly.

Future fix: patch cfdemSolverPiso.C to force runTime iteration / explicit
field write at coupling intervals, then rebuild the solver.
