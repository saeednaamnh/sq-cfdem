# LIGGGHTS 3.8.0 side — VERIFIED against source: NO PATCH NEEDED

Inspection of LIGGGHTS-PUBLIC src/fix_cfd_coupling_force.cpp confirms the
entire DEM-side superquadric exchange already exists:

- push to CFD (when superquadric mode is active):
  volume (scalar-atom), area (scalar-atom), shape (vector-atom),
  blockiness (vector2D-atom), quaternion (quaternion-atom),
  plus x, v, radius; omega is pushed when transfer_torque = yes.
- pull from CFD: dragforce AND hdtorque (vector-atom); the fix adds
  hdtorque[i] into torque[i] each step (vectorAdd3D in post_force),
  i.e. hydrodynamic torque is applied alongside contact torque with
  no integrator modification required.
- use_torque_ defaults to true; keyword `transfer_torque yes|no`.
- superquadric mode auto-enables torque transfer (use_superquadric_
  sets use_torque_ = true).
- USE_SUPERQUADRICS = "ON" is present in MAKE/Makefile.user_default
  (line ~25) — copy to Makefile.user and keep ON.

Coupled-run DEM input script (template):
  fix cfd  all couple/cfd couple_every 100 mpi
  fix fdrg all couple/cfd/force transfer_torque yes

Remaining LIGGGHTS-side work: NONE beyond rebuilding as a library
(`make auto`) from the same source tree with superquadrics ON.
