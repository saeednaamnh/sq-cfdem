"""Unit tests of sqfoam_tool/sqfoam/case_checks.py on a minimal synthetic DEM template."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "sqfoam_tool", "sqfoam"))
import case_checks as C
T = """variable    a           equal   0.2371
variable    b           equal   0.2371
variable    c           equal   0.7114
variable    wall_vel    equal   1.0
fix  integrate layerGroup nve/superquadric integration_scheme 1
fix  move_up upperGroup move linear ${wall_vel} 0 0
fix  move_lo lowerGroup move linear -${wall_vel} 0 0
run  5000000
run 0
"""
fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print("  FAIL:", msg)
check(C.integrator(T)[1] == "layerGroup", "integrator group")
check(C.wall_velocity(T) == 1.0 and C.wall_motion_signs(T) == (1, -1), "wall motion")
check(C.dem_gravity(T) is None and C.active_runs(T) == ["run  5000000"], "gravity / runs")
new, g = C.inject_coupling(T.replace("layerGroup", "free_particles"))
check(g == "free_particles" and C.coupling_group(new) == "free_particles", "coupling on the integrator's group")
try:
    C.inject_coupling(T.replace("nve/superquadric", "nve/sphere")); check(False, "no error without integrator")
except C.CaseError:
    pass
check(C.dem_gravity("fix g all gravity 9.81 vector 0 0 -1\n") == (0.0, 0.0, -9.81), "gravity vector")
U = "wallBottom { type fixedValue; value uniform (-1 0 0); } wallTop { type fixedValue; value uniform (1 0 0); }"
G = "value           ( 0 0 0 );"
CP = ("forceModels ( HolzerSommerfeldDrag spheroidRotationTorque );\n"
      "HolzerSommerfeldDragProps { voidageCorrection on; }\ncouplingFlag cfdCoupled;")
dem, flag = C.add_coupling_flag(new.replace("run  5000000", "# run disabled"))
check(flag == "cfdCoupled", "coupling flag added")
blk = dem[dem.index("property/atom cfdCoupled"):dem.index("couple/cfd")]
check("run             0" in blk and blk.count("property/atom cfdCoupled 0.0") == 2, "flag block: run 0 and both walls set to 0")
R = dict((i, s) for i, s, _ in C.audit(dem, G, U, CP, expected_axes=(0.2371, 0.2371, 0.7114)))
check(all(s == "OK" for s in R.values()), "clean case audits OK: %s" % R)
R = dict((i, s) for i, s, _ in C.audit(T, "value ( 0 0 -9.81 );", U.replace("(1 0 0)", "(2 0 0)"),
                                        CP.replace("voidageCorrection on;", ""), expected_axes=(0.2, 0.2, 1.0)))
check([R[k] for k in ("A1", "A2", "A3", "A4", "A5", "A8")] == ["FAIL"] * 6, "broken case: %s" % R)
R = dict((i, s) for i, s, _ in C.audit(dem, G, U, CP.replace("couplingFlag cfdCoupled;", ""), blockmesh_txt=None))
check(R["A10"] == "FAIL", "A10 must fail without couplingFlag")
R = dict((i, s) for i, s, _ in C.audit(new, G, U, CP))
check(R["A10"] == "FAIL", "A10 must fail when the walls have no flag")
R = dict((i, s) for i, s, _ in C.audit(dem.replace("run             0\n", "", 1), G, U, CP))
check(R["A10"] == "FAIL", "A10 must fail without run 0 between the fix and the set commands")
BM = """vertices ( (0 0 0)(15 0 0)(15 10 0)(0 10 0)(0 0 10)(15 0 10)(15 10 10)(0 10 10) );
boundary
( cycXa { type cyclic; neighbourPatch cycXb; faces ((0 4 7 3)); }
  cycXb { type cyclic; neighbourPatch cycXa; faces ((1 2 6 5)); }
  cycYa { type cyclic; neighbourPatch cycYb; faces ((0 1 5 4)); }
  cycYb { type cyclic; neighbourPatch cycYa; faces ((3 7 6 2)); }
  wallBottom { type wall; faces ((0 3 2 1)); }
  wallTop    { type wall; faces ((4 5 6 7)); } );"""
check(C.mesh_periodic(BM) == (1, 1, 0), "periodic directions from the mesh geometry")
check(C.mesh_periodic(BM.replace("type cyclic; neighbourPatch cycYb", "type wall; neighbourPatch cycYb")
                        .replace("type cyclic; neighbourPatch cycYa", "type wall; neighbourPatch cycYa")) == (1, 0, 0), "x only")
VF = lambda v: CP + "\ndividedSuperquadricProps { alphaMin 0.1; %s }" % v
for txt, want in ((VF("periodicDirections (1 1 0);"), "OK"), (VF("periodicDirections (1 0 0);"), "FAIL"), (VF(""), "FAIL")):
    R = dict((i, s) for i, s, _ in C.audit(dem, G, U, txt, blockmesh_txt=BM))
    check(R.get("A9") == want, "A9 expected %s, got %s" % (want, R.get("A9")))
print("== test_case_checks: %s" % ("all checks passed" if not fails else "%d FAILED" % fails))
sys.exit(1 if fails else 0)
