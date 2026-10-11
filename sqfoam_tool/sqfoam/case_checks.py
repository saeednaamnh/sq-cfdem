"""
case_checks.py -- read a LIGGGHTS shear template or a finished CFD-DEM case and
report configuration faults. Pure functions, no side effects; used by
`sqfoam new` (to build cases correctly or refuse) and `sqfoam check` (to audit).

Every check corresponds to a fault found in earlier cases:
  A1 coupling fixes present, on the integrator's own group  (silent decoupling)
  A2 fluid gravity equal to the DEM gravity                  (one-sided buoyancy)
  A3 fluid wall velocities equal to the DEM wall velocities  (shear mismatch)
  A4 no active 'run N' in the DEM script                     (blocking run -> no coupling)
  A5 drag voidage correction stated explicitly               (silent default)
  A6 hydrodynamic torque model present                       (no fluid torque)
  A7 lift model only with explicit 'experimental true'       (uncalibrated forces)
  A8 template particle semi-axes equal to the requested shape (wrong template)
"""
import re

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


class CaseError(Exception):
    """A template or case that cannot be built or run correctly."""


# --------------------------------------------------------------------- DEM side
def integrator(txt):
    """(line, group) of the 'fix <id> <group> nve/superquadric ...' line, or None."""
    m = re.search(r"(?m)^(fix\s+\S+\s+(\S+)\s+nve/superquadric\b.*)$", txt)
    return (m.group(1), m.group(2)) if m else None


def coupling_group(txt):
    """Group of the 'couple/cfd/force/implicit' fix, or None if absent."""
    m = re.search(r"(?m)^fix\s+\S+\s+(\S+)\s+couple/cfd/force/implicit\b", txt)
    return m.group(1) if m else None


def wall_velocity(txt):
    """Value of 'variable wall_vel equal <number>', or None."""
    m = re.search(r"(?m)^variable\s+wall_vel\s+equal\s+(" + NUM + r")\s*$", txt)
    return float(m.group(1)) if m else None


def wall_motion_signs(txt):
    """(upper, lower) sign of the x-velocity in 'fix ... upperGroup/lowerGroup move linear'."""
    def sign(group):
        m = re.search(r"(?m)^fix\s+\S+\s+" + group + r"\s+move\s+linear\s+(-?)\$\{wall_vel\}\s+0\s+0", txt)
        return None if not m else (-1 if m.group(1) else 1)
    return sign("upperGroup"), sign("lowerGroup")


def dem_gravity(txt):
    """Gravity vector from 'fix ... gravity <mag> vector x y z' (None if no gravity fix).
    Raises CaseError if the gravity fix uses variables that cannot be evaluated here."""
    m = re.search(r"(?m)^fix\s+\S+\s+\S+\s+gravity\s+(\S+)\s+vector\s+(\S+)\s+(\S+)\s+(\S+)", txt)
    if not m:
        return None
    try:
        g, x, y, z = (float(v) for v in m.groups())
    except ValueError:
        raise CaseError("the DEM gravity fix uses variables (%s); write numbers so the"
                        " fluid gravity can be matched" % m.group(0).strip())
    n = (x * x + y * y + z * z) ** 0.5
    return (g * x / n, g * y / n, g * z / n)


def active_runs(txt):
    """Lines 'run N' with N > 0 that are not commented out."""
    out = []
    for m in re.finditer(r"(?m)^run\s+(\S+)", txt):
        if m.group(1) not in ("0",):
            out.append(m.group(0).strip())
    return out


def semi_axes(txt):
    """(a, b, c) from 'variable a/b/c equal <number>' in the template, or None."""
    v = {}
    for k in "abc":
        m = re.search(r"(?m)^variable\s+" + k + r"\s+equal\s+(" + NUM + r")\s*$", txt)
        if not m:
            return None
        v[k] = float(m.group(1))
    return v["a"], v["b"], v["c"]


def shape_of(txt, shapes, tol=1e-6):
    """Name of the entry of `shapes` whose semi-axes match the DEM script's a, b, c (or None)."""
    ax = semi_axes(txt)
    if ax is None:
        return None
    for name, s in shapes.items():
        if max(abs(x - y) for x, y in zip(ax, (s["a"], s["b"], s["c"]))) <= tol:
            return name
    return None


def wall_groups(txt):
    """Groups moved as walls ('fix <id> <group> move linear ...')."""
    return re.findall(r"(?m)^fix\s+\S+\s+(\S+)\s+move\s+linear\b", txt)


FLAG = "cfdCoupled"


def add_coupling_flag(txt, flag=FLAG):
    """Give every atom a per-atom property 'flag' (1 = coupled to the fluid) and set it to 0
    for the wall groups, so the fluid side can leave the walls out (couplingFlag). Inserted
    before the coupling fixes, after the contact model, so its 'run 0' (required by LIGGGHTS
    between 'fix property/atom' and 'set') is a plain DEM initialisation.
    Returns (new_txt, flag or None if there are no wall groups)."""
    walls = wall_groups(txt)
    if not walls:
        return txt, None
    m = re.search(r"(?m)^fix\s+\S+\s+\S+\s+couple/cfd\b.*$", txt) or \
        re.search(r"(?m)^fix\s+\S+\s+\S+\s+nve/superquadric\b.*$", txt)
    if not m:
        raise CaseError("no coupling fix or nve/superquadric integrator to insert the coupling flag before")
    block = ("# coupling flag: 1 = coupled to the fluid, 0 = not (the wall spheres);\n"
             "# read by the fluid side through couplingFlag in couplingProperties\n"
             "fix             %s all property/atom %s scalar no no no 1.0\n"
             "run             0\n" % (flag, flag) +
             "".join("set             group %s property/atom %s 0.0\n" % (g, flag) for g in walls) + "\n")
    return txt[:m.start()] + block + txt[m.start():], flag


def coupling_flag(cp_txt):
    """Value of the top-level 'couplingFlag <word>;' entry, or None if absent."""
    m = re.search(r"(?m)^\s*couplingFlag\s+(\S+?)\s*;", cp_txt)
    return m.group(1) if m else None


def inject_coupling(txt, couple_every=100):
    """Insert the CFDEM coupling fixes above the nve/superquadric integrator, on the
    integrator's own group. Raises CaseError instead of returning an uncoupled script."""
    found = integrator(txt)
    if not found:
        raise CaseError("no 'fix <id> <group> nve/superquadric' line in the DEM template:"
                        " refusing to write a case that would run without coupling")
    line, group = found
    if "couple/cfd" in txt:
        g = coupling_group(txt)
        if g != group:
            raise CaseError("template already has coupling fixes, on group '%s' but the"
                            " integrator acts on '%s'" % (g, group))
        return txt, group
    block = ("fix  cfd  all       couple/cfd couple_every %d mpi\n"
             "fix  cfd2 %s couple/cfd/force/implicit transfer_torque yes\n\n%s"
             % (couple_every, group, line))
    new = txt.replace(line, block, 1)
    if coupling_group(new) != group:
        raise CaseError("coupling insertion failed (internal error)")
    return new, group


# --------------------------------------------------------------------- CFD side
def foam_vector(txt, entry):
    """First 'entry ... ( x y z )' vector after a keyword, e.g. the value in constant/g."""
    m = re.search(entry + r"[^(;]*\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)", txt)
    return tuple(float(v) for v in m.groups()) if m else None


def patch_value(u_txt, patch):
    """Fixed velocity vector of a patch in 0/U, or None."""
    m = re.search(patch + r"\s*\{[^}]*value\s+uniform\s*\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)", u_txt)
    return tuple(float(v) for v in m.groups()) if m else None


def force_models(cp_txt):
    m = re.search(r"forceModels\s*\(([^)]*)\)", cp_txt)
    return m.group(1).split() if m else []


def props_switch(cp_txt, model, key):
    """Value of 'key' inside '<model>Props { ... }', or None if not stated."""
    m = re.search(model + r"Props\s*\{([^}]*)\}", cp_txt)
    if not m:
        return None
    k = re.search(r"\b" + key + r"\s+(\S+?)\s*;", m.group(1))
    return k.group(1) if k else None


def mesh_periodic(bm_txt):
    """(px, py, pz): 1 for each direction normal to a cyclic patch of a blockMeshDict,
    found from the patch face geometry (not from patch names)."""
    vs = bm_txt.split("vertices")[1].split(");")[0]
    V = [tuple(float(v) for v in m.groups()) for m in
         re.finditer(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)", vs)]
    per = [0, 0, 0]
    for m in re.finditer(r"\w+\s*\{([^}]*)\}", bm_txt.split("boundary")[-1] if "boundary" in bm_txt else bm_txt):
        body = m.group(1)
        if not re.search(r"type\s+cyclic\s*;", body):
            continue
        for f in re.finditer(r"\(\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*\)", body):
            a, b, c = (V[int(i)] for i in f.groups()[:3])
            u = [b[i] - a[i] for i in range(3)]; w = [c[i] - a[i] for i in range(3)]
            n = [u[1]*w[2] - u[2]*w[1], u[2]*w[0] - u[0]*w[2], u[0]*w[1] - u[1]*w[0]]
            per[max(range(3), key=lambda i: abs(n[i]))] = 1
    return tuple(per)


def props_vector(cp_txt, model, key):
    """Vector value '(x y z)' of 'key' inside '<model>Props { ... }', or None."""
    m = re.search(model + r"Props\s*\{([^}]*)\}", cp_txt)
    if not m:
        return None
    k = re.search(r"\b" + key + r"\s*\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)\s*;", m.group(1))
    return tuple(int(float(v)) for v in k.groups()) if k else None


# --------------------------------------------------------------------- audit
def audit(dem_txt, g_txt, u_txt, cp_txt, expected_axes=None, tol=1e-9, blockmesh_txt=None):
    """Static audit of a case. Returns a list of (id, 'OK'|'FAIL'|'WARN', message)."""
    R = []
    integ = integrator(dem_txt)
    cg = coupling_group(dem_txt)
    if not integ:
        R.append(("A1", "FAIL", "no nve/superquadric integrator in the DEM script"))
    elif cg is None:
        R.append(("A1", "FAIL", "no couple/cfd/force/implicit fix: the case would run UNCOUPLED"))
    elif cg != integ[1]:
        R.append(("A1", "FAIL", "coupling fix on group '%s' but integrator on '%s'" % (cg, integ[1])))
    else:
        R.append(("A1", "OK", "coupling on the integrator's group '%s'" % cg))

    try:
        gd = dem_gravity(dem_txt) or (0.0, 0.0, 0.0)
        gf = foam_vector(g_txt, "value")
        if gf is None:
            R.append(("A2", "FAIL", "cannot read the fluid gravity in constant/g"))
        elif max(abs(x - y) for x, y in zip(gd, gf)) > 1e-12:
            R.append(("A2", "FAIL", "fluid gravity %s but DEM gravity %s: one-sided buoyancy" % (gf, gd)))
        else:
            R.append(("A2", "OK", "fluid and DEM gravity both %s" % (gf,)))
    except CaseError as e:
        R.append(("A2", "FAIL", str(e)))

    uw = wall_velocity(dem_txt)
    up, lo = wall_motion_signs(dem_txt)
    top, bot = patch_value(u_txt, "wallTop"), patch_value(u_txt, "wallBottom")
    if uw is None or up is None or lo is None:
        R.append(("A3", "FAIL", "cannot read wall_vel and the wall 'move linear' fixes in the DEM script"))
    elif top is None or bot is None:
        R.append(("A3", "FAIL", "cannot read the wallTop/wallBottom velocities in 0/U"))
    elif abs(top[0] - up * uw) > tol or abs(bot[0] - lo * uw) > tol or any(abs(v) > tol for v in top[1:] + bot[1:]):
        R.append(("A3", "FAIL", "fluid walls move at %s (top) and %s (bottom), DEM walls at %+g and %+g in x"
                  % (top, bot, up * uw, lo * uw)))
    else:
        R.append(("A3", "OK", "fluid and DEM walls both at %+g (top) and %+g (bottom)" % (up * uw, lo * uw)))

    runs = active_runs(dem_txt)
    R.append(("A4", "FAIL", "active %s in the DEM script blocks the coupled time loop" % runs) if runs
             else ("A4", "OK", "no active 'run N' in the DEM script"))

    fm = force_models(cp_txt)
    if "HolzerSommerfeldDrag" in fm:
        v = props_switch(cp_txt, "HolzerSommerfeldDrag", "voidageCorrection")
        R.append(("A5", "FAIL", "voidageCorrection not stated in HolzerSommerfeldDragProps") if v is None
                 else ("A5", "OK", "voidageCorrection %s" % v))
    else:
        R.append(("A5", "WARN", "HolzerSommerfeldDrag not among the force models"))
    R.append(("A6", "OK", "spheroidRotationTorque present") if "spheroidRotationTorque" in fm
             else ("A6", "WARN", "no spheroidRotationTorque: the fluid exerts no torque on the particles"))
    if "spheroidLiftTorque" in fm:
        e = props_switch(cp_txt, "spheroidLiftTorque", "experimental")
        R.append(("A7", "OK" if e in ("true", "yes", "on") else "FAIL",
                  "spheroidLiftTorque enabled (experimental %s)" % e))
    else:
        R.append(("A7", "OK", "no lift model"))

    if expected_axes:
        ax = semi_axes(dem_txt)
        if ax is None:
            R.append(("A8", "WARN", "cannot read the semi-axes a, b, c in the DEM script"))
        elif max(abs(x - y) for x, y in zip(ax, expected_axes)) > 1e-6:
            R.append(("A8", "FAIL", "DEM script semi-axes %s differ from the requested shape %s" % (ax, expected_axes)))
        else:
            R.append(("A8", "OK", "semi-axes %s" % (ax,)))
    cf = coupling_flag(cp_txt)
    walls = wall_groups(dem_txt)
    if cf is None:
        R.append(("A10", "FAIL", "couplingFlag not stated in couplingProperties"))
    elif walls:
        problems = []
        fx = re.search(r"(?m)^fix\s+(\S+)\s+all\s+property/atom\s+" + re.escape(cf) + r"\s+scalar\s+\S+\s+\S+\s+\S+\s+(" + NUM + r")", dem_txt)
        if not fx or float(fx.group(2)) != 1.0:
            problems.append("no 'fix ... all property/atom %s scalar ... 1.0' in the DEM script" % cf)
        else:
            sets = [re.search(r"(?m)^set\s+group\s+" + re.escape(g) + r"\s+property/atom\s+" + re.escape(cf) + r"\s+(" + NUM + r")", dem_txt) for g in walls]
            for g, m in zip(walls, sets):
                if not m or float(m.group(1)) != 0.0:
                    problems.append("wall group %s not set to %s 0" % (g, cf))
            first_set = min((m.start() for m in sets if m), default=None)
            if first_set is not None and not re.search(r"(?m)^run\s+0\b", dem_txt[fx.end():first_set]):
                problems.append("no 'run 0' between the property/atom fix and the set commands (LIGGGHTS requires it)")
        R.append(("A10", "FAIL", "; ".join(problems)) if problems
                 else ("A10", "OK", "wall groups %s not coupled to the fluid (couplingFlag %s)" % (walls, cf)))
    else:
        R.append(("A10", "OK", "no wall groups; couplingFlag %s" % cf))

    if blockmesh_txt is not None and "dividedSuperquadric" in cp_txt:
        want = mesh_periodic(blockmesh_txt)
        have = props_vector(cp_txt, "dividedSuperquadric", "periodicDirections")
        if have is None:
            R.append(("A9", "FAIL", "periodicDirections not stated in dividedSuperquadricProps"
                      " (the mesh is periodic in %s)" % (want,)))
        elif have != want:
            R.append(("A9", "FAIL", "periodicDirections %s but the mesh's cyclic patches are %s" % (have, want)))
        else:
            R.append(("A9", "OK", "periodicDirections %s, as the mesh's cyclic patches" % (have,)))
    return R
