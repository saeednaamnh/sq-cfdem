import os, json
from pathlib import Path

CONF = Path.home()/".sqfoam.conf"

# per-shape geometry + dry shear rate (from Table I)
SHAPES = {
    "sphere":    dict(a=0.3420,b=0.3420,c=0.3420,N=3724, gd=0.06984, ca=1.0),
    "prolate_2": dict(a=0.2714,b=0.2714,c=0.5429,N=3780, gd=0.06758, ca=2.0),
    "prolate_3": dict(a=0.2371,b=0.2371,c=0.7114,N=3804, gd=0.06971, ca=3.0),
    "prolate_5": dict(a=0.2000,b=0.2000,c=1.0000,N=3830, gd=0.02225, ca=5.0),
}
RHO_P, KN = 1.0, 2e5

DEFAULTS = {
    "CFDEM_LIGGGHTS_SRC_DIR": os.environ.get("CFDEM_LIGGGHTS_SRC_DIR",
        str(Path.home()/"LIGGGHTS/LIGGGHTS-PUBLIC/src")),
    "TUTORIAL_CFD": str(Path.home()/
        "CFDEM/CFDEMcoupling-PUBLIC/tutorials/cfdemSolverPiso/ErgunTestMPI/CFD"),
    "DRY_RESTARTS": str(Path.home()/"Desktop/aspect_ratio_sweep/cases"),
    "EXTRACT_PIPELINE": str(Path.home()/"Desktop/aspect_ratio_sweep"),
    "POSTPROCESS_PY": str(Path.home()/
        "Desktop/aspect_ratio_sweep/postprocess/postprocess_shear.py"),
    "NCORES": "4",
    "SOLVER": "cfdemSolverPisoSQ",
}

def load():
    if CONF.is_file():
        return {**DEFAULTS, **json.load(open(CONF))}
    return dict(DEFAULTS)

def save(cfg):
    json.dump(cfg, open(CONF,"w"), indent=2)

import math
def contact_dt(s):
    V=(4/3)*math.pi*s["a"]*s["b"]*s["c"]; mij=RHO_P*V/2
    e=0.98
    eta0=math.sqrt(KN/mij*(math.log(e))**2/(math.pi**2+(math.log(e))**2))
    om0=math.sqrt(KN/mij-eta0**2)
    return (math.pi/om0)/50.0

def dp_of(s):
    return 2.0*(s["a"]*s["b"]*s["c"])**(1.0/3.0)

def fluid_props(s, St, rho_f=1e-3):
    dp=dp_of(s); tau_p=St/s["gd"]
    mu_f=RHO_P*dp*dp/(18.0*tau_p)
    return rho_f, mu_f/rho_f, mu_f, tau_p
