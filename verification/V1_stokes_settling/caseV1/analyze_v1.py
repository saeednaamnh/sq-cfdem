#!/usr/bin/env python3
"""
analyze_v1.py - V1 verification: Stokes drag on a fixed prolate spheroid.

Compares the measured hydrodynamic force (from DEM/post/forces_v1.txt)
against TWO references:

  1. Oberbeck/Perrin (exact Stokes solution, unbounded fluid)
       -> physics truth. Deviation includes wall blockage (~3-5% here,
          slip walls, a_eq/L ~ 0.025) and correlation error.
  2. Holzer-Sommerfeld correlation at the same Re and orientation
       -> what OUR force model is supposed to reproduce. Measured vs H-S
          isolates IMPLEMENTATION error (target <= 1-2%).

Stokes linearity: F(alpha) = F_par cos^2(a) + F_perp sin^2(a).

Usage:
  python3 analyze_v1.py --forces DEM/post/forces_v1.txt --alpha 90 \
      [--a 0.6785e-3 --c 1.357e-3 --U 0.05 --nu 1.5e-3 --rho 10]
  For the sphere control run: add --sphere (uses R = a).
"""
import argparse, math, sys

def perrin_factors(lam):
    """Perrin translation friction factors K_par, K_perp for a prolate
    spheroid (lam = c/a > 1), normalized by the equal-VOLUME sphere.
    Sphere limit lam->1 gives (1,1)."""
    if abs(lam - 1.0) < 1e-9:
        return 1.0, 1.0
    e = math.sqrt(1.0 - 1.0/lam**2)
    L = math.log((1.0 + e)/(1.0 - e))
    # Perrin (1934); normalized to sphere of equal volume: R_eq = a*lam^(1/3)
    # X_A, Y_A are normalized by the MAJOR semi-axis c (Kim & Karrila);
    # multiply by lam^(2/3) = c/R_eq to renormalize by equal-volume sphere.
    X_A = (8.0/3.0)*e**3 / (-2.0*e + (1.0 + e**2)*L)
    Y_A = (16.0/3.0)*e**3 / ( 2.0*e + (3.0*e**2 - 1.0)*L)
    s   = lam**(2.0/3.0)
    return s*X_A, s*Y_A

def holzer_sommerfeld_F(rho, nu, U, a, c, alpha_rad):
    """Force from the H-S correlation exactly as implemented in
    spheroidGeometry.H (mirrored here in python)."""
    dEq = 2.0*(a*a*c)**(1.0/3.0)
    Re  = dEq*U/nu
    cosA = abs(math.cos(alpha_rad))
    # surface area (prolate, exact)
    if c <= a*(1+1e-12):
        A = 4.0*math.pi*a*a
    else:
        e = math.sqrt(1.0 - (a*a)/(c*c))
        A = 2.0*math.pi*a*a*(1.0 + (c/(a*e))*math.asin(e))
    Phi   = math.pi*dEq*dEq/A
    Aperp = math.pi*a*math.sqrt(a*a*cosA*cosA + c*c*(1.0-cosA*cosA))
    PhiP  = 0.25*math.pi*dEq*dEq/Aperp
    Apar_denom = max(0.5*A - Aperp, 1e-30)
    PhiL  = 0.25*math.pi*dEq*dEq/Apar_denom
    Cd = ( 8.0/(Re*math.sqrt(PhiL)) + 16.0/(Re*math.sqrt(Phi))
         + 3.0/(math.sqrt(Re)*Phi**0.75)
         + 0.42*10.0**(0.4*max(0.0,-math.log10(Phi))**0.2)/PhiP )
    return 0.5*rho*Cd*Aperp*U*U, Re

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--forces", required=True)
    p.add_argument("--alpha", type=float, required=True,
                   help="angle between symmetry axis and flow, degrees")
    p.add_argument("--a", type=float, default=0.6785e-3)
    p.add_argument("--c", type=float, default=1.357e-3)
    p.add_argument("--U", type=float, default=0.05)
    p.add_argument("--nu", type=float, default=1.5e-3)
    p.add_argument("--rho", type=float, default=10.0)
    p.add_argument("--sphere", action="store_true")
    p.add_argument("--tail", type=int, default=5,
                   help="average the last N force samples")
    args = p.parse_args()

    a, c = (args.a, args.a) if args.sphere else (args.a, args.c)
    mu = args.rho*args.nu
    alpha = math.radians(args.alpha)

    # measured force: average of last N samples of Fz
    rows = []
    with open(args.forces) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) >= 4:
                rows.append([float(x) for x in parts[:4]])
    if not rows:
        sys.exit("no force samples found")
    tail = rows[-args.tail:]
    Fz = sum(r[3] for r in tail)/len(tail)
    Fx = sum(r[1] for r in tail)/len(tail)
    Fy = sum(r[2] for r in tail)/len(tail)
    Fmeas = abs(Fz)

    # Oberbeck/Perrin exact
    lam = c/a
    Req = a*lam**(1.0/3.0)
    F_sph = 6.0*math.pi*mu*Req*args.U
    Kpar, Kperp = perrin_factors(lam)
    K_alpha = Kpar*math.cos(alpha)**2 + Kperp*math.sin(alpha)**2
    F_exact = F_sph*K_alpha

    # H-S (our model's target)
    F_hs, Re = holzer_sommerfeld_F(args.rho, args.nu, args.U, a, c, alpha)

    print(f"shape: a={a*1e3:.4f} mm  c={c*1e3:.4f} mm  lam={lam:.3f}"
          f"  alpha={args.alpha:.1f} deg  Re={Re:.4f}")
    print(f"measured   F = {Fmeas:.6e} N   (Fx={Fx:.2e}, Fy={Fy:.2e})")
    print(f"H-S model  F = {F_hs:.6e} N   -> measured/H-S    = "
          f"{Fmeas/F_hs:.4f}  ({(Fmeas/F_hs-1)*100:+.2f}%)"
          f"   [implementation check, target |dev| <= 2%]")
    print(f"Oberbeck   F = {F_exact:.6e} N   -> measured/exact  = "
          f"{Fmeas/F_exact:.4f}  ({(Fmeas/F_exact-1)*100:+.2f}%)"
          f"   [physics check; includes ~3-5% blockage +"
          f" correlation error]")
    print(f"  (Perrin factors: K_par={Kpar:.4f}  K_perp={Kperp:.4f}"
          f"  K(alpha)={K_alpha:.4f};  sphere ref F_sph={F_sph:.4e} N)")

if __name__ == "__main__":
    main()
