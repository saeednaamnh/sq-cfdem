#!/usr/bin/env python3
"""
analyze_v3.py - V3: Jeffery orbit period of a free prolate spheroid in shear.

Reads the multi-frame LIGGGHTS dump DEM/post/quat_v3.txt (id quat1..4),
reconstructs the symmetry axis k(t), tracks the in-plane angle
theta(t) = atan2(kx, kz) in the shear (x-z) plane, unwraps it, and
measures the tumbling period from the total swept angle (robust for
Jeffery's strongly non-uniform rotation rate).

Pass gate: |T_measured - T_Jeffery| / T_Jeffery <= 5%,
T_Jeffery = (2 pi / gdot) (lam + 1/lam).

Usage: python3 analyze_v3.py [--dump DEM/post/quat_v3.txt]
       [--gdot 2.0] [--lam 2.0] [--dt-per-step 1e-5] [--discard 1.0]
"""
import argparse, math, sys

def frames(path):
    ts, quat = None, None
    out = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            ts = int(lines[i+1]); i += 2
        elif lines[i].startswith("ITEM: ATOMS"):
            j = i + 1
            while j < len(lines) and not lines[j].startswith("ITEM"):
                p = lines[j].split()
                if len(p) == 5:
                    out.append((ts, [float(v) for v in p[1:]]))
                j += 1
            i = j
        else:
            i += 1
    return out

def axis(q):
    w, x, y, z = q
    return (2*(x*z + w*y), 2*(y*z - w*x), 1 - 2*(x*x + y*y))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default="DEM/post/quat_v3.txt")
    ap.add_argument("--gdot", type=float, default=2.0)
    ap.add_argument("--lam", type=float, default=2.0)
    ap.add_argument("--dt-per-step", type=float, default=1e-5)
    ap.add_argument("--discard", type=float, default=1.0,
                    help="seconds to discard while Couette flow develops")
    a = ap.parse_args()

    fr = frames(a.dump)
    if len(fr) < 10:
        sys.exit(f"only {len(fr)} frames found - run too short?")

    t0 = None
    ts_list, th_list = [], []
    prev, offset = None, 0.0
    for step, q in fr:
        t = step*a.dt_per_step
        kx, ky, kz = axis(q)
        th = math.atan2(kx, kz)
        if prev is not None:
            d = th - prev
            if d >  math.pi: offset -= 2*math.pi
            if d < -math.pi: offset += 2*math.pi
        prev = th
        ts_list.append(t)
        th_list.append(th + offset)

    # discard start-up, fit total swept angle
    pairs = [(t, th) for t, th in zip(ts_list, th_list) if t >= a.discard]
    if len(pairs) < 10:
        sys.exit("not enough post-discard samples")
    t1, th1 = pairs[0]
    t2, th2 = pairs[-1]
    swept = abs(th2 - th1)                      # radians (pi per half-tumble)
    n_half = swept/math.pi
    if n_half < 1.0:
        sys.exit(f"particle swept only {math.degrees(swept):.1f} deg after "
                 f"t={a.discard}s - not a full half-tumble; check torque path "
                 f"(is the axis drifting to the vorticity axis? ky trend) or "
                 f"extend endTime")
    T_meas = (t2 - t1)/(n_half/2.0)             # 2 half-tumbles per period
    T_j = (2*math.pi/a.gdot)*(a.lam + 1.0/a.lam)
    dev = (T_meas/T_j - 1)*100

    # log-spin check: axis should stay in shear plane (ky ~ const)
    ky0 = axis(fr[0][1])[1]
    kyE = axis(fr[-1][1])[1]

    print(f"frames: {len(fr)}   window: {t1:.2f}-{t2:.2f} s"
          f"   swept: {math.degrees(swept):.1f} deg ({n_half:.2f} half-tumbles)")
    print(f"measured period  T = {T_meas:.4f} s")
    print(f"Jeffery period   T = {T_j:.4f} s   (gdot={a.gdot}, lam={a.lam})")
    print(f"deviation = {dev:+.2f}%   [pass gate: |dev| <= 5%]")
    print(f"axis ky: start {ky0:+.3f} -> end {kyE:+.3f} "
          f"(should stay ~0: in-plane tumbling)")
    print("VERDICT:", "PASS" if abs(dev) <= 5.0 else "FAIL")

if __name__ == "__main__":
    main()
