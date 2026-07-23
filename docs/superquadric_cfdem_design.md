# Superquadric CFDEM Extension — Design & Verification Plan

Target stack: LIGGGHTS-PUBLIC 3.8.0 (superquadric DEM, already in use) + OpenFOAM 5.x + CFDEMcoupling-PUBLIC (fork).
Scope: unresolved (volume-averaged) CFD-DEM for prolate spheroids, orientation-dependent drag/lift/torque, satellite-point void fraction. GPL fork published on GitHub.

---

## 1. Where the code goes

OpenFOAM is untouched except for solver-level smoothing utilities. LIGGGHTS is nearly untouched
(one small patch, see 2.1). All new classes live in the CFDEMcoupling fork:

```
src/lagrangian/cfdemParticle/
├── cfdemCloud/            → cfdemCloudSQ            (new, subclass)
├── subModels/
│   ├── dataExchangeModel/ → twoWayMPI_SQ            (extend field push/pull)
│   ├── voidFractionModel/ → SQsatelliteVoidFraction  (new)
│   ├── forceModel/
│   │   ├── HolzerSommerfeldDrag/                    (new)
│   │   ├── spheroidLiftTorque/                      (new; Zastawny/Sanjeevi-type)
│   │   └── DiFeliceVoidageCorrection/               (reuse/wrap existing)
│   ├── averagingModel/    → reuse (dense/dilute)
│   ├── locateModel/       → reuse (engineSearch), extended for satellites
│   └── smoothingModel/    → reuse diffusion-based smoother (constDiffSmoothing)
```

### 1.1 `cfdemCloudSQ` — extra per-particle state

Additional arrays beyond the spherical cloud (N = local particle count):

| array          | size  | source                    |
|----------------|-------|---------------------------|
| `shape_`       | N × 3 | semi-axes (a, b, c)       |
| `blockiness_`  | N × 2 | (n1, n2) — = 2 for ellipsoids |
| `quat_`        | N × 4 | quaternion (w, x, y, z)   |
| `omega_`       | N × 3 | angular velocity (world frame) |
| `torqueExp_`   | N × 3 | hydrodynamic torque → pushed back to DEM |

Derived cached quantities (computed once for monodisperse shape, per particle otherwise):
volume-equivalent diameter d_eq, sphericity Φ, body-frame satellite point set (Sec. 3).

### 1.2 LIGGGHTS-side patch (small)

`fix couple/cfd` communicates per-atom vectors by name. Superquadric atom style already stores
`shape`, `blockiness`, `quaternion`, `omega`, `torque` as per-atom properties; the patch registers
these names in the push/pull lists of `fix_cfd_coupling_force` so `twoWayMPI` can transfer them.
Pull to CFD: shape, blockiness, quat, omega, v, x. Push to DEM: dragforce **and** `hdtorque`
(new per-atom vector added to the coupling fix, applied in the superquadric integrator alongside
contact torque).

---

## 2. Coupling loop (per CFD time step, coupling interval n_c DEM steps)

```
1. pull x, v, omega, quat, shape from LIGGGHTS
2. locate centroid cell for each particle (engineSearch, cached seed)
3. SQsatelliteVoidFraction:
     - map satellite points → cells, accumulate solid volume per cell
     - store per-particle {cellID, weight} lists for step 5
4. interpolate fluid velocity (and ∇p, ∇·τ if needed) to particle centroids
     - cellPoint interpolation; for particles longer than a cell,
       average u_f over the particle's own weighted cell list instead (flag)
5. force models:
     - HolzerSommerfeldDrag  → F_d per particle
     - spheroidLiftTorque    → F_L, T_pitch, T_spin per particle
     - buoyancy/∇p term      → F_∇p
     - distribute −(F_d+F_L) and torque reaction to cells using weights from 3
6. push F, T to LIGGGHTS; run n_c DEM steps
7. smooth voidage + momentum source fields (constDiffSmoothing, length ℓ_s)
8. solve fluid momentum (cfdemSolverPiso with Ksl implicit drag treatment)
```

Stability notes: keep coupling interval so that particle relaxation time τ_p = ρ_p d_eq²/(18 μ_f)
is resolved by ≥ 10 CFD steps at the lowest St studied; use implicit (Ksl) drag coupling to avoid
stiffness at low St.

---

## 3. Satellite-point void fraction — algorithm

**Precompute (once per shape, body frame):**

```
input: semi-axes (a,b,c), target count Ns (default 29–100, convergence-tested)
option A (deterministic, preferred for spheroids):
    concentric-shell lattice: scale a template sphere lattice by (a,b,c)
option B (general superquadric):
    rejection-sample uniform points in bounding box [-a,a]×[-b,b]×[-c,c]
    keep x if F(x) ≤ 0, where
    F(x) = ( (|x/a|^(2/n2) + |y/b|^(2/n2))^(n2/n1) + |z/c|^(2/n1) ) − 1
assign each satellite sub-volume Vp/Ns
```

**Per coupling step, per particle p:**

```
R = rotationMatrix(quat_p)                    # body → world
for s in satellites:
    xw = x_p + R · xs                          # world position
    cell = locate(xw, seed = centroidCell)     # walk-search from centroid
    if cell is local:
        alphaSolid[cell]   += Vp/Ns
        weight_p[cell]     += 1/Ns             # for force distribution
    else:
        stash (xw, Vp/Ns, particleID) → parallel transfer buffer
parallel: exchange stashed satellites to owning processors, accumulate there
voidfraction = 1 − alphaSolid/Vcell            # clip at alphaMin (e.g. 0.10)
```

**Parallel correctness (the real work):** a c/a = 5 particle near a processor boundary can have
satellites on 2–4 ranks. Two safe designs: (i) satellite scatter-gather as above (simplest,
message size ~ Ns × boundary particles); (ii) rely on LIGGGHTS ghost atoms + halo width ≥ c so each
rank computes its own cells from ghosts (no messages, larger halos). Implement (i) first; it is
independent of DEM halo settings.

**Smoothing:** at ν ≈ 0.46 with 2c comparable to Δx_CFD, raw voidage is noisy and locally < alphaMin.
Apply implicit diffusion smoothing to voidage and momentum source with smoothing length
ℓ_s ≈ (1–2)·2c, and report sensitivity to ℓ_s in the paper (this is a known modeling knob, not a bug).

---

## 4. Force and torque models

Definitions: slip u_r = u_f@p − v_p, Re = ρ_f |u_r| d_eq / μ_f, symmetry axis k = R·ẑ,
incidence angle α = angle(k, u_r) ∈ [0, π/2].

### 4.1 Projected & sphericity quantities (analytic for spheroids)

```
A_perp(α) = π a sqrt(a² cos²α + c² sin²α)        # projected area ⊥ flow
Φ   = sphericity = A_sphere(d_eq)/A_spheroid      # closed form for prolate
Φ⊥  = crosswise sphericity = (π d_eq²/4)/A_perp(α)
Φ∥  = lengthwise sphericity (Hölzer–Sommerfeld definition, closed form)
```

### 4.2 Drag — Hölzer & Sommerfeld (2008), orientation-aware default

```
C_D = 8/(Re √Φ∥) + 16/(Re √Φ) + 3/(√Re · Φ^{3/4})
      + 0.42 · 10^{0.4(−log Φ)^0.2} / Φ⊥
F_d = ½ ρ_f C_D A_perp(α) |u_r| u_r · f_ε
```

Voidage correction f_ε = ε^{−β}, Di Felice: β = 3.7 − 0.65 exp[−(1.5 − log₁₀Re)²/2].
Optional alternative model: sin²-interpolation form C_D(α) = C_D,0 + (C_D,90 − C_D,0) sin²α with
endpoint values from Sanjeevi & Padding DNS fits — run both, report sensitivity.

### 4.3 Lift and pitching torque (the term that moves your orientation physics)

```
C_L(α, Re): Zastawny/Sanjeevi form, ∝ sinα cosα at leading order
F_L = ½ ρ_f C_L A_ref |u_r|² · ê_L,  ê_L = (k×u_r)×u_r / |...|
C_T(α, Re): pitching-torque fit, ∝ sinα cosα
T_pitch = ½ ρ_f C_T A_ref (d_eq/2) |u_r|² · ê_T,  ê_T = k×u_r/|k×u_r|
T_spin  = rotational damping: C_R(Re_ω) form vs. relative spin (u_f vorticity/2 − ω)
```

Caveat to state honestly in both papers: published DNS coefficient fits exist only for specific
aspect ratios (e.g. 2.5, 4); interpolate in c/a and verify against ladder cases V2–V3. If the
science demands it, generate your own fits later with PR-DNS (MFiX IBM-superDEM route) — that is a
separate paper.

### 4.4 Fluid-side reaction

Distribute −F and the torque-equivalent momentum couple to cells using the satellite weight list
(consistent with where the solid actually is), not the centroid cell only. Torque reaction on the
fluid may be neglected at first (standard in unresolved CFD-DEM) — flag as `torqueReactionOnFluid off`
and document.

---

## 5. Verification ladder (each = one test case in the repo, CI-runnable)

| # | Case | Reference | Pass tolerance |
|---|------|-----------|----------------|
| V0 | Sphere regression: single settling sphere; Ergun packed bed | stock CFDEMcoupling results | ≤ 1–2 % (must not break spheres) |
| V1 | Stokes settling, fixed prolate spheroid, α = 0° and 90°, Re < 0.1 | Oberbeck/Perrin analytical friction factors | ≤ 3 % terminal velocity |
| V2 | Fixed spheroid drag & lift vs. α at Re = 1, 10, 100 | Hölzer–Sommerfeld / Sanjeevi–Padding DNS | ≤ 10 % C_D, C_L (within correlation scatter) |
| V3 | Freely rotating spheroid in simple shear, low Re | Jeffery orbit period T = (2π/γ̇)(λ + 1/λ) | ≤ 5 % on period; orbit shape qualitative |
| V4 | Packed bed of rods, pressure drop vs. superficial velocity | Nemec–Levec / non-spherical Ergun | ≤ 15 % |
| V5 | Fluidized bed of elongated particles | Vollmari et al. (superquadrics), Mahajan et al. (spherocylinders, pseudo-2D) | U_mf ≤ 10–15 %; bed expansion & orientation distribution qualitative–semiquantitative |
| V6 | **Dry-limit regression of your shear cell** (fluid forces off or St → ∞) | your Table I (φ, ψ, κ_w, i₁, T) | within published error bars |
| V7 | Satellite-count & smoothing-length convergence study | self | reported, monotone convergence; Ns and ℓ_s chosen on plateau |

V6 is the bridge between the software paper and the physics paper: it proves the coupled code
reduces exactly to your published dry results.

---

## 6. Phased plan

**Phase 1 (weeks):** cfdemCloudSQ + data exchange patch; centroid/equivalent-sphere voidage;
Hölzer–Sommerfeld drag + pitching torque. Run V0–V3. → enough to *start* shear-cell physics runs
at moderate-to-high St.

**Phase 2 (the months):** SQsatelliteVoidFraction with parallel scatter-gather; smoothing;
V4–V5–V7. → release version; software paper (SoftwareX / Comput. Phys. Commun. / Powder Technol.).

**Phase 3 (physics paper):** wall-driven shear cell, c/a = {2, 3, 5}, St sweep (≥ 4 values per
shape + dry limit), re-extract (φ, ψ, κ_w) with the existing least-squares pipeline. Headline
questions: does φ ≈ 1.2 relax toward Bretherton B(λ) as St ↓; ψ(St) vs. Folgar–Tucker C_I.

**Repo hygiene:** GPLv3 (inherited), pinned versions (OpenFOAM 5.x, LIGGGHTS 3.8.0 + patch as a
diff), Allrun scripts per verification case, CITATION.cff, Zenodo DOI on release.
