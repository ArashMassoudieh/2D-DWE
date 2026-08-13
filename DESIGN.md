# 2D-DWE — design note

A two-dimensional overland-flow solver for the Coyote Canyon (NHA WS3 cluster C6) spill
problem, and for parameter sweeps that HEC-RAS makes impractical.

Status: **design only.** No code written yet. Decisions marked **[OPEN]** need settling
before implementation starts.

---

## 1. What this is for

The C6 study established, with a 1-D steady HEC-RAS model, that the 100-year water surface
in Figueredo Wash stands **3.3 ft above the controlling divide** separating the wash from the
community, over **786 ft of the 2,620 ft community frontage**, at 38 of 80 valid transects.

What that model cannot produce — because its cross sections do not span the flow and its
domain stops ~2,300 ft short of the nearest structure — is the volume crossing the divide,
the depth and velocity on the community bench, and how far the spill travels. Those are the
quantities a mitigation decision needs.

This solver exists to produce them, and to do it cheaply enough that the questions that
actually matter can be asked hundreds of times:

- what berm height on the north-west subdivision boundary is required, and where must it tie in;
- how sensitive is that to Manning's *n*, to the saddle elevation, and to the design discharge;
- which of the 84 structures are exposed, and at what depth.

**It is not intended to be the sealed deliverable.** HEC-RAS 7.0.1 is the model of record for
the report; a reviewer can audit it and cannot audit this. This is a second opinion and a
sweep engine, validated against RAS on the same domain.

---

## 2. Why diffusive wave is defensible here

Diffusive wave drops the inertial terms and asserts that the friction slope equals the
water-surface slope. That is invalid at a control section, where flow passes through
critical. The question was whether the divide acts as one.

Measured on the 1 m LiDAR along the flow direction at all 38 exceeded transects:

| quantity | value |
|---|---|
| crest width within 1 ft of the peak | median 27 ft |
| **crest width within 3 ft of the peak** | **median 304 ft** (p10 73, p90 614) |
| local slope through the crest | median 0.66 %, p90 2.85 % |
| Froude number, 3.3 ft depth, *n* = 0.05 | **0.52 — subcritical** |

Head-to-crest-length is 3.3 / 304 ≈ **0.011**. A broad-crested weir needs roughly h/L > 0.08
before it controls. This is a long, gentle rise carrying gradually varied subcritical flow.
There is no crest in the hydraulic sense, and no weir link is needed.

**The caveat, recorded so it is not forgotten:** at the steepest ~10 % of exceeded stations
the local slope reaches 2.85 %, giving Fr ≈ 1.08. Diffusive wave will bias the rating there.
Scheme A below retains local acceleration and handles those stations better, which is the
main reason it is preferred.

---

## 3. Governing equations

Conservation of mass on a cell, with `H = z + h` the water-surface elevation:

```
∂h/∂t + ∇·q = s
```

Two momentum closures are specified. **They share every other part of the code** — the face
logic, wetting and drying, mass accounting, boundaries and I/O are identical. Only the
`update_fluxes()` kernel differs, so both can be built and compared.

### Scheme A — local inertial (recommended)

Retains the local acceleration term (Bates, Horritt & Fewtrell, 2010, *J. Hydrol.* 387).
Per face, explicit:

```
q^{n+1} = ( q^n  −  g·h_f·Δt·∂H/∂x )
          ─────────────────────────────────
          ( 1 + g·Δt·n²·|q^n| / h_f^{7/3} )
```

then

```
h_i^{n+1} = h_i^n − (Δt/A)·Σ_f q_f·Δy  +  Δt·s_i
```

**Why this is preferred:**

- **Explicit.** No linear system, no preconditioner, no PETSc. Every cell updates
  independently — the natural GPU kernel.
- **Courant-limited timestep**, `Δt ≤ α·Δx/√(g·h_max)` with α ≈ 0.7. At Δx = 1 m and
  h = 1 m that is ≈ 0.22 s: workable.
- **More accurate than pure DW** near Fr ≈ 1, which is exactly the 10 % of stations flagged
  in §2.
- Well validated for floodplain inundation; it is the basis of LISFLOOD-FP.

### Scheme B — semi-implicit diffusive wave

Strict DW. Substituting Manning into the momentum closure gives a nonlinear parabolic
equation in `H`:

```
∂h/∂t = ∇·( K·∇H ) + s ,        K = h^{5/3} / ( n·|∇H|^{1/2} )
```

**Explicit time-stepping is not viable for this form.** Stability requires
`Δt ≤ Δx²/(4K)`, and `K → ∞` as `|∇H| → 0`. Even with the gradient floor of §4.2
(S_min = 1e-5), at h = 1 m and n = 0.05 the conductance reaches ≈ 6.3 × 10³ m²/s, giving
Δt ≈ 4 × 10⁻⁵ s at 1 m cells. Hopeless. Scheme B must therefore be semi-implicit:

```
(A/Δt)·h_i^{n+1} + Σ_f (K_f/Δx)·( h_i^{n+1} − h_j^{n+1} )
        = (A/Δt)·h_i^n − Σ_f (K_f/Δx)·( z_i − z_j ) + A·s_i
```

with `K_f` lagged at level *n*. The operator is a **symmetric positive definite 5-point
stencil**: diagonal `A/Δt + Σ_f K_f/Δx`, off-diagonals `−K_f/Δx`.

**[OPEN]** Build Scheme A first. Implement B only if the comparison against RAS suggests the
local-acceleration term is doing something unexpected.

---

## 4. Discretisation

Cell-centred finite volume on a uniform structured grid, Δx = Δy. Cell `(i,j)` stores `z`,
`h`, `n`; faces store `q` (Scheme A) or `K` (Scheme B).

### 4.1 The face-depth rule — the single most important line in the code

Never average depths across a face. Use the depth over the **higher** of the two bed
elevations:

```cpp
const double zf = std::max(z[i], z[j]);          // higher bed of the pair
const double hf = std::max(H[i], H[j]) - zf;     // depth over that bed
if (hf <= H_DRY) { q[f] = 0.0; continue; }       // face closed
```

This is what stops water climbing a rise. An averaged face depth leaks water onto dry high
ground, and the leak is invisible in the output — it looks like a plausible flood.

### 4.2 Gradient floor (Scheme B only)

`K ∝ |∇H|^{−1/2}` is singular on flat water:

```cpp
const double S = std::max(std::fabs(H[i]-H[j])/dx, S_FLOOR);   // S_FLOOR = 1e-5
K[f] = (dy/n_f) * std::pow(hf, 5.0/3.0) / std::sqrt(S);
```

### 4.3 Grid resolution

The C6 rectangle (`outputs_2D/boundary.shp`, 4,351 × 2,050 m, 8.92 km²):

| Δx | cells | state @ ~10 doubles | crest under-representation (median / p90 / worst) |
|---|---|---|---|
| 5 m | 357 k | 29 MB | 0.66 / 1.88 / 3.02 ft |
| 3 m | 991 k | 79 MB | 0.32 / 0.87 / 1.72 ft |
| 2 m | 2.23 M | 178 MB | 0.26 / 0.74 / 1.18 ft |
| **1 m** | **8.92 M** | **714 MB** | **none — it is the DEM** |

The margin under test is 3.3 ft. Storing one elevation per cell means a face that crosses
the rise sees it low by the amounts above — a one-directional bias, always toward more spill.

**Develop at 3 m; produce results at 1 m.** At native resolution the discretisation
compromise disappears entirely, and with it the need for breaklines, sub-grid face profiles,
or any argument about smoothing the DEM.

**Do not pre-smooth the terrain.** A moving average is a low-pass filter: it lowers ridges
and raises channel inverts. Measured on this crest — 5×5 m lowers it 0.26 ft median /
1.00 ft worst; 15×15 m lowers it 0.95 / 4.61 ft. There is also nothing to fix: no cell in
this DEM departs from its own 5×5 median by more than 1.8 ft.

---

## 5. Wetting and drying

Where these codes fail. Five rules, in order of importance.

1. **The face-depth rule of §4.1.** Faces, not cells, carry the wet/dry decision.

2. **Positivity limiter — non-negotiable.** After computing fluxes, before applying them:

   ```cpp
   double out = 0.0;
   for (f : faces(i)) if (q[f] > 0) out += q[f]*dy;     // leaving cell i
   const double avail = h[i]*A;
   if (out*dt > avail) {
       const double a = avail/(out*dt);
       for (f : faces(i)) if (q[f] > 0) q[f] *= a;      // scale BOTH sides of the face
   }
   ```

   Guarantees `h ≥ 0` and conserves mass exactly. Without it, negative depths appear at the
   wetting front and the run either dies or — worse — does not.

3. **Hysteresis.** Wet at 5 mm, dry at 2 mm. Equal thresholds make cells flicker wet/dry
   every step and the solution never settles.

4. **Account for the discarded film.** When a cell drops below `H_DRY` its residual depth is
   zeroed. That mass goes into a `discarded` accumulator, never into thin air.

5. **Scheme B only — no active set needed.** Keeping the `A/Δt` diagonal means dry cells
   decouple and solve trivially as `h^{n+1} = h^n = 0`. Mask them from the Krylov solve as an
   optimisation, not for correctness.

---

## 6. Boundary conditions

| boundary | treatment | C6 value |
|---|---|---|
| inflow | discharge distributed over the face, hydrograph in time | west edge, on the wash: **32,628 cfs** (6-hr) / **29,884 cfs** (24-hr) |
| outflow | normal depth: `q = (1/n)·h^{5/3}·√S_bed` | full east edge — the low side of the domain, 6,130–6,138 ft |
| walls | closed faces, `q = 0` | remaining edges |

The reach has **one** inflow and no lateral inflows, so the discharge is the same wherever
the inflow face is placed. Hydrographs come from the converged HEC-HMS run
(`Run_6hr_100yr.dss`, `Run_24hr_100yr.dss`), 3-minute interval.

**Run the 6-hour storm as primary.** It governs at every one of the 121 RAS cross sections,
and its control window is 9 h against 36 h for the 24-hour event.

**Watch for ponding against the east boundary.** If water banks up there, the domain is too
short and the rectangle must be extended east.

---

## 7. Mass balance

Not optional, and not a debugging aid — it is the acceptance test.

```
inflow − outflow − Δstorage − discarded  <  0.1 % of cumulative inflow
```

Track all four cumulatively; report every N steps and at the end of the run. **If this does
not close, no result derived from the run means anything.** A solver that produces a
plausible-looking depth map with a 5 % mass error is worse than one that crashes.

---

## 8. Linear solver (Scheme B only)

Scheme A needs none. For B:

- The operator is **SPD, 5-point, structured** → conjugate gradient.
- **Matrix-free.** `K_f` changes every step because it depends on depth; assembling a sparse
  matrix each step — and on GPU, transferring it — would dominate the runtime. Store four
  face conductances per cell and apply the operator as a stencil loop.
- **Jacobi preconditioning** is sufficient given the strong `A/Δt` diagonal. Seed each solve
  with the previous step's `h`; expect 10–30 iterations.
- **Measure the iteration count before choosing hardware.** If it stays in the tens,
  matrix-free Jacobi-CG is right and PETSc buys little. If it reaches the hundreds, the
  system needs algebraic multigrid — which requires an assembled matrix, and the whole
  calculus changes.

On libraries: Armadillo is excellent for dense work but its sparse path is `spsolve` →
SuperLU, a *direct* factorisation whose fill-in will exhaust memory at 10⁶ unknowns. Use
`arma::vec` for the vector algebra (dots, axpys, norms) and write the operator apply
yourself. Eigen's `ConjugateGradient` + `IncompleteCholesky` is the alternative if a
ready-made iterative solver is wanted.

---

## 9. Architecture

The one rule that matters: **the CPU and GPU paths share a single operator definition**, so
porting is not a rewrite.

```
2D-DWE/
  DESIGN.md
  src/
    grid.hpp          index maps, geometry, neighbour access
    fields.hpp        z, h, n, q — SoA layout, no AoS
    faces.hpp         face-depth rule, conductance, closure  ← shared, hot
    scheme_inertial.* Scheme A flux update
    scheme_dw.*       Scheme B assembly + matrix-free CG
    wetdry.hpp        thresholds, hysteresis, positivity limiter
    massbalance.hpp   accumulators and the closure report
    bc.*              inflow hydrograph, normal depth, walls
    io_gdal.*         GeoTIFF in/out, CSV hydrograph
    main.cpp
  tests/
  cases/c6/           config for the Coyote Canyon domain
```

Guidance:

- **Structure of arrays**, flat `std::vector<double>`, index `k = j*nx + i`. No cell objects,
  no virtual dispatch in the inner loops.
- Keep the face kernel free of branches where possible; the wet/dry test is one unavoidable
  branch and it is well predicted.
- CPU parallelism: `#pragma omp parallel for` over rows. GPU later: the same loop body
  becomes the kernel body.
- **Never touch GDAL inside the timestep.** Read once, write on an output cadence.

**[OPEN]** Build system — CMake assumed. GDAL is the only hard dependency for I/O.

---

## 10. Verification and validation

In order. Do not skip to the last one.

1. **Mass conservation** on every run (§7).
2. **Steady uniform flow on a tilted plane.** Constant inflow on a constant slope must
   converge to the Manning normal depth analytically. Tests the friction term and the face
   conductance.
3. **Rainfall on an impermeable plane.** Tests the source term and the rising limb against
   the kinematic analytical solution.
4. **Wetting front on a dry adverse slope.** No analytical answer, but the front must stop
   where the head runs out, and mass must close. Tests the wet/dry machinery — the part most
   likely to be wrong.
5. **A closed depression.** Water must fill it to the spill elevation and no further. Directly
   analogous to the divide problem and it will expose a leaky face rule immediately.
6. **HEC-RAS on the C6 domain**, matched resolution, same terrain, same *n*, same hydrograph.
   Compare: peak depth field, spill volume across the divide, arrival time at the community,
   and depth at each of the 84 structures. Differences are expected — RAS 2-D uses sub-grid
   face profiles — and should shrink as Δx → 1 m.

---

## 11. Expected performance

Scheme A, explicit, at Δt ≈ 0.2 s (Courant, 1 m cells) over the 9-hour 6-hour-storm window
is ≈ 1.6 × 10⁵ steps. At 8.92 M cells that is ≈ 1.5 × 10¹² cell-updates.

- CPU, OpenMP, ~10⁸ updates/s/core over 8 cores → several hours per run.
- GPU, ~10¹⁰ updates/s → minutes per run.

At 3 m (991 k cells, Δt ≈ 0.4 s) the CPU run is minutes, which is why development happens
there.

**The GPU case rests on resolution, not impatience.** 1 m removes the discretisation
compromise that dominates §4.3; the speed is what makes 1 m affordable for a sweep.

---

## 12. Open decisions

| | |
|---|---|
| **[OPEN]** | Scheme A only, or both A and B? Recommendation: A first. |
| **[OPEN]** | Build system and dependency policy (CMake + GDAL assumed). |
| **[OPEN]** | Manning's *n* by NLCD class — screening values exist in `outputs_2D/manning_n_2d.csv`; confirm against NHA Manual Tables 3.1–3.4 before results are relied on. |
| **[OPEN]** | Output cadence and format — GeoTIFF depth/velocity snapshots plus a per-structure CSV time series is the assumption. |
| **[OPEN]** | Whether the domain needs extending east (§6). |
| **[OPEN]** | Single or double precision. Start double; single is worth testing on GPU only after mass balance is proven. |

---

## 13. Inputs already prepared

Generated by `NHA/PythonScripts/ras2d/C6-2d-inputs.py` into
`WS3_GIS/NM15-097/outputs_2D/`:

| file | contents |
|---|---|
| `boundary.shp` | domain rectangle, 4,351 × 2,050 m |
| `terrain_2d.tif` | 1 m LiDAR, 4,751 × 2,450, EPSG:26912 |
| `landcover_2d.tif` + `manning_n_2d.csv` | roughness by NLCD class |
| `for_ras/bc_lines.shp` | inflow (west) and outflow (east) faces |
| `for_ras/breaklines.shp` | divide crest, channel, two berm alignments |

Breaklines are a RAS concept and are **not** needed by this solver at 1 m — they are the
mechanism by which RAS compensates for cells coarser than the terrain. They remain useful
here as the berm alignments for the mitigation sweep.

---

## 14. References

- Bates, P.D., Horritt, M.S. & Fewtrell, T.J. (2010). A simple inertial formulation of the
  shallow water equations for efficient two-dimensional flood inundation modelling.
  *Journal of Hydrology*, 387, 33–45.
- de Almeida, G.A.M. & Bates, P. (2013). Applicability of the local inertial approximation of
  the shallow water equations to flood modeling. *Water Resources Research*, 49, 4833–4844.
- USACE HEC. *HEC-RAS 2D Modeling User's Manual*.
- NHA. *Drainage Design Manual* (final).
- C6 study: `WS3_GIS/NM15-097/Report_for_C6/main.tex` and its `HANDOFF.md`.

---

## Appendix A — implementation log

**2026-08-12 — Scheme A core implemented and verified.**

Built with `qmake6` (plain `qmake` on this machine is a broken Qt5 symlink).
Solver core has no third-party dependency; GDAL and VTK are opt-in via
`CONFIG+=with_gdal` / `CONFIG+=with_vtk`.

### Verification results

| test | criterion | result |
|---|---|---|
| tilted plane | Manning normal depth, analytical | **0.2276 m expected, 0.2276 m obtained (+0.00 %)** |
| tilted plane | mass closure | 7.8e-15 |
| closed depression | no water over the 10.0 m rim | **0.00000 m — the face rule holds** |
| closed depression | drains toward the 8.0 m notch | 8.43 → 8.14 → 8.04 m at t = 2/10/40 ks |
| closed depression | mass closure | 5.0e-15 |

Mass closes to machine precision in both, which is what the Kahan summation in
`Fields::volume()` is for — a naive sum loses enough at 10⁷ cells to muddy a
0.1 % test.

### Two things learned from the first runs

**The pond test must be pre-filled, not fed.** The first version injected water
at the west domain edge, which sits *outside* the bowl on the rim plateau. It
reported a rim leak that was really water correctly sitting where it had been
put. Pre-filling the bowl and letting it drain removes the boundary forcing and
isolates the face rule, which is what the test is for.

**The pond drains asymptotically and never reaches the spill level in finite
time** — notch discharge tends to zero with the head. The acceptance criterion is
therefore that the level approaches the notch *from above* and never settles
below it; settling below would mean the notch was leaking. A fixed tolerance at
a fixed time is the wrong test and was replaced.

**A consequence worth carrying forward:** in the first (wrong) pond run, 32 % of
the inflow was lost to sub-`h_dry` film discarding as a thin sheet spread over a
flat plateau. That is visible only because the discarded volume is accounted for
rather than dropped. At C6 depths it will be negligible, but on any thin-sheet
problem `h_dry` is a real mass sink and the accumulator must be watched.

### Still to write

`io_gdal.*` (terrain and roughness rasters, the C6 case), `viewer_vtk.*`
(interactive rendering; results already write as `.vti` and open in ParaView),
and the remaining rungs of the validation ladder — rainfall on a plane, the
wetting front on an adverse slope, and the HEC-RAS comparison.
