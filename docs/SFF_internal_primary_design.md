# Design: internally-computed 1D primary field for the Secondary Field Formulation (SFF)

**Status:** draft design (no code written yet)
**Author:** A. Kelbert, with Claude (Opus 4.8)
**Scope:** spherical-coordinate 3D MT ModEM; integrate the ModEM-global1d analytic
1D spherical solver as an in-memory primary-field generator for SFF, replacing the
current read-primary-field-from-file workflow of the `-E` job.

---

## 1. Motivation and goals

The Secondary Field Formulation (SFF) solves for a *scattered* field `dE` driven by an
interior source `−iωμ₀(σ − σ₁ᴅ)·E₁ᴅ`, with the total field reported as `E = E₁ᴅ + dE`.
Today the primary field `E₁ᴅ` must be precomputed by a separate program
(ModEM-global1d) and written to a `.esoln` file, then read back by ModEM's `-E` job
(`read_solnVectorMTX` into `eAllPrimary`). This design **removes that file round-trip**:
ModEM computes `E₁ᴅ` internally, per period and per mode, using the vendored global1d
solver, driven by a spherical-harmonic source `.prm` file.

Goals:

1. **Near-term:** extend the existing `-E` (SECONDARY_FIELD) job so the primary field is
   computed internally from a source `.prm` given on the command line, bypassing primary
   `E₁ᴅ` file writing/reading. Behavior otherwise identical to today's `-E`.
2. **Long-term:** promote `SFF = .true.` (namelist) to a general modifier available to
   *all* jobs, including `-F` (FORWARD) and `-I` (INVERSE), with the source `.prm`
   optional (defaulting to the standard spherical MT source when omitted).
3. Use ModEM's existing MPI parallelization over **periods × modes**.
4. Write the total-field output **independently per period and per mode**, as global1d
   does.
5. Behave correctly and predictably under **Cartesian** builds (where internal spherical
   primary computation is not meaningful): allow a user-supplied primary field, but
   refuse to internally compute a spherical source, with a clear diagnostic.

Non-goals (this phase): changing the secondary solve, the data functionals, the inverse
algorithms, or the physics of the interior-source construction.

---

## 2. Background: how `-E`/SFF works today

Traced across the current tree:

- **Job dispatch.** `-E` → `SECONDARY_FIELD = 'E'` ([UserCtrl.f90:10](../f90/UserCtrl.f90#L10)).
  Argument list ([UserCtrl.f90:386-443](../f90/UserCtrl.f90#L386-L443)):
  `-E rFile_Model rFile_Model1D rFile_EMsoln1D rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl rFile_EMrhs]`.
  `SECONDARY_FIELD` shares the `FORWARD` case block in the main driver
  ([Mod3DMT.f90:192-204](../f90/Mod3DMT.f90#L192-L204)).
- **Primary field ingestion.** When `PRIMARY_E_FROM_FILE` is set, both the MPI and serial
  paths read the whole primary solution set into `eAllPrimary` (a `solnVectorMTX_t`) and
  the 1D model into `sigmaPrimary` (a `modelParam_t`), on **every node**
  ([Mod3DMT.f90:142-160](../f90/Mod3DMT.f90#L142-L160)).
- **Transmitter typing.** A transmitter's `Tx_type` (`'MT'`, `'CSEM'`, `'SFF'`, `'TIDE'`,
  `'GLOBAL'`) is set from the **data file** transmitter-type block
  ([DataIO_ASCII.f90:683](../f90/3D_MT/DataIO_ASCII.f90#L683),
  [transmitters.f90:84-88, 339-340](../f90/3D_MT/DICT/transmitters.f90#L339-L340)).
- **Interior source.** For `Tx_type=='SFF'`, `fwdSetup` builds the edge anomaly
  `condAnomaly = ModelParamToEdge(σ) − ModelParamToEdge(σ_primary)`, then
  `b0 = −ISIGN·iωμ₀·condAnomaly·E_P` with `E_P = eAllPrimary%solns(iTx)%pol(iMode)`;
  boundary conditions are zeroed
  ([ForwardSolver.f90:274-284, 734-747](../f90/3D_MT/ForwardSolver.f90#L734-L747)).
- **Secondary solve + total field.** `fwdSolve` zeroes `e0`, solves `dE`, then adds the
  primary back: `e0%pol(iMode) = E_P + dE`
  ([ForwardSolver.f90:849-867](../f90/3D_MT/ForwardSolver.f90#L849-L867)).
- **MPI granularity.** Jobs are distributed per `(transmitter, polarization)` pair, i.e.
  periods × modes (`Sub_MPI.f90`, job count `Σ solns(iTx)%nPol`). Mode is carried by
  `e0%Pol_index`, not array order, and `b0` uses the same fake index.
- **Dormant scaffolding.** `ctrl%SFF`, `primary_field`, `primary_field_file` namelist
  keys already exist ([UserCtrl.f90:88-92, 1047-1086](../f90/UserCtrl.f90#L1047-L1086)),
  currently unused.

**Key equivalence used throughout this design:** the in-memory `eAllPrimary` we intend to
build is, by construction, exactly what `read_solnVectorMTX` would load from the `.esoln`
that global1d already knows how to write. This is both the port target and the primary
validation anchor.

## 3. Coordinate system: compile-time, and the guard

The coordinate system is a **compile-time parameter** `gridCoords`, from
[math_constants.f90:70-71](../f90/UTILS/math_constants.f90#L70-L71)
(`CARTESIAN='Cartesian'`, `SPHERICAL='Spherical'`), fixed by which grid module is
compiled: `GridCalcS.f90` sets `gridCoords = SPHERICAL`
([GridCalcS.f90:38](../f90/3D_MT/GridCalcS.f90#L38)), `GridCalc.f90` sets
`gridCoords = CARTESIAN` ([GridCalc.f90:31](../f90/3D_MT/GridCalc.f90#L31)). The
`BUILD_SPHERICAL` CMake flag selects the spherical build and appends `_SPH` to the
executable name ([CMakeLists.txt:12-16](../f90/CMakeLists.txt#L12-L16)).

Consequence: **internal primary computation is available only in the `_SPH` build**, and
the guard is a single compile-time-resolved runtime check on `gridCoords`. See §7.11.

---

## 4. Requirements

| # | Requirement | Phase |
|---|---|---|
| R1 | `-E` computes `E₁ᴅ` internally from a source `.prm`, bypassing primary-field file I/O | near-term |
| R2 | Primary conductivity `σ₁ᴅ`: user-supplied model **or** default = 3D model averaged at each depth | near-term |
| R3 | Source `.prm`: user-supplied **or** default = standard spherical MT source (degree-1, two polarizations) | near-term (default may be phase 2) |
| R4 | Parallelize primary computation over periods × modes using existing MPI | near-term |
| R5 | Output total field independently per (period, mode), global1d-style | near-term |
| R6 | `SFF=.true.` namelist enables SFF for `-F` and `-I` as well | long-term |
| R7 | Cartesian build: allow user-supplied primary field; refuse internal spherical computation with a clear error | near-term (guard), full CSEM-style primary later |

---

## 5. Architecture overview

```
                         data file (periods, Tx_type='SFF', modes)
                                        │
   source .prm ───────┐                │
   (or default MT)    ▼                ▼
   σ₁ᴅ model ──► ┌──────────────┐   txDict ──► MPI distribute (Tx × Pol)
   (or depth-    │  Primary1D   │                    │
    averaged)    │  module      │        ┌───────────┴───────────┐
                 │              │        ▼                       ▼
   ModEM grid ──►│ compute_     │     worker k               worker m
                 │ primary1D()  │     (its Tx,Pol subset)    (its subset)
                 └──────┬───────┘        │                       │
                        │  fills eAllPrimary%solns(iTx)%pol(iMode) in memory
                        ▼                ▼                       ▼
              existing SFF fwdSetup/fwdSolve  (condAnomaly, b0, dE, E=E₁ᴅ+dE)
                        │
                        ▼
              per-(period,mode) total-field output
```

The only structural change to the solve path is the *provenance* of `eAllPrimary`:
computed in-memory instead of read from `rFile_EMsoln1D`. `condAnomaly`, `b0`, the
secondary solve, and `E = E₁ᴅ + dE` are unchanged.

---

## 6. Detailed design

### 6.1 New module `Primary1D` (`f90/3D_MT/PRIMARY1D/Primary1D.f90`)

Wraps the vendored global1d physics and exposes a single high-level entry point. Proposed
interface (names indicative):

```fortran
module Primary1D
  ! Depends on: ModEM grid_t, modelParam_t, solnVectorMTX_t, txDict,
  !             and vendored global1d: field1d / field1d_s2 / output_convention
  implicit none
  public :: compute_primary1D, set_primary_source, primary1D_available

contains

  !> True only in the spherical build (gridCoords==SPHERICAL).
  logical function primary1D_available()

  !> Fill eAllPrimary in memory for the transmitters/modes required.
  !!  grid        : ModEM computational grid (defines the edges we solve on)
  !!  sigma1D     : primary (1D) conductivity model  [modelParam_t]
  !!  srcFile     : source .prm path ('' => use built-in default MT source)
  !!  iTxList     : which transmitters this rank must fill (MPI subset); empty => all
  !!  eAllPrimary : (out) solnVectorMTX_t, edge E-field per (iTx, iMode)
  subroutine compute_primary1D(grid, sigma1D, srcFile, iTxList, eAllPrimary)
end module
```

Internally, `compute_primary1D`:

1. Translates `grid` → global1d grid object (§6.2).
2. Translates `sigma1D` → global1d layered `earth` (§6.3).
3. Reads/loads the source `.prm` into per-(period,mode) coefficient blocks (§6.4),
   matching them against `txDict` periods and `nPol`.
4. For each required `(iTx, iMode)`: calls `sourceField1d[_s2](earth, lmax, coeff_block,
   period, g1d_grid, h1d, e1d)`, applies the chosen output convention (§6.5), and maps the
   resulting edge cvector into `eAllPrimary%solns(iTx)%pol(iMode)` (§6.6).

### 6.2 Grid translation (translate, do not interpolate)

Build the global1d grid **from ModEM's `grid`** so the analytic field is evaluated on
exactly ModEM's edges. This eliminates interpolation error and makes the cvector →
solnVector mapping a fixed axis relabel (§6.6). Mapping of quantities:

- ModEM `grid` spherical axes (x=latitude/colatitude, y=longitude, z=depth/radius) →
  global1d `grid%th`, `grid%ph`, `grid%r`, `grid%dr`, air-layer count.
- Requires the air layers and radial discretization to be **identical** on both sides
  (they are, since we build one grid from the other).

Open sub-question: whether to construct a full global1d `grid` derived type or refactor
`sourceField1d` to accept ModEM's `grid_t` directly. Preference: a thin adapter that
populates the global1d grid type, to keep the vendored physics untouched.

### 6.3 Primary conductivity `σ₁ᴅ` (supplied or depth-averaged)

Two sources, in priority order:

1. **User-supplied** (`rFile_Model1D` today; a namelist key long-term): read as
   `modelParam_t` → global1d `earth` layered profile (log-conductivity per layer). This is
   also the `sigmaPrimary` used for `condAnomaly`, guaranteeing consistency.
2. **Default = lateral average of the 3D model at each depth.** When no 1D model is given,
   derive `σ₁ᴅ(z)` by averaging the 3D `σ` over all `(x,y)` cells within each earth layer.
   Design decisions to fix:
   - **Averaging space:** arithmetic in log₁₀σ (i.e. geometric mean of σ) is the natural
     choice for a log-parametrized model and is what "background halfspace/layered" usually
     means; flagged for confirmation.
   - **Weighting:** area/volume-weighted vs. plain cell average. Area-weighted is more
     defensible on a stretched grid.
   - Produce a `modelParam_t` that is laterally uniform (so `sigmaPrimary` and the global1d
     `earth` see the same profile), then proceed exactly as case 1.

Both cases converge to: one `sigmaPrimary` model feeding **both** the global1d earth and
the SFF anomaly.

### 6.4 Source `.prm` (supplied or default)

- **Supplied:** the path given on the command line (near-term) / namelist
  `primary_field_file` (long-term). Parsed with global1d's `read_modelParam`, which now
  supports multi-mode blocks and optional inversion columns. Its `periods × modes` layout
  must match the SFF transmitters and `nPol` in `txDict`; validated at setup with a clear
  error naming any mismatch (reuse global1d's period-group consistency check).
- **Default (no file):** the standard **spherical MT source** — degree-1 external field,
  two polarizations (the global1d Mode1/Mode2 = zonal `l=1,m=0` and non-zonal `l=1,m=±1`
  blocks). Emitted for every period in `txDict`. Implementation options: (a) ship a
  canonical `.prm` and load it; (b) synthesize the coefficient blocks in code. Preference:
  (b), so the default needs no on-disk asset and always matches `txDict`'s periods.

### 6.5 Field convention / normalization

The interior source `−iωμ₀(σ−σ₁ᴅ)E₁ᴅ` and the total `E = E₁ᴅ + dE` require the **genuine
physical primary E-field** in ModEM's internal `e^{−iωt}`, raw SI units — i.e. global1d's
**native** convention (`SUNEGBERT2012`), **not** the `EGBERTKELBERT2012_MODEM` comparison
rescaling (÷ iωμ₀G), which is a bookkeeping transform for matching ModEM's boundary-E
*source normalization* and would double-count here.

The source amplitude ("unit external multipole") sets the absolute scale of the total
field. Since impedances `Z = E/H` are scale-invariant this does not affect `-E`/`-F`/`-I`
data functionals, but it **does** affect the reported total-field magnitudes (R5). The
chosen normalization must therefore be documented in the output header (reuse global1d's
convention-metadata header), so downstream tools interpret the field scale correctly.

### 6.6 Filling `eAllPrimary` (cvector ↔ solnVector)

global1d produces H,E as edge cvectors in its own `(phi, theta, r)` staggering. ModEM's
`solnVector%pol(iMode)` is an edge cvector in `(x=lat, y=lon, z=depth)`. The mapping is the
**same axis relabel already implemented** in global1d's MODEM `.esoln` writer:
`%x(phi/east) → Ey`, `%y(theta/north) → Ex`, with the theta/r index handling that
`apply_output_convention` performs. Because we solve on ModEM's own grid (§6.2), this is a
pure index/label remap with no resampling. Factor this remap into a small shared routine so
the in-memory fill and the existing file writer cannot drift apart.

### 6.7 Command-line interface — near-term (`-E`)

Extend the `-E` argument list so the third slot may be **either** a primary `.esoln`
(current behavior) **or** a source `.prm` (new). Two candidate encodings:

- **(A) File-type sniff:** keep the positional list; detect whether
  `rFile_EMsoln1D`/that slot is a `.prm` source vs. an `.esoln` solution (by extension
  and/or header magic). `.prm` ⇒ internal computation; `.esoln` ⇒ current file read.
- **(B) Explicit flag:** add the `SFF=.true.` namelist as the switch that reinterprets the
  same slot as a source `.prm` and triggers internal computation.

Preference: **(B)**, because it is exactly the long-term mechanism (R6) and avoids
magic-based dispatch; (A) can be a convenience on top. Under (B), the near-term `-E` call
becomes:

```
# internal primary (spherical build), source from .prm, SFF=.true. in namelist:
-E rFile_Model [rFile_Model1D] source.prm rFile_Data wFile_Data ...
#   rFile_Model1D optional: omit => depth-averaged default (R2)
```

`σ₁ᴅ` supplied-vs-default (R2) is resolved by whether the 1D-model slot is present.

### 6.8 `SFF = .true.` namelist — long-term (`-F`, `-I`, …)

Promote `ctrl%SFF` to a general modifier read from the optional namelist
([UserCtrl.f90:1047-1086](../f90/UserCtrl.f90#L1047-L1086), already parsed). When set:

- Any job that runs a forward solve (`-F`, `-I`, `-J`, …) treats its transmitters as SFF:
  compute the internal primary once per `(period, mode)`, build the interior source, solve
  the secondary, report the total field / derived data.
- `primary_field_file` (namelist) supplies the source `.prm`; absent ⇒ default MT source
  (§6.4).
- The transmitter typing question: today `Tx_type='SFF'` comes from the data file. For a
  global `SFF=.true.` modifier we must decide whether the flag (a) forces all `MT`
  transmitters to be treated as `SFF`, or (b) still requires the data file to mark them.
  Preference: `SFF=.true.` makes plane-wave/`MT` transmitters use the SFF path globally,
  so existing MT data files work unchanged. **Decision needed.**

This is deliberately staged **after** the `-E` path is proven, to avoid perturbing the
inverse solvers until the primary generator is trusted.

### 6.9 MPI strategy

The 1D analytic field is cheap and embarrassingly parallel, and each worker already knows
its `(iTx, iMode)` assignment. **Each worker computes its own primary field** for its
assigned transmitters/modes (mirroring today's code, which already reads the full primary
file redundantly on every node) — avoiding large-E-field MPI messaging. Concretely:

- Master distributes `txDict`, `grid`, `sigma1D`, and the (small) source coefficient set —
  the first three are already broadcast; the source `.prm` (or default) is small.
- Replace the "read primary on each node" block ([Mod3DMT.f90:142-160](../f90/Mod3DMT.f90#L142-L160))
  with `call compute_primary1D(...)` on each node, restricted to that node's `iTxList`.
- Serial build: same call, `iTxList = all`.

### 6.10 Output per period/mode

Keep `eAll` (total field) in memory as `solnVectorMTX_t`; write **one file per
(period, mode)**. Lightest path: route SFF total-field output through a per-`(iTx,iMode)`
writer rather than one multi-Tx `.esoln`. Reuse global1d's naming/metadata
(`…T01.mode01.efield`, convention header) and/or ModEM's per-solution NCI writer
(`storeSolnsInFile`) — **decision needed** on which container format to standardize on.

### 6.11 Cartesian behavior and the guard

Under the Cartesian build (`gridCoords == CARTESIAN`):

- **Refuse internal spherical primary computation.** If SFF-with-internal-primary is
  requested (source `.prm` given / `SFF=.true.` without a primary file), stop early with a
  precise message, e.g.:
  *"Internal 1D primary-field computation from a spherical-harmonic source is only
  available in the spherical build (Mod3DMT_*_SPH). This is a Cartesian build. Supply a
  precomputed primary field via `-E … rFile_EMsoln1D …`, or rebuild with
  `-DBUILD_SPHERICAL=ON`."*
- **Still allow a user-supplied primary field** (the existing CSEM/SFF-from-file path):
  a general primary field on the Cartesian grid remains valid input; only the *internal
  spherical source* is disallowed.

Because `gridCoords` is a compile-time parameter, `primary1D_available()` is effectively a
compile-time constant; the guard collapses to a trivial branch and cannot be mis-set at
runtime. The guard must fire **before** any grid/source translation, at argument/namelist
validation time, so Cartesian users get an immediate, actionable error.

---

## 7. Build integration

Vendor the minimal global1d source set into ModEM's build (proposed
`f90/3D_MT/PRIMARY1D/`): `field1d.f90`, `field1d_s2.f90`, `output_convention.f90`, plus the
source-`.prm` reader. **Module-collision risk** — both trees define `GridDef`,
`sg_vector`, `math_constants`, and a `read_modelParam`/`ModelSpace`. Resolution policy:

- **Reuse ModEM's** `grid`, `sg_vector`, `math_constants` (this is what makes §6.2's
  "translate, don't interpolate" clean). The vendored physics must be adapted to ModEM's
  `cvector`/`grid_t` at its boundary, or wrapped by the `Primary1D` adapter.
- **Namespace** anything genuinely global1d-specific that would clash (e.g. its layered
  `read_modelParam` for the source `.prm`) under a `primary1d_` prefix or a dedicated
  module, so it does not shadow ModEM's model reader.
- Add the new directory to both the CMake build and the classic Makefile builds; gate the
  physics compilation on the spherical build if the vendored code assumes spherical
  geometry.

---

## 8. Error handling and edge cases

1. **Cartesian + internal source** → hard stop, message as in §6.11.
2. **Source `.prm` periods/modes ≠ `txDict`** → hard stop naming the first mismatch.
3. **Grid mismatch** (air layers / radial discretization inconsistent between supplied
   primary and 3D grid) → cannot occur when we build the grid from ModEM's (§6.2); if a
   *file* primary is used (Cartesian or explicit), keep the existing checks.
4. **`nPol` per SFF transmitter** must equal the number of modes in the source block for
   that period.
5. **Depth-averaging with mixed air/earth layers** → average earth layers only; air is
   fixed/insulating.
6. **Regression safety:** with `SFF=.false.` and no source `.prm`, every existing job path
   must be byte-for-byte unchanged.

---

## 9. Validation plan (staged)

1. **Equivalence anchor (unit).** For a fixed model/grid/source, assert the in-memory
   `eAllPrimary` equals `read_solnVectorMTX` of the `.esoln` global1d writes standalone,
   component-by-component to machine precision. Isolates the port from physics.
2. **Grid-translation check.** global1d on the ModEM-derived grid == global1d on its own
   `.grd` at shared edges.
3. **End-to-end vs. current `-E`.** File-based `-E` (global1d `.esoln` → `-E`) vs. the new
   internal-primary `-E` on the same problem: total-field output and predicted data match.
4. **Depth-averaged default.** For a laterally-uniform 3D model, the depth-averaged `σ₁ᴅ`
   reproduces the explicit 1D model exactly; anomaly → 0 ⇒ total field == pure primary.
5. **MPI vs. serial** identical; scaling over periods × modes.
6. **Physics cross-check.** Reuse `testing/test_vs_modem_1D` halfspace/Cartesian
   benchmarks; SFF total field over a 1D model reduces to the primary and matches the
   standalone global1d / analytic values (the `c=iωμ₀G` framework already characterizes
   the relationship).
7. **Cartesian guard.** Cartesian build + internal-source request → clean, specific error;
   Cartesian build + file primary → works as today.
8. **Regression.** `-F`/`-E`/`-I` unchanged when `SFF=.false.`.

---

## 10. Milestones / phasing

**Phase 0 — decisions.** Resolve the open questions in §11.

**Phase 1 — vendor + compile.** Land global1d sources under `f90/3D_MT/PRIMARY1D/`; resolve
module collisions; spherical build compiles with the new module linked but unused.

**Phase 2 — `Primary1D` core (serial).** Grid + model translation, source loading,
`eAllPrimary` fill. Pass validation #1–2.

**Phase 3 — wire into `-E`.** `SFF=.true.` reinterprets the source slot; bypass the file
read; guard Cartesian. Pass #3, #7.

**Phase 4 — defaults.** Depth-averaged `σ₁ᴅ` (#4) and built-in default MT source.

**Phase 5 — MPI.** Worker-local compute; pass #5.

**Phase 6 — output + cross-checks.** Per-(period,mode) writer (#R5); physics checks #6;
regression #8.

**Phase 7 — long-term `SFF` modifier.** Extend to `-F`/`-I` (§6.8) once the `-E` path is
trusted.

---

## 11. Open questions / decisions needed

1. **Command switch (near-term):** confirm `SFF=.true.` namelist as the trigger that
   reinterprets the `-E` source slot (preferred), vs. file-extension sniffing.
2. **Depth-averaging** space (log vs linear) and weighting (area/volume vs plain).
3. **Default source** synthesized in code vs. shipped `.prm`.
4. **Output container** for per-(period,mode) total field: global1d `.efield`/`.esoln`
   naming vs. ModEM NCI per-solution files.
5. **`SFF=.true.` transmitter semantics (long-term):** does the flag force `MT`
   transmitters onto the SFF path, or must the data file still mark `Tx_type='SFF'`?
6. **Grid adapter shape:** populate a global1d grid type vs. teach `sourceField1d` to read
   ModEM's `grid_t` directly.
7. **S1 vs S2 solver** for the primary (both validated; S2 is natively in the paper's
   convention). Default choice + whether to expose it as an option.

---

## 12. Reference map (verified locations)

- Job letters / `-E` usage: `f90/UserCtrl.f90:9-10, 237-245, 386-443`
- Primary-from-file ingestion: `f90/Mod3DMT.f90:142-160`
- SFF interior source: `f90/3D_MT/ForwardSolver.f90:274-284, 734-747`
- SFF secondary solve + total field: `f90/3D_MT/ForwardSolver.f90:849-867`
- Transmitter typing from data file: `f90/3D_MT/DataIO_ASCII.f90:683`;
  `f90/3D_MT/DICT/transmitters.f90:84-88, 339-340`
- MPI per-(Tx,Pol) distribution: `f90/3D_MT/Sub_MPI.f90`
- Coordinate parameter: `f90/UTILS/math_constants.f90:70-71`;
  `f90/3D_MT/GridCalcS.f90:38`; `f90/3D_MT/GridCalc.f90:31`
- Spherical build flag: `f90/CMakeLists.txt:12-16`
- Dormant SFF namelist scaffolding: `f90/UserCtrl.f90:88-92, 1047-1086`
- global1d entry points: `ModEM-global1d/FWD1D.f90:400-455`;
  `EARTH/FWD/field1d.f90` (`sourceField1d`), `field1d_s2.f90` (`sourceField1d_s2`),
  `EARTH/FWD/output_convention.f90`
