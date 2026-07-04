# Architecture

## Governing equations

Quasi-1D Euler equations, conservative form, area-weighted:

```
∂(A·U)/∂t + ∂(A·F(U))/∂x = S(U, x)
U = [ρ, ρu, ρE]
F = [ρu, ρu² + p, u(ρE + p)]
S = area-change pressure term + wall friction + wall heat transfer
```

`A(x)` is duct cross-sectional area (piecewise or continuous per duct).
Conservative variables are the stored state in every cell; primitives
(`p, u, T, a`) are always derived, never stored as primary state.
The finite-volume update advances area-weighted conservative quantities
internally (`A·U`) and divides back to cell-local `U` after each step.
For variable-area ducts, face fluxes use `A_face·F`, and the pressure
area source `p·dA/dx` is included inside each explicit predictor/corrector
so a static uniform-pressure duct remains well balanced.

## Interior solver interface

The duct interior advance is a selectable solver model. Every interior
solver must use the same contract:

- Conservative cell states `U = [ρ, ρu, ρE]` are the only primary state.
- Boundary objects supply ghost-cell states; solver implementations do
  not special-case boundary types.
- The model supplies one global explicit `Δt` for the whole step.
- The solver returns the same `StepReport` diagnostics for clipping and
  fallback behavior.
- Source terms, artificial dissipation, and positivity enforcement are
  applied through shared hooks unless a future decision explicitly
  documents why a solver needs different handling.
- Constant-area and variable-area ducts go through the same solver
  selection path. A solver is not Venturi-ready until it passes static
  variable-area well-balance and area-weighted mass conservation tests.

The GUI and validation viewers select a solver through model/run
configuration. They must not fork into separate Lax–Wendroff and
MacCormack viewers.

## Baseline solver: two-step Lax–Wendroff (Richtmyer)

Per duct, per step:

1. **Predictor** (half-step, cell-face midpoints):
   `U*_{i+1/2} = ½(U_i + U_{i+1}) − (Δt/2Δx)(F_{i+1} − F_i)`
2. **Corrector** (full step, cell centers):
   `U_i^{n+1} = U_i^n − (Δt/Δx)(F*_{i+1/2} − F*_{i-1/2})`

Friction and heat-transfer source terms are applied via operator
splitting: advance the interior-solver step, then integrate sources over
`Δt` for each cell. Keep source integration explicit (forward Euler) for
v1 — revisit only if it causes stability problems at Δt set by the CFL
condition below. The geometry pressure source is part of the explicit
area-weighted predictor/corrector itself, not a post-step split source,
because it must balance `A·p` flux gradients at rest.

## Solver: MacCormack predictor-corrector

MacCormack is the approved second interior solver model. It must match
the same milestone validation goals as the Lax–Wendroff baseline to
remain exposed as a normal user-selectable model.

Per duct, per step for the homogeneous conservative update:

1. **Predictor** (forward difference):
   `Ū_i = U_i^n − (Δt/Δx)(F_{i+1}^n − F_i^n)`
2. **Corrector** (backward difference and average):
   `U_i^{n+1} = ½(U_i^n + Ū_i − (Δt/Δx)(F̂_i − F̂_{i-1}))`
   where `F̂_i = F(Ū_i)`.

Use the same conservative variables, area-weighted flux/source handling,
ghost-cell boundary interface, global timestep, source splitting,
artificial dissipation, and positivity/fallback policy as the
Lax–Wendroff solver. Do not add limiters, Riemann solvers, implicit
stepping, or shock detectors as part of the MacCormack branch; those
remain out of scope unless a separate decision changes the product scope.

## Boundaries (v1 — simple, not MoC)

All boundary types are implemented as **pluggable objects attached to
duct ends**, each responsible for supplying the ghost-cell state(s) the
interior stencil needs. This interface must not change when MoC
boundaries are added later — only the internals of each boundary object
change.

- **Closed end**: reflect velocity (u → −u) in the ghost cell, mirror
  ρ and p. Zero mass flux at the wall.
- **Open end (to ambient)**: fix ambient static pressure in the ghost
  cell; extrapolate density/velocity from the interior (zeroth-order
  extrapolation for v1). Accept that this is a simplification — it
  will slightly over/under-reflect versus a true non-reflecting or
  characteristic boundary. Acceptable for v1 per PRODUCT_SPEC.
- **Multi-pipe junction**: constant-pressure junction model. All pipe
  ends at the junction share one instantaneous pressure; mass and energy
  are balanced across the junction each step. No pressure-loss
  coefficients in v1 (add later if junction losses matter).
- **Valve/orifice**: quasi-steady compressible orifice equation using a
  discharge coefficient `Cd` (constant or simple lookup table vs. flow
  area). Given upstream/downstream stagnation states and an
  instantaneous flow area, compute mass flow and impose it as the
  boundary flux for that duct end — see DECISIONS.md.

Ghost cells are the mechanism for all boundary types: each boundary
object writes one or two ghost cells per step before the interior
predictor/corrector runs.

## Artificial dissipation (simplest option)

Add a **Lapidus-type or simple second-difference artificial viscosity**
term to the flux or directly to `U` after the corrector step:

```
U_i ← U_i + ε · (U_{i+1} − 2U_i + U_{i-1})
```

with a single global scalar `ε` (start around 0.05–0.2 of the CFL-limited
value, tune empirically). This is the minimum viable stabilizer for
dispersive ringing in LW and MacCormack — no flux limiter, no TVD switch,
no shock detector.
Wall heat transfer (below) also contributes physical damping and reduces
how much artificial dissipation is needed — tune the two together.

## Time stepping

- **Global explicit Δt**, shared by every duct in the model — no local
  time stepping (this is a time-accurate problem; local time stepping
  would corrupt wave timing).
- Each step: compute `Δt_i = CFL · Δx_i / (|u_i| + a_i)` in every cell of
  every duct, take the minimum, apply a **10% safety margin**:
  `Δt = 0.9 · min(Δt_i)`.
- Recompute every step (wave speeds change as the flow evolves).

## Gas properties

- Single effective gas for v1 (no species transport yet).
- **Temperature-dependent γ(T) and cp(T)**: use a polynomial or table fit
  (e.g. NASA-polynomial-style or a simple piecewise-linear fit) valid
  over the intake-to-exhaust temperature range. Compute `a = sqrt(γ(T)·R·T)`
  and pressure from the ideal gas law each step from the derived
  temperature.
- Structure the gas-property module as a swappable interface
  (`get_gamma(T)`, `get_cp(T)`, `get_R(composition)`) so that later
  swapping in a mixture-averaged multi-species version doesn't touch
  the interior solver — species mass fractions would just become
  additional advected scalars feeding `get_R`/`get_cp`.

## Wall heat transfer

Simple correlation-based heat loss per cell per step, based on local gas
temperature, wall temperature (fixed or slowly-varying), duct diameter,
and a standard pipe-flow Nusselt-number correlation (e.g. Reynolds-based
Colburn/Dittus-Boelter style). Implemented as a source term in `S(U,x)`
(see Governing equations). Two purposes: physical realism, and damping
of numerical ringing left over after artificial dissipation — tune both
together, not one against the other in isolation.

## Positivity and robustness

- After every corrector step, clip density and pressure to a small
  positive floor; log a warning (not silent) if a floor is hit.
- Before the corrector consumes predictor output, sanity-check predictor
  states for negative density/pressure; if hit, fall back to a
  first-order (upwind) flux for that face for that step only.
- These are safety nets, not a substitute for correct dissipation tuning
  — repeated floor-hits indicate a tuning or boundary bug, not a
  case to silently absorb.

## Module layout

```
src/
  gas_properties.*     get_gamma(T), get_cp(T), get_R(), ideal-gas helpers
  duct.*                cell arrays (U, A(x)), boundary attachment,
                        calls selected interior solver
  solvers/
    mod.*               SolverKind / solver interface
    lax_wendroff.*      two-step LW implementation
    mac_cormack.*       MacCormack implementation
  boundaries/
    closed_end.*
    open_end.*
    junction.*
    valve_orifice.*
  timestep.*             global CFL Δt computation
  model.*                assembles ducts + boundaries, runs the loop
tests/
  test_organ_pipe.*      Milestone 1
  test_boundaries.*       Milestone 2
  test_junction.*         Milestone 3
  test_orifice.*          Milestone 4
  test_sod.*              Milestone 5
  test_solver_parity.*    shared full-slice cross-solver comparisons
```

0D components (cylinders, plenums) are out of scope for this project
phase — see `PRODUCT_SPEC.md`. If added later, they should sit above
`model.*` and consume the boundary interfaces, not be threaded into
`duct.*`.

## Testing

Testing is a specific goal, tracked per milestone, not an incidental
byproduct of development:

- Unit tests live next to the code they test (gas property functions,
  each boundary type, the dissipation term, the positivity clip).
- Validation/regression tests live in `tests/`, one file per milestone,
  each asserting against an independently computed reference value
  (analytic frequency, conservation identity, hand-computed orifice
  flow, Sod's exact solution) with an explicit numeric tolerance.
- Full-slice validation tests should use a shared scenario harness
  parameterized by `SolverKind` where possible. The same physical setup,
  boundary configuration, probe definitions, and reference calculations
  should run for both Lax–Wendroff and MacCormack; only tolerances should
  differ when the numerical method justifies it.
- Cross-solver parity tests compare accepted solver models against each
  other after each model has passed the independent reference for the
  scenario. Parity tests are regression evidence, not the primary
  validation reference.
- A milestone's test must pass before starting the next milestone's
  implementation — see `AGENTS.md` § Build order.

## Upgrade seams (do not design these away)

- Boundary objects: interior scheme must never special-case a boundary
  type; this is required for the MoC swap-in described in DECISIONS.md.
- Solver interface: model, validation, and GUI code select an interior
  solver by configuration and consume the same snapshots/reports. Do not
  fork the GUI or duplicate full-slice case setup per solver.
- Gas properties module: must be swappable for multi-species without
  touching `duct.*`.
- Source-term integration: currently forward Euler; if stiffness becomes
  a problem (fine grids, high heat-transfer coefficients), this is the
  one place an implicit sub-step could be introduced without changing
  the explicit interior wave solver.
