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

## Interior scheme: two-step Lax–Wendroff (Richtmyer)

Per duct, per step:

1. **Predictor** (half-step, cell-face midpoints):
   `U*_{i+1/2} = ½(U_i + U_{i+1}) − (Δt/2Δx)(F_{i+1} − F_i)`
2. **Corrector** (full step, cell centers):
   `U_i^{n+1} = U_i^n − (Δt/Δx)(F*_{i+1/2} − F*_{i-1/2})`

Source terms (area change, friction, heat transfer) are applied via
operator splitting: advance the homogeneous LW step, then integrate
sources over `Δt` for each cell. Keep source integration explicit
(forward Euler) for v1 — revisit only if it causes stability problems at
Δt set by the CFL condition below.

## Boundaries (v1 — simple, not MoC)

All boundary types are implemented as **pluggable objects attached to
duct ends**, each responsible for supplying the ghost-cell state(s) the
interior LW stencil needs. This interface must not change when MoC
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
value, tune empirically). This is the minimum viable stabilizer for LW's
dispersive ringing — no flux limiter, no TVD switch, no shock detector.
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
  duct.*                cell arrays (U, A(x)), interior LW step, dissipation
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
- A milestone's test must pass before starting the next milestone's
  implementation — see `AGENTS.md` § Build order.

## Upgrade seams (do not design these away)

- Boundary objects: interior scheme must never special-case a boundary
  type; this is required for the MoC swap-in described in DECISIONS.md.
- Gas properties module: must be swappable for multi-species without
  touching `duct.*`.
- Source-term integration: currently forward Euler; if stiffness becomes
  a problem (fine grids, high heat-transfer coefficients), this is the
  one place an implicit sub-step could be introduced without changing
  the explicit interior wave solver.
