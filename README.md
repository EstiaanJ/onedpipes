# 1D Duct Gas-Dynamics Solver

A quasi-1D unsteady compressible flow solver for ducts (pipes) of varying
cross-sectional area. This is the current, and only, scope of work.

This is an engineering-prediction tool, not an acoustics or
shock-capturing tool. Sound quality, exact shock structure, and
discontinuity preservation are explicitly out of scope (see
`docs/PRODUCT_SPEC.md`).

This library is intended for engine-simulation

## Current implementation

The checked-in Rust solver currently covers milestones 1 through 5 from
`docs/PRODUCT_SPEC.md`:

- Solves the quasi-1D Euler equations in each duct using a selectable
  explicit interior solver: two-step Lax–Wendroff (Richtmyer) or
  MacCormack predictor-corrector.
- Provides pluggable closed-end and open-end boundary models with pulse
  reflection validation.
- Provides a constant-pressure multi-pipe junction model with coupled
  model stepping and mass/energy conservation validation.
- Provides a quasi-steady valve/orifice model with coupled model
  stepping and hand-computed compressible discharge validation.
- Includes a Sod shock tube validation against an independent exact
  Riemann solution, with tolerances chosen for the expected smearing of
  the unlimited LW scheme.
- Tracks temperature-dependent gas properties (γ(T), cp(T)) for a single
  effective gas for now, structured so per-species tracking can be added
  later without restructuring the solver.
- Includes a minimal artificial-dissipation term for numerical stability.
- Supports variable-area duct profiles in the interior solver path, with
  area-weighted updates and well-balance tests for Venturi-shaped
  geometry.

All v1 validation milestones are implemented for the current models.
Steady Venturi-effect validation still needs pressure/stagnation-state
or mass-flow inlet/outlet boundaries; Helmholtz validation still needs
connected side-branch/cavity support.
Wall heat transfer remains planned milestone work.

## Milestone TODO

- [x] Milestone 1: bare closed-closed duct wave propagation validated
  against analytic organ-pipe frequency.
- [x] Milestone 2: closed-end and open-end boundary reflection behavior
  validated independently.
- [x] Milestone 3: multi-pipe junction conserves mass and energy across
  2-3 connected pipes.
- [x] Milestone 3 subtask: constant-pressure junction core solves shared
  pressure and balanced port flux diagnostics.
- [x] Milestone 3 subtask: attach junction coupling to multi-duct model
  stepping.
- [x] Milestone 3 subtask: validate 2-pipe, 3-pipe, and model-coupled
  junction conservation in `tests/`.
- [x] Milestone 4: valve/orifice boundary matches hand-computed
  compressible orifice discharge.
- [x] Milestone 4 subtask: validate subcritical and choked pressure
  ratios in `tests/`.
- [x] Milestone 4 subtask: attach orifice flux coupling to multi-duct
  model stepping.
- [x] Milestone 5: Sod shock tube wave structure and speeds match the
  analytic reference within expected LW smearing.
- [x] Milestone 5 subtask: validate density, velocity, pressure, and
  shock position against an exact Riemann solution in `tests/`.

## Solver TODO

- [x] Introduce a shared interior-solver selection API while preserving
  current Lax–Wendroff behavior.
- [x] Add MacCormack predictor-corrector as a peer solver using the same
  conservative state, ghost-cell boundaries, global timestep,
  artificial-dissipation hook, and positivity diagnostics.
- [x] Parameterize milestone 1 and 2 full-slice validation cases over
  solver method.
- [ ] Parameterize future full-slice validation cases over solver method
  where possible; both solvers must pass independent references before
  cross-solver comparisons are treated as proof/regression evidence.
- [x] Add solver selection to the shared GUI/viewer instead of creating a
  separate MacCormack viewer.
- [ ] Keep future branch work coordinated: Lax–Wendroff and MacCormack
  milestone work should continue against the same public model/run
  concepts.

## What it doesn't do (yet)

- No Method of Characteristics (MoC) boundaries — planned upgrade, see
  `docs/DECISIONS.md`.
- No multi-species transport yet — single effective gas for now.
- No shock-capturing / TVD limiters — not needed for this application
  (see `docs/PRODUCT_SPEC.md` for why).
- No built-in 0D engine components (cylinders, plenums, combustion);
  external coupling points are exposed for a host engine simulator.

## Testing

Testing is a first-class goal, not an afterthought. Every boundary type
and gas-property model ships with a unit test, and every milestone in
`docs/PRODUCT_SPEC.md` has a corresponding validation case checked into
`tests/` before the next milestone starts (e.g. organ-pipe resonant
frequency, Sod shock tube, junction mass/energy conservation, orifice
discharge flow). See `AGENTS.md` § Build order.

## Repository map

```
README.md                 You are here
AGENTS.md                 Instructions for AI coding agents working on this repo
docs/PRODUCT_SPEC.md       What is being built and why
docs/ARCHITECTURE.md      How it is technically structured
docs/DECISIONS.md          Key decisions, rationale, and future upgrade paths
docs/PUBLIC_API.md         Public library API for engine-simulator integration
docs/TODO.md               Branch-level solver and validation TODOs
src/                       Solver source (created during implementation)
tests/                     Validation cases (organ-pipe resonance, Sod tube, etc.)
```

## Getting started

See `docs/ARCHITECTURE.md` for the module layout and `docs/PRODUCT_SPEC.md`
for the milestone plan. See `docs/PUBLIC_API.md` for the integration API.
The completed build order is: bare-pipe wave propagation → closed/open
end boundaries → junctions → valve/orifice boundary → Sod shock tube.
