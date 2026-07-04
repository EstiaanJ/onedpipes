# Project TODO

This file tracks planning work that cuts across the milestone docs. The
canonical product goals remain in `PRODUCT_SPEC.md`; this TODO exists to
coordinate branch-level work.

## Active Lax–Wendroff Milestone Track

- [x] Milestone 1: bare closed-closed duct wave propagation validated
  against analytic organ-pipe frequency.
- [x] Milestone 2: closed-end and open-end boundary reflection behavior
  validated independently.
- [ ] Milestone 3: attach junction coupling to multi-duct model stepping
  and validate 2-3 pipe mass/energy conservation.
- [ ] Milestone 4: valve/orifice boundary matches hand-computed
  compressible orifice discharge.
- [ ] Milestone 5: Sod shock tube wave structure and speeds match the
  analytic reference within expected unlimited-scheme smearing.

## Parallel MacCormack Track

- [x] Keep MacCormack development on its own branch until it passes the
  same milestone gates as the Lax–Wendroff baseline for the implemented
  scope.
- [x] Add an interior-solver selection API (`SolverKind` or equivalent)
  without changing the conservative state, gas-property, boundary, or
  global-timestep contracts.
- [x] Implement MacCormack predictor-corrector with the same source
  splitting, artificial dissipation, positivity floors, and diagnostic
  reports used by Lax–Wendroff.
- [x] Promote MacCormack to a user-selectable model for milestones 1 and
  2 after its validation cases pass against independent references.
- [ ] Keep MacCormack gated per future milestone as junction,
  valve/orifice, and Sod coverage are added.

## Shared Full-Slice Tests And GUI

- [x] Parameterize milestone 1 and 2 full-slice validation cases over
  solver method where the physical setup is solver-independent.
- [x] Reuse case setup, boundary configuration, probes, reference
  calculations, and GUI snapshots for both solvers; allow
  method-specific tolerances only where the numerical behavior justifies
  it.
- [x] Add an organ-pipe cross-solver parity check after both solvers pass
  the independent frequency reference.
- [x] Add current-scope pulse parity checks for centroid position over
  time, probe arrival time, pressure/velocity profile similarity,
  acoustic pressure-velocity relation, closed-pipe mass drift, and a
  pulse parameter sweep over grid size, pulse width, amplitude, and
  artificial viscosity.
- [x] Add an organ-pipe reference sweep over length, pressure/density
  state, grid size, perturbation amplitude, and CFL.
- [ ] Extend cross-solver parity checks to future milestones after both
  solvers pass each case's independent reference.
- [x] Keep one shared GUI/viewer with solver selection in the run/model
  configuration; do not fork the GUI by solver method.

## Variable-Area And Venturi Readiness

- [x] Add variable-area duct configuration with cell and face area
  profiles.
- [x] Advance area-weighted conservative quantities internally for both
  Lax–Wendroff and MacCormack.
- [x] Include the `p·dA/dx` geometry source inside each explicit
  predictor/corrector so a static Venturi-shaped duct remains
  well-balanced.
- [x] Add variable-area tests for static well-balance, area-weighted mass
  conservation, and resolved Venturi throat geometry.
- [ ] Add pressure/stagnation-state or prescribed mass-flow inlet/outlet
  boundaries so steady Venturi effect validation can be run against
  Bernoulli/continuity hand calculations.

## Deferred Validation Cases Requiring New Model Features

- [ ] Helmholtz resonator: requires either a validated side-branch/cavity
  representation or a junction/cavity boundary model. Reference should be
  the standard Helmholtz frequency formula with end-correction terms.
- [ ] Venturi tube/effect: variable-area geometry is available; the
  remaining blocker is pressure/stagnation-state or prescribed mass-flow
  inlet/outlet boundaries for a steady reference case.
- [ ] Steady prescribed mass-flow tests: require inlet/outlet boundary
  conditions that impose mass flow or stagnation state, then validate
  conservation and thermodynamic state against hand calculations.
- [ ] Unsteady forced mass-flow/pulse train tests: require a time-varying
  inlet boundary, with sweeps over forcing frequency, pulse width,
  temperature, and duct diameter/area.
