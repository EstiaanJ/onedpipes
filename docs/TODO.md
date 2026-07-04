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

- [ ] Keep MacCormack development on its own branch until it passes the
  same milestone gates as the Lax–Wendroff baseline for the implemented
  scope.
- [ ] Add an interior-solver selection API (`SolverKind` or equivalent)
  without changing the conservative state, gas-property, boundary, or
  global-timestep contracts.
- [ ] Implement MacCormack predictor-corrector with the same source
  splitting, artificial dissipation, positivity floors, and diagnostic
  reports used by Lax–Wendroff.
- [ ] Promote MacCormack to a user-selectable model only after its
  milestone validation cases pass against independent references.

## Shared Full-Slice Tests And GUI

- [ ] Parameterize full-slice validation cases over solver method wherever
  the physical setup is solver-independent.
- [ ] Reuse case setup, boundary configuration, probes, reference
  calculations, and GUI snapshots for both solvers; allow
  method-specific tolerances only where the numerical behavior justifies
  it.
- [ ] Add cross-solver parity checks after both solvers pass a case's
  independent reference.
- [ ] Keep one shared GUI/viewer with solver selection in the run/model
  configuration; do not fork the GUI by solver method.
