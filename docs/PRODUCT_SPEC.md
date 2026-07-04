# Product Spec

## Goal

Build and validate a standalone 1D (quasi-1D) unsteady compressible
duct-flow solver with selectable explicit interior solver models:

- **Wave propagation accuracy** — organ-pipe resonant frequency matches
  the analytic value within a few percent.
- **Boundary correctness** — open end, closed end, junction, and
  valve/orifice boundaries each produce physically correct behavior
  (reflection sign/magnitude, mass/energy conservation, discharge flow).
- **Solver-method parity** — Lax–Wendroff and MacCormack run the same
  physical cases through the same model, boundary, and GUI surfaces, and
  agree within method-appropriate tolerances after each has passed the
  independent reference validation for that case.
- **Sod shock tube** — run as a standard verification case for the
  interior solver's conservation and wave-speed behavior. This checks
  that the scheme is *correct*, not that it captures shocks sharply: see
  Explicit non-goals below. Do not expect crisp shock resolution from an
  unlimited second-order method without limiters — some smearing/ringing
  at the discontinuity is expected and acceptable.


## Explicit non-goals

These are deliberately excluded because they don't affect the goals
above and would add cost/complexity for no benefit:

- Exhaust/intake **sound quality or audio synthesis**
- **Shock-capturing accuracy** — a diffusive scheme that smears a shock
  is acceptable here; TVD/WENO limiters and Riemann solvers are out of
  scope
- **Exact preservation of discontinuities/contact surfaces**
- Multi-dimensional flow effects (swirl, in-cylinder CFD)



## Physical scope (v1)

- Quasi-1D compressible flow in ducts of varying cross-sectional area
- User-selectable explicit interior solvers:
  - two-step Lax–Wendroff (Richtmyer), the current baseline
  - MacCormack predictor-corrector, a peer model with the same validation
    goals
- Single effective gas with temperature-dependent γ(T), cp(T)
- Simple wall heat transfer (damping mechanism + physical realism)
- Boundary types: open end, closed end, multi-pipe junction, valve/orifice
  flow restriction

## Explicitly deferred (planned, not v1)

| Feature | Deferred because | Tracked in |
|---|---|---|
| MoC boundary treatment | Simple boundary models are enough to get a working model; MoC is an accuracy upgrade, not a functionality requirement | DECISIONS.md |
| Multi-species transport (air, fuel vapor, inert/EGR, combustion products) | Adds a transport equation per species; not needed for the current goals | DECISIONS.md |

## Testing strategy

Testing is a specific, tracked goal, not incidental:

- Every gas-property function and boundary type ships with a unit test
  alongside the code that introduces it.
- Every milestone below has one or more validation cases checked into
  `tests/` with an explicit numeric tolerance, and the next milestone
  does not start until its predecessor's validation case passes.
- Validation cases compare against an analytic or hand-computed
  reference (resonant frequency formula, known reflection behavior,
  conservation identities, compressible orifice flow equation, Sod's
  published exact solution) — not against the solver's own prior output.
- Full-slice validation cases should be shared and parameterized over
  solver method wherever the physical setup is the same. Solver-specific
  tolerances are allowed for dispersive/diffusive behavior, but the case
  setup, boundary configuration, probes, and reference calculation should
  be common.
- Cross-solver comparison is an additional regression/proof signal, not
  a replacement for independent references. Lax–Wendroff and MacCormack
  should both pass the reference first, then be compared against each
  other for frequency, conserved quantities, wave speeds, and probe
  histories where applicable.

## Milestones (must pass in order per solver)

The Lax–Wendroff branch remains the baseline milestone path. The
MacCormack branch follows the same milestone order; it may remain
user-selectable for milestone scopes it has passed, and future milestone
behavior stays gated until the matching validation passes. Shared
full-slice tests should cover both methods as soon as the MacCormack
implementation reaches the relevant milestone.

1. [x] **Bare duct wave propagation** — closed-closed pipe reproduces the
   analytic organ-pipe resonant frequency within a few percent.
2. [x] **Boundary models** — open end, closed end each independently produce
   physically correct reflection sign/magnitude on a test pulse.
3. [ ] **Junction** — 2–3 pipe junction conserves mass and energy to numerical
   tolerance.
4. [ ] **Valve/orifice boundary** — steady-state mass flow through the
   boundary matches hand-computed compressible orifice flow for a range
   of pressure ratios.
5. [ ] **Sod shock tube** — density/velocity/pressure profiles match the
   analytic Sod solution qualitatively (correct wave structure and
   speeds), within the smearing expected of an unlimited second-order
   explicit scheme.

## Success criteria for v1

- Milestones 1–5 all pass their validation cases in `tests/`.
- Resonance behavior is quantitatively correct without requiring
  shock-capturing or MoC.
- Lax–Wendroff and MacCormack are selectable solver models behind the
  same public model API and shared GUI, with common full-slice validation
  coverage wherever the setup is solver-independent.
- Codebase structured so MoC boundaries and species transport can be
  added later without rewriting the interior solver (see
  `ARCHITECTURE.md`).
