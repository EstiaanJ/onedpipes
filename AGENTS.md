# Agent Operating Instructions

Read `docs/PRODUCT_SPEC.md` and `docs/ARCHITECTURE.md` before writing code.
`docs/DECISIONS.md` explains *why* choices were made — check it before
"improving" something that looks suboptimal; it may be a deliberate
simplicity/accuracy tradeoff with a documented upgrade path.

## Priorities, in order

1. **Get a fully working end-to-end model first.** Simplicity beats
   accuracy at this stage. Do not add sophistication (limiters, MoC,
   multi-species, implicit solvers) unless it's already in the spec.
2. **Numerical stability and positivity** over accuracy. Never let density,
   pressure, or temperature go non-physical silently — clip and log, per
   `ARCHITECTURE.md` § Positivity.
3. **Testing is mandatory, not incidental.** No boundary type, gas
   model, or numerical feature is done until it has a unit test and,
   where it corresponds to a milestone, a validation case in `tests/`
   checked against an independent reference (analytic, hand-computed, or
   a published solution) — see `PRODUCT_SPEC.md` § Testing strategy.
4. **Validate incrementally.** Do not move to the next milestone until
   the current one passes its validation case (see Milestones in
   `PRODUCT_SPEC.md`).

## Build order

1. Single duct, two-step Lax–Wendroff, closed ends both sides →
   validate: standing wave / organ-pipe resonant frequency matches analytic
   value within a few percent.
2. Add open-end and closed-end boundary models independently → validate
   with known reflection behavior.
3. Add junction boundary (2–3 pipes) → validate mass/energy conservation
   across the junction.
4. Add valve/port orifice boundary → validate steady-flow discharge against
   a hand-computed compressible orifice flow.

Do not proceed to the next step until the current one has a passing
validation case checked into `tests/`.

## Coding conventions

- Conservative variables `U = [ρ, ρu, ρE]` (times duct area, see
  `ARCHITECTURE.md`) are the source of truth in every duct cell. Primitive
  variables (`p, u, T, a`) are derived on demand, never stored as the
  primary state.
- One global `Δt` per solver step across the entire model — no local time
  stepping. See `ARCHITECTURE.md` § Time Stepping.
- Boundary conditions are pluggable objects/functions attached to duct
  ends, not special-cased branches inside the interior solver loop. This
  is required to swap in MoC boundaries later without touching the
  interior scheme.
- Every new boundary type or gas model must ship with a unit test and a
  one-paragraph note in `docs/DECISIONS.md` if it introduces a
  simplification.
- Do not introduce a second numerical scheme (e.g. MacCormack, upwind,
  MoC) as an alternative interior solver without discussing it as a
  decision first — see `DECISIONS.md`.

## When something looks wrong

- Oscillations near a valve event or blowdown → check artificial
  dissipation coefficient and wall heat transfer before assuming the
  interior scheme is broken.
- Drifting mean pressure over long runs → check boundary mass conservation
  first, not the interior flux.
- NaN/negative density or pressure → check the positivity clipping is
  actually being hit before the corrector step consumes bad predictor
  data.

## Out of scope — do not implement unless asked

- Shock-capturing (Riemann solvers, TVD/WENO limiters)
- Implicit time integration for the interior wave solver
- Multi-dimensional (2D/3D) effects
- 0D engine simulation (cylinders, plenums, combustion, valve timing,
  VE/scavenging metrics) — this project is a 1D duct gas-dynamics
  library only, at this stage
