# 1D Duct Gas-Dynamics Solver

A quasi-1D unsteady compressible flow solver for ducts (pipes) of varying
cross-sectional area. This is the current, and only, scope of work.

This is an engineering-prediction tool, not an acoustics or
shock-capturing tool. Sound quality, exact shock structure, and
discontinuity preservation are explicitly out of scope (see
`docs/PRODUCT_SPEC.md`).

This library is inteded for engine-simulation

## Current implementation

The checked-in Rust solver currently covers milestone 1 from
`docs/PRODUCT_SPEC.md`:

- Solves the quasi-1D Euler equations in each duct using a two-step
  Lax–Wendroff (Richtmyer) scheme.
- Provides pluggable closed-end and open-end boundary models with pulse
  reflection validation.
- Tracks temperature-dependent gas properties (γ(T), cp(T)) for a single
  effective gas for now, structured so per-species tracking can be added
  later without restructuring the solver.
- Includes a minimal artificial-dissipation term for numerical stability.

Junctions, valve/orifice boundaries, and wall heat transfer remain
planned milestone work.

## What it doesn't do (yet)

- No Method of Characteristics (MoC) boundaries — planned upgrade, see
  `docs/DECISIONS.md`.
- No multi-species transport yet — single effective gas for now.
- No shock-capturing / TVD limiters — not needed for this application
  (see `docs/PRODUCT_SPEC.md` for why).
- No 0D engine coupling (cylinders, plenums, combustion) — out of scope
  for this phase of the project.

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
src/                       Solver source (created during implementation)
tests/                     Validation cases (organ-pipe resonance, Sod tube, etc.)
```

## Getting started

See `docs/ARCHITECTURE.md` for the module layout and `docs/PRODUCT_SPEC.md`
for the milestone plan. The build order is: bare-pipe wave propagation →
closed/open end boundaries → junctions → valve/orifice boundary.
