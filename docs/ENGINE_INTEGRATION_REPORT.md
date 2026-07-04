# Engine Integration Findings

The cylinder temperature is not the primary cause of the GUI frame-rate issue.
The same runaway temperatures occur in non-GUI fixed-speed engine tests, and
the GUI only displays telemetry. Low frame rate follows from the simulation
running on the UI thread: each GUI frame catches up wall-clock time with many
engine steps, and unstable pipe coupling forces many 1D substeps plus repeated
positivity recovery.

## Debugging And Testing Needs

- Add external-boundary stress tests that impose engine-like pulsed mass/energy
  flow, especially exhaust blowdown from a hot/high-pressure 0D chamber.
- Add long-run coupled tests that assert bounded pressure, temperature, mass,
  and energy over many engine cycles.
- Treat repeated positivity clipping or fallback faces as test failures above a
  small threshold; current recovery can hide runaway until the host model fails.
- External `Flow` boundaries should never panic on an end state the duct solver
  just produced. Clamped/fallible primitive conversion is needed on all
  robustness paths.
- Provide per-step diagnostics for clipped cells, fallback faces, boundary mass
  limits, and energy imbalance, with pipe/end identifiers.

## Species And Lambda

The engine needs species carry-through for oxygen, fuel vapor, inert gas, and
combustion products. Without advected species, exhaust lambda and residual gas
composition are unavailable or wrong. Add passive scalar transport per duct
cell, expose end-cell species fractions in `ExternalPort`, and include species
fluxes in boundary diagnostics. The host engine needs enough exhaust-port
composition to reconstruct lambda from oxygen, fuel, and products, not only
bulk air properties.

## API Recommendations

The current `ExternalBoundaryControl::Flow` is too blunt for valve coupling:
the host computes one fixed flow from stale pipe/chamber states, then the pipe
solver applies it through its substeps. Prefer one of:

- a callback/substep API where the host recomputes boundary flow from updated
  pipe-end and 0D states each internal substep;
- a first-class 0D reservoir/valve boundary that takes chamber pressure,
  temperature, volume, species, valve area, and discharge coefficient;
- a bounded external-flow control with explicit per-step mass/energy limits and
  returned accepted flow.

For violent exhaust blowdown, add validation cases with hot chamber-to-pipe
discharge, choked flow, flow reversal, and acoustic reflection. The API should
return accepted mass, energy, and species transfer so the engine can conserve
its 0D chamber exactly.
