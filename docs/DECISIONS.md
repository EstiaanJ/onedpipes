# Decisions



## Interior scheme: two-step Lax–Wendroff (Richtmyer)

**Decision**: use two-step LW for all duct interiors


## Boundary conditions: simple models now, MoC later

**Decision**: v1 uses simple, non-characteristic boundary models (mirror
for closed end, fixed-pressure/extrapolation for open end,
constant-pressure for junctions, quasi-steady orifice for valves) instead
of characteristic-based (MoC) boundaries.


**Upgrade path**: boundary conditions are implemented as pluggable
objects (see ARCHITECTURE.md) specifically so that swapping in
MoC-based boundaries later is a localized change — replace the object
behind each duct-end interface, leave the interior LW solver and the
rest of the model untouched. Do this per-boundary-type as needed (the
valve/orifice boundary is the highest-value target for an eventual MoC
or characteristic-based upgrade, since its accuracy dominates any
downstream use of this boundary).

## Junction boundary: constant pressure with linear acoustic port update

**Decision**: the first junction implementation solves one instantaneous
shared pressure from all connected port states using a linear acoustic
relation at each duct end, then reports per-port mass and energy fluxes.

**Why**: this matches the v1 constant-pressure junction model without
introducing Method of Characteristics machinery or pressure-loss
coefficients before the basic 2-3 pipe conservation milestone is proven.

**Tradeoff accepted**: the linear port update enforces mass balance
directly. Energy conservation is exact for hand-balanced cases with equal
total enthalpy across the ports; more general mixing and loss behavior is
left for a later junction upgrade after the model-level coupling is in
place and validated.

## Artificial dissipation: simplest possible (basic 2nd-difference / Lapidus-style)

**Decision**: a single global scalar artificial-viscosity coefficient
applied as a second-difference smoothing term, not a flux limiter, TVD
scheme, or shock detector.

**Why**: explicitly requested as "the simplest way possible." 

**Tradeoff accepted**:
too high damps genuine tuning/resonance amplitude, too low leaves
ringing. It must be tuned jointly with wall heat transfer (see below),
which provides physical damping for the same numerical symptom.

## Wall heat transfer as a damping mechanism

**Decision**: implement a simple wall heat-loss correlation in ducts, in
part specifically to help damp numerical ringing, not purely for thermal
accuracy.

**Why**: explicitly requested



## Gas properties: temperature-dependent γ(T), cp(T), single species

**Decision**: model γ and cp as functions of temperature for a single
effective gas, but do not implement multi-species transport yet.

**Why**: explicitly flagged as important for accuracy


**Upgrade path**: gas properties are behind a swappable interface
(`get_gamma(T)`, `get_cp(T)`, `get_R()`) specifically so species mass
fractions can later become additional advected scalars that feed these
functions, without changing the interior LW solver. When implemented,
each species is a passively-advected scalar riding in the same
predictor/corrector update as ρ, ρu, ρE.

## Time stepping: global explicit Δt, CFL with 10% margin

**Decision**: one shared `Δt` for the whole model each step, computed as
`0.9 × CFL-limit`, no local time stepping.

**Why**: this is a time-accurate unsteady problem (wave timing matters
for resonance prediction) — local time stepping is a steady-state
convergence trick that would corrupt the physical timing of reflected
waves. The 10% margin (0.9 factor) gives stability headroom against the
nonlinear CFL limit (LW's linear stability limit is CFL = 1) without
being so conservative it needlessly slows the run.


## Positivity/robustness: floors + predictor sanity check + local fallback

**Decision**: clip negative density/pressure post-corrector (with
logging), and check predictor states before the corrector consumes them,
falling back to first-order flux locally if a predictor state is
unphysical.

**Why**: strong blowdown events can push LW predictor states negative
even when the true solution is smooth; a hard crash on this is worse
than a logged, localized correction for a first working model.

**Tradeoff accepted**: repeated fallback/clipping in the same region is a
sign of a tuning problem (too little dissipation, boundary error, Δx too
coarse) — it must not be treated as an acceptable steady-state; the
build-order milestones in PRODUCT_SPEC exist partly to catch this early,
in isolated test cases, rather than in the full duct/junction/orifice
model.
