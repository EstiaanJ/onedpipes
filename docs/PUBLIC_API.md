# Public API Specification

This document describes the public Rust API intended for embedding
`onedpipes` as the 1D duct-flow module inside a larger engine simulator.
The crate does not implement 0D engine components; instead, it exposes
pipe-end state and accepts boundary controls so an external cylinder,
plenum, or controller can couple to the 1D network.

## Conventions

- SI units are assumed throughout.
- `State` stores conservative variables per unit duct area:
  `rho`, `rho*u`, `rho*E`.
- `Primitive` values (`rho`, `u`, `p`, `temperature`, `sound_speed`) are
  derived from `State`.
- One global time step is used per `Model` step.
- `DuctEnd::Left` and `DuctEnd::Right` identify the negative-x and
  positive-x ends of a pipe.
- For `ExternalBoundaryControl::Flow`, positive `mass_flow_out` and
  `energy_flow_out` mean flow leaves the 1D pipe and enters the external
  0D component.

## Typical Engine-Simulator Loop

```rust
use onedpipes::{
    DuctConfig, ExternalBoundaryControl, ExternalBoundaryId, Model, ModelBoundary,
    State, TemperatureDependentAir,
};

let gas = TemperatureDependentAir::new();
let initial = State::from_primitive(1.2, 0.0, 101_325.0, gas);
let mut model = Model::new(0.5);

let pipe = model.add_uniform_duct(
    gas,
    DuctConfig::new(0.6, 80, 3.0e-4),
    initial,
    ModelBoundary::external(0),
    ModelBoundary::open(101_325.0),
);

for port in model.external_ports() {
    let primitive = port.state.primitive(gas);
    let _pressure_seen_by_0d = primitive.p;
    model.set_external_boundary_control(
        port.external_id,
        ExternalBoundaryControl::Flow {
            mass_flow_out: 0.0,
            energy_flow_out: 0.0,
        },
    );
}

let report = model.step();
let cells = model.pipe_primitive_cells(pipe);
```

## Assembly Types

### `Model<G>`

Owns a network of 1D pipes using gas model `G: GasProperties`.

Public functions:

- `Model::new(cfl: f64) -> Model<G>`  
  Creates an empty model. `cfl` must be positive.
- `add_duct(duct) -> PipeId`  
  Adds a preconstructed duct with model-level boundary types.
- `add_uniform_duct(gas, config, initial_state, left_boundary, right_boundary) -> PipeId`  
  Convenience constructor for a pipe initialized to one state.
- `add_uniform_duct_with_species(gas, config, initial_state, initial_species, left_boundary, right_boundary) -> PipeId`  
  Convenience constructor for a pipe initialized to one state and one
  passive species composition.
- `ducts()`, `ducts_mut()`  
  Low-level access to the stored ducts in insertion order.
- `pipe(pipe_id)`, `pipe_mut(pipe_id)`  
  Access one pipe by stable `PipeId`.
- `pipe_cells(pipe_id) -> &[State]`  
  Returns conservative cell states.
- `pipe_species_cells(pipe_id) -> &[SpeciesFractions]`  
  Returns passive species fractions in each cell.
- `pipe_primitive_cells(pipe_id) -> Vec<Primitive>`  
  Returns derived primitive cell states.
- `pipe_end_state(PipeEnd) -> State`  
  Returns the conservative state at a pipe end.
- `pipe_end_species(PipeEnd) -> SpeciesFractions`  
  Returns passive species fractions at a pipe end.
- `pipe_total_mass(pipe_id) -> f64`, `pipe_total_energy(pipe_id) -> f64`  
  Return area-integrated inventories for one pipe.
- `time() -> f64`  
  Returns current model time.
- `step() -> StepReport`  
  Advances by the global CFL time step.
- `step_with_dt(dt) -> StepReport`  
  Advances by an explicit time step.
- `step_with_dt_and_external_callback(dt, substeps, callback) -> StepReport`  
  Advances by smaller pipe substeps. Before each substep, the callback
  receives fresh `ExternalPort` snapshots and returns boundary controls.
- `run_until(end_time) -> StepReport`  
  Advances until the requested model time.
- `junction_diagnostics() -> Vec<JunctionDiagnostic>`  
  Solves current junctions and reports conservation diagnostics without
  advancing the model.
- `orifice_diagnostics() -> Vec<OrificeDiagnostic>`  
  Solves current orifices and reports discharge diagnostics without
  advancing the model.
- `external_ports() -> Vec<ExternalPort>`  
  Returns current pipe-end states exposed to external 0D components.
- `set_external_boundary_control(external_id, control)`  
  Supplies the boundary input required before stepping a model with
  external boundaries.
- `clear_external_boundary_controls()`  
  Clears stored external boundary controls.

### `ModelBoundary`

Boundary specification used when adding pipes to `Model`.

Variants and constructors:

- `ModelBoundary::Closed`
- `ModelBoundary::open(ambient_pressure)`
- `ModelBoundary::junction(junction_id)`
- `ModelBoundary::orifice(orifice_id, ValveOrifice)`
- `ModelBoundary::external(external_id)`

Junctions may connect two or more pipe ends with the same `junction_id`.
Orifices must connect exactly two pipe ends with the same `orifice_id`
and the same `ValveOrifice` parameters. External boundaries are controlled
by the embedding simulator.

### Identifiers and Port Snapshots

- `PipeId(pub usize)`  
  Stable identifier returned when a pipe is added.
- `ExternalBoundaryId(pub usize)`  
  Stable identifier for a 0D coupling boundary.
- `PipeEnd { pipe_id, end }`  
  Identifies one pipe end.
- `ExternalPort { external_id, pipe_id, end, area, state, species }`  
  Snapshot passed to external 0D code, including passive composition for
  lambda/residual diagnostics.
- `ExternalBoundaryControl`  
  `GhostState(State)`, `Flow { mass_flow_out, energy_flow_out }`, or
  `BoundedFlow { mass_flow_out, energy_flow_out, max_mass_transfer,
  max_energy_transfer, inflow_species }`.

## Pipe and State Types

### `DuctConfig`

Pipe discretization and robustness settings.

- `DuctConfig::new(length, cells, area)`  
  Requires positive `length` and `area`, and at least four cells.
- Fields: `length`, `cells`, `area`, `artificial_viscosity`,
  `density_floor`, `pressure_floor`.
- `dx() -> f64` returns cell width.

### `Duct<G, L, R>`

Low-level single-pipe solver. Engine integrations should normally use
`Model`, but `Duct` remains public for isolated validation and advanced
use.

Public functions include constructors (`new`, `from_initializer`), state
accessors (`cells`, `primitive_cells`, `end_state`, `set_cell`),
inventory accessors (`total_mass`, `total_energy`), and stepping methods
(`step`, `step_with_boundary_overrides`,
`step_with_boundary_controls`).

### `State`

Conservative state.

- `State::from_primitive(rho, u, p, gas)` creates a conservative state.
- `try_primitive(gas) -> Result<Primitive, PrimitiveError>` derives
  primitive values without silently accepting nonphysical states.
- `primitive(gas) -> Primitive` derives primitive values and panics if
  the state is nonphysical.
- `primitive_clamped(gas) -> Primitive` is for robustness/fallback paths.
- `flux(gas)` and `flux_clamped(gas)` return Euler fluxes.

### `Primitive`

Derived state: `rho`, `u`, `p`, `temperature`, `sound_speed`.

### `StepReport`

Per-step robustness counters:

- `clipped_cells`
- `fallback_faces`
- `clipped_cell_indices`
- `fallback_face_indices`
- `pipe_diagnostics`
- `external_boundary_diagnostics`

Repeated nonzero values are a sign of a tuning, boundary, or grid issue.
External diagnostics report requested and accepted flow, integrated
mass/energy transfer, species transfer, and whether bounded-flow limits
were active.

### `SpeciesFractions`

Passive composition for oxygen, fuel vapor, inert gas, and combustion
products. `lambda(stoich_fuel_oxygen_ratio)` returns a host-side lambda
estimate when fuel vapor is present. These species are transported for
diagnostics and coupling; they do not yet alter the gas-property model.

## Boundary Models

### Closed and Open Ends

- `ClosedEnd` reflects velocity and preserves pressure/density.
- `OpenEnd::new(ambient_pressure)` fixes ghost-cell static pressure.

### Junction

- `MultiPipeJunction::solve(&[JunctionPort], gas) -> JunctionSolution`
- `JunctionPort::new(state, end, area)`
- `JunctionSolution` contains shared pressure, per-port boundary states,
  per-port flows, `mass_residual()`, and `energy_residual()`.

### Valve/Orifice

- `ValveOrifice::new(discharge_coefficient, flow_area)`
- `mass_flow(upstream, downstream, gas) -> OrificeFlow`
- `OrificeFlow` reports mass flow, energy flow, pressure ratio, critical
  pressure ratio, choking state, and upstream stagnation values.

## Gas Model Interface

`GasProperties` is the swappable gas-property trait:

- `r()`
- `cp(temperature)`
- `cv(temperature)`
- `gamma(temperature)`
- `internal_energy_from_temperature(temperature)`
- `temperature_from_internal_energy(internal_energy)`

`TemperatureDependentAir::new()` is the provided single-effective-gas
model.

## Validation and Viewer Helpers

The `validation` exports (`OrganPipeConfig`, `OrganPipeRun`, `RunReport`,
`ScalarField`, `Snapshot`) support the included organ-pipe viewer and
validation workflow. They are public but not required for engine
integration.
