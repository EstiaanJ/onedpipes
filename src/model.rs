use std::collections::BTreeMap;

use crate::{
    boundaries::{
        BoundaryCondition, ClosedEnd, DuctEnd, JunctionPort, JunctionSolution, MultiPipeJunction,
        OpenEnd, OrificeFlow, ValveOrifice,
    },
    duct::{
        BoundaryOverride, Duct, DuctConfig, ExternalBoundaryStepDiagnostic, PipeStepDiagnostic,
        StepReport,
    },
    gas_properties::GasProperties,
    species::{SpeciesFractions, SpeciesMass},
    state::State,
};

/// Stable identifier for a pipe stored in a [`Model`].
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PipeId(pub usize);

/// Stable identifier for a model boundary controlled by an external 0D component.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ExternalBoundaryId(pub usize);

/// A specific end of a pipe.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PipeEnd {
    pub pipe_id: PipeId,
    pub end: DuctEnd,
}

#[derive(Clone, Debug)]
pub enum ModelBoundary {
    Closed,
    Open {
        ambient_pressure: f64,
    },
    Junction {
        junction_id: usize,
    },
    Orifice {
        orifice_id: usize,
        valve: ValveOrifice,
    },
    External {
        external_id: ExternalBoundaryId,
    },
}

/// Snapshot of a pipe end exposed to an external 0D component.
#[derive(Clone, Copy, Debug)]
pub struct ExternalPort {
    pub external_id: ExternalBoundaryId,
    pub pipe_id: PipeId,
    pub end: DuctEnd,
    pub area: f64,
    pub state: State,
    pub species: SpeciesFractions,
}

/// Boundary input supplied by an external 0D component.
#[derive(Clone, Copy, Debug)]
pub enum ExternalBoundaryControl {
    /// Supply a ghost-cell state directly.
    GhostState(State),
    /// Supply an integrated boundary flow.
    ///
    /// Positive `mass_flow_out` and `energy_flow_out` mean flow leaving
    /// the 1D pipe and entering the external 0D component.
    Flow {
        mass_flow_out: f64,
        energy_flow_out: f64,
    },
    /// Supply a boundary flow with per-step transfer limits.
    ///
    /// Limits are integrated amounts over the next model step, not rates.
    /// `inflow_species` is used when accepted flow enters the 1D pipe.
    BoundedFlow {
        mass_flow_out: f64,
        energy_flow_out: f64,
        max_mass_transfer: f64,
        max_energy_transfer: f64,
        inflow_species: SpeciesFractions,
    },
}

#[derive(Clone, Debug)]
pub struct JunctionDiagnostic {
    pub junction_id: usize,
    pub solution: JunctionSolution,
}

#[derive(Clone, Debug)]
pub struct OrificeDiagnostic {
    pub orifice_id: usize,
    pub flow: OrificeFlow,
}

#[derive(Clone, Debug)]
struct SolvedJunction {
    junction_id: usize,
    connections: Vec<(usize, DuctEnd)>,
    solution: JunctionSolution,
}

#[derive(Clone, Debug)]
struct SolvedOrifice {
    orifice_id: usize,
    connections: [(usize, DuctEnd); 2],
    flow: OrificeFlow,
}

impl ModelBoundary {
    pub fn open(ambient_pressure: f64) -> Self {
        Self::Open { ambient_pressure }
    }

    pub fn junction(junction_id: usize) -> Self {
        Self::Junction { junction_id }
    }

    pub fn orifice(orifice_id: usize, valve: ValveOrifice) -> Self {
        Self::Orifice { orifice_id, valve }
    }

    pub fn external(external_id: usize) -> Self {
        Self::External {
            external_id: ExternalBoundaryId(external_id),
        }
    }
}

impl<G: GasProperties> BoundaryCondition<G> for ModelBoundary {
    fn ghost_state(&self, interior: State, end: DuctEnd, gas: G) -> State {
        match *self {
            Self::Closed => ClosedEnd.ghost_state(interior, end, gas),
            Self::Open { ambient_pressure } => {
                OpenEnd::new(ambient_pressure).ghost_state(interior, end, gas)
            }
            Self::Junction { .. } => {
                panic!("junction boundaries require Model::step_with_dt coupling")
            }
            Self::Orifice { .. } => {
                panic!("orifice boundaries require Model::step_with_dt coupling")
            }
            Self::External { .. } => {
                panic!("external boundaries require Model::step_with_dt coupling")
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct Model<G>
where
    G: GasProperties,
{
    ducts: Vec<Duct<G, ModelBoundary, ModelBoundary>>,
    external_controls: BTreeMap<ExternalBoundaryId, ExternalBoundaryControl>,
    cfl: f64,
    time: f64,
}

impl<G> Model<G>
where
    G: GasProperties,
{
    pub fn new(cfl: f64) -> Self {
        assert!(cfl > 0.0);
        Self {
            ducts: Vec::new(),
            external_controls: BTreeMap::new(),
            cfl,
            time: 0.0,
        }
    }

    pub fn add_duct(&mut self, duct: Duct<G, ModelBoundary, ModelBoundary>) -> PipeId {
        let pipe_id = PipeId(self.ducts.len());
        self.ducts.push(duct);
        pipe_id
    }

    pub fn add_uniform_duct(
        &mut self,
        gas: G,
        config: DuctConfig,
        initial_state: State,
        left_boundary: ModelBoundary,
        right_boundary: ModelBoundary,
    ) -> PipeId {
        self.add_duct(Duct::new(
            gas,
            config,
            initial_state,
            left_boundary,
            right_boundary,
        ))
    }

    pub fn add_uniform_duct_with_species(
        &mut self,
        gas: G,
        config: DuctConfig,
        initial_state: State,
        initial_species: SpeciesFractions,
        left_boundary: ModelBoundary,
        right_boundary: ModelBoundary,
    ) -> PipeId {
        self.add_duct(Duct::new_with_species(
            gas,
            config,
            initial_state,
            initial_species,
            left_boundary,
            right_boundary,
        ))
    }

    pub fn ducts(&self) -> &[Duct<G, ModelBoundary, ModelBoundary>] {
        &self.ducts
    }

    pub fn ducts_mut(&mut self) -> &mut [Duct<G, ModelBoundary, ModelBoundary>] {
        &mut self.ducts
    }

    pub fn pipe(&self, pipe_id: PipeId) -> &Duct<G, ModelBoundary, ModelBoundary> {
        &self.ducts[pipe_id.0]
    }

    pub fn pipe_mut(&mut self, pipe_id: PipeId) -> &mut Duct<G, ModelBoundary, ModelBoundary> {
        &mut self.ducts[pipe_id.0]
    }

    pub fn pipe_cells(&self, pipe_id: PipeId) -> &[State] {
        self.pipe(pipe_id).cells()
    }

    pub fn pipe_species_cells(&self, pipe_id: PipeId) -> &[SpeciesFractions] {
        self.pipe(pipe_id).species_cells()
    }

    pub fn pipe_primitive_cells(&self, pipe_id: PipeId) -> Vec<crate::Primitive> {
        self.pipe(pipe_id).primitive_cells()
    }

    pub fn pipe_end_state(&self, pipe_end: PipeEnd) -> State {
        self.pipe(pipe_end.pipe_id).end_state(pipe_end.end)
    }

    pub fn pipe_end_species(&self, pipe_end: PipeEnd) -> SpeciesFractions {
        self.pipe(pipe_end.pipe_id).end_species(pipe_end.end)
    }

    pub fn pipe_total_mass(&self, pipe_id: PipeId) -> f64 {
        self.pipe(pipe_id).total_mass()
    }

    pub fn pipe_total_energy(&self, pipe_id: PipeId) -> f64 {
        self.pipe(pipe_id).total_energy()
    }

    pub fn time(&self) -> f64 {
        self.time
    }

    pub fn junction_diagnostics(&self) -> Vec<JunctionDiagnostic> {
        self.solve_junctions()
            .into_iter()
            .map(|solved| JunctionDiagnostic {
                junction_id: solved.junction_id,
                solution: solved.solution,
            })
            .collect()
    }

    pub fn orifice_diagnostics(&self) -> Vec<OrificeDiagnostic> {
        self.solve_orifices()
            .into_iter()
            .map(|solved| OrificeDiagnostic {
                orifice_id: solved.orifice_id,
                flow: solved.flow,
            })
            .collect()
    }

    pub fn external_ports(&self) -> Vec<ExternalPort> {
        let mut ports = Vec::new();
        for (duct_index, duct) in self.ducts.iter().enumerate() {
            if let ModelBoundary::External { external_id } = duct.left_boundary() {
                ports.push(ExternalPort {
                    external_id: *external_id,
                    pipe_id: PipeId(duct_index),
                    end: DuctEnd::Left,
                    area: duct.config().area,
                    state: duct.end_state(DuctEnd::Left),
                    species: duct.end_species(DuctEnd::Left),
                });
            }
            if let ModelBoundary::External { external_id } = duct.right_boundary() {
                ports.push(ExternalPort {
                    external_id: *external_id,
                    pipe_id: PipeId(duct_index),
                    end: DuctEnd::Right,
                    area: duct.config().area,
                    state: duct.end_state(DuctEnd::Right),
                    species: duct.end_species(DuctEnd::Right),
                });
            }
        }
        ports
    }

    pub fn set_external_boundary_control(
        &mut self,
        external_id: ExternalBoundaryId,
        control: ExternalBoundaryControl,
    ) {
        self.external_controls.insert(external_id, control);
    }

    pub fn clear_external_boundary_controls(&mut self) {
        self.external_controls.clear();
    }

    pub fn step(&mut self) -> StepReport {
        let dt = self.global_timestep();
        self.step_with_dt(dt)
    }

    pub fn step_with_dt(&mut self, dt: f64) -> StepReport {
        let (boundary_overrides, mut total) = self.boundary_overrides(dt);
        for (duct_index, (duct, (left_override, right_override))) in self
            .ducts
            .iter_mut()
            .zip(boundary_overrides.into_iter())
            .enumerate()
        {
            let mut report = duct.step_with_boundary_controls(dt, left_override, right_override);
            if report.clipped_cells > 0 || report.fallback_faces > 0 {
                report.pipe_diagnostics.push(PipeStepDiagnostic {
                    pipe_index: duct_index,
                    clipped_cell_indices: report.clipped_cell_indices.clone(),
                    fallback_face_indices: report.fallback_face_indices.clone(),
                });
            }
            total.absorb(report);
        }
        self.time += dt;
        total
    }

    pub fn step_with_dt_and_external_callback<F>(
        &mut self,
        dt: f64,
        substeps: usize,
        mut callback: F,
    ) -> StepReport
    where
        F: FnMut(ExternalPort) -> ExternalBoundaryControl,
    {
        assert!(substeps > 0);
        let sub_dt = dt / substeps as f64;
        let mut total = StepReport::default();
        for _ in 0..substeps {
            self.clear_external_boundary_controls();
            for port in self.external_ports() {
                let control = callback(port);
                self.set_external_boundary_control(port.external_id, control);
            }
            total.absorb(self.step_with_dt(sub_dt));
        }
        total
    }

    pub fn run_until(&mut self, end_time: f64) -> StepReport {
        let mut total = StepReport::default();
        while self.time < end_time {
            let mut dt = self.global_timestep();
            if self.time + dt > end_time {
                dt = end_time - self.time;
            }
            let report = self.step_with_dt(dt);
            total.clipped_cells += report.clipped_cells;
            total.fallback_faces += report.fallback_faces;
        }
        total
    }

    fn global_timestep(&self) -> f64 {
        assert!(!self.ducts.is_empty());
        let min_dt = self
            .ducts
            .iter()
            .map(|duct| duct.config().dx() / duct.max_signal_speed())
            .fold(f64::INFINITY, f64::min);
        0.9 * self.cfl * min_dt
    }

    fn boundary_overrides(
        &self,
        dt: f64,
    ) -> (Vec<(BoundaryOverride, BoundaryOverride)>, StepReport) {
        let mut overrides =
            vec![(BoundaryOverride::default(), BoundaryOverride::default()); self.ducts.len()];
        let mut report = StepReport::default();
        for solved in self.solve_junctions() {
            for ((duct_index, end), boundary_state) in solved
                .connections
                .into_iter()
                .zip(solved.solution.boundary_states)
            {
                match end {
                    DuctEnd::Left => {
                        assert!(overrides[duct_index].0.ghost_state.is_none());
                        overrides[duct_index].0 = BoundaryOverride::ghost_with_species(
                            boundary_state,
                            self.ducts[duct_index].end_species(end),
                        );
                    }
                    DuctEnd::Right => {
                        assert!(overrides[duct_index].1.ghost_state.is_none());
                        overrides[duct_index].1 = BoundaryOverride::ghost_with_species(
                            boundary_state,
                            self.ducts[duct_index].end_species(end),
                        );
                    }
                }
            }
        }
        for solved in self.solve_orifices() {
            let [first, second] = solved.connections;
            let first_state = self.ducts[first.0].end_state(first.1);
            let second_state = self.ducts[second.0].end_state(second.1);
            let first_area = self.ducts[first.0].config().area;
            let second_area = self.ducts[second.0].config().area;
            let gas = self.ducts[first.0].gas();
            let first_flux = boundary_flux_from_outflow(
                first_state,
                first.1,
                first_area,
                solved.flow.mass_flow,
                solved.flow.energy_flow,
                gas,
            );
            let first_species_flux = species_flux_from_outflow(
                first.1,
                first_area,
                solved.flow.mass_flow,
                self.ducts[first.0].end_species(first.1),
                self.ducts[second.0].end_species(second.1),
            );
            let second_flux = boundary_flux_from_outflow(
                second_state,
                second.1,
                second_area,
                -solved.flow.mass_flow,
                -solved.flow.energy_flow,
                gas,
            );
            let second_species_flux = species_flux_from_outflow(
                second.1,
                second_area,
                -solved.flow.mass_flow,
                self.ducts[second.0].end_species(second.1),
                self.ducts[first.0].end_species(first.1),
            );

            set_flux_override(
                &mut overrides[first.0],
                first.1,
                first_state,
                self.ducts[first.0].end_species(first.1),
                first_flux,
                first_species_flux,
            );
            set_flux_override(
                &mut overrides[second.0],
                second.1,
                second_state,
                self.ducts[second.0].end_species(second.1),
                second_flux,
                second_species_flux,
            );
        }
        for (duct_index, end, external_id) in self.external_connections() {
            let duct = &self.ducts[duct_index];
            let state = duct.end_state(end);
            let species = duct.end_species(end);
            let control = self.external_controls.get(&external_id).unwrap_or_else(|| {
                panic!(
                    "external boundary {:?} requires control before Model::step_with_dt",
                    external_id
                )
            });
            let boundary_override = match *control {
                ExternalBoundaryControl::GhostState(ghost_state) => {
                    BoundaryOverride::ghost_with_species(ghost_state, species)
                }
                ExternalBoundaryControl::Flow {
                    mass_flow_out,
                    energy_flow_out,
                } => {
                    let face_flux = boundary_flux_from_outflow(
                        state,
                        end,
                        duct.config().area,
                        mass_flow_out,
                        energy_flow_out,
                        duct.gas(),
                    );
                    let species_flux = species_flux_from_outflow(
                        end,
                        duct.config().area,
                        mass_flow_out,
                        species,
                        SpeciesFractions::AIR,
                    );
                    report
                        .external_boundary_diagnostics
                        .push(external_diagnostic(
                            external_id,
                            PipeId(duct_index),
                            end,
                            mass_flow_out,
                            mass_flow_out,
                            energy_flow_out,
                            energy_flow_out,
                            dt,
                            species_transfer_from_outflow(
                                mass_flow_out,
                                species,
                                SpeciesFractions::AIR,
                                dt,
                            ),
                            false,
                        ));
                    BoundaryOverride::flux_with_species(state, species, face_flux, species_flux)
                }
                ExternalBoundaryControl::BoundedFlow {
                    mass_flow_out,
                    energy_flow_out,
                    max_mass_transfer,
                    max_energy_transfer,
                    inflow_species,
                } => {
                    let (accepted_mass_flow_out, mass_limited) =
                        limit_flow_rate(mass_flow_out, max_mass_transfer, dt);
                    let (accepted_energy_flow_out, energy_limited) =
                        limit_flow_rate(energy_flow_out, max_energy_transfer, dt);
                    let face_flux = boundary_flux_from_outflow(
                        state,
                        end,
                        duct.config().area,
                        accepted_mass_flow_out,
                        accepted_energy_flow_out,
                        duct.gas(),
                    );
                    let species_flux = species_flux_from_outflow(
                        end,
                        duct.config().area,
                        accepted_mass_flow_out,
                        species,
                        inflow_species,
                    );
                    report
                        .external_boundary_diagnostics
                        .push(external_diagnostic(
                            external_id,
                            PipeId(duct_index),
                            end,
                            mass_flow_out,
                            accepted_mass_flow_out,
                            energy_flow_out,
                            accepted_energy_flow_out,
                            dt,
                            species_transfer_from_outflow(
                                accepted_mass_flow_out,
                                species,
                                inflow_species,
                                dt,
                            ),
                            mass_limited || energy_limited,
                        ));
                    BoundaryOverride::flux_with_species(state, species, face_flux, species_flux)
                }
            };
            match end {
                DuctEnd::Left => {
                    assert!(overrides[duct_index].0.ghost_state.is_none());
                    overrides[duct_index].0 = boundary_override;
                }
                DuctEnd::Right => {
                    assert!(overrides[duct_index].1.ghost_state.is_none());
                    overrides[duct_index].1 = boundary_override;
                }
            }
        }
        (overrides, report)
    }

    fn solve_junctions(&self) -> Vec<SolvedJunction> {
        let mut groups: BTreeMap<usize, Vec<(usize, DuctEnd, JunctionPort)>> = BTreeMap::new();
        for (duct_index, duct) in self.ducts.iter().enumerate() {
            if let ModelBoundary::Junction { junction_id } = duct.left_boundary() {
                groups.entry(*junction_id).or_default().push((
                    duct_index,
                    DuctEnd::Left,
                    self.port_for(duct_index, DuctEnd::Left),
                ));
            }
            if let ModelBoundary::Junction { junction_id } = duct.right_boundary() {
                groups.entry(*junction_id).or_default().push((
                    duct_index,
                    DuctEnd::Right,
                    self.port_for(duct_index, DuctEnd::Right),
                ));
            }
        }

        let mut solved_junctions = Vec::with_capacity(groups.len());
        for (junction_id, ports) in groups {
            assert!(
                ports.len() >= 2,
                "junction {junction_id} must connect at least two duct ends"
            );
            let port_states: Vec<_> = ports.iter().map(|(_, _, port)| *port).collect();
            let gas = self.ducts[ports[0].0].gas();
            let solution = MultiPipeJunction.solve(&port_states, gas);
            assert!(
                solution.mass_residual().abs() < 1.0e-8,
                "junction {junction_id} mass residual = {}",
                solution.mass_residual()
            );
            assert!(
                solution.energy_residual().abs() < 1.0e-3,
                "junction {junction_id} energy residual = {}",
                solution.energy_residual()
            );

            solved_junctions.push(SolvedJunction {
                junction_id,
                connections: ports
                    .into_iter()
                    .map(|(duct_index, end, _)| (duct_index, end))
                    .collect(),
                solution,
            });
        }
        solved_junctions
    }

    fn solve_orifices(&self) -> Vec<SolvedOrifice> {
        let mut groups: BTreeMap<usize, Vec<(usize, DuctEnd, ValveOrifice)>> = BTreeMap::new();
        for (duct_index, duct) in self.ducts.iter().enumerate() {
            if let ModelBoundary::Orifice { orifice_id, valve } = duct.left_boundary() {
                groups
                    .entry(*orifice_id)
                    .or_default()
                    .push((duct_index, DuctEnd::Left, *valve));
            }
            if let ModelBoundary::Orifice { orifice_id, valve } = duct.right_boundary() {
                groups
                    .entry(*orifice_id)
                    .or_default()
                    .push((duct_index, DuctEnd::Right, *valve));
            }
        }

        let mut solved_orifices = Vec::with_capacity(groups.len());
        for (orifice_id, ports) in groups {
            assert_eq!(
                ports.len(),
                2,
                "orifice {orifice_id} must connect exactly two duct ends"
            );
            let first = ports[0];
            let second = ports[1];
            assert_eq!(
                first.2, second.2,
                "orifice {orifice_id} must use the same ValveOrifice at both ends"
            );
            let gas = self.ducts[first.0].gas();
            let first_state = self.ducts[first.0].end_state(first.1);
            let second_state = self.ducts[second.0].end_state(second.1);
            let flow = first.2.mass_flow(first_state, second_state, gas);
            solved_orifices.push(SolvedOrifice {
                orifice_id,
                connections: [(first.0, first.1), (second.0, second.1)],
                flow,
            });
        }
        solved_orifices
    }

    fn port_for(&self, duct_index: usize, end: DuctEnd) -> JunctionPort {
        let duct = &self.ducts[duct_index];
        JunctionPort::new(duct.end_state(end), end, duct.config().area)
    }

    fn external_connections(&self) -> Vec<(usize, DuctEnd, ExternalBoundaryId)> {
        let mut connections = Vec::new();
        for (duct_index, duct) in self.ducts.iter().enumerate() {
            if let ModelBoundary::External { external_id } = duct.left_boundary() {
                connections.push((duct_index, DuctEnd::Left, *external_id));
            }
            if let ModelBoundary::External { external_id } = duct.right_boundary() {
                connections.push((duct_index, DuctEnd::Right, *external_id));
            }
        }
        connections
    }
}

fn set_flux_override(
    overrides: &mut (BoundaryOverride, BoundaryOverride),
    end: DuctEnd,
    ghost_state: State,
    ghost_species: SpeciesFractions,
    face_flux: State,
    face_species_flux: SpeciesMass,
) {
    match end {
        DuctEnd::Left => {
            assert!(overrides.0.face_flux.is_none());
            overrides.0 = BoundaryOverride::flux_with_species(
                ghost_state,
                ghost_species,
                face_flux,
                face_species_flux,
            );
        }
        DuctEnd::Right => {
            assert!(overrides.1.face_flux.is_none());
            overrides.1 = BoundaryOverride::flux_with_species(
                ghost_state,
                ghost_species,
                face_flux,
                face_species_flux,
            );
        }
    }
}

fn boundary_flux_from_outflow<G: GasProperties>(
    port_state: State,
    end: DuctEnd,
    area: f64,
    mass_flow_out: f64,
    energy_flow_out: f64,
    gas: G,
) -> State {
    let prim = port_state.primitive_clamped(gas);
    let coordinate_sign = match end {
        DuctEnd::Left => -1.0,
        DuctEnd::Right => 1.0,
    };
    let mass_flux = coordinate_sign * mass_flow_out / area;
    let velocity = mass_flux / prim.rho;
    State {
        rho: mass_flux,
        momentum: mass_flux * velocity + prim.p,
        rho_total_energy: coordinate_sign * energy_flow_out / area,
    }
}

fn limit_flow_rate(flow_rate: f64, max_transfer: f64, dt: f64) -> (f64, bool) {
    assert!(dt > 0.0);
    let max_rate = max_transfer.abs() / dt;
    if !max_rate.is_finite() {
        return (flow_rate, false);
    }
    let accepted = flow_rate.clamp(-max_rate, max_rate);
    (accepted, accepted != flow_rate)
}

fn species_flux_from_outflow(
    end: DuctEnd,
    area: f64,
    mass_flow_out: f64,
    pipe_species: SpeciesFractions,
    inflow_species: SpeciesFractions,
) -> SpeciesMass {
    let coordinate_sign = match end {
        DuctEnd::Left => -1.0,
        DuctEnd::Right => 1.0,
    };
    let mass_flux = coordinate_sign * mass_flow_out / area;
    if mass_flow_out >= 0.0 {
        pipe_species.scale(mass_flux)
    } else {
        inflow_species.scale(mass_flux)
    }
}

fn species_transfer_from_outflow(
    mass_flow_out: f64,
    pipe_species: SpeciesFractions,
    inflow_species: SpeciesFractions,
    dt: f64,
) -> SpeciesMass {
    if mass_flow_out >= 0.0 {
        pipe_species.scale(mass_flow_out * dt)
    } else {
        inflow_species.scale(mass_flow_out * dt)
    }
}

#[allow(clippy::too_many_arguments)]
fn external_diagnostic(
    external_id: ExternalBoundaryId,
    pipe_id: PipeId,
    end: DuctEnd,
    requested_mass_flow_out: f64,
    accepted_mass_flow_out: f64,
    requested_energy_flow_out: f64,
    accepted_energy_flow_out: f64,
    dt: f64,
    species_transferred_out: SpeciesMass,
    limited: bool,
) -> ExternalBoundaryStepDiagnostic {
    ExternalBoundaryStepDiagnostic {
        external_id: external_id.0,
        pipe_index: pipe_id.0,
        end,
        requested_mass_flow_out,
        accepted_mass_flow_out,
        requested_energy_flow_out,
        accepted_energy_flow_out,
        mass_transferred_out: accepted_mass_flow_out * dt,
        energy_transferred_out: accepted_energy_flow_out * dt,
        species_transferred_out,
        limited,
    }
}

#[cfg(test)]
mod tests {
    use super::{Model, ModelBoundary};
    use crate::{DuctEnd, GasProperties, SpeciesFractions, ValveOrifice};
    use crate::{
        duct::DuctConfig,
        gas_properties::TemperatureDependentAir,
        model::{ExternalBoundaryControl, ExternalBoundaryId, PipeEnd},
        state::State,
    };

    #[test]
    fn model_supports_mixed_boundary_types() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::Closed,
            ModelBoundary::junction(0),
        );
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::junction(0),
            ModelBoundary::open(101_325.0),
        );

        assert_eq!(model.ducts().len(), 2);
    }

    #[test]
    fn model_solves_junction_before_advancing_ducts() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::Closed,
            ModelBoundary::junction(0),
        );
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::junction(0),
            ModelBoundary::open(101_325.0),
        );

        let report = model.step_with_dt(1.0e-6);

        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
        assert!(model.time() > 0.0);
    }

    #[test]
    fn model_reports_junction_conservation_diagnostics() {
        let gas = TemperatureDependentAir::new();
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            State::from_primitive(1.2, 24.0, 140_000.0, gas),
            ModelBoundary::Closed,
            ModelBoundary::junction(0),
        );
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 0.6),
            State::from_primitive(0.9, 0.0, 90_000.0, gas),
            ModelBoundary::junction(0),
            ModelBoundary::open(101_325.0),
        );
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 0.4),
            State::from_primitive(1.4, 0.0, 95_000.0, gas),
            ModelBoundary::junction(0),
            ModelBoundary::open(101_325.0),
        );

        let diagnostics = model.junction_diagnostics();

        assert_eq!(diagnostics.len(), 1);
        assert_eq!(diagnostics[0].junction_id, 0);
        assert!(diagnostics[0].solution.mass_residual().abs() < 1.0e-10);
        assert!(diagnostics[0].solution.energy_residual().abs() < 1.0e-5);
    }

    #[test]
    fn model_reports_orifice_discharge_diagnostics() {
        let gas = TemperatureDependentAir::new();
        let valve = ValveOrifice::new(0.8, 1.0e-4);
        let upstream = State::from_primitive(200_000.0 / (gas.r() * 300.0), 0.0, 200_000.0, gas);
        let downstream = State::from_primitive(140_000.0 / (gas.r() * 300.0), 0.0, 140_000.0, gas);
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            upstream,
            ModelBoundary::Closed,
            ModelBoundary::orifice(0, valve),
        );
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            downstream,
            ModelBoundary::orifice(0, valve),
            ModelBoundary::open(101_325.0),
        );

        let diagnostics = model.orifice_diagnostics();

        assert_eq!(diagnostics.len(), 1);
        assert_eq!(diagnostics[0].orifice_id, 0);
        assert!(diagnostics[0].flow.mass_flow > 0.0);
    }

    #[test]
    fn model_steps_with_orifice_flux_boundary() {
        let gas = TemperatureDependentAir::new();
        let valve = ValveOrifice::new(0.8, 1.0e-4);
        let upstream = State::from_primitive(200_000.0 / (gas.r() * 300.0), 0.0, 200_000.0, gas);
        let downstream = State::from_primitive(140_000.0 / (gas.r() * 300.0), 0.0, 140_000.0, gas);
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            upstream,
            ModelBoundary::Closed,
            ModelBoundary::orifice(0, valve),
        );
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            downstream,
            ModelBoundary::orifice(0, valve),
            ModelBoundary::open(101_325.0),
        );

        let report = model.step_with_dt(1.0e-7);

        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
    }

    #[test]
    fn model_api_exposes_pipe_ids_and_state_queries() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        let pipe_id = model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::Closed,
            ModelBoundary::Open {
                ambient_pressure: 101_325.0,
            },
        );

        assert_eq!(pipe_id.0, 0);
        assert_eq!(model.pipe_cells(pipe_id).len(), 8);
        assert_eq!(
            model.pipe_end_state(PipeEnd {
                pipe_id,
                end: DuctEnd::Left,
            }),
            state
        );
        assert!(model.pipe_total_mass(pipe_id) > 0.0);
        assert!(model.pipe_total_energy(pipe_id) > 0.0);
    }

    #[test]
    fn external_ports_expose_pipe_end_state_for_zero_d_coupling() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        let pipe_id = model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::Closed,
            ModelBoundary::external(7),
        );

        let ports = model.external_ports();

        assert_eq!(ports.len(), 1);
        assert_eq!(ports[0].external_id, ExternalBoundaryId(7));
        assert_eq!(ports[0].pipe_id, pipe_id);
        assert_eq!(ports[0].end, DuctEnd::Right);
        assert_eq!(ports[0].state, state);
    }

    #[test]
    fn model_steps_with_external_ghost_state_control() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::Closed,
            ModelBoundary::external(0),
        );
        model.set_external_boundary_control(
            ExternalBoundaryId(0),
            ExternalBoundaryControl::GhostState(state),
        );

        let report = model.step_with_dt(1.0e-6);

        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
    }

    #[test]
    fn model_steps_with_external_flow_control() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::external(0),
            ModelBoundary::Closed,
        );
        model.set_external_boundary_control(
            ExternalBoundaryId(0),
            ExternalBoundaryControl::Flow {
                mass_flow_out: 0.0,
                energy_flow_out: 0.0,
            },
        );

        let report = model.step_with_dt(1.0e-6);

        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
    }

    #[test]
    fn bounded_external_flow_reports_accepted_transfer_and_species() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        let pipe_id = model.add_uniform_duct_with_species(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            SpeciesFractions::AIR,
            ModelBoundary::external(0),
            ModelBoundary::Closed,
        );
        let exhaust = SpeciesFractions::EXHAUST;
        model.set_external_boundary_control(
            ExternalBoundaryId(0),
            ExternalBoundaryControl::BoundedFlow {
                mass_flow_out: -0.20,
                energy_flow_out: -40_000.0,
                max_mass_transfer: 1.0e-8,
                max_energy_transfer: 2.0e-3,
                inflow_species: exhaust,
            },
        );

        let report = model.step_with_dt(1.0e-6);

        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
        assert_eq!(report.external_boundary_diagnostics.len(), 1);
        let diagnostic = &report.external_boundary_diagnostics[0];
        assert!(diagnostic.limited);
        assert!((diagnostic.accepted_mass_flow_out + 0.01).abs() < 1.0e-12);
        assert!((diagnostic.mass_transferred_out + 1.0e-8).abs() < 1.0e-14);
        assert!(diagnostic.species_transferred_out.products < 0.0);
        assert!(model.pipe_species_cells(pipe_id)[0].products > SpeciesFractions::AIR.products);
    }

    #[test]
    fn external_callback_recomputes_controls_each_substep() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut model = Model::new(0.5);
        model.add_uniform_duct(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ModelBoundary::external(0),
            ModelBoundary::Closed,
        );
        let mut calls = 0;

        let report = model.step_with_dt_and_external_callback(4.0e-7, 4, |port| {
            calls += 1;
            let prim = port.state.primitive(gas);
            ExternalBoundaryControl::Flow {
                mass_flow_out: -1.0e-6 * (prim.p / 101_325.0),
                energy_flow_out: 0.0,
            }
        });

        assert_eq!(calls, 4);
        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
        assert_eq!(report.external_boundary_diagnostics.len(), 4);
    }
}
