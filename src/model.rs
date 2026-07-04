use std::collections::BTreeMap;

use crate::{
    boundaries::{
        BoundaryCondition, ClosedEnd, DuctEnd, JunctionPort, JunctionSolution, MultiPipeJunction,
        OpenEnd, OrificeFlow, ValveOrifice,
    },
    duct::{BoundaryOverride, Duct, DuctConfig, StepReport},
    gas_properties::GasProperties,
    state::State,
};

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
        }
    }
}

#[derive(Clone, Debug)]
pub struct Model<G>
where
    G: GasProperties,
{
    ducts: Vec<Duct<G, ModelBoundary, ModelBoundary>>,
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
            cfl,
            time: 0.0,
        }
    }

    pub fn add_duct(&mut self, duct: Duct<G, ModelBoundary, ModelBoundary>) {
        self.ducts.push(duct);
    }

    pub fn add_uniform_duct(
        &mut self,
        gas: G,
        config: DuctConfig,
        initial_state: State,
        left_boundary: ModelBoundary,
        right_boundary: ModelBoundary,
    ) {
        self.add_duct(Duct::new(
            gas,
            config,
            initial_state,
            left_boundary,
            right_boundary,
        ));
    }

    pub fn ducts(&self) -> &[Duct<G, ModelBoundary, ModelBoundary>] {
        &self.ducts
    }

    pub fn ducts_mut(&mut self) -> &mut [Duct<G, ModelBoundary, ModelBoundary>] {
        &mut self.ducts
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

    pub fn step(&mut self) -> StepReport {
        let dt = self.global_timestep();
        self.step_with_dt(dt)
    }

    pub fn step_with_dt(&mut self, dt: f64) -> StepReport {
        let boundary_overrides = self.boundary_overrides();
        let mut total = StepReport::default();
        for (duct, (left_override, right_override)) in
            self.ducts.iter_mut().zip(boundary_overrides.into_iter())
        {
            let report = duct.step_with_boundary_controls(dt, left_override, right_override);
            total.clipped_cells += report.clipped_cells;
            total.fallback_faces += report.fallback_faces;
        }
        self.time += dt;
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

    fn boundary_overrides(&self) -> Vec<(BoundaryOverride, BoundaryOverride)> {
        let mut overrides =
            vec![(BoundaryOverride::default(), BoundaryOverride::default()); self.ducts.len()];
        for solved in self.solve_junctions() {
            for ((duct_index, end), boundary_state) in solved
                .connections
                .into_iter()
                .zip(solved.solution.boundary_states)
            {
                match end {
                    DuctEnd::Left => {
                        assert!(overrides[duct_index].0.ghost_state.is_none());
                        overrides[duct_index].0 = BoundaryOverride::ghost(boundary_state);
                    }
                    DuctEnd::Right => {
                        assert!(overrides[duct_index].1.ghost_state.is_none());
                        overrides[duct_index].1 = BoundaryOverride::ghost(boundary_state);
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
            let second_flux = boundary_flux_from_outflow(
                second_state,
                second.1,
                second_area,
                -solved.flow.mass_flow,
                -solved.flow.energy_flow,
                gas,
            );

            set_flux_override(&mut overrides[first.0], first.1, first_state, first_flux);
            set_flux_override(
                &mut overrides[second.0],
                second.1,
                second_state,
                second_flux,
            );
        }
        overrides
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
}

fn set_flux_override(
    overrides: &mut (BoundaryOverride, BoundaryOverride),
    end: DuctEnd,
    ghost_state: State,
    face_flux: State,
) {
    match end {
        DuctEnd::Left => {
            assert!(overrides.0.face_flux.is_none());
            overrides.0 = BoundaryOverride::flux(ghost_state, face_flux);
        }
        DuctEnd::Right => {
            assert!(overrides.1.face_flux.is_none());
            overrides.1 = BoundaryOverride::flux(ghost_state, face_flux);
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
    let prim = port_state.primitive(gas);
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

#[cfg(test)]
mod tests {
    use super::{Model, ModelBoundary};
    use crate::{GasProperties, ValveOrifice};
    use crate::{duct::DuctConfig, gas_properties::TemperatureDependentAir, state::State};

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
}
