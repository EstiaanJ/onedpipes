use crate::{
    boundaries::{BoundaryCondition, DuctEnd},
    gas_properties::GasProperties,
    solvers::{self, SolverKind},
    state::{Primitive, State},
};

#[derive(Clone, Debug)]
pub struct DuctConfig {
    pub length: f64,
    pub cells: usize,
    pub area: f64,
    pub cell_areas: Vec<f64>,
    pub face_areas: Vec<f64>,
    pub artificial_viscosity: f64,
    pub density_floor: f64,
    pub pressure_floor: f64,
    pub solver: SolverKind,
}

impl DuctConfig {
    pub fn new(length: f64, cells: usize, area: f64) -> Self {
        assert!(length > 0.0);
        assert!(cells >= 4);
        assert!(area > 0.0);
        Self {
            length,
            cells,
            area,
            cell_areas: vec![area; cells],
            face_areas: vec![area; cells + 1],
            artificial_viscosity: 0.02,
            density_floor: 1.0e-8,
            pressure_floor: 1.0,
            solver: SolverKind::default(),
        }
    }

    pub fn dx(&self) -> f64 {
        self.length / self.cells as f64
    }

    pub fn with_area_profile<F>(length: f64, cells: usize, mut area_at_x: F) -> Self
    where
        F: FnMut(f64) -> f64,
    {
        let mut config = Self::new(length, cells, area_at_x(0.5 * length / cells as f64));
        let dx = config.dx();
        config.cell_areas = (0..cells)
            .map(|i| area_at_x((i as f64 + 0.5) * dx))
            .collect();
        config.face_areas = (0..=cells).map(|i| area_at_x(i as f64 * dx)).collect();
        config.validate_areas();
        config.area = config.cell_areas.iter().sum::<f64>() / cells as f64;
        config
    }

    pub fn cell_area(&self, index: usize) -> f64 {
        self.cell_areas[index]
    }

    pub fn face_area(&self, index: usize) -> f64 {
        self.face_areas[index]
    }

    pub fn area_gradient(&self, index: usize) -> f64 {
        (self.face_areas[index + 1] - self.face_areas[index]) / self.dx()
    }

    pub fn has_variable_area(&self) -> bool {
        self.cell_areas
            .iter()
            .any(|area| (*area - self.area).abs() > 1.0e-12 * self.area)
            || self
                .face_areas
                .iter()
                .any(|area| (*area - self.area).abs() > 1.0e-12 * self.area)
    }

    fn validate_areas(&self) {
        assert_eq!(self.cell_areas.len(), self.cells);
        assert_eq!(self.face_areas.len(), self.cells + 1);
        assert!(
            self.cell_areas
                .iter()
                .all(|area| area.is_finite() && *area > 0.0)
        );
        assert!(
            self.face_areas
                .iter()
                .all(|area| area.is_finite() && *area > 0.0)
        );
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct StepReport {
    pub clipped_cells: usize,
    pub fallback_faces: usize,
}

#[derive(Clone, Debug)]
pub struct Duct<G, L, R>
where
    G: GasProperties,
    L: BoundaryCondition<G>,
    R: BoundaryCondition<G>,
{
    gas: G,
    config: DuctConfig,
    cells: Vec<State>,
    left_boundary: L,
    right_boundary: R,
}

impl<G, L, R> Duct<G, L, R>
where
    G: GasProperties,
    L: BoundaryCondition<G>,
    R: BoundaryCondition<G>,
{
    pub fn new(
        gas: G,
        config: DuctConfig,
        initial_state: State,
        left_boundary: L,
        right_boundary: R,
    ) -> Self {
        Self {
            gas,
            cells: vec![initial_state; config.cells],
            config,
            left_boundary,
            right_boundary,
        }
    }

    pub fn from_initializer<F>(
        gas: G,
        config: DuctConfig,
        left_boundary: L,
        right_boundary: R,
        mut initializer: F,
    ) -> Self
    where
        F: FnMut(f64) -> State,
    {
        let dx = config.dx();
        let cells = (0..config.cells)
            .map(|i| initializer((i as f64 + 0.5) * dx))
            .collect();
        Self {
            gas,
            config,
            cells,
            left_boundary,
            right_boundary,
        }
    }

    pub fn gas(&self) -> G {
        self.gas
    }

    pub fn config(&self) -> DuctConfig {
        self.config.clone()
    }

    pub fn cells(&self) -> &[State] {
        &self.cells
    }

    pub fn primitive_cells(&self) -> Vec<Primitive> {
        self.cells
            .iter()
            .map(|state| state.primitive(self.gas))
            .collect()
    }

    pub fn set_cell(&mut self, index: usize, state: State) {
        self.cells[index] = state;
    }

    pub fn max_signal_speed(&self) -> f64 {
        self.cells
            .iter()
            .map(|state| {
                let prim = state.primitive(self.gas);
                prim.u.abs() + prim.sound_speed
            })
            .fold(0.0, f64::max)
    }

    pub fn step(&mut self, dt: f64) -> StepReport {
        let dx = self.config.dx();
        let lambda = dt / dx;
        let extended = self.extended_states();
        let extended_areas = self.extended_areas();
        let old_cells = self.cells.clone();
        let solver_output = match self.config.solver {
            SolverKind::LaxWendroff => solvers::advance_lax_wendroff(
                &old_cells,
                &extended,
                &self.config.cell_areas,
                &extended_areas,
                &self.config.face_areas,
                lambda,
                self.gas,
                self.config.density_floor,
                self.config.pressure_floor,
            ),
            SolverKind::MacCormack => solvers::advance_mac_cormack(
                &old_cells,
                &extended,
                &self.config.cell_areas,
                &extended_areas,
                &self.config.face_areas,
                lambda,
                self.gas,
                self.config.density_floor,
                self.config.pressure_floor,
                |predicted| self.extended_states_for(predicted),
            ),
        };

        let mut next_cells = solver_output.cells;
        let mut report = StepReport {
            clipped_cells: 0,
            fallback_faces: solver_output.fallback_faces,
        };

        if self.config.artificial_viscosity > 0.0 {
            for i in 1..next_cells.len() - 1 {
                let laplacian = old_cells[i + 1]
                    .minus(old_cells[i].scale(2.0))
                    .plus(old_cells[i - 1]);
                next_cells[i] =
                    next_cells[i].add_scaled(laplacian, self.config.artificial_viscosity);
            }
        }

        for state in &mut next_cells {
            if self.enforce_positivity(state) {
                report.clipped_cells += 1;
            }
        }

        self.cells = next_cells;
        report
    }

    fn extended_states(&self) -> Vec<State> {
        self.extended_states_for(&self.cells)
    }

    fn extended_states_for(&self, cells: &[State]) -> Vec<State> {
        debug_assert_eq!(cells.len(), self.cells.len());
        let mut extended = Vec::with_capacity(self.cells.len() + 2);
        extended.push(
            self.left_boundary
                .ghost_state(cells[0], DuctEnd::Left, self.gas),
        );
        extended.extend_from_slice(cells);
        extended.push(self.right_boundary.ghost_state(
            cells[cells.len() - 1],
            DuctEnd::Right,
            self.gas,
        ));
        extended
    }

    fn extended_areas(&self) -> Vec<f64> {
        let mut extended = Vec::with_capacity(self.cells.len() + 2);
        extended.push(self.config.cell_areas[0]);
        extended.extend_from_slice(&self.config.cell_areas);
        extended.push(self.config.cell_areas[self.config.cell_areas.len() - 1]);
        extended
    }

    fn enforce_positivity(&self, state: &mut State) -> bool {
        let original_rho = state.rho;
        let original_internal_energy = state.specific_internal_energy();
        let prim = state.primitive(self.gas);
        let clipped_rho = state.rho.max(self.config.density_floor);
        let clipped_p = prim.p.max(self.config.pressure_floor);
        if clipped_rho != original_rho || clipped_p != prim.p || original_internal_energy <= 0.0 {
            *state = State::from_primitive(clipped_rho, prim.u, clipped_p, self.gas);
            eprintln!(
                "warning: positivity floor applied: rho={:.6e}, p={:.6e}",
                clipped_rho, clipped_p
            );
            true
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Duct, DuctConfig};
    use crate::{
        boundaries::ClosedEnd, gas_properties::TemperatureDependentAir, solvers::SolverKind,
        state::State,
    };

    fn assert_uniform_closed_duct_remains_uniform(solver: SolverKind) {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let config = DuctConfig {
            solver,
            ..DuctConfig::new(1.0, 32, 1.0)
        };
        let mut duct = Duct::new(gas, config, state, ClosedEnd, ClosedEnd);
        let dt = 0.4 * duct.config().dx() / duct.max_signal_speed();
        let report = duct.step(dt);
        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
        for prim in duct.primitive_cells() {
            assert!((prim.p - 101_325.0).abs() < 1.0e-8);
            assert!(prim.u.abs() < 1.0e-10);
        }
    }

    #[test]
    fn lax_wendroff_uniform_closed_duct_remains_uniform() {
        assert_uniform_closed_duct_remains_uniform(SolverKind::LaxWendroff);
    }

    #[test]
    fn mac_cormack_uniform_closed_duct_remains_uniform() {
        assert_uniform_closed_duct_remains_uniform(SolverKind::MacCormack);
    }

    fn assert_unphysical_predictor_uses_fallback_flux(solver: SolverKind) {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let config = DuctConfig {
            solver,
            ..DuctConfig::new(1.0, 8, 1.0)
        };
        let mut duct = Duct::new(gas, config, state, ClosedEnd, ClosedEnd);
        let bad_state = State {
            rho: -1.2,
            momentum: 0.0,
            rho_total_energy: -10.0,
        };
        duct.set_cell(4, bad_state);
        duct.set_cell(5, bad_state);
        let report = duct.step(1.0e-6);
        assert!(report.fallback_faces > 0);
        assert!(report.clipped_cells > 0);
    }

    #[test]
    fn lax_wendroff_unphysical_predictor_uses_fallback_flux() {
        assert_unphysical_predictor_uses_fallback_flux(SolverKind::LaxWendroff);
    }

    #[test]
    fn mac_cormack_unphysical_predictor_uses_fallback_flux() {
        assert_unphysical_predictor_uses_fallback_flux(SolverKind::MacCormack);
    }
}
