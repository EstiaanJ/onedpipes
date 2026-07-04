use crate::{
    boundaries::{BoundaryCondition, DuctEnd},
    gas_properties::GasProperties,
    state::{Primitive, State},
};

#[derive(Clone, Copy, Debug)]
pub struct DuctConfig {
    pub length: f64,
    pub cells: usize,
    pub area: f64,
    pub artificial_viscosity: f64,
    pub density_floor: f64,
    pub pressure_floor: f64,
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
            artificial_viscosity: 0.02,
            density_floor: 1.0e-8,
            pressure_floor: 1.0,
        }
    }

    pub fn dx(&self) -> f64 {
        self.length / self.cells as f64
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct StepReport {
    pub clipped_cells: usize,
    pub fallback_faces: usize,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct BoundaryOverride {
    pub ghost_state: Option<State>,
    pub face_flux: Option<State>,
}

impl BoundaryOverride {
    pub fn ghost(ghost_state: State) -> Self {
        Self {
            ghost_state: Some(ghost_state),
            face_flux: None,
        }
    }

    pub fn flux(ghost_state: State, face_flux: State) -> Self {
        Self {
            ghost_state: Some(ghost_state),
            face_flux: Some(face_flux),
        }
    }
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
        self.config
    }

    pub fn cells(&self) -> &[State] {
        &self.cells
    }

    pub fn left_boundary(&self) -> &L {
        &self.left_boundary
    }

    pub fn right_boundary(&self) -> &R {
        &self.right_boundary
    }

    pub fn primitive_cells(&self) -> Vec<Primitive> {
        self.cells
            .iter()
            .map(|state| state.primitive(self.gas))
            .collect()
    }

    pub fn end_state(&self, end: DuctEnd) -> State {
        match end {
            DuctEnd::Left => self.cells[0],
            DuctEnd::Right => self.cells[self.cells.len() - 1],
        }
    }

    pub fn set_cell(&mut self, index: usize, state: State) {
        self.cells[index] = state;
    }

    pub fn total_mass(&self) -> f64 {
        self.cells
            .iter()
            .map(|state| state.rho * self.config.area * self.config.dx())
            .sum()
    }

    pub fn total_energy(&self) -> f64 {
        self.cells
            .iter()
            .map(|state| state.rho_total_energy * self.config.area * self.config.dx())
            .sum()
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
        self.step_with_boundary_overrides(dt, None, None)
    }

    pub fn step_with_boundary_overrides(
        &mut self,
        dt: f64,
        left_ghost: Option<State>,
        right_ghost: Option<State>,
    ) -> StepReport {
        self.step_with_boundary_controls(
            dt,
            BoundaryOverride {
                ghost_state: left_ghost,
                face_flux: None,
            },
            BoundaryOverride {
                ghost_state: right_ghost,
                face_flux: None,
            },
        )
    }

    pub fn step_with_boundary_controls(
        &mut self,
        dt: f64,
        left_boundary: BoundaryOverride,
        right_boundary: BoundaryOverride,
    ) -> StepReport {
        let dx = self.config.dx();
        let lambda = dt / dx;
        let extended = self.extended_states(left_boundary.ghost_state, right_boundary.ghost_state);
        let mut face_fluxes = Vec::with_capacity(self.cells.len() + 1);
        let mut report = StepReport::default();

        for face in 0..=self.cells.len() {
            let left = extended[face];
            let right = extended[face + 1];
            if self.is_physical(left) && self.is_physical(right) {
                let predicted = left.plus(right).scale(0.5).add_scaled(
                    right.flux(self.gas).minus(left.flux(self.gas)),
                    -0.5 * lambda,
                );
                if !self.is_physical(predicted) {
                    report.fallback_faces += 1;
                    face_fluxes.push(self.rusanov_flux(left, right));
                    continue;
                }
                face_fluxes.push(predicted.flux(self.gas));
            } else {
                report.fallback_faces += 1;
                face_fluxes.push(self.rusanov_flux(left, right));
            }
        }
        if let Some(face_flux) = left_boundary.face_flux {
            face_fluxes[0] = face_flux;
        }
        if let Some(face_flux) = right_boundary.face_flux {
            let right_face = self.cells.len();
            face_fluxes[right_face] = face_flux;
        }

        let old_cells = self.cells.clone();
        let mut next_cells = Vec::with_capacity(self.cells.len());
        for i in 0..self.cells.len() {
            let area_weighted_flux_delta = face_fluxes[i + 1]
                .minus(face_fluxes[i])
                .scale(self.config.area);
            let next =
                self.cells[i].add_scaled(area_weighted_flux_delta, -lambda / self.config.area);
            next_cells.push(next);
        }

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

    fn extended_states(&self, left_ghost: Option<State>, right_ghost: Option<State>) -> Vec<State> {
        let mut extended = Vec::with_capacity(self.cells.len() + 2);
        extended.push(left_ghost.unwrap_or_else(|| {
            self.left_boundary
                .ghost_state(self.cells[0], DuctEnd::Left, self.gas)
        }));
        extended.extend_from_slice(&self.cells);
        extended.push(right_ghost.unwrap_or_else(|| {
            self.right_boundary.ghost_state(
                self.cells[self.cells.len() - 1],
                DuctEnd::Right,
                self.gas,
            )
        }));
        extended
    }

    fn is_physical(&self, state: State) -> bool {
        if !state.rho.is_finite() || state.rho <= self.config.density_floor {
            return false;
        }
        let Ok(prim) = state.try_primitive(self.gas) else {
            return false;
        };
        prim.p.is_finite() && prim.p > self.config.pressure_floor
    }

    fn enforce_positivity(&self, state: &mut State) -> bool {
        let original_rho = state.rho;
        let original_internal_energy = state.specific_internal_energy();
        let prim = state.primitive_clamped(self.gas);
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

    fn rusanov_flux(&self, left: State, right: State) -> State {
        let left_prim = left.primitive_clamped(self.gas);
        let right_prim = right.primitive_clamped(self.gas);
        let wave_speed = (left_prim.u.abs() + left_prim.sound_speed)
            .max(right_prim.u.abs() + right_prim.sound_speed);
        left.flux_clamped(self.gas)
            .plus(right.flux_clamped(self.gas))
            .scale(0.5)
            .add_scaled(right.minus(left), -0.5 * wave_speed)
    }
}

#[cfg(test)]
mod tests {
    use super::{BoundaryOverride, Duct, DuctConfig};
    use crate::{boundaries::ClosedEnd, gas_properties::TemperatureDependentAir, state::State};

    #[test]
    fn uniform_closed_duct_remains_uniform() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut duct = Duct::new(
            gas,
            DuctConfig::new(1.0, 32, 1.0),
            state,
            ClosedEnd,
            ClosedEnd,
        );
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
    fn conserved_inventory_scales_with_duct_area() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let small = Duct::new(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ClosedEnd,
            ClosedEnd,
        );
        let large = Duct::new(
            gas,
            DuctConfig::new(1.0, 8, 2.5),
            state,
            ClosedEnd,
            ClosedEnd,
        );

        assert!((large.total_mass() / small.total_mass() - 2.5).abs() < 1.0e-12);
        assert!((large.total_energy() / small.total_energy() - 2.5).abs() < 1.0e-12);
    }

    #[test]
    fn step_accepts_precomputed_boundary_ghosts_for_coupled_boundaries() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut duct = Duct::new(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ClosedEnd,
            ClosedEnd,
        );
        let ghost = State::from_primitive(1.2, 0.0, 101_325.0, gas);

        let report = duct.step_with_boundary_overrides(1.0e-6, Some(ghost), Some(ghost));

        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
    }

    #[test]
    fn step_accepts_boundary_flux_overrides() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut duct = Duct::new(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ClosedEnd,
            ClosedEnd,
        );
        let ghost = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let left_flux = State {
            rho: 0.0,
            momentum: 101_325.0,
            rho_total_energy: 0.0,
        };
        let right_flux = left_flux;

        let report = duct.step_with_boundary_controls(
            1.0e-6,
            BoundaryOverride::flux(ghost, left_flux),
            BoundaryOverride::flux(ghost, right_flux),
        );

        assert_eq!(report.clipped_cells, 0);
        assert_eq!(report.fallback_faces, 0);
    }

    #[test]
    fn unphysical_predictor_uses_fallback_flux() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
        let mut duct = Duct::new(
            gas,
            DuctConfig::new(1.0, 8, 1.0),
            state,
            ClosedEnd,
            ClosedEnd,
        );
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
}
