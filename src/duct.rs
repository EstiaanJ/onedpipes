use crate::{
    boundaries::{BoundaryCondition, DuctEnd},
    gas_properties::GasProperties,
    solvers::{self, SolverKind},
    species::{SpeciesFractions, SpeciesMass},
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

#[derive(Clone, Debug, Default)]
pub struct StepReport {
    pub clipped_cells: usize,
    pub fallback_faces: usize,
    pub clipped_cell_indices: Vec<usize>,
    pub fallback_face_indices: Vec<usize>,
    pub pipe_diagnostics: Vec<PipeStepDiagnostic>,
    pub external_boundary_diagnostics: Vec<ExternalBoundaryStepDiagnostic>,
}

impl StepReport {
    pub fn absorb(&mut self, other: Self) {
        self.clipped_cells += other.clipped_cells;
        self.fallback_faces += other.fallback_faces;
        self.clipped_cell_indices.extend(other.clipped_cell_indices);
        self.fallback_face_indices
            .extend(other.fallback_face_indices);
        self.pipe_diagnostics.extend(other.pipe_diagnostics);
        self.external_boundary_diagnostics
            .extend(other.external_boundary_diagnostics);
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PipeStepDiagnostic {
    pub pipe_index: usize,
    pub clipped_cell_indices: Vec<usize>,
    pub fallback_face_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ExternalBoundaryStepDiagnostic {
    pub external_id: usize,
    pub pipe_index: usize,
    pub end: DuctEnd,
    pub requested_mass_flow_out: f64,
    pub accepted_mass_flow_out: f64,
    pub requested_energy_flow_out: f64,
    pub accepted_energy_flow_out: f64,
    pub mass_transferred_out: f64,
    pub energy_transferred_out: f64,
    pub species_transferred_out: SpeciesMass,
    pub limited: bool,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct BoundaryOverride {
    pub ghost_state: Option<State>,
    pub ghost_species: Option<SpeciesFractions>,
    pub face_flux: Option<State>,
    pub face_species_flux: Option<SpeciesMass>,
}

impl BoundaryOverride {
    pub fn ghost(ghost_state: State) -> Self {
        Self {
            ghost_state: Some(ghost_state),
            ghost_species: None,
            face_flux: None,
            face_species_flux: None,
        }
    }

    pub fn ghost_with_species(ghost_state: State, ghost_species: SpeciesFractions) -> Self {
        Self {
            ghost_state: Some(ghost_state),
            ghost_species: Some(ghost_species),
            face_flux: None,
            face_species_flux: None,
        }
    }

    pub fn flux(ghost_state: State, face_flux: State) -> Self {
        Self {
            ghost_state: Some(ghost_state),
            ghost_species: None,
            face_flux: Some(face_flux),
            face_species_flux: None,
        }
    }

    pub fn flux_with_species(
        ghost_state: State,
        ghost_species: SpeciesFractions,
        face_flux: State,
        face_species_flux: SpeciesMass,
    ) -> Self {
        Self {
            ghost_state: Some(ghost_state),
            ghost_species: Some(ghost_species),
            face_flux: Some(face_flux),
            face_species_flux: Some(face_species_flux),
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
    species: Vec<SpeciesFractions>,
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
            species: vec![SpeciesFractions::AIR; config.cells],
            config,
            left_boundary,
            right_boundary,
        }
    }

    pub fn new_with_species(
        gas: G,
        config: DuctConfig,
        initial_state: State,
        initial_species: SpeciesFractions,
        left_boundary: L,
        right_boundary: R,
    ) -> Self {
        Self {
            gas,
            cells: vec![initial_state; config.cells],
            species: vec![initial_species.normalized(); config.cells],
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
            species: vec![SpeciesFractions::AIR; config.cells],
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

    pub fn species_cells(&self) -> &[SpeciesFractions] {
        &self.species
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

    pub fn end_species(&self, end: DuctEnd) -> SpeciesFractions {
        match end {
            DuctEnd::Left => self.species[0],
            DuctEnd::Right => self.species[self.species.len() - 1],
        }
    }

    pub fn set_cell(&mut self, index: usize, state: State) {
        self.cells[index] = state;
    }

    pub fn set_species(&mut self, index: usize, species: SpeciesFractions) {
        self.species[index] = species.normalized();
    }

    pub fn total_mass(&self) -> f64 {
        self.cells
            .iter()
            .enumerate()
            .map(|(i, state)| state.rho * self.config.cell_area(i) * self.config.dx())
            .sum()
    }

    pub fn total_energy(&self) -> f64 {
        self.cells
            .iter()
            .enumerate()
            .map(|(i, state)| state.rho_total_energy * self.config.cell_area(i) * self.config.dx())
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
                ghost_species: None,
                face_flux: None,
                face_species_flux: None,
            },
            BoundaryOverride {
                ghost_state: right_ghost,
                ghost_species: None,
                face_flux: None,
                face_species_flux: None,
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
        let extended_species =
            self.extended_species(left_boundary.ghost_species, right_boundary.ghost_species);
        let old_cells = self.cells.clone();
        let extended_areas = self.extended_areas();
        let mut solver_output = match self.config.solver {
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
                |predicted| {
                    self.extended_states_for(
                        predicted,
                        left_boundary.ghost_state,
                        right_boundary.ghost_state,
                    )
                },
            ),
        };
        let mut report = StepReport {
            clipped_cells: 0,
            fallback_faces: solver_output.fallback_faces,
            clipped_cell_indices: Vec::new(),
            fallback_face_indices: solver_output.fallback_face_indices,
            pipe_diagnostics: Vec::new(),
            external_boundary_diagnostics: Vec::new(),
        };

        if let Some(face_flux) = left_boundary.face_flux {
            solver_output.face_fluxes[0] = face_flux.scale(self.config.face_area(0));
        }
        if let Some(face_flux) = right_boundary.face_flux {
            let right_face = self.cells.len();
            solver_output.face_fluxes[right_face] =
                face_flux.scale(self.config.face_area(right_face));
        }

        let face_species_fluxes: Vec<_> = (0..=self.cells.len())
            .map(|face| {
                let mass_flux = solver_output.face_fluxes[face].rho;
                if face == 0 {
                    if let Some(face_species_flux) = left_boundary.face_species_flux {
                        return face_species_flux.scale(self.config.face_area(face));
                    }
                } else if face == self.cells.len() {
                    if let Some(face_species_flux) = right_boundary.face_species_flux {
                        return face_species_flux.scale(self.config.face_area(face));
                    }
                }
                upwind_species_flux(
                    mass_flux,
                    extended_species[face],
                    extended_species[face + 1],
                )
            })
            .collect();

        let old_species_mass: Vec<_> = old_cells
            .iter()
            .zip(self.species.iter().copied())
            .map(|(state, species)| SpeciesMass::from_density(state.rho, species))
            .collect();
        let mut next_species_mass = Vec::with_capacity(self.cells.len());
        for i in 0..self.cells.len() {
            next_species_mass.push(
                old_species_mass[i]
                    .add_scaled(face_species_fluxes[i + 1], -lambda)
                    .add_scaled(face_species_fluxes[i], lambda),
            );
        }
        let mut next_cells = solver_output.cells;

        if self.config.artificial_viscosity > 0.0 {
            for i in 1..next_cells.len() - 1 {
                let laplacian = old_cells[i + 1]
                    .minus(old_cells[i].scale(2.0))
                    .plus(old_cells[i - 1]);
                next_cells[i] =
                    next_cells[i].add_scaled(laplacian, self.config.artificial_viscosity);
            }
        }

        for (index, state) in next_cells.iter_mut().enumerate() {
            if self.enforce_positivity(state) {
                report.clipped_cells += 1;
                report.clipped_cell_indices.push(index);
            }
        }

        self.species = next_species_mass
            .into_iter()
            .map(SpeciesMass::fractions)
            .collect();
        self.cells = next_cells;
        report
    }

    fn extended_states(&self, left_ghost: Option<State>, right_ghost: Option<State>) -> Vec<State> {
        self.extended_states_for(&self.cells, left_ghost, right_ghost)
    }

    fn extended_states_for(
        &self,
        cells: &[State],
        left_ghost: Option<State>,
        right_ghost: Option<State>,
    ) -> Vec<State> {
        debug_assert_eq!(cells.len(), self.cells.len());
        let mut extended = Vec::with_capacity(self.cells.len() + 2);
        extended.push(left_ghost.unwrap_or_else(|| {
            self.left_boundary
                .ghost_state(cells[0], DuctEnd::Left, self.gas)
        }));
        extended.extend_from_slice(cells);
        extended.push(right_ghost.unwrap_or_else(|| {
            self.right_boundary
                .ghost_state(cells[cells.len() - 1], DuctEnd::Right, self.gas)
        }));
        extended
    }

    fn extended_species(
        &self,
        left_ghost: Option<SpeciesFractions>,
        right_ghost: Option<SpeciesFractions>,
    ) -> Vec<SpeciesFractions> {
        let mut extended = Vec::with_capacity(self.species.len() + 2);
        extended.push(left_ghost.unwrap_or(self.species[0]));
        extended.extend_from_slice(&self.species);
        extended.push(right_ghost.unwrap_or(self.species[self.species.len() - 1]));
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
}

fn upwind_species_flux(
    mass_flux: f64,
    left: SpeciesFractions,
    right: SpeciesFractions,
) -> SpeciesMass {
    if mass_flux >= 0.0 {
        left.scale(mass_flux)
    } else {
        right.scale(mass_flux)
    }
}

#[cfg(test)]
mod tests {
    use super::{BoundaryOverride, Duct, DuctConfig};
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
