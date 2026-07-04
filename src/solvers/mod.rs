mod lax_wendroff;
mod mac_cormack;

use crate::{gas_properties::GasProperties, state::State};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SolverKind {
    #[default]
    LaxWendroff,
    MacCormack,
}

impl SolverKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::LaxWendroff => "Lax-Wendroff",
            Self::MacCormack => "MacCormack",
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct SolverOutput {
    pub cells: Vec<State>,
    pub fallback_faces: usize,
}

pub(crate) fn advance_lax_wendroff<G: GasProperties>(
    cells: &[State],
    extended: &[State],
    cell_areas: &[f64],
    extended_areas: &[f64],
    face_areas: &[f64],
    lambda: f64,
    gas: G,
    density_floor: f64,
    pressure_floor: f64,
) -> SolverOutput {
    lax_wendroff::advance(
        cells,
        extended,
        cell_areas,
        extended_areas,
        face_areas,
        lambda,
        gas,
        density_floor,
        pressure_floor,
    )
}

pub(crate) fn advance_mac_cormack<G, F>(
    cells: &[State],
    extended: &[State],
    cell_areas: &[f64],
    extended_areas: &[f64],
    face_areas: &[f64],
    lambda: f64,
    gas: G,
    density_floor: f64,
    pressure_floor: f64,
    extend_predicted: F,
) -> SolverOutput
where
    G: GasProperties,
    F: FnOnce(&[State]) -> Vec<State>,
{
    mac_cormack::advance(
        cells,
        extended,
        cell_areas,
        extended_areas,
        face_areas,
        lambda,
        gas,
        density_floor,
        pressure_floor,
        extend_predicted,
    )
}

pub(crate) fn area_weighted_state(state: State, area: f64) -> State {
    state.scale(area)
}

pub(crate) fn unweighted_state(area_weighted: State, area: f64) -> State {
    area_weighted.scale(1.0 / area)
}

pub(crate) fn area_weighted_flux<G: GasProperties>(state: State, area: f64, gas: G) -> State {
    state.flux(gas).scale(area)
}

pub(crate) fn is_physical<G: GasProperties>(
    state: State,
    gas: G,
    density_floor: f64,
    pressure_floor: f64,
) -> bool {
    if !state.rho.is_finite() || state.rho <= density_floor {
        return false;
    }
    let internal_energy = state.specific_internal_energy();
    if !internal_energy.is_finite() || internal_energy <= 0.0 {
        return false;
    }
    let prim = state.primitive(gas);
    prim.p.is_finite() && prim.p > pressure_floor
}

pub(crate) fn rusanov_flux<G: GasProperties>(left: State, right: State, gas: G) -> State {
    let left_prim = left.primitive(gas);
    let right_prim = right.primitive(gas);
    let wave_speed = (left_prim.u.abs() + left_prim.sound_speed)
        .max(right_prim.u.abs() + right_prim.sound_speed);
    left.flux(gas)
        .plus(right.flux(gas))
        .scale(0.5)
        .add_scaled(right.minus(left), -0.5 * wave_speed)
}

pub(crate) fn rusanov_cell_update<G: GasProperties>(
    cell: State,
    left_ghost_or_neighbor: State,
    right_ghost_or_neighbor: State,
    cell_area: f64,
    left_face_area: f64,
    right_face_area: f64,
    lambda: f64,
    gas: G,
) -> State {
    let left_flux = rusanov_flux(left_ghost_or_neighbor, cell, gas).scale(left_face_area);
    let right_flux = rusanov_flux(cell, right_ghost_or_neighbor, gas).scale(right_face_area);
    let next_area_weighted =
        area_weighted_state(cell, cell_area).add_scaled(right_flux.minus(left_flux), -lambda);
    unweighted_state(next_area_weighted, cell_area)
}
