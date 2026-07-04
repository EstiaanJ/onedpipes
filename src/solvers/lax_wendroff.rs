use crate::{
    gas_properties::GasProperties,
    solvers::{
        SolverOutput, area_weighted_flux, area_weighted_state, is_physical, rusanov_flux,
        unweighted_state,
    },
    state::State,
};

pub(crate) fn advance<G: GasProperties>(
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
    debug_assert_eq!(extended.len(), cells.len() + 2);
    debug_assert_eq!(cell_areas.len(), cells.len());
    debug_assert_eq!(extended_areas.len(), extended.len());
    debug_assert_eq!(face_areas.len(), cells.len() + 1);

    let mut face_fluxes = Vec::with_capacity(cells.len() + 1);
    let mut fallback_faces = 0;

    for face in 0..=cells.len() {
        let left = extended[face];
        let right = extended[face + 1];
        let left_area = extended_areas[face];
        let right_area = extended_areas[face + 1];
        let face_area = face_areas[face];
        let predictor_area = 0.5 * (left_area + right_area);
        let mut predicted_area_weighted = area_weighted_state(left, left_area)
            .plus(area_weighted_state(right, right_area))
            .scale(0.5)
            .add_scaled(
                area_weighted_flux(right, right_area, gas)
                    .minus(area_weighted_flux(left, left_area, gas)),
                -0.5 * lambda,
            );
        let left_pressure = left.primitive(gas).p;
        let right_pressure = right.primitive(gas).p;
        predicted_area_weighted.momentum +=
            0.5 * lambda * 0.5 * (left_pressure + right_pressure) * (right_area - left_area);
        let predicted = unweighted_state(predicted_area_weighted, predictor_area);

        if is_physical(predicted, gas, density_floor, pressure_floor) {
            face_fluxes.push(area_weighted_flux(predicted, face_area, gas));
        } else {
            fallback_faces += 1;
            face_fluxes.push(rusanov_flux(left, right, gas).scale(face_area));
        }
    }

    let mut next_cells = Vec::with_capacity(cells.len());
    for i in 0..cells.len() {
        let mut next_area_weighted = area_weighted_state(cells[i], cell_areas[i])
            .add_scaled(face_fluxes[i + 1].minus(face_fluxes[i]), -lambda);
        next_area_weighted.momentum +=
            lambda * cells[i].primitive(gas).p * (face_areas[i + 1] - face_areas[i]);
        next_cells.push(unweighted_state(next_area_weighted, cell_areas[i]));
    }

    SolverOutput {
        cells: next_cells,
        fallback_faces,
    }
}
