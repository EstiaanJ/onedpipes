use crate::{
    gas_properties::GasProperties,
    solvers::{
        SolverOutput, area_weighted_flux, area_weighted_state, is_physical, rusanov_cell_update,
        unweighted_state,
    },
    state::State,
};

pub(crate) fn advance<G, F>(
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
    debug_assert_eq!(extended.len(), cells.len() + 2);
    debug_assert_eq!(cell_areas.len(), cells.len());
    debug_assert_eq!(extended_areas.len(), extended.len());
    debug_assert_eq!(face_areas.len(), cells.len() + 1);

    let mut fallback_faces = 0;
    let mut predictor_is_physical = Vec::with_capacity(cells.len());
    let mut predicted_area_weighted_cells = Vec::with_capacity(cells.len());
    let mut predicted_cells = Vec::with_capacity(cells.len());

    for i in 0..cells.len() {
        let cell = extended[i + 1];
        let right = extended[i + 2];
        let cell_area = extended_areas[i + 1];
        let right_area = extended_areas[i + 2];
        let mut predicted_area_weighted = area_weighted_state(cell, cell_area).add_scaled(
            area_weighted_flux(right, right_area, gas)
                .minus(area_weighted_flux(cell, cell_area, gas)),
            -lambda,
        );
        predicted_area_weighted.momentum +=
            lambda * cell.primitive(gas).p * (right_area - cell_area);
        let predicted = unweighted_state(predicted_area_weighted, cell_areas[i]);
        let physical = is_physical(predicted, gas, density_floor, pressure_floor);
        predictor_is_physical.push(physical);
        if physical {
            predicted_area_weighted_cells.push(predicted_area_weighted);
            predicted_cells.push(predicted);
        } else {
            fallback_faces += 1;
            predicted_area_weighted_cells.push(area_weighted_state(cell, cell_area));
            predicted_cells.push(cell);
        }
    }

    let predicted_extended = extend_predicted(&predicted_cells);
    debug_assert_eq!(predicted_extended.len(), cells.len() + 2);

    let mut next_cells = Vec::with_capacity(cells.len());
    for i in 0..cells.len() {
        let predicted_left = predicted_extended[i];
        let predicted_cell = predicted_extended[i + 1];
        let can_correct = predictor_is_physical[i]
            && is_physical(predicted_left, gas, density_floor, pressure_floor)
            && is_physical(predicted_cell, gas, density_floor, pressure_floor);

        if can_correct {
            let predicted_left_area = extended_areas[i];
            let predicted_cell_area = cell_areas[i];
            let mut corrected_area_weighted = area_weighted_state(cells[i], cell_areas[i])
                .plus(predicted_area_weighted_cells[i])
                .add_scaled(
                    area_weighted_flux(predicted_cell, predicted_cell_area, gas)
                        .minus(area_weighted_flux(predicted_left, predicted_left_area, gas)),
                    -lambda,
                )
                .scale(0.5);
            corrected_area_weighted.momentum += 0.5
                * lambda
                * predicted_cell.primitive(gas).p
                * (predicted_cell_area - predicted_left_area);
            let corrected = unweighted_state(corrected_area_weighted, cell_areas[i]);

            if is_physical(corrected, gas, density_floor, pressure_floor) {
                next_cells.push(corrected);
                continue;
            }
        }

        fallback_faces += 1;
        next_cells.push(rusanov_cell_update(
            cells[i],
            extended[i],
            extended[i + 2],
            cell_areas[i],
            face_areas[i],
            face_areas[i + 1],
            lambda,
            gas,
        ));
    }

    SolverOutput {
        cells: next_cells,
        fallback_faces,
    }
}
