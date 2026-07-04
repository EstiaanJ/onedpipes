use onedpipes::{ClosedEnd, Duct, DuctConfig, SolverKind, State, TemperatureDependentAir};

fn venturi_area(x: f64) -> f64 {
    let throat = (-((x - 0.5) / 0.16).powi(2)).exp();
    1.0 - 0.45 * throat
}

fn variable_area_config(solver: SolverKind) -> DuctConfig {
    DuctConfig {
        solver,
        artificial_viscosity: 0.0,
        ..DuctConfig::with_area_profile(1.0, 160, venturi_area)
    }
}

fn total_area_weighted_mass(duct: &Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>) -> f64 {
    let dx = duct.config().dx();
    duct.cells()
        .iter()
        .enumerate()
        .map(|(i, state)| state.rho * duct.config().cell_area(i) * dx)
        .sum()
}

fn assert_static_variable_area_duct_remains_uniform(solver: SolverKind) {
    let gas = TemperatureDependentAir::new();
    let state = State::from_primitive(1.2, 0.0, 101_325.0, gas);
    let mut duct = Duct::new(
        gas,
        variable_area_config(solver),
        state,
        ClosedEnd,
        ClosedEnd,
    );
    let dt = 0.45 * duct.config().dx() / duct.max_signal_speed();

    for _ in 0..40 {
        let report = duct.step(dt);
        assert_eq!(report.clipped_cells, 0);
    }

    for prim in duct.primitive_cells() {
        assert!(
            (prim.p - 101_325.0).abs() < 1.0e-6,
            "{} static variable-area duct pressure drifted to {}",
            solver.label(),
            prim.p
        );
        assert!(
            prim.u.abs() < 1.0e-9,
            "{} static variable-area duct developed velocity {}",
            solver.label(),
            prim.u
        );
    }
}

#[test]
fn lax_wendroff_static_venturi_shape_remains_well_balanced() {
    assert_static_variable_area_duct_remains_uniform(SolverKind::LaxWendroff);
}

#[test]
fn mac_cormack_static_venturi_shape_remains_well_balanced() {
    assert_static_variable_area_duct_remains_uniform(SolverKind::MacCormack);
}

fn assert_variable_area_closed_pipe_conserves_area_weighted_mass(solver: SolverKind) {
    let gas = TemperatureDependentAir::new();
    let mut duct = Duct::from_initializer(
        gas,
        variable_area_config(solver),
        ClosedEnd,
        ClosedEnd,
        |x| {
            let shape = (-((x - 0.25) / 0.04).powi(2)).exp();
            State::from_primitive(1.2 + 1.0e-3 * shape, 0.0, 101_325.0, gas)
        },
    );
    let initial_mass = total_area_weighted_mass(&duct);
    let dt = 0.35 * duct.config().dx() / duct.max_signal_speed();

    for _ in 0..120 {
        let report = duct.step(dt);
        assert_eq!(report.clipped_cells, 0);
    }

    let final_mass = total_area_weighted_mass(&duct);
    let relative_drift = ((final_mass - initial_mass) / initial_mass).abs();
    assert!(
        relative_drift < 2.0e-10,
        "{} variable-area mass drift = {relative_drift:.6e}",
        solver.label()
    );
}

#[test]
fn lax_wendroff_variable_area_closed_pipe_conserves_mass() {
    assert_variable_area_closed_pipe_conserves_area_weighted_mass(SolverKind::LaxWendroff);
}

#[test]
fn mac_cormack_variable_area_closed_pipe_conserves_mass() {
    assert_variable_area_closed_pipe_conserves_area_weighted_mass(SolverKind::MacCormack);
}

#[test]
fn venturi_area_profile_has_a_resolved_throat() {
    let config = variable_area_config(SolverKind::MacCormack);
    let throat_area = (0..config.cells)
        .map(|i| config.cell_area(i))
        .fold(f64::INFINITY, f64::min);
    let inlet_area = config.face_area(0);
    let outlet_area = config.face_area(config.cells);

    assert!(throat_area < 0.65 * inlet_area);
    assert!((inlet_area - outlet_area).abs() < 1.0e-12);
    assert!(config.has_variable_area());
}
