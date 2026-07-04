use onedpipes::{
    DuctConfig, DuctEnd, JunctionPort, Model, ModelBoundary, MultiPipeJunction, State,
    TemperatureDependentAir,
};

#[test]
fn two_pipe_junction_conserves_asymmetric_mass_and_energy() {
    let gas = TemperatureDependentAir::new();
    let ports = [
        JunctionPort::new(
            State::from_primitive(1.2, 18.0, 130_000.0, gas),
            DuctEnd::Right,
            1.0,
        ),
        JunctionPort::new(
            State::from_primitive(0.95, 0.0, 92_000.0, gas),
            DuctEnd::Left,
            1.0,
        ),
    ];

    let solution = MultiPipeJunction.solve(&ports, gas);

    assert!(
        solution.mass_residual().abs() < 1.0e-10,
        "mass residual = {}",
        solution.mass_residual()
    );
    assert!(
        solution.energy_residual().abs() < 1.0e-5,
        "energy residual = {}",
        solution.energy_residual()
    );
}

#[test]
fn three_pipe_junction_conserves_hand_balanced_mass_and_energy() {
    let gas = TemperatureDependentAir::new();
    let rho = 1.2;
    let p = 101_325.0;
    let speed = 5.0;
    let ports = [
        JunctionPort::new(
            State::from_primitive(rho, speed, p, gas),
            DuctEnd::Right,
            2.0,
        ),
        JunctionPort::new(
            State::from_primitive(rho, speed, p, gas),
            DuctEnd::Left,
            1.0,
        ),
        JunctionPort::new(
            State::from_primitive(rho, speed, p, gas),
            DuctEnd::Left,
            1.0,
        ),
    ];

    let solution = MultiPipeJunction.solve(&ports, gas);

    assert!(
        (solution.pressure - p).abs() < 1.0e-8,
        "balanced equal-pressure reference should keep the junction pressure unchanged"
    );
    assert!(
        solution.mass_residual().abs() < 1.0e-10,
        "mass residual = {}",
        solution.mass_residual()
    );
    assert!(
        solution.energy_residual().abs() < 1.0e-5,
        "energy residual = {}",
        solution.energy_residual()
    );
}

#[test]
fn three_pipe_junction_conserves_asymmetric_mass_and_energy() {
    let gas = TemperatureDependentAir::new();
    let ports = [
        JunctionPort::new(
            State::from_primitive(1.2, 24.0, 140_000.0, gas),
            DuctEnd::Right,
            1.0,
        ),
        JunctionPort::new(
            State::from_primitive(0.9, 0.0, 90_000.0, gas),
            DuctEnd::Left,
            0.6,
        ),
        JunctionPort::new(
            State::from_primitive(1.4, 0.0, 95_000.0, gas),
            DuctEnd::Left,
            0.4,
        ),
    ];

    let solution = MultiPipeJunction.solve(&ports, gas);

    for state in &solution.boundary_states {
        assert!((state.primitive(gas).p - solution.pressure).abs() < 1.0e-8);
    }
    assert!(
        solution.mass_residual().abs() < 1.0e-10,
        "mass residual = {}",
        solution.mass_residual()
    );
    assert!(
        solution.energy_residual().abs() < 1.0e-5,
        "energy residual = {}",
        solution.energy_residual()
    );
}

#[test]
fn model_coupled_three_pipe_junction_reports_conservation() {
    let gas = TemperatureDependentAir::new();
    let mut model = Model::new(0.5);
    model.add_uniform_duct(
        gas,
        DuctConfig::new(1.0, 12, 1.0),
        State::from_primitive(1.2, 24.0, 140_000.0, gas),
        ModelBoundary::Closed,
        ModelBoundary::junction(0),
    );
    model.add_uniform_duct(
        gas,
        DuctConfig::new(1.0, 12, 0.6),
        State::from_primitive(0.9, 0.0, 90_000.0, gas),
        ModelBoundary::junction(0),
        ModelBoundary::open(101_325.0),
    );
    model.add_uniform_duct(
        gas,
        DuctConfig::new(1.0, 12, 0.4),
        State::from_primitive(1.4, 0.0, 95_000.0, gas),
        ModelBoundary::junction(0),
        ModelBoundary::open(101_325.0),
    );

    let diagnostics = model.junction_diagnostics();

    assert_eq!(diagnostics.len(), 1);
    assert_eq!(diagnostics[0].junction_id, 0);
    assert!(
        diagnostics[0].solution.mass_residual().abs() < 1.0e-10,
        "mass residual = {}",
        diagnostics[0].solution.mass_residual()
    );
    assert!(
        diagnostics[0].solution.energy_residual().abs() < 1.0e-5,
        "energy residual = {}",
        diagnostics[0].solution.energy_residual()
    );

    let report = model.step_with_dt(1.0e-6);
    assert_eq!(report.clipped_cells, 0);
    assert_eq!(report.fallback_faces, 0);
}
