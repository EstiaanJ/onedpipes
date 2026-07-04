use onedpipes::{
    DuctConfig, ExternalBoundaryControl, ExternalBoundaryId, GasProperties, Model, ModelBoundary,
    SpeciesFractions, State, TemperatureDependentAir,
};

fn state_at(gas: TemperatureDependentAir, pressure: f64, temperature: f64) -> State {
    State::from_primitive(pressure / (gas.r() * temperature), 0.0, pressure, gas)
}

fn enthalpy(gas: TemperatureDependentAir, temperature: f64) -> f64 {
    gas.internal_energy_from_temperature(temperature) + gas.r() * temperature
}

#[test]
fn external_hot_blowdown_pulse_stays_bounded_and_diagnostic() {
    let gas = TemperatureDependentAir::new();
    let initial = state_at(gas, 101_325.0, 300.0);
    let mut model = Model::new(0.35);
    let pipe = model.add_uniform_duct_with_species(
        gas,
        DuctConfig {
            artificial_viscosity: 0.035,
            ..DuctConfig::new(0.8, 48, 3.0e-4)
        },
        initial,
        SpeciesFractions::AIR,
        ModelBoundary::external(0),
        ModelBoundary::open(101_325.0),
    );
    let exhaust = SpeciesFractions::EXHAUST;
    let hot_enthalpy = enthalpy(gas, 1200.0);
    let mut total_clipped = 0;
    let mut total_fallback = 0;
    let mut limited_steps = 0;

    for step in 0..240 {
        let pulse = if step < 80 { -0.012 } else { 0.003 };
        model.set_external_boundary_control(
            ExternalBoundaryId(0),
            ExternalBoundaryControl::BoundedFlow {
                mass_flow_out: pulse,
                energy_flow_out: pulse * hot_enthalpy,
                max_mass_transfer: 2.0e-10,
                max_energy_transfer: 0.30,
                inflow_species: exhaust,
            },
        );
        let report = model.step_with_dt(5.0e-8);
        total_clipped += report.clipped_cells;
        total_fallback += report.fallback_faces;
        limited_steps += report
            .external_boundary_diagnostics
            .iter()
            .filter(|diagnostic| diagnostic.limited)
            .count();
    }

    assert_eq!(total_clipped, 0);
    assert_eq!(total_fallback, 0);
    assert!(limited_steps > 0);

    for primitive in model.pipe_primitive_cells(pipe) {
        assert!(primitive.p.is_finite());
        assert!(primitive.temperature.is_finite());
        assert!(primitive.p > 50_000.0 && primitive.p < 300_000.0);
        assert!(primitive.temperature > 250.0 && primitive.temperature < 1500.0);
    }
    assert!(model.pipe_species_cells(pipe)[0].products > 0.0);
}

#[test]
fn pipe_end_species_and_lambda_are_exposed_for_host_engine() {
    let gas = TemperatureDependentAir::new();
    let initial = state_at(gas, 101_325.0, 300.0);
    let rich = SpeciesFractions::new(0.12, 0.04, 0.74, 0.10);
    let mut model = Model::new(0.5);
    let pipe = model.add_uniform_duct_with_species(
        gas,
        DuctConfig::new(0.3, 8, 1.0e-3),
        initial,
        rich,
        ModelBoundary::external(3),
        ModelBoundary::Closed,
    );

    let port = model.external_ports().remove(0);

    assert_eq!(port.pipe_id, pipe);
    assert_eq!(port.species, rich);
    assert!((port.species.lambda(3.4).unwrap() - 0.12 / (3.4 * 0.04)).abs() < 1.0e-12);
}
