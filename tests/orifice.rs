use onedpipes::{
    DuctConfig, GasProperties, Model, ModelBoundary, State, TemperatureDependentAir, ValveOrifice,
};

fn state_at(gas: TemperatureDependentAir, pressure: f64, temperature: f64) -> State {
    State::from_primitive(pressure / (gas.r() * temperature), 0.0, pressure, gas)
}

fn hand_computed_orifice_flow(
    gas: TemperatureDependentAir,
    discharge_coefficient: f64,
    area: f64,
    upstream_pressure: f64,
    upstream_temperature: f64,
    downstream_pressure: f64,
) -> f64 {
    let gamma = gas.gamma(upstream_temperature);
    let pressure_ratio = downstream_pressure / upstream_pressure;
    let critical_pressure_ratio = (2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0));
    let flow_factor = if pressure_ratio <= critical_pressure_ratio {
        (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    } else {
        let pressure_term =
            pressure_ratio.powf(2.0 / gamma) - pressure_ratio.powf((gamma + 1.0) / gamma);
        (2.0 / (gamma - 1.0) * pressure_term).sqrt()
    };

    discharge_coefficient
        * area
        * upstream_pressure
        * (gamma / (gas.r() * upstream_temperature)).sqrt()
        * flow_factor
}

#[test]
fn orifice_matches_hand_computed_discharge_over_pressure_ratios() {
    let gas = TemperatureDependentAir::new();
    let discharge_coefficient = 0.82;
    let area = 1.4e-4;
    let valve = ValveOrifice::new(discharge_coefficient, area);
    let upstream_pressure = 210_000.0;
    let upstream_temperature = 320.0;
    let upstream = state_at(gas, upstream_pressure, upstream_temperature);

    for pressure_ratio in [0.90, 0.75, 0.60, 0.50, 0.35] {
        let downstream_pressure = upstream_pressure * pressure_ratio;
        let downstream = state_at(gas, downstream_pressure, upstream_temperature);
        let expected = hand_computed_orifice_flow(
            gas,
            discharge_coefficient,
            area,
            upstream_pressure,
            upstream_temperature,
            downstream_pressure,
        );

        let flow = valve.mass_flow(upstream, downstream, gas);
        let relative_error = ((flow.mass_flow - expected) / expected).abs();

        assert!(
            relative_error < 1.0e-12,
            "pressure_ratio={pressure_ratio}, measured={}, expected={}, relative_error={relative_error}",
            flow.mass_flow,
            expected
        );
    }
}

#[test]
fn model_coupled_orifice_reports_hand_computed_discharge() {
    let gas = TemperatureDependentAir::new();
    let discharge_coefficient = 0.82;
    let area = 1.4e-4;
    let valve = ValveOrifice::new(discharge_coefficient, area);
    let upstream_pressure = 210_000.0;
    let downstream_pressure = 126_000.0;
    let temperature = 320.0;
    let upstream = state_at(gas, upstream_pressure, temperature);
    let downstream = state_at(gas, downstream_pressure, temperature);
    let expected = hand_computed_orifice_flow(
        gas,
        discharge_coefficient,
        area,
        upstream_pressure,
        temperature,
        downstream_pressure,
    );
    let mut model = Model::new(0.5);
    model.add_uniform_duct(
        gas,
        DuctConfig::new(1.0, 12, 1.0),
        upstream,
        ModelBoundary::Closed,
        ModelBoundary::orifice(0, valve),
    );
    model.add_uniform_duct(
        gas,
        DuctConfig::new(1.0, 12, 1.0),
        downstream,
        ModelBoundary::orifice(0, valve),
        ModelBoundary::open(101_325.0),
    );

    let diagnostics = model.orifice_diagnostics();

    assert_eq!(diagnostics.len(), 1);
    let relative_error = ((diagnostics[0].flow.mass_flow - expected) / expected).abs();
    assert!(
        relative_error < 1.0e-12,
        "measured={}, expected={}, relative_error={relative_error}",
        diagnostics[0].flow.mass_flow,
        expected
    );

    let report = model.step_with_dt(1.0e-7);
    assert_eq!(report.clipped_cells, 0);
    assert_eq!(report.fallback_faces, 0);
}
