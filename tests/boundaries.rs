use onedpipes::{
    BoundaryCondition, ClosedEnd, Duct, DuctConfig, OpenEnd, SolverKind, State,
    TemperatureDependentAir,
};

fn right_going_pulse_state(
    gas: TemperatureDependentAir,
    x: f64,
    center: f64,
    width: f64,
    base_rho: f64,
    base_p: f64,
    amplitude: f64,
) -> State {
    let base = State::from_primitive(base_rho, 0.0, base_p, gas).primitive(gas);
    let shape = (-((x - center) / width).powi(2)).exp();
    let dp = amplitude * shape;
    let rho = base_rho + dp / (base.sound_speed * base.sound_speed);
    let u = dp / (base_rho * base.sound_speed);
    State::from_primitive(rho, u, base_p + dp, gas)
}

fn reflected_pressure_peak<R>(solver: SolverKind, right_boundary: R) -> (f64, f64)
where
    R: BoundaryCondition<TemperatureDependentAir> + Copy,
{
    let gas = TemperatureDependentAir::new();
    let length = 1.0;
    let cells = 220;
    let base_rho = 1.2;
    let base_p = 101_325.0;
    let pulse_center = 0.28;
    let pulse_width = 0.035;
    let pulse_amplitude = 120.0;
    let probe_x = 0.72;
    let probe_cell = ((probe_x / length) * cells as f64).floor() as usize;
    let config = DuctConfig {
        artificial_viscosity: 0.003,
        solver,
        ..DuctConfig::new(length, cells, 1.0)
    };
    let mut duct = Duct::from_initializer(gas, config, ClosedEnd, right_boundary, |x| {
        right_going_pulse_state(
            gas,
            x,
            pulse_center,
            pulse_width,
            base_rho,
            base_p,
            pulse_amplitude,
        )
    });

    let sound_speed = State::from_primitive(base_rho, 0.0, base_p, gas)
        .primitive(gas)
        .sound_speed;
    let incident_time = (probe_x - pulse_center) / sound_speed;
    let reflected_time = (2.0 * length - pulse_center - probe_x) / sound_speed;
    let end_time = reflected_time + 0.35 * (length / sound_speed);
    let mut time = 0.0;
    let mut incident_peak = 0.0_f64;
    let mut reflected_positive_peak = 0.0_f64;
    let mut reflected_negative_peak = 0.0_f64;

    while time < end_time {
        let mut dt = 0.9 * 0.5 * duct.config().dx() / duct.max_signal_speed();
        if time + dt > end_time {
            dt = end_time - time;
        }
        let report = duct.step(dt);
        assert_eq!(report.clipped_cells, 0);
        time += dt;

        let pressure_deviation = duct.cells()[probe_cell].primitive(gas).p - base_p;
        if (time - incident_time).abs() < 0.22 * (length / sound_speed) {
            incident_peak = incident_peak.max(pressure_deviation);
        }
        if (time - reflected_time).abs() < 0.22 * (length / sound_speed) {
            reflected_positive_peak = reflected_positive_peak.max(pressure_deviation);
            reflected_negative_peak = reflected_negative_peak.min(pressure_deviation);
        }
    }

    assert!(
        incident_peak > 0.25 * pulse_amplitude,
        "incident peak too small: {incident_peak}"
    );

    (reflected_positive_peak, reflected_negative_peak)
}

#[test]
fn lax_wendroff_closed_end_reflects_pressure_pulse_with_same_sign() {
    let (positive_peak, negative_peak) =
        reflected_pressure_peak(SolverKind::LaxWendroff, ClosedEnd);
    assert!(
        positive_peak > 20.0,
        "closed end should return a positive pressure pulse, got {positive_peak}"
    );
    assert!(
        positive_peak > negative_peak.abs(),
        "closed-end reflected pulse should be dominated by same-sign pressure"
    );
}

#[test]
fn mac_cormack_closed_end_reflects_pressure_pulse_with_same_sign() {
    let (positive_peak, negative_peak) = reflected_pressure_peak(SolverKind::MacCormack, ClosedEnd);
    assert!(
        positive_peak > 20.0,
        "closed end should return a positive pressure pulse, got {positive_peak}"
    );
    assert!(
        positive_peak > negative_peak.abs(),
        "closed-end reflected pulse should be dominated by same-sign pressure"
    );
}

#[test]
fn lax_wendroff_open_end_reflects_pressure_pulse_with_opposite_sign() {
    let (positive_peak, negative_peak) =
        reflected_pressure_peak(SolverKind::LaxWendroff, OpenEnd::new(101_325.0));
    assert!(
        negative_peak < -20.0,
        "open end should return a negative pressure pulse, got {negative_peak}"
    );
    assert!(
        negative_peak.abs() > positive_peak,
        "open-end reflected pulse should be dominated by inverted pressure"
    );
}

#[test]
fn mac_cormack_open_end_reflects_pressure_pulse_with_opposite_sign() {
    let (positive_peak, negative_peak) =
        reflected_pressure_peak(SolverKind::MacCormack, OpenEnd::new(101_325.0));
    assert!(
        negative_peak < -20.0,
        "open end should return a negative pressure pulse, got {negative_peak}"
    );
    assert!(
        negative_peak.abs() > positive_peak,
        "open-end reflected pulse should be dominated by inverted pressure"
    );
}
