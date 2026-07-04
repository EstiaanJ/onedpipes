use onedpipes::{ClosedEnd, Duct, DuctConfig, State, TemperatureDependentAir};

fn first_positive_peak_time(samples: &[(f64, f64)], min_time: f64) -> Option<f64> {
    samples.windows(3).find_map(|window| {
        let (t, p) = window[1];
        if t > min_time && p > window[0].1 && p >= window[2].1 {
            Some(t)
        } else {
            None
        }
    })
}

#[test]
fn closed_closed_pipe_matches_fundamental_resonance_frequency() {
    let gas = TemperatureDependentAir::new();
    let length = 1.0;
    let cells = 96;
    let base_rho = 1.2;
    let base_p = 101_325.0;
    let amplitude = 1.0e-3;
    let config = DuctConfig {
        artificial_viscosity: 0.01,
        ..DuctConfig::new(length, cells, 1.0)
    };
    let mut duct = Duct::from_initializer(gas, config, ClosedEnd, ClosedEnd, |x| {
        let mode = (std::f64::consts::PI * x / length).cos();
        State::from_primitive(base_rho, 0.0, base_p * (1.0 + amplitude * mode), gas)
    });

    let sound_speed = State::from_primitive(base_rho, 0.0, base_p, gas)
        .primitive(gas)
        .sound_speed;
    let expected_frequency = sound_speed / (2.0 * length);
    let expected_period = 1.0 / expected_frequency;
    let end_time = 3.0 * expected_period;
    let probe_cell = 0;
    let mut time = 0.0;
    let mut samples = Vec::new();

    while time < end_time {
        let signal_speed = duct.max_signal_speed();
        let mut dt = 0.9 * 0.55 * duct.config().dx() / signal_speed;
        if time + dt > end_time {
            dt = end_time - time;
        }
        let report = duct.step(dt);
        assert_eq!(report.clipped_cells, 0);
        time += dt;
        samples.push((time, duct.cells()[probe_cell].primitive(gas).p));
    }

    let first_peak = first_positive_peak_time(&samples, 0.25 * expected_period).unwrap();
    let second_peak =
        first_positive_peak_time(&samples, first_peak + 0.5 * expected_period).unwrap();
    let measured_frequency = 1.0 / (second_peak - first_peak);
    let relative_error = ((measured_frequency - expected_frequency) / expected_frequency).abs();

    assert!(
        relative_error < 0.03,
        "measured={measured_frequency:.4}, expected={expected_frequency:.4}, relative_error={relative_error:.4}"
    );
}
