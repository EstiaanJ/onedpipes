use onedpipes::{
    ClosedEnd, Duct, DuctConfig, OrganPipeConfig, OrganPipeRun, SolverKind, State,
    TemperatureDependentAir,
};

const LENGTH: f64 = 1.0;
const AREA: f64 = 1.0;
const BASE_RHO: f64 = 1.2;
const BASE_P: f64 = 101_325.0;

#[derive(Clone, Copy, Debug)]
struct PulseCase {
    cells: usize,
    center: f64,
    width: f64,
    amplitude: f64,
    artificial_viscosity: f64,
    cfl: f64,
}

impl Default for PulseCase {
    fn default() -> Self {
        Self {
            cells: 240,
            center: 0.20,
            width: 0.035,
            amplitude: 120.0,
            artificial_viscosity: 0.003,
            cfl: 0.45,
        }
    }
}

fn base_sound_speed(gas: TemperatureDependentAir) -> f64 {
    State::from_primitive(BASE_RHO, 0.0, BASE_P, gas)
        .primitive(gas)
        .sound_speed
}

fn right_going_pulse_state(
    gas: TemperatureDependentAir,
    x: f64,
    center: f64,
    width: f64,
    amplitude: f64,
) -> State {
    let sound_speed = base_sound_speed(gas);
    let shape = (-((x - center) / width).powi(2)).exp();
    let dp = amplitude * shape;
    let rho = BASE_RHO + dp / (sound_speed * sound_speed);
    let u = dp / (BASE_RHO * sound_speed);
    State::from_primitive(rho, u, BASE_P + dp, gas)
}

fn pulse_duct(
    solver: SolverKind,
    case: PulseCase,
) -> Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd> {
    let gas = TemperatureDependentAir::new();
    let config = DuctConfig {
        artificial_viscosity: case.artificial_viscosity,
        solver,
        ..DuctConfig::new(LENGTH, case.cells, AREA)
    };
    Duct::from_initializer(gas, config, ClosedEnd, ClosedEnd, |x| {
        right_going_pulse_state(gas, x, case.center, case.width, case.amplitude)
    })
}

fn fixed_dt(case: PulseCase) -> f64 {
    let gas = TemperatureDependentAir::new();
    case.cfl * (LENGTH / case.cells as f64) / base_sound_speed(gas)
}

fn step_to(
    duct: &mut Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>,
    time: &mut f64,
    target_time: f64,
    dt: f64,
) {
    while *time < target_time {
        let step_dt = (target_time - *time).min(dt);
        let report = duct.step(step_dt);
        assert_eq!(report.clipped_cells, 0);
        *time += step_dt;
    }
}

fn pressure_centroid(
    duct: &Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>,
    gas: TemperatureDependentAir,
) -> f64 {
    let dx = duct.config().dx();
    let (weighted_x, weight) =
        duct.cells()
            .iter()
            .enumerate()
            .fold((0.0, 0.0), |(weighted_x, weight), (i, state)| {
                let x = (i as f64 + 0.5) * dx;
                let dp = (state.primitive(gas).p - BASE_P).max(0.0);
                (weighted_x + x * dp, weight + dp)
            });
    weighted_x / weight
}

fn probe_peak_time(solver: SolverKind, case: PulseCase, probe_x: f64) -> f64 {
    let gas = TemperatureDependentAir::new();
    let mut duct = pulse_duct(solver, case);
    let dt = fixed_dt(case);
    let probe_cell = ((probe_x / LENGTH) * case.cells as f64).floor() as usize;
    let end_time = (probe_x - case.center) / base_sound_speed(gas) + 5.0 * dt;
    let mut time = 0.0;
    let mut samples = Vec::new();

    while time < end_time {
        let target_time = (time + dt).min(end_time);
        step_to(&mut duct, &mut time, target_time, dt);
        let pressure = duct.cells()[probe_cell].primitive(gas).p;
        samples.push((time, pressure));
    }

    samples
        .windows(3)
        .filter_map(|window| {
            let (t, p) = window[1];
            (p > window[0].1 && p >= window[2].1).then_some((t, p))
        })
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(t, _)| t)
        .expect("right-going pulse should produce a probe peak before reflection")
}

fn total_mass(duct: &Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>) -> f64 {
    let dx = duct.config().dx();
    duct.cells().iter().map(|state| state.rho * dx * AREA).sum()
}

fn rms_acoustic_relation_error(
    duct: &Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>,
    gas: TemperatureDependentAir,
) -> f64 {
    let sound_speed = base_sound_speed(gas);
    let mut squared_error = 0.0;
    let mut count = 0;
    for state in duct.cells() {
        let prim = state.primitive(gas);
        let dp = prim.p - BASE_P;
        if dp > 1.0 {
            let expected_u = dp / (BASE_RHO * sound_speed);
            squared_error += ((prim.u - expected_u) / expected_u).powi(2);
            count += 1;
        }
    }
    (squared_error / count as f64).sqrt()
}

fn measured_organ_pipe_frequency(config: OrganPipeConfig) -> (f64, f64) {
    let mut run = OrganPipeRun::new(config);
    let end_time = 3.0 / run.expected_frequency();

    while run.time() < end_time {
        run.step();
    }

    assert_eq!(
        run.report().clipped_cells,
        0,
        "{} organ-pipe run should not hit positivity floors",
        config.solver.label()
    );
    let measured = run.measured_frequency().unwrap_or_else(|| {
        panic!(
            "{} organ-pipe run did not produce two peaks",
            config.solver.label()
        )
    });
    (measured, run.expected_frequency())
}

fn normalized_profile_differences(
    left: &Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>,
    right: &Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>,
    gas: TemperatureDependentAir,
) -> (f64, f64) {
    assert_eq!(left.cells().len(), right.cells().len());

    let mut pressure_error = 0.0;
    let mut pressure_scale = 0.0;
    let mut velocity_error = 0.0;
    let mut velocity_scale = 0.0;

    for (left_state, right_state) in left.cells().iter().zip(right.cells()) {
        let left_prim = left_state.primitive(gas);
        let right_prim = right_state.primitive(gas);
        let left_dp = left_prim.p - BASE_P;
        let right_dp = right_prim.p - BASE_P;

        pressure_error += (left_dp - right_dp).powi(2);
        pressure_scale += 0.5 * (left_dp.powi(2) + right_dp.powi(2));
        velocity_error += (left_prim.u - right_prim.u).powi(2);
        velocity_scale += 0.5 * (left_prim.u.powi(2) + right_prim.u.powi(2));
    }

    (
        (pressure_error / pressure_scale).sqrt(),
        (velocity_error / velocity_scale).sqrt(),
    )
}

#[test]
fn organ_pipe_frequency_matches_between_lax_wendroff_and_mac_cormack() {
    let (lax_wendroff, _) = measured_organ_pipe_frequency(OrganPipeConfig {
        solver: SolverKind::LaxWendroff,
        ..OrganPipeConfig::default()
    });
    let (mac_cormack, _) = measured_organ_pipe_frequency(OrganPipeConfig {
        solver: SolverKind::MacCormack,
        ..OrganPipeConfig::default()
    });
    let relative_difference =
        (lax_wendroff - mac_cormack).abs() / (0.5 * (lax_wendroff + mac_cormack));

    assert!(
        relative_difference < 0.02,
        "Lax-Wendroff={lax_wendroff:.4}, MacCormack={mac_cormack:.4}, relative_difference={relative_difference:.4}"
    );
}

#[test]
fn organ_pipe_frequency_sweep_matches_reference_for_both_solvers() {
    let cases = [
        OrganPipeConfig::default(),
        OrganPipeConfig {
            length: 0.75,
            cells: 128,
            base_density: 1.0,
            base_pressure: 90_000.0,
            perturbation_amplitude: 8.0e-4,
            cfl: 0.50,
            ..OrganPipeConfig::default()
        },
        OrganPipeConfig {
            length: 1.30,
            cells: 160,
            base_density: 1.5,
            base_pressure: 140_000.0,
            perturbation_amplitude: 5.0e-4,
            cfl: 0.45,
            ..OrganPipeConfig::default()
        },
    ];

    for base_config in cases {
        let (lw, expected) = measured_organ_pipe_frequency(OrganPipeConfig {
            solver: SolverKind::LaxWendroff,
            ..base_config
        });
        let (mc, _) = measured_organ_pipe_frequency(OrganPipeConfig {
            solver: SolverKind::MacCormack,
            ..base_config
        });
        let lw_reference_error = ((lw - expected) / expected).abs();
        let mc_reference_error = ((mc - expected) / expected).abs();
        let cross_solver_error = (lw - mc).abs() / (0.5 * (lw + mc));

        assert!(
            lw_reference_error < 0.04,
            "LW organ-pipe sweep case={base_config:?}: measured={lw:.4}, expected={expected:.4}, error={lw_reference_error:.4}"
        );
        assert!(
            mc_reference_error < 0.04,
            "MacCormack organ-pipe sweep case={base_config:?}: measured={mc:.4}, expected={expected:.4}, error={mc_reference_error:.4}"
        );
        assert!(
            cross_solver_error < 0.025,
            "organ-pipe sweep case={base_config:?}: LW={lw:.4}, MacCormack={mc:.4}, cross_error={cross_solver_error:.4}"
        );
    }
}

#[test]
fn right_going_pulse_centroid_tracks_same_position_over_time() {
    let gas = TemperatureDependentAir::new();
    let case = PulseCase::default();
    let sound_speed = base_sound_speed(gas);
    let targets = [0.00020, 0.00045, 0.00075];
    let dt = fixed_dt(case);
    let mut lw = pulse_duct(SolverKind::LaxWendroff, case);
    let mut mc = pulse_duct(SolverKind::MacCormack, case);
    let mut lw_time = 0.0;
    let mut mc_time = 0.0;

    for target in targets {
        step_to(&mut lw, &mut lw_time, target, dt);
        step_to(&mut mc, &mut mc_time, target, dt);
        let lw_x = pressure_centroid(&lw, gas);
        let mc_x = pressure_centroid(&mc, gas);
        let expected_x = case.center + sound_speed * target;

        assert!(
            (lw_x - expected_x).abs() < 0.025,
            "LW centroid target={target:.6}: measured={lw_x:.5}, expected={expected_x:.5}"
        );
        assert!(
            (mc_x - expected_x).abs() < 0.025,
            "MacCormack centroid target={target:.6}: measured={mc_x:.5}, expected={expected_x:.5}"
        );
        assert!(
            (lw_x - mc_x).abs() < 0.012,
            "centroids diverged at t={target:.6}: LW={lw_x:.5}, MacCormack={mc_x:.5}"
        );
    }
}

#[test]
fn pressure_and_velocity_profiles_stay_close_after_pulse_propagation() {
    let gas = TemperatureDependentAir::new();
    let case = PulseCase {
        cells: 300,
        amplitude: 150.0,
        width: 0.030,
        artificial_viscosity: 0.004,
        ..PulseCase::default()
    };
    let target_time = 0.00070;
    let dt = fixed_dt(case);
    let mut lw = pulse_duct(SolverKind::LaxWendroff, case);
    let mut mc = pulse_duct(SolverKind::MacCormack, case);
    let mut lw_time = 0.0;
    let mut mc_time = 0.0;

    step_to(&mut lw, &mut lw_time, target_time, dt);
    step_to(&mut mc, &mut mc_time, target_time, dt);

    let (pressure_difference, velocity_difference) = normalized_profile_differences(&lw, &mc, gas);

    assert!(
        pressure_difference < 0.20,
        "normalized pressure-profile difference = {pressure_difference}"
    );
    assert!(
        velocity_difference < 0.20,
        "normalized velocity-profile difference = {velocity_difference}"
    );
}

#[test]
fn right_going_pulse_probe_arrival_times_match() {
    let gas = TemperatureDependentAir::new();
    let case = PulseCase {
        cells: 320,
        ..PulseCase::default()
    };
    let probes = [0.42, 0.58, 0.72];

    for probe_x in probes {
        let expected = (probe_x - case.center) / base_sound_speed(gas);
        let lw = probe_peak_time(SolverKind::LaxWendroff, case, probe_x);
        let mc = probe_peak_time(SolverKind::MacCormack, case, probe_x);

        assert!(
            (lw - expected).abs() < 3.5e-5,
            "LW probe={probe_x:.2}: peak={lw:.7}, expected={expected:.7}"
        );
        assert!(
            (mc - expected).abs() < 3.5e-5,
            "MacCormack probe={probe_x:.2}: peak={mc:.7}, expected={expected:.7}"
        );
        assert!(
            (lw - mc).abs() < 2.0e-5,
            "probe={probe_x:.2}: LW={lw:.7}, MacCormack={mc:.7}"
        );
    }
}

#[test]
fn right_going_pulse_preserves_acoustic_pressure_velocity_relation() {
    let gas = TemperatureDependentAir::new();
    let case = PulseCase {
        amplitude: 80.0,
        ..PulseCase::default()
    };
    let target_time = 0.00050;
    let dt = fixed_dt(case);
    let mut lw = pulse_duct(SolverKind::LaxWendroff, case);
    let mut mc = pulse_duct(SolverKind::MacCormack, case);
    let mut lw_time = 0.0;
    let mut mc_time = 0.0;

    step_to(&mut lw, &mut lw_time, target_time, dt);
    step_to(&mut mc, &mut mc_time, target_time, dt);

    let lw_error = rms_acoustic_relation_error(&lw, gas);
    let mc_error = rms_acoustic_relation_error(&mc, gas);

    assert!(
        lw_error < 0.18,
        "LW acoustic relation RMS error = {lw_error}"
    );
    assert!(
        mc_error < 0.18,
        "MacCormack acoustic relation RMS error = {mc_error}"
    );
    assert!(
        (lw_error - mc_error).abs() < 0.08,
        "acoustic relation errors diverged: LW={lw_error}, MacCormack={mc_error}"
    );
}

#[test]
fn closed_pipe_mass_drift_matches_for_traveling_pulse() {
    let case = PulseCase {
        center: 0.35,
        amplitude: 180.0,
        artificial_viscosity: 0.0,
        ..PulseCase::default()
    };
    let dt = fixed_dt(case);
    let end_time = 0.0012;
    let mut lw = pulse_duct(SolverKind::LaxWendroff, case);
    let mut mc = pulse_duct(SolverKind::MacCormack, case);
    let lw_mass_initial = total_mass(&lw);
    let mc_mass_initial = total_mass(&mc);
    let mut lw_time = 0.0;
    let mut mc_time = 0.0;

    step_to(&mut lw, &mut lw_time, end_time, dt);
    step_to(&mut mc, &mut mc_time, end_time, dt);

    let lw_relative_drift = ((total_mass(&lw) - lw_mass_initial) / lw_mass_initial).abs();
    let mc_relative_drift = ((total_mass(&mc) - mc_mass_initial) / mc_mass_initial).abs();

    assert!(
        lw_relative_drift < 5.0e-9,
        "LW relative mass drift = {lw_relative_drift:.6e}"
    );
    assert!(
        mc_relative_drift < 5.0e-9,
        "MacCormack relative mass drift = {mc_relative_drift:.6e}"
    );
    assert!(
        (lw_relative_drift - mc_relative_drift).abs() < 5.0e-9,
        "mass drift mismatch: LW={lw_relative_drift:.6e}, MacCormack={mc_relative_drift:.6e}"
    );
}

#[test]
fn pulse_parameter_sweep_keeps_solver_positions_close() {
    let gas = TemperatureDependentAir::new();
    let cases = [
        PulseCase {
            cells: 160,
            width: 0.030,
            amplitude: 60.0,
            ..PulseCase::default()
        },
        PulseCase {
            cells: 220,
            width: 0.045,
            amplitude: 160.0,
            artificial_viscosity: 0.006,
            ..PulseCase::default()
        },
        PulseCase {
            cells: 320,
            center: 0.16,
            width: 0.025,
            amplitude: 240.0,
            cfl: 0.35,
            ..PulseCase::default()
        },
    ];

    for case in cases {
        let target_time = 0.00065;
        let dt = fixed_dt(case);
        let mut lw = pulse_duct(SolverKind::LaxWendroff, case);
        let mut mc = pulse_duct(SolverKind::MacCormack, case);
        let mut lw_time = 0.0;
        let mut mc_time = 0.0;

        step_to(&mut lw, &mut lw_time, target_time, dt);
        step_to(&mut mc, &mut mc_time, target_time, dt);
        let lw_x = pressure_centroid(&lw, gas);
        let mc_x = pressure_centroid(&mc, gas);

        assert!(
            (lw_x - mc_x).abs() < 0.018,
            "case={case:?}: LW centroid={lw_x:.5}, MacCormack centroid={mc_x:.5}"
        );
    }
}
