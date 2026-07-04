use onedpipes::{DuctEnd, JunctionPort, MultiPipeJunction, State, TemperatureDependentAir};

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
