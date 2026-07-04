use onedpipes::{ClosedEnd, Duct, DuctConfig, GasProperties, State, TemperatureDependentAir};

fn main() {
    let gas = TemperatureDependentAir::new();
    let ambient = State::from_primitive(1.2, 0.0, 101_325.0, gas);
    let duct = Duct::new(
        gas,
        DuctConfig::new(1.0, 80, 1.0),
        ambient,
        ClosedEnd,
        ClosedEnd,
    );
    let prim = duct.cells()[0].primitive(gas);
    println!(
        "onedpipes: closed-closed duct ready, cells={}, dx={:.5} m, a0={:.2} m/s, gamma0={:.4}",
        duct.config().cells,
        duct.config().dx(),
        prim.sound_speed,
        gas.gamma(prim.temperature)
    );
}
