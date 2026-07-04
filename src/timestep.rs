use crate::{boundaries::BoundaryCondition, duct::Duct, gas_properties::GasProperties};

pub fn global_timestep<G, L, R>(ducts: &[Duct<G, L, R>], cfl: f64) -> f64
where
    G: GasProperties,
    L: BoundaryCondition<G>,
    R: BoundaryCondition<G>,
{
    assert!(!ducts.is_empty());
    assert!(cfl > 0.0);
    let min_dt = ducts
        .iter()
        .map(|duct| duct.config().dx() / duct.max_signal_speed())
        .fold(f64::INFINITY, f64::min);
    0.9 * cfl * min_dt
}
