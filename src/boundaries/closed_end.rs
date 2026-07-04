use crate::{
    boundaries::{BoundaryCondition, DuctEnd},
    gas_properties::GasProperties,
    state::State,
};

#[derive(Clone, Copy, Debug, Default)]
pub struct ClosedEnd;

impl<G: GasProperties> BoundaryCondition<G> for ClosedEnd {
    fn ghost_state(&self, interior: State, _end: DuctEnd, gas: G) -> State {
        let prim = interior.primitive(gas);
        State::from_primitive(prim.rho, -prim.u, prim.p, gas)
    }
}

#[cfg(test)]
mod tests {
    use super::ClosedEnd;
    use crate::{
        boundaries::{BoundaryCondition, DuctEnd},
        gas_properties::TemperatureDependentAir,
        state::State,
    };

    #[test]
    fn closed_end_reflects_velocity_and_preserves_pressure_density() {
        let gas = TemperatureDependentAir::new();
        let interior = State::from_primitive(1.2, 15.0, 101_325.0, gas);
        let ghost = ClosedEnd.ghost_state(interior, DuctEnd::Left, gas);
        let interior_prim = interior.primitive(gas);
        let ghost_prim = ghost.primitive(gas);
        assert!((ghost_prim.rho - interior_prim.rho).abs() < 1.0e-12);
        assert!((ghost_prim.p - interior_prim.p).abs() < 1.0e-8);
        assert!((ghost_prim.u + interior_prim.u).abs() < 1.0e-12);
    }
}
