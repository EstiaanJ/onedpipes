use crate::{
    boundaries::{BoundaryCondition, DuctEnd},
    gas_properties::GasProperties,
    state::State,
};

#[derive(Clone, Copy, Debug)]
pub struct OpenEnd {
    ambient_pressure: f64,
}

impl OpenEnd {
    pub fn new(ambient_pressure: f64) -> Self {
        assert!(ambient_pressure > 0.0);
        Self { ambient_pressure }
    }

    pub fn ambient_pressure(&self) -> f64 {
        self.ambient_pressure
    }
}

impl<G: GasProperties> BoundaryCondition<G> for OpenEnd {
    fn ghost_state(&self, interior: State, _end: DuctEnd, gas: G) -> State {
        let prim = interior.primitive(gas);
        State::from_primitive(prim.rho, prim.u, self.ambient_pressure, gas)
    }
}

#[cfg(test)]
mod tests {
    use super::OpenEnd;
    use crate::{
        boundaries::{BoundaryCondition, DuctEnd},
        gas_properties::TemperatureDependentAir,
        state::State,
    };

    #[test]
    fn open_end_fixes_ambient_pressure_and_extrapolates_density_velocity() {
        let gas = TemperatureDependentAir::new();
        let boundary = OpenEnd::new(100_000.0);
        let interior = State::from_primitive(1.2, 15.0, 101_325.0, gas);
        let ghost = boundary.ghost_state(interior, DuctEnd::Right, gas);
        let interior_prim = interior.primitive(gas);
        let ghost_prim = ghost.primitive(gas);
        assert!((ghost_prim.rho - interior_prim.rho).abs() < 1.0e-12);
        assert!((ghost_prim.u - interior_prim.u).abs() < 1.0e-12);
        assert!((ghost_prim.p - boundary.ambient_pressure()).abs() < 1.0e-8);
    }
}
