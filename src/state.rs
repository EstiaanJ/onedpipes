use crate::gas_properties::GasProperties;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct State {
    pub rho: f64,
    pub momentum: f64,
    pub rho_total_energy: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Primitive {
    pub rho: f64,
    pub u: f64,
    pub p: f64,
    pub temperature: f64,
    pub sound_speed: f64,
}

impl State {
    pub fn from_primitive<G: GasProperties>(rho: f64, u: f64, p: f64, gas: G) -> Self {
        let temperature = p / (rho * gas.r());
        let internal_energy = gas.internal_energy_from_temperature(temperature);
        Self {
            rho,
            momentum: rho * u,
            rho_total_energy: rho * (internal_energy + 0.5 * u * u),
        }
    }

    pub fn primitive<G: GasProperties>(&self, gas: G) -> Primitive {
        let rho = self.rho.max(1.0e-12);
        let u = self.momentum / rho;
        let internal_energy = self.specific_internal_energy().max(1.0);
        let temperature = gas.temperature_from_internal_energy(internal_energy);
        let p = (rho * gas.r() * temperature).max(1.0);
        let sound_speed = (gas.gamma(temperature) * gas.r() * temperature).sqrt();
        Primitive {
            rho,
            u,
            p,
            temperature,
            sound_speed,
        }
    }

    pub fn specific_internal_energy(&self) -> f64 {
        let rho = self.rho.max(1.0e-12);
        let u = self.momentum / rho;
        self.rho_total_energy / rho - 0.5 * u * u
    }

    pub fn flux<G: GasProperties>(&self, gas: G) -> Self {
        let prim = self.primitive(gas);
        Self {
            rho: self.momentum,
            momentum: self.momentum * prim.u + prim.p,
            rho_total_energy: prim.u * (self.rho_total_energy + prim.p),
        }
    }

    pub fn add_scaled(self, other: Self, scale: f64) -> Self {
        Self {
            rho: self.rho + scale * other.rho,
            momentum: self.momentum + scale * other.momentum,
            rho_total_energy: self.rho_total_energy + scale * other.rho_total_energy,
        }
    }

    pub fn scale(self, scale: f64) -> Self {
        Self {
            rho: self.rho * scale,
            momentum: self.momentum * scale,
            rho_total_energy: self.rho_total_energy * scale,
        }
    }

    pub fn plus(self, other: Self) -> Self {
        Self {
            rho: self.rho + other.rho,
            momentum: self.momentum + other.momentum,
            rho_total_energy: self.rho_total_energy + other.rho_total_energy,
        }
    }

    pub fn minus(self, other: Self) -> Self {
        Self {
            rho: self.rho - other.rho,
            momentum: self.momentum - other.momentum,
            rho_total_energy: self.rho_total_energy - other.rho_total_energy,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::State;
    use crate::gas_properties::TemperatureDependentAir;

    #[test]
    fn primitive_round_trip_preserves_inputs() {
        let gas = TemperatureDependentAir::new();
        let state = State::from_primitive(1.18, 12.0, 101_325.0, gas);
        let prim = state.primitive(gas);
        assert!((prim.rho - 1.18).abs() < 1.0e-12);
        assert!((prim.u - 12.0).abs() < 1.0e-12);
        assert!((prim.p - 101_325.0).abs() < 1.0e-8);
    }
}
