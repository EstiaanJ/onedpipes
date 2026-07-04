use crate::{gas_properties::GasProperties, state::State};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveOrifice {
    pub discharge_coefficient: f64,
    pub flow_area: f64,
}

impl ValveOrifice {
    pub fn new(discharge_coefficient: f64, flow_area: f64) -> Self {
        assert!(discharge_coefficient > 0.0);
        assert!(flow_area >= 0.0);
        Self {
            discharge_coefficient,
            flow_area,
        }
    }

    pub fn mass_flow<G: GasProperties>(
        &self,
        upstream: State,
        downstream: State,
        gas: G,
    ) -> OrificeFlow {
        let upstream_prim = upstream.primitive(gas);
        let downstream_prim = downstream.primitive(gas);
        let upstream_total = stagnation_state(upstream, gas);
        let downstream_total = stagnation_state(downstream, gas);

        if downstream_total.pressure > upstream_total.pressure {
            let reverse = self.mass_flow(downstream, upstream, gas);
            return OrificeFlow {
                mass_flow: -reverse.mass_flow,
                energy_flow: -reverse.energy_flow,
                pressure_ratio: reverse.pressure_ratio,
                critical_pressure_ratio: reverse.critical_pressure_ratio,
                choked: reverse.choked,
                upstream_stagnation_pressure: reverse.upstream_stagnation_pressure,
                upstream_stagnation_temperature: reverse.upstream_stagnation_temperature,
            };
        }

        let gamma = gas.gamma(upstream_prim.temperature);
        let critical_pressure_ratio = (2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0));
        if self.flow_area == 0.0 || downstream_prim.p >= upstream_total.pressure {
            return OrificeFlow {
                mass_flow: 0.0,
                energy_flow: 0.0,
                pressure_ratio: 1.0,
                critical_pressure_ratio,
                choked: false,
                upstream_stagnation_pressure: upstream_total.pressure,
                upstream_stagnation_temperature: upstream_total.temperature,
            };
        }

        let pressure_ratio = (downstream_prim.p / upstream_total.pressure).clamp(1.0e-12, 1.0);
        let choked = pressure_ratio <= critical_pressure_ratio;
        let flow_factor = if choked {
            (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        } else {
            let pressure_term =
                pressure_ratio.powf(2.0 / gamma) - pressure_ratio.powf((gamma + 1.0) / gamma);
            (2.0 / (gamma - 1.0) * pressure_term).max(0.0).sqrt()
        };
        let mass_flow = self.discharge_coefficient
            * self.flow_area
            * upstream_total.pressure
            * (gamma / (gas.r() * upstream_total.temperature)).sqrt()
            * flow_factor;
        let energy_flow = mass_flow * upstream_total.enthalpy;

        OrificeFlow {
            mass_flow,
            energy_flow,
            pressure_ratio,
            critical_pressure_ratio,
            choked,
            upstream_stagnation_pressure: upstream_total.pressure,
            upstream_stagnation_temperature: upstream_total.temperature,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrificeFlow {
    pub mass_flow: f64,
    pub energy_flow: f64,
    pub pressure_ratio: f64,
    pub critical_pressure_ratio: f64,
    pub choked: bool,
    pub upstream_stagnation_pressure: f64,
    pub upstream_stagnation_temperature: f64,
}

#[derive(Clone, Copy, Debug)]
struct StagnationState {
    pressure: f64,
    temperature: f64,
    enthalpy: f64,
}

fn stagnation_state<G: GasProperties>(state: State, gas: G) -> StagnationState {
    let prim = state.primitive(gas);
    let gamma = gas.gamma(prim.temperature);
    let mach = prim.u / prim.sound_speed;
    let temperature = prim.temperature * (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
    let pressure = prim.p * (temperature / prim.temperature).powf(gamma / (gamma - 1.0));
    let enthalpy = gas.internal_energy_from_temperature(prim.temperature)
        + gas.r() * prim.temperature
        + 0.5 * prim.u * prim.u;
    StagnationState {
        pressure,
        temperature,
        enthalpy,
    }
}

#[cfg(test)]
mod tests {
    use super::ValveOrifice;
    use crate::{State, gas_properties::GasProperties, gas_properties::TemperatureDependentAir};

    fn state_at(gas: TemperatureDependentAir, pressure: f64, temperature: f64) -> State {
        State::from_primitive(pressure / (gas.r() * temperature), 0.0, pressure, gas)
    }

    #[test]
    fn equal_stagnation_pressure_has_zero_mass_flow() {
        let gas = TemperatureDependentAir::new();
        let valve = ValveOrifice::new(0.8, 1.0e-4);
        let state = state_at(gas, 100_000.0, 300.0);

        let flow = valve.mass_flow(state, state, gas);

        assert_eq!(flow.mass_flow, 0.0);
        assert_eq!(flow.energy_flow, 0.0);
    }

    #[test]
    fn low_downstream_pressure_chokes_flow() {
        let gas = TemperatureDependentAir::new();
        let valve = ValveOrifice::new(0.8, 1.0e-4);
        let upstream = state_at(gas, 200_000.0, 300.0);
        let downstream = state_at(gas, 70_000.0, 300.0);

        let flow = valve.mass_flow(upstream, downstream, gas);

        assert!(flow.choked);
        assert!(flow.mass_flow > 0.0);
    }

    #[test]
    fn higher_right_pressure_reverses_flow_sign() {
        let gas = TemperatureDependentAir::new();
        let valve = ValveOrifice::new(0.8, 1.0e-4);
        let left = state_at(gas, 90_000.0, 300.0);
        let right = state_at(gas, 130_000.0, 300.0);

        let flow = valve.mass_flow(left, right, gas);

        assert!(flow.mass_flow < 0.0);
        assert!(flow.energy_flow < 0.0);
    }
}
