use crate::{boundaries::DuctEnd, gas_properties::GasProperties, state::State};

#[derive(Clone, Copy, Debug)]
pub struct JunctionPort {
    pub state: State,
    pub end: DuctEnd,
    pub area: f64,
}

impl JunctionPort {
    pub fn new(state: State, end: DuctEnd, area: f64) -> Self {
        assert!(area > 0.0);
        Self { state, end, area }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortFlow {
    pub mass_flow_out: f64,
    pub energy_flow_out: f64,
}

#[derive(Clone, Debug)]
pub struct JunctionSolution {
    pub pressure: f64,
    pub boundary_states: Vec<State>,
    pub port_flows: Vec<PortFlow>,
}

impl JunctionSolution {
    pub fn mass_residual(&self) -> f64 {
        self.port_flows.iter().map(|flow| flow.mass_flow_out).sum()
    }

    pub fn energy_residual(&self) -> f64 {
        self.port_flows
            .iter()
            .map(|flow| flow.energy_flow_out)
            .sum()
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct MultiPipeJunction;

impl MultiPipeJunction {
    pub fn solve<G: GasProperties>(&self, ports: &[JunctionPort], gas: G) -> JunctionSolution {
        assert!(ports.len() >= 2);

        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for port in ports {
            let prim = port.state.primitive(gas);
            let outward_velocity = outward_velocity(prim.u, port.end);
            numerator +=
                port.area * prim.p / prim.sound_speed + port.area * prim.rho * outward_velocity;
            denominator += port.area / prim.sound_speed;
        }

        let pressure = numerator / denominator;
        let mut boundary_outward_velocities = Vec::with_capacity(ports.len());
        let mut mass_flows_out = Vec::with_capacity(ports.len());
        let mut inflow_energy = 0.0;
        let mut inflow_mass = 0.0;

        for port in ports {
            let prim = port.state.primitive(gas);
            let interior_outward_velocity = outward_velocity(prim.u, port.end);
            let boundary_outward_velocity =
                interior_outward_velocity + (prim.p - pressure) / (prim.rho * prim.sound_speed);
            let mass_flow_out = port.area * prim.rho * boundary_outward_velocity;
            let duct_velocity = duct_velocity(boundary_outward_velocity, port.end);
            let boundary_state = State::from_primitive(prim.rho, duct_velocity, pressure, gas);
            let total_enthalpy = total_enthalpy(boundary_state, pressure);

            if mass_flow_out > 0.0 {
                inflow_mass += mass_flow_out;
                inflow_energy += mass_flow_out * total_enthalpy;
            }
            boundary_outward_velocities.push(boundary_outward_velocity);
            mass_flows_out.push(mass_flow_out);
        }

        let mixed_total_enthalpy = if inflow_mass > 0.0 {
            inflow_energy / inflow_mass
        } else {
            0.0
        };
        let mut boundary_states = Vec::with_capacity(ports.len());
        let mut port_flows = Vec::with_capacity(ports.len());

        for ((port, boundary_outward_velocity), mass_flow_out) in ports
            .iter()
            .zip(boundary_outward_velocities)
            .zip(mass_flows_out)
        {
            let boundary_velocity = duct_velocity(boundary_outward_velocity, port.end);
            let boundary_state = if mass_flow_out < 0.0 && inflow_mass > 0.0 {
                state_from_pressure_velocity_total_enthalpy(
                    pressure,
                    boundary_velocity,
                    mixed_total_enthalpy,
                    gas,
                )
            } else {
                let prim = port.state.primitive(gas);
                State::from_primitive(prim.rho, boundary_velocity, pressure, gas)
            };
            let energy_flow_out = mass_flow_out * total_enthalpy(boundary_state, pressure);

            boundary_states.push(boundary_state);
            port_flows.push(PortFlow {
                mass_flow_out,
                energy_flow_out,
            });
        }

        JunctionSolution {
            pressure,
            boundary_states,
            port_flows,
        }
    }
}

fn total_enthalpy(state: State, pressure: f64) -> f64 {
    (state.rho_total_energy + pressure) / state.rho
}

fn state_from_pressure_velocity_total_enthalpy<G: GasProperties>(
    pressure: f64,
    velocity: f64,
    total_enthalpy: f64,
    gas: G,
) -> State {
    let static_enthalpy = (total_enthalpy - 0.5 * velocity * velocity).max(gas.r());
    let mut low = 1.0;
    let mut high = 10_000.0;
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        let h_mid = gas.internal_energy_from_temperature(mid) + gas.r() * mid;
        if h_mid < static_enthalpy {
            low = mid;
        } else {
            high = mid;
        }
    }
    let temperature = 0.5 * (low + high);
    let rho = pressure / (gas.r() * temperature);
    State::from_primitive(rho, velocity, pressure, gas)
}

fn outward_velocity(duct_velocity: f64, end: DuctEnd) -> f64 {
    match end {
        DuctEnd::Left => -duct_velocity,
        DuctEnd::Right => duct_velocity,
    }
}

fn duct_velocity(outward_velocity: f64, end: DuctEnd) -> f64 {
    match end {
        DuctEnd::Left => -outward_velocity,
        DuctEnd::Right => outward_velocity,
    }
}

#[cfg(test)]
mod tests {
    use super::{JunctionPort, MultiPipeJunction};
    use crate::{boundaries::DuctEnd, gas_properties::TemperatureDependentAir, state::State};

    #[test]
    fn junction_returns_one_shared_pressure_for_all_ports() {
        let gas = TemperatureDependentAir::new();
        let p = 101_325.0;
        let ports = [
            JunctionPort::new(State::from_primitive(1.2, 8.0, p, gas), DuctEnd::Right, 1.0),
            JunctionPort::new(State::from_primitive(1.2, 4.0, p, gas), DuctEnd::Left, 1.0),
        ];

        let solution = MultiPipeJunction.solve(&ports, gas);

        for state in &solution.boundary_states {
            let prim = state.primitive(gas);
            assert!((prim.p - solution.pressure).abs() < 1.0e-8);
        }
    }

    #[test]
    fn junction_balances_mass_flux() {
        let gas = TemperatureDependentAir::new();
        let p = 101_325.0;
        let ports = [
            JunctionPort::new(
                State::from_primitive(1.2, 18.0, p, gas),
                DuctEnd::Right,
                1.0,
            ),
            JunctionPort::new(State::from_primitive(1.2, 0.0, p, gas), DuctEnd::Left, 1.0),
            JunctionPort::new(State::from_primitive(1.2, 0.0, p, gas), DuctEnd::Left, 1.0),
        ];

        let solution = MultiPipeJunction.solve(&ports, gas);

        assert!(solution.mass_residual().abs() < 1.0e-10);
    }

    #[test]
    fn junction_balances_energy_for_asymmetric_total_enthalpy() {
        let gas = TemperatureDependentAir::new();
        let ports = [
            JunctionPort::new(
                State::from_primitive(1.2, 24.0, 140_000.0, gas),
                DuctEnd::Right,
                1.0,
            ),
            JunctionPort::new(
                State::from_primitive(0.9, 0.0, 90_000.0, gas),
                DuctEnd::Left,
                0.6,
            ),
            JunctionPort::new(
                State::from_primitive(1.4, 0.0, 95_000.0, gas),
                DuctEnd::Left,
                0.4,
            ),
        ];

        let solution = MultiPipeJunction.solve(&ports, gas);

        assert!(solution.mass_residual().abs() < 1.0e-10);
        assert!(solution.energy_residual().abs() < 1.0e-5);
    }
}
