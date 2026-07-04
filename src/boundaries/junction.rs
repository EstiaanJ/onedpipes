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
        let mut boundary_states = Vec::with_capacity(ports.len());
        let mut port_flows = Vec::with_capacity(ports.len());

        for port in ports {
            let prim = port.state.primitive(gas);
            let interior_outward_velocity = outward_velocity(prim.u, port.end);
            let boundary_outward_velocity =
                interior_outward_velocity + (prim.p - pressure) / (prim.rho * prim.sound_speed);
            let boundary_velocity = duct_velocity(boundary_outward_velocity, port.end);
            let boundary_state = State::from_primitive(prim.rho, boundary_velocity, pressure, gas);
            let mass_flow_out = port.area * prim.rho * boundary_outward_velocity;
            let total_enthalpy = (boundary_state.rho_total_energy + pressure) / prim.rho;
            let energy_flow_out = mass_flow_out * total_enthalpy;

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
}
