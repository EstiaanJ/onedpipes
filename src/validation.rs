use crate::{
    boundaries::ClosedEnd,
    duct::Duct,
    gas_properties::{GasProperties, TemperatureDependentAir},
    solvers::SolverKind,
    state::State,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ScalarField {
    Density,
    Pressure,
    Temperature,
}

impl ScalarField {
    pub fn label(self) -> &'static str {
        match self {
            Self::Density => "Density",
            Self::Pressure => "Pressure",
            Self::Temperature => "Temperature",
        }
    }

    pub fn units(self) -> &'static str {
        match self {
            Self::Density => "kg/m^3",
            Self::Pressure => "Pa",
            Self::Temperature => "K",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Snapshot {
    pub time: f64,
    pub x: Vec<f64>,
    pub density: Vec<f64>,
    pub pressure: Vec<f64>,
    pub temperature: Vec<f64>,
    pub c_plus: Vec<f64>,
    pub c_minus: Vec<f64>,
}

impl Snapshot {
    pub fn values(&self, field: ScalarField) -> &[f64] {
        match field {
            ScalarField::Density => &self.density,
            ScalarField::Pressure => &self.pressure,
            ScalarField::Temperature => &self.temperature,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RunReport {
    pub clipped_cells: usize,
    pub fallback_faces: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct OrganPipeConfig {
    pub length: f64,
    pub cells: usize,
    pub area: f64,
    pub cfl: f64,
    pub base_density: f64,
    pub base_pressure: f64,
    pub perturbation_amplitude: f64,
    pub artificial_viscosity: f64,
    pub snapshot_interval: f64,
    pub max_history: usize,
    pub solver: SolverKind,
}

impl Default for OrganPipeConfig {
    fn default() -> Self {
        Self {
            length: 1.0,
            cells: 96,
            area: 1.0,
            cfl: 0.55,
            base_density: 1.2,
            base_pressure: 101_325.0,
            perturbation_amplitude: 1.0e-3,
            artificial_viscosity: 0.01,
            snapshot_interval: 2.0e-5,
            max_history: 360,
            solver: SolverKind::default(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct OrganPipeRun {
    gas: TemperatureDependentAir,
    duct: Duct<TemperatureDependentAir, ClosedEnd, ClosedEnd>,
    config: OrganPipeConfig,
    time: f64,
    next_snapshot_time: f64,
    history: Vec<Snapshot>,
    probe_pressure: Vec<(f64, f64)>,
    clipped_cells: usize,
    fallback_faces: usize,
    expected_frequency: f64,
}

impl OrganPipeRun {
    pub fn new(config: OrganPipeConfig) -> Self {
        let gas = TemperatureDependentAir::new();
        let duct_config = crate::duct::DuctConfig {
            artificial_viscosity: config.artificial_viscosity,
            solver: config.solver,
            ..crate::duct::DuctConfig::new(config.length, config.cells, config.area)
        };
        let duct = Duct::from_initializer(gas, duct_config, ClosedEnd, ClosedEnd, |x| {
            let mode = (std::f64::consts::PI * x / config.length).cos();
            State::from_primitive(
                config.base_density,
                0.0,
                config.base_pressure * (1.0 + config.perturbation_amplitude * mode),
                gas,
            )
        });
        let sound_speed =
            State::from_primitive(config.base_density, 0.0, config.base_pressure, gas)
                .primitive(gas)
                .sound_speed;
        let expected_frequency = sound_speed / (2.0 * config.length);
        let mut run = Self {
            gas,
            duct,
            config,
            time: 0.0,
            next_snapshot_time: 0.0,
            history: Vec::new(),
            probe_pressure: Vec::new(),
            clipped_cells: 0,
            fallback_faces: 0,
            expected_frequency,
        };
        run.record_snapshot();
        run
    }

    pub fn config(&self) -> OrganPipeConfig {
        self.config
    }

    pub fn time(&self) -> f64 {
        self.time
    }

    pub fn history(&self) -> &[Snapshot] {
        &self.history
    }

    pub fn latest_snapshot(&self) -> &Snapshot {
        self.history
            .last()
            .expect("organ-pipe run always records an initial snapshot")
    }

    pub fn probe_pressure(&self) -> &[(f64, f64)] {
        &self.probe_pressure
    }

    pub fn expected_frequency(&self) -> f64 {
        self.expected_frequency
    }

    pub fn measured_frequency(&self) -> Option<f64> {
        let expected_period = 1.0 / self.expected_frequency;
        let first = first_positive_peak_time(&self.probe_pressure, 0.25 * expected_period)?;
        let second = first_positive_peak_time(&self.probe_pressure, first + 0.5 * expected_period)?;
        Some(1.0 / (second - first))
    }

    pub fn report(&self) -> RunReport {
        RunReport {
            clipped_cells: self.clipped_cells,
            fallback_faces: self.fallback_faces,
        }
    }

    pub fn step(&mut self) -> RunReport {
        let mut dt = 0.9 * self.config.cfl * self.duct.config().dx() / self.duct.max_signal_speed();
        if self.time + dt > self.next_snapshot_time && self.next_snapshot_time > self.time {
            dt = self.next_snapshot_time - self.time;
        }
        let report = self.duct.step(dt);
        self.time += dt;
        self.clipped_cells += report.clipped_cells;
        self.fallback_faces += report.fallback_faces;
        self.probe_pressure
            .push((self.time, self.duct.cells()[0].primitive(self.gas).p));
        if self.time >= self.next_snapshot_time {
            self.record_snapshot();
            self.next_snapshot_time = self.time + self.config.snapshot_interval;
        }
        self.report()
    }

    fn record_snapshot(&mut self) {
        let dx = self.duct.config().dx();
        let mut snapshot = Snapshot {
            time: self.time,
            x: Vec::with_capacity(self.duct.cells().len()),
            density: Vec::with_capacity(self.duct.cells().len()),
            pressure: Vec::with_capacity(self.duct.cells().len()),
            temperature: Vec::with_capacity(self.duct.cells().len()),
            c_plus: Vec::with_capacity(self.duct.cells().len()),
            c_minus: Vec::with_capacity(self.duct.cells().len()),
        };
        for (i, state) in self.duct.cells().iter().enumerate() {
            let prim = state.primitive(self.gas);
            let gamma = self.gas.gamma(prim.temperature);
            let acoustic = 2.0 * prim.sound_speed / (gamma - 1.0);
            snapshot.x.push((i as f64 + 0.5) * dx);
            snapshot.density.push(prim.rho);
            snapshot.pressure.push(prim.p);
            snapshot.temperature.push(prim.temperature);
            snapshot.c_plus.push(prim.u + acoustic);
            snapshot.c_minus.push(prim.u - acoustic);
        }
        self.history.push(snapshot);
        if self.history.len() > self.config.max_history {
            let excess = self.history.len() - self.config.max_history;
            self.history.drain(0..excess);
        }
    }
}

fn first_positive_peak_time(samples: &[(f64, f64)], min_time: f64) -> Option<f64> {
    samples.windows(3).find_map(|window| {
        let (t, p) = window[1];
        if t > min_time && p > window[0].1 && p >= window[2].1 {
            Some(t)
        } else {
            None
        }
    })
}

#[cfg(test)]
mod tests {
    use super::{OrganPipeConfig, OrganPipeRun};
    use crate::solvers::SolverKind;

    fn assert_organ_pipe_viewer_run_records_snapshots_and_probe_data(solver: SolverKind) {
        let mut run = OrganPipeRun::new(OrganPipeConfig {
            cells: 32,
            snapshot_interval: 1.0e-4,
            max_history: 16,
            solver,
            ..OrganPipeConfig::default()
        });

        for _ in 0..20 {
            run.step();
        }

        assert!(!run.history().is_empty());
        assert!(!run.probe_pressure().is_empty());
        assert_eq!(run.latest_snapshot().x.len(), 32);
        assert_eq!(run.report().clipped_cells, 0);
    }

    #[test]
    fn lax_wendroff_organ_pipe_viewer_run_records_snapshots_and_probe_data() {
        assert_organ_pipe_viewer_run_records_snapshots_and_probe_data(SolverKind::LaxWendroff);
    }

    #[test]
    fn mac_cormack_organ_pipe_viewer_run_records_snapshots_and_probe_data() {
        assert_organ_pipe_viewer_run_records_snapshots_and_probe_data(SolverKind::MacCormack);
    }
}
