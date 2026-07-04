use crate::{
    boundaries::BoundaryCondition,
    duct::{Duct, StepReport},
    gas_properties::GasProperties,
    timestep::global_timestep,
};

#[derive(Clone, Debug)]
pub struct Model<G, L, R>
where
    G: GasProperties,
    L: BoundaryCondition<G>,
    R: BoundaryCondition<G>,
{
    ducts: Vec<Duct<G, L, R>>,
    cfl: f64,
    time: f64,
}

impl<G, L, R> Model<G, L, R>
where
    G: GasProperties,
    L: BoundaryCondition<G>,
    R: BoundaryCondition<G>,
{
    pub fn new(cfl: f64) -> Self {
        assert!(cfl > 0.0);
        Self {
            ducts: Vec::new(),
            cfl,
            time: 0.0,
        }
    }

    pub fn add_duct(&mut self, duct: Duct<G, L, R>) {
        self.ducts.push(duct);
    }

    pub fn ducts(&self) -> &[Duct<G, L, R>] {
        &self.ducts
    }

    pub fn ducts_mut(&mut self) -> &mut [Duct<G, L, R>] {
        &mut self.ducts
    }

    pub fn time(&self) -> f64 {
        self.time
    }

    pub fn step(&mut self) -> StepReport {
        let dt = global_timestep(&self.ducts, self.cfl);
        self.step_with_dt(dt)
    }

    pub fn step_with_dt(&mut self, dt: f64) -> StepReport {
        let mut total = StepReport::default();
        for duct in &mut self.ducts {
            let report = duct.step(dt);
            total.clipped_cells += report.clipped_cells;
            total.fallback_faces += report.fallback_faces;
        }
        self.time += dt;
        total
    }

    pub fn run_until(&mut self, end_time: f64) -> StepReport {
        let mut total = StepReport::default();
        while self.time < end_time {
            let mut dt = global_timestep(&self.ducts, self.cfl);
            if self.time + dt > end_time {
                dt = end_time - self.time;
            }
            let report = self.step_with_dt(dt);
            total.clipped_cells += report.clipped_cells;
            total.fallback_faces += report.fallback_faces;
        }
        total
    }
}
