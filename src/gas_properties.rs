pub trait GasProperties: Copy {
    fn r(&self) -> f64;
    fn cp(&self, temperature: f64) -> f64;

    fn cv(&self, temperature: f64) -> f64 {
        self.cp(temperature) - self.r()
    }

    fn gamma(&self, temperature: f64) -> f64 {
        self.cp(temperature) / self.cv(temperature)
    }

    fn internal_energy_from_temperature(&self, temperature: f64) -> f64;
    fn temperature_from_internal_energy(&self, internal_energy: f64) -> f64;
}

#[derive(Clone, Copy, Debug)]
pub struct TemperatureDependentAir {
    r: f64,
    cp_ref: f64,
    cp_slope: f64,
    t_ref: f64,
    min_temperature: f64,
    max_temperature: f64,
}

impl Default for TemperatureDependentAir {
    fn default() -> Self {
        Self {
            r: 287.05,
            cp_ref: 1005.0,
            cp_slope: 0.10,
            t_ref: 300.0,
            min_temperature: 100.0,
            max_temperature: 4000.0,
        }
    }
}

impl TemperatureDependentAir {
    pub fn new() -> Self {
        Self::default()
    }

    fn clamp_temperature(&self, temperature: f64) -> f64 {
        temperature.clamp(self.min_temperature, self.max_temperature)
    }
}

impl GasProperties for TemperatureDependentAir {
    fn r(&self) -> f64 {
        self.r
    }

    fn cp(&self, temperature: f64) -> f64 {
        let t = self.clamp_temperature(temperature);
        self.cp_ref + self.cp_slope * (t - self.t_ref)
    }

    fn internal_energy_from_temperature(&self, temperature: f64) -> f64 {
        let t = self.clamp_temperature(temperature);
        let cv_intercept = self.cp_ref - self.r - self.cp_slope * self.t_ref;
        cv_intercept * t + 0.5 * self.cp_slope * t * t
    }

    fn temperature_from_internal_energy(&self, internal_energy: f64) -> f64 {
        let cv_intercept = self.cp_ref - self.r - self.cp_slope * self.t_ref;
        if self.cp_slope.abs() < f64::EPSILON {
            return (internal_energy / cv_intercept)
                .clamp(self.min_temperature, self.max_temperature);
        }

        let discriminant = cv_intercept * cv_intercept + 2.0 * self.cp_slope * internal_energy;
        if !discriminant.is_finite() || discriminant <= 0.0 {
            return self.min_temperature;
        }

        let t = (-cv_intercept + discriminant.sqrt()) / self.cp_slope;
        self.clamp_temperature(t)
    }
}

#[cfg(test)]
mod tests {
    use super::{GasProperties, TemperatureDependentAir};

    #[test]
    fn temperature_energy_round_trip_is_stable() {
        let gas = TemperatureDependentAir::new();
        for temperature in [220.0, 300.0, 800.0, 1500.0] {
            let energy = gas.internal_energy_from_temperature(temperature);
            let recovered = gas.temperature_from_internal_energy(energy);
            assert!((recovered - temperature).abs() < 1.0e-9);
        }
    }

    #[test]
    fn gamma_decreases_as_cp_rises() {
        let gas = TemperatureDependentAir::new();
        assert!(gas.gamma(1200.0) < gas.gamma(300.0));
        assert!(gas.gamma(300.0) > 1.0);
    }
}
