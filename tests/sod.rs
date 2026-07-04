use onedpipes::{ClosedEnd, Duct, DuctConfig, GasProperties, State};

#[derive(Clone, Copy, Debug)]
struct ConstantGammaGas {
    gamma: f64,
    r: f64,
}

impl ConstantGammaGas {
    fn new(gamma: f64, r: f64) -> Self {
        Self { gamma, r }
    }
}

impl GasProperties for ConstantGammaGas {
    fn r(&self) -> f64 {
        self.r
    }

    fn cp(&self, _temperature: f64) -> f64 {
        self.gamma * self.r / (self.gamma - 1.0)
    }

    fn internal_energy_from_temperature(&self, temperature: f64) -> f64 {
        self.r * temperature / (self.gamma - 1.0)
    }

    fn temperature_from_internal_energy(&self, internal_energy: f64) -> f64 {
        internal_energy * (self.gamma - 1.0) / self.r
    }
}

#[derive(Clone, Copy, Debug)]
struct ExactState {
    rho: f64,
    u: f64,
    p: f64,
}

impl ExactState {
    fn sound_speed(self, gamma: f64) -> f64 {
        (gamma * self.p / self.rho).sqrt()
    }
}

#[derive(Clone, Copy, Debug)]
struct SodExact {
    gamma: f64,
    x0: f64,
    left: ExactState,
    right: ExactState,
    p_star: f64,
    u_star: f64,
}

impl SodExact {
    fn new(gamma: f64, x0: f64, left: ExactState, right: ExactState) -> Self {
        let p_star = solve_star_pressure(gamma, left, right);
        let u_star = 0.5
            * (left.u + right.u + pressure_function(gamma, p_star, right).0
                - pressure_function(gamma, p_star, left).0);
        Self {
            gamma,
            x0,
            left,
            right,
            p_star,
            u_star,
        }
    }

    fn sample(self, x: f64, t: f64) -> ExactState {
        let s = (x - self.x0) / t;
        if s <= self.u_star {
            self.sample_left(s)
        } else {
            self.sample_right(s)
        }
    }

    fn wave_positions(self, t: f64) -> SodWavePositions {
        let a_left = self.left.sound_speed(self.gamma);
        let a_right = self.right.sound_speed(self.gamma);
        let a_star_left =
            a_left * (self.p_star / self.left.p).powf((self.gamma - 1.0) / (2.0 * self.gamma));
        let shock_speed = self.right.u
            + a_right
                * ((self.gamma + 1.0) / (2.0 * self.gamma) * self.p_star / self.right.p
                    + (self.gamma - 1.0) / (2.0 * self.gamma))
                    .sqrt();

        SodWavePositions {
            rarefaction_head: self.x0 + (self.left.u - a_left) * t,
            rarefaction_tail: self.x0 + (self.u_star - a_star_left) * t,
            contact: self.x0 + self.u_star * t,
            shock: self.x0 + shock_speed * t,
        }
    }

    fn sample_left(self, s: f64) -> ExactState {
        let gamma = self.gamma;
        let left = self.left;
        let a_left = left.sound_speed(gamma);
        let a_star = a_left * (self.p_star / left.p).powf((gamma - 1.0) / (2.0 * gamma));
        let head_speed = left.u - a_left;
        let tail_speed = self.u_star - a_star;

        if s <= head_speed {
            left
        } else if s >= tail_speed {
            ExactState {
                rho: left.rho * (self.p_star / left.p).powf(1.0 / gamma),
                u: self.u_star,
                p: self.p_star,
            }
        } else {
            let u = 2.0 / (gamma + 1.0) * (a_left + 0.5 * (gamma - 1.0) * left.u + s);
            let a = 2.0 / (gamma + 1.0) * (a_left + 0.5 * (gamma - 1.0) * (left.u - s));
            ExactState {
                rho: left.rho * (a / a_left).powf(2.0 / (gamma - 1.0)),
                u,
                p: left.p * (a / a_left).powf(2.0 * gamma / (gamma - 1.0)),
            }
        }
    }

    fn sample_right(self, s: f64) -> ExactState {
        let gamma = self.gamma;
        let right = self.right;
        let a_right = right.sound_speed(gamma);
        let shock_speed = right.u
            + a_right
                * ((gamma + 1.0) / (2.0 * gamma) * self.p_star / right.p
                    + (gamma - 1.0) / (2.0 * gamma))
                    .sqrt();

        if s >= shock_speed {
            right
        } else {
            ExactState {
                rho: right.rho
                    * ((self.p_star / right.p + (gamma - 1.0) / (gamma + 1.0))
                        / ((gamma - 1.0) / (gamma + 1.0) * self.p_star / right.p + 1.0)),
                u: self.u_star,
                p: self.p_star,
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct SodWavePositions {
    rarefaction_head: f64,
    rarefaction_tail: f64,
    contact: f64,
    shock: f64,
}

fn solve_star_pressure(gamma: f64, left: ExactState, right: ExactState) -> f64 {
    let a_left = left.sound_speed(gamma);
    let a_right = right.sound_speed(gamma);
    let mut p = (0.5 * (left.p + right.p)
        - 0.125 * (right.u - left.u) * (left.rho + right.rho) * (a_left + a_right))
        .max(1.0e-8);

    for _ in 0..30 {
        let (f_left, df_left) = pressure_function(gamma, p, left);
        let (f_right, df_right) = pressure_function(gamma, p, right);
        let correction = (f_left + f_right + right.u - left.u) / (df_left + df_right);
        p = (p - correction).max(1.0e-8);
        if correction.abs() / p < 1.0e-12 {
            break;
        }
    }
    p
}

fn pressure_function(gamma: f64, p: f64, state: ExactState) -> (f64, f64) {
    if p > state.p {
        let a = 2.0 / ((gamma + 1.0) * state.rho);
        let b = (gamma - 1.0) / (gamma + 1.0) * state.p;
        let root = (a / (p + b)).sqrt();
        (
            (p - state.p) * root,
            root * (1.0 - 0.5 * (p - state.p) / (p + b)),
        )
    } else {
        let sound_speed = state.sound_speed(gamma);
        let exponent = (gamma - 1.0) / (2.0 * gamma);
        (
            2.0 * sound_speed / (gamma - 1.0) * ((p / state.p).powf(exponent) - 1.0),
            (1.0 / (state.rho * sound_speed)) * (p / state.p).powf(-(gamma + 1.0) / (2.0 * gamma)),
        )
    }
}

#[test]
fn sod_shock_tube_matches_exact_wave_structure_and_speeds() {
    let gamma = 1.4;
    let gas = ConstantGammaGas::new(gamma, 1.0);
    let left = ExactState {
        rho: 1.0,
        u: 0.0,
        p: 1.0,
    };
    let right = ExactState {
        rho: 0.125,
        u: 0.0,
        p: 0.1,
    };
    let exact = SodExact::new(gamma, 0.5, left, right);
    let length = 1.0;
    let cells = 500;
    let end_time = 0.20;
    let config = DuctConfig {
        artificial_viscosity: 0.015,
        pressure_floor: 1.0e-8,
        ..DuctConfig::new(length, cells, 1.0)
    };
    let mut duct = Duct::from_initializer(gas, config, ClosedEnd, ClosedEnd, |x| {
        let state = if x < exact.x0 { left } else { right };
        State::from_primitive(state.rho, state.u, state.p, gas)
    });
    let mut time = 0.0;
    let mut total_clipped_cells = 0;
    let mut total_fallback_faces = 0;

    while time < end_time {
        let mut dt = 0.9 * 0.45 * duct.config().dx() / duct.max_signal_speed();
        if time + dt > end_time {
            dt = end_time - time;
        }
        let report = duct.step(dt);
        total_clipped_cells += report.clipped_cells;
        total_fallback_faces += report.fallback_faces;
        time += dt;
    }

    assert_eq!(total_clipped_cells, 0);

    let dx = duct.config().dx();
    let mut rho_l1 = 0.0;
    let mut velocity_l1 = 0.0;
    let mut pressure_l1 = 0.0;
    let mut exact_rho_l1 = 0.0;
    let mut exact_pressure_l1 = 0.0;
    for (i, state) in duct.cells().iter().enumerate() {
        let x = (i as f64 + 0.5) * dx;
        let numerical = state.primitive(gas);
        let reference = exact.sample(x, end_time);
        rho_l1 += (numerical.rho - reference.rho).abs();
        velocity_l1 += (numerical.u - reference.u).abs();
        pressure_l1 += (numerical.p - reference.p).abs();
        exact_rho_l1 += reference.rho.abs();
        exact_pressure_l1 += reference.p.abs();
    }

    let cells = duct.cells().len() as f64;
    assert!(
        rho_l1 / exact_rho_l1 < 0.30,
        "relative density L1 error = {}",
        rho_l1 / exact_rho_l1
    );
    assert!(
        pressure_l1 / exact_pressure_l1 < 0.35,
        "relative pressure L1 error = {}",
        pressure_l1 / exact_pressure_l1
    );
    assert!(
        velocity_l1 / cells < 0.20,
        "mean absolute velocity error = {}, fallback_faces={total_fallback_faces}",
        velocity_l1 / cells
    );

    let wave_positions = exact.wave_positions(end_time);
    let shock_index = duct
        .cells()
        .iter()
        .enumerate()
        .skip_while(|(i, _)| (*i as f64 + 0.5) * dx < wave_positions.contact)
        .find_map(|(i, state)| {
            let x = (i as f64 + 0.5) * dx;
            let p = state.primitive(gas).p;
            if p < 0.5 * (exact.p_star + right.p) {
                Some((i, x))
            } else {
                None
            }
        })
        .map(|(_, x)| x)
        .expect("shock transition should be present");

    assert!(
        (shock_index - wave_positions.shock).abs() < 0.05,
        "shock position measured={shock_index}, exact={}",
        wave_positions.shock
    );
    assert!(wave_positions.rarefaction_head < wave_positions.rarefaction_tail);
    assert!(wave_positions.rarefaction_tail < wave_positions.contact);
    assert!(wave_positions.contact < wave_positions.shock);
}
