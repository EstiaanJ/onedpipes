#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeciesFractions {
    pub oxygen: f64,
    pub fuel_vapor: f64,
    pub inert: f64,
    pub products: f64,
}

impl SpeciesFractions {
    pub const AIR: Self = Self {
        oxygen: 0.21,
        fuel_vapor: 0.0,
        inert: 0.79,
        products: 0.0,
    };

    pub const EXHAUST: Self = Self {
        oxygen: 0.02,
        fuel_vapor: 0.0,
        inert: 0.73,
        products: 0.25,
    };

    pub fn new(oxygen: f64, fuel_vapor: f64, inert: f64, products: f64) -> Self {
        Self {
            oxygen,
            fuel_vapor,
            inert,
            products,
        }
        .normalized()
    }

    pub fn normalized(self) -> Self {
        let oxygen = finite_nonnegative(self.oxygen);
        let fuel_vapor = finite_nonnegative(self.fuel_vapor);
        let inert = finite_nonnegative(self.inert);
        let products = finite_nonnegative(self.products);
        let sum = oxygen + fuel_vapor + inert + products;
        if sum <= 0.0 {
            return Self::AIR;
        }
        Self {
            oxygen: oxygen / sum,
            fuel_vapor: fuel_vapor / sum,
            inert: inert / sum,
            products: products / sum,
        }
    }

    pub fn scale(self, scale: f64) -> SpeciesMass {
        SpeciesMass {
            oxygen: self.oxygen * scale,
            fuel_vapor: self.fuel_vapor * scale,
            inert: self.inert * scale,
            products: self.products * scale,
        }
    }

    pub fn lambda(self, stoich_fuel_oxygen_ratio: f64) -> Option<f64> {
        if self.fuel_vapor <= 0.0 || stoich_fuel_oxygen_ratio <= 0.0 {
            return None;
        }
        Some(self.oxygen / (stoich_fuel_oxygen_ratio * self.fuel_vapor))
    }
}

impl Default for SpeciesFractions {
    fn default() -> Self {
        Self::AIR
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpeciesMass {
    pub oxygen: f64,
    pub fuel_vapor: f64,
    pub inert: f64,
    pub products: f64,
}

impl SpeciesMass {
    pub fn from_density(density: f64, fractions: SpeciesFractions) -> Self {
        fractions.scale(density)
    }

    pub fn fractions(self) -> SpeciesFractions {
        SpeciesFractions::new(self.oxygen, self.fuel_vapor, self.inert, self.products)
    }

    pub fn add_scaled(self, other: Self, scale: f64) -> Self {
        Self {
            oxygen: self.oxygen + scale * other.oxygen,
            fuel_vapor: self.fuel_vapor + scale * other.fuel_vapor,
            inert: self.inert + scale * other.inert,
            products: self.products + scale * other.products,
        }
    }

    pub fn scale(self, scale: f64) -> Self {
        Self {
            oxygen: self.oxygen * scale,
            fuel_vapor: self.fuel_vapor * scale,
            inert: self.inert * scale,
            products: self.products * scale,
        }
    }
}

fn finite_nonnegative(value: f64) -> f64 {
    if value.is_finite() {
        value.max(0.0)
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::SpeciesFractions;

    #[test]
    fn species_fractions_normalize_inputs() {
        let species = SpeciesFractions::new(2.0, 1.0, -5.0, 1.0);

        assert!((species.oxygen - 0.5).abs() < 1.0e-12);
        assert!((species.fuel_vapor - 0.25).abs() < 1.0e-12);
        assert_eq!(species.inert, 0.0);
        assert!((species.products - 0.25).abs() < 1.0e-12);
    }

    #[test]
    fn lambda_uses_available_oxygen_and_fuel() {
        let species = SpeciesFractions::new(0.20, 0.05, 0.70, 0.05);

        assert!((species.lambda(3.4).unwrap() - 0.20 / (3.4 * 0.05)).abs() < 1.0e-12);
        assert_eq!(SpeciesFractions::AIR.lambda(3.4), None);
    }
}
