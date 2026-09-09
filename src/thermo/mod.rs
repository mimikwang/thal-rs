mod loader;
pub mod params;

const MIN_ENTROPY: f64 = -3224.0;

/// Holds enthalpy and entropy data
#[derive(Default, Debug, Copy, Clone, PartialEq)]
pub struct Thermo {
    /// Entropy
    pub ds: f64,
    /// Enthalpy
    pub dh: f64,
}

impl Thermo {
    /// Default constructor
    pub fn new() -> Self {
        Self::default()
    }

    /// Infinite
    pub fn with_inf() -> Self {
        Self {
            ds: 0.0,
            dh: f64::INFINITY,
        }
    }

    /// Initial values for base pairs
    pub fn init_base_pairs() -> Self {
        Self {
            ds: MIN_ENTROPY,
            dh: 0.0,
        }
    }

    /// Initial values for non base pairs
    pub fn init_not_base_pairs() -> Self {
        Self {
            ds: -1.0,
            dh: f64::INFINITY,
        }
    }

    /// Duplex
    pub fn init_duplex() -> Self {
        Self {
            ds: -5.7,
            dh: 200.0,
        }
    }

    /// Construct with values
    pub fn with_values(ds: f64, dh: f64) -> Self {
        Self { ds, dh }
    }

    /// Calculate the gibbs free energy at room temperature
    pub fn dg_with_temp(&self, temp: f64) -> f64 {
        self.dh - temp * self.ds
    }

    /// Calculate the gibbs free energy
    pub fn dg(&self) -> f64 {
        self.dg_with_temp(310.15)
    }

    /// Add if the rhs has finite enthalpy
    pub fn add_finite(&mut self, rhs: Option<Self>) {
        if let Some(t) = rhs
            && f64::is_finite(t.dh)
        {
            *self += t;
        }
    }
}

impl std::ops::Add for Thermo {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            ds: self.ds + rhs.ds,
            dh: self.dh + rhs.dh,
        }
    }
}

impl std::ops::AddAssign for Thermo {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl std::ops::Sub for Thermo {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            ds: self.ds - rhs.ds,
            dh: self.dh - rhs.dh,
        }
    }
}

impl std::ops::SubAssign for Thermo {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

// Pick the thermo with the lowest dg
pub fn lowest_dg(t1: Thermo, t2: Thermo) -> Thermo {
    if t1.dg() < t2.dg() {
        return t1;
    }
    t2
}
