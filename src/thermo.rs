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

    /// Construct with values
    pub fn with_values(ds: f64, dh: f64) -> Self {
        Self { ds, dh }
    }

    /// Calculate the gibbs free energy
    pub fn dg(&self, temp: f64) -> f64 {
        self.dh - temp * self.ds
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
