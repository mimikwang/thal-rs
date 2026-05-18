use crate::errors::Result;
use crate::thermo::Thermo;

/// A 2D matrix represented by a 1D vector
#[derive(Debug, PartialEq)]
pub struct Matrix {
    values: Vec<Thermo>,
    width: usize,
}

impl Matrix {
    /// Constructor
    pub fn new(height: usize, width: usize) -> Self {
        Self {
            values: vec![Thermo::default(); width * height],
            width,
        }
    }

    /// Get element at row and col
    pub fn get(&self, row: usize, col: usize) -> Result<Thermo> {
        let thermo = self.values.get(self.ind(col, row)).ok_or("out of bounds")?;
        Ok(*thermo)
    }

    /// Set a value
    pub fn set(&mut self, row: usize, col: usize, value: Thermo) -> Result<()> {
        let ind = self.ind(col, row);
        let thermo = self.values.get_mut(ind).ok_or("out of bounds")?;
        *thermo = value;
        Ok(())
    }

    fn ind(&self, row: usize, col: usize) -> usize {
        col * self.width + row
    }
}
