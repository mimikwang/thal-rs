use crate::errors::Result;

/// A 2D matrix represented by a 1D vector
#[derive(Debug, PartialEq)]
pub struct Matrix<T> {
    values: Vec<T>,
    width: usize,
}

impl<T: Default + Clone + Copy> Matrix<T> {
    /// Constructor
    pub fn new(height: usize, width: usize) -> Self {
        Self {
            values: vec![T::default(); width * height],
            width,
        }
    }

    pub fn with_value(height: usize, width: usize, value: T) -> Self {
        Self {
            values: vec![value; width * height],
            width,
        }
    }

    /// Get element at row and col
    pub fn get(&self, row: usize, col: usize) -> Result<T> {
        let val = self.values.get(self.ind(col, row)).ok_or("out of bounds")?;
        Ok(*val)
    }

    /// Set a value
    pub fn set(&mut self, row: usize, col: usize, value: T) -> Result<()> {
        let ind = self.ind(col, row);
        let val = self.values.get_mut(ind).ok_or("out of bounds")?;
        *val = value;
        Ok(())
    }

    fn ind(&self, row: usize, col: usize) -> usize {
        col * self.width + row
    }
}
