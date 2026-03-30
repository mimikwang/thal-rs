use crate::{matrix::Matrix, seq::is_base_pair, thermo::{Thermo, params::ThermoParams}};

/// Minimum free energy calculator based on zuker's algorithm
///
/// Found in this paper: https://pmc.ncbi.nlm.nih.gov/articles/instance/326673/pdf/nar00394-0137.pdf
#[derive(Debug)]
pub struct Mfe {
    /// W contains the minimum free energy of all possible admissible structures formed from the
    /// subsequence S_ij.
    w: Matrix,
    /// V contains the minimum free energy of all possible admissible structures formed from S_ij
    /// in which S_i and S_j base pair with each other. If S_i and S_j cannot base pair, then V_ij
    /// is infinity.
    v: Matrix,
    /// Current row position - i.e. i in the paper
    current_row: usize,
    /// Current column position - i.e. j in the paper
    current_col: usize,
    /// The first sequence in 5' to 3' - corresponding to the i index
    seq1: Vec<u8>,
    seq1_length: usize,
    /// The second sequence in 3' to 5' - corresponding to the j index
    seq2: Vec<u8>,
    seq2_length: usize,

    /// Finished
    is_done: bool,

    // The thermo parameters
    thermo_params: ThermoParams,
}

impl Mfe {
    /// Construct a new MFE
    pub fn new(seq1: Vec<u8>, seq2: Vec<u8>, thermo_params: ThermoParams) -> Self {
        let seq1_length = seq1.len();
        let seq2_length = seq2.len();

        Self {
            w: Matrix::new(seq1_length, seq2_length),
            v: Matrix::new(seq1_length, seq2_length),
            current_row: 0,
            current_col: 0,
            seq1,
            seq1_length,
            seq2,
            seq2_length,
            is_done: false,
            thermo_params,
        }
    }

    /// Run through algorithm
    pub fn calculate(&mut self) -> Result<(), &'static str> {
        if self.is_done {
            return Err("already calculated");
        }

        let mut i = 0;
        let mut j = 0;
        loop {
            Self::fill_v(&mut self.v, i, j, &self.seq1, &self.seq2)?;
            Self::fill_w(&mut self.w, i, j, &self.seq1, &self.seq2)?;
            (i, j) = Self::next_move(&self.v, &self.w, i, j);

            if i >= self.seq1_length || j >= self.seq2_length {
                break;
            }
        }
        self.is_done = true;
        Ok(())
    }

    /// Return the minimum free energy
    pub fn mfe(&self) -> Option<Thermo> {
        if !self.is_done {
            return None;
        }

        Some(Thermo::new())
    }

    /// Retrieve the base for a sequence at a specific position
    fn get_base(seq: &[u8], pos: usize) -> Result<&u8, &'static str> {
        seq.get(pos).ok_or("out of bounds")
    }

    /// For V, check to see if S_i and S_j can base pair. If not, then it's set to infinity
    fn fill_v(
        v: &mut Matrix,
        row: usize,
        col: usize,
        seq1: &[u8],
        seq2: &[u8],
    ) -> Result<(), &'static str> {
        let b1 = Self::get_base(seq1, row)?;
        let b2 = Self::get_base(seq2, col)?;

        if is_base_pair(b1, b2) {
            // TODO: Add actual thermo for base pairs
            v.set(row, col, Thermo::new())?;
            return Ok(());
        }

        v.set(row, col, Thermo::with_inf())?;
        Ok(())
    }

    /// Fill w with values - either the hairpin, stacking region, buldge loop, or interior loop.
    /// We'll largely ignore bifurcation loops.
    fn fill_w(
        w: &mut Matrix,
        row: usize,
        col: usize,
        seq1: &[u8],
        seq2: &[u8],
    ) -> Result<(), &'static str> {
        Ok(())
    }

    /// Decide on the next move based on the matrices
    fn next_move(v: &Matrix, w: &Matrix, row: usize, col: usize) -> (usize, usize) {
        (row + 1, col + 1)
    }
}
