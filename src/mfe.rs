use crate::{
    errors::Result,
    matrix::Matrix,
    seq::is_base_pair,
    thermo::{Thermo, lowest_dg, params::ThermoParams},
};

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
    pub fn calculate(&mut self) -> Result<()> {
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
    fn get_base(seq: &[u8], pos: usize) -> Result<&u8> {
        seq.get(pos).ok_or("out of bounds")
    }

    /// For V, check to see if S_i and S_j can base pair. If not, then it's set to infinity
    fn fill_v(v: &mut Matrix, row: usize, col: usize, seq1: &[u8], seq2: &[u8]) -> Result<()> {
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
    fn fill_w(w: &mut Matrix, row: usize, col: usize, seq1: &[u8], seq2: &[u8]) -> Result<()> {
        Ok(())
    }

    /// Decide on the next move based on the matrices
    fn next_move(v: &Matrix, w: &Matrix, row: usize, col: usize) -> (usize, usize) {
        (row + 1, col + 1)
    }
}

/// Given 2 pairs of bases, calculate the lowest end thermo
///
/// bases1 is always oriented 5' -> 3' and bases2 is always oriented 3' -> 5'
///
/// The best choices can be that they stack:
///
/// 5' -> 3'
/// b1_0 b1_1
/// b2_0 b2_1
/// 3' -> 5'
///
/// This is only possible if b1_1 and b2_1 are WC pairs
///
/// Another choice is that they're dangling
///
/// b1_0
///     \
///      b1_1
///      b2_1
///     /
/// b2_0
fn calc_lowest_end_thermo(
    bases1: [u8; 2],
    bases2: [u8; 2],
    thermo_params: &ThermoParams,
) -> Result<Thermo> {
    let b1_0 = bases1[0];
    let b1_1 = bases1[1];
    let b2_0 = bases2[0];
    let b2_1 = bases2[1];

    let thermo_inf = Thermo::with_inf();
    let base = thermo_params.at_penalty(&b1_0, &b2_0);

    if is_base_pair(&b1_1, &b2_1) {
        let mut stacked = base;
        if let Some(ts) = thermo_params.get_tstack(&bases1, &bases2)? {
            stacked += ts;
        }
        return Ok(lowest_dg(thermo_inf, stacked));
    }

    let mut dangling = base;
    if let Some(d1) = thermo_params.get_dangle(&bases1, &[b2_1])? {
        dangling += d1;
    }

    if let Some(d2) = thermo_params.get_dangle(&[b1_1], &bases2)? {
        dangling += d2;
    }

    Ok(lowest_dg(thermo_inf, dangling))
}
