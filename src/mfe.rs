use crate::{matrix::Matrix, thermo::Thermo};

/// Minimum free energy calculator based on zuker's algorithm
///
/// Found in this paper: https://pmc.ncbi.nlm.nih.gov/articles/instance/326673/pdf/nar00394-0137.pdf
#[derive(Debug)]
pub struct MFE {
    /// W contains the minimum free energy of all possible admissible structures formed from the
    /// subsequence S_ij.
    W: Matrix,
    /// V contains the minimum free energy of all possible admissible structures formed from S_ij
    /// in which S_i and S_j base pair with each other. If S_i and S_j cannot base pair, then V_ij
    /// is infinity.
    V: Matrix,
    /// Current row position - i.e. i in the paper
    current_row: usize,
    /// Current column position - i.e. j in the paper
    current_col: usize,
    /// The first sequence in 5' to 3' - corresponding to the i index
    seq1: Vec<u8>,
    seq1_length: usize,
    /// The second sequence in 5' to 3' - corresponding to the j index
    seq2: Vec<u8>,
    seq2_length: usize,

    /// Finished
    is_done: bool,
}

impl MFE {
    /// Construct a new MFE
    pub fn new(seq1: Vec<u8>, seq2: Vec<u8>) -> Self {
        let seq1_length = seq1.len();
        let seq2_length = seq2.len();

        Self {
            W: Matrix::new(seq1_length, seq2_length),
            V: Matrix::new(seq1_length, seq2_length),
            current_row: 0,
            current_col: 0,
            seq1,
            seq1_length,
            seq2,
            seq2_length,
            is_done: false,
        }
    }

    /// Run through algorithm
    pub fn calculate(&mut self) -> Result<(), &str> {
        if self.is_done {
            return Err("already calculated");
        }
        // TODO: Calculation here
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
}
