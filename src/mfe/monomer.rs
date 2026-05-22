use crate::{
    errors::Result,
    matrix::Matrix,
    mfe::{calc::ThermoCalc, common},
    seq::is_base_pair,
    thermo::{Thermo, params::ThermoParams},
};

const MIN_HAIRPIN_LOOP: usize = 3;
const MAX_HAIRPIN_LOOP: usize = 30;

/// Used to calculate thermos for monomers
#[derive(Debug)]
pub struct Monomer<'a> {
    mat: Matrix<Thermo>,
    traceback: Matrix<(i8, i8)>,
    seq: Vec<u8>,
    seq_r: Vec<u8>,
    thermo_calc: ThermoCalc<'a>,
    params: &'a ThermoParams,
}

impl<'a> Monomer<'a> {
    pub fn new(seq: &[u8], params: &'a ThermoParams) -> Result<Self> {
        let mut seq = seq.to_vec();
        common::pad_seq(&mut seq);
        let mut seq_r = seq.clone();
        seq_r.reverse();

        let mat = Self::init_matrix(&seq, &seq_r)?;

        Ok(Self {
            mat,
            traceback: Self::init_traceback(seq.len()),
            seq,
            seq_r,
            thermo_calc: ThermoCalc::new(params),
            params,
        })
    }

    fn hairpin(&self, i: usize, j: usize, traceback: bool) -> Result<Thermo> {
        let mut i = i;
        let mut j = j;
        let loop_size = j - i - 1;
        let mut thermo = Thermo::with_inf();

        if loop_size < 3 {
            return Ok(thermo);
        }

        let seq_len = self.seq.len();
        if i <= seq_len && seq_len < j {
            return Ok(thermo);
        } else if i > seq_len {
            i -= seq_len;
            j -= seq_len;
        }

        thermo = Thermo::new();
        if loop_size <= 30 {
            thermo.add_finite(self.params.get_hairpin(loop_size - 1));
        } else {
            thermo.add_finite(self.params.get_hairpin(29));
        }

        if loop_size > 3 {
            thermo.add_finite(
                self.params
                    .get_tstack(&self.seq[i..=i + 1], &[self.seq[j], self.seq[j - 1]])?,
            );
        } else if loop_size == 3 {
            thermo += ThermoParams::at_penalty(&self.seq[i], &self.seq[j]);
            thermo.add_finite(self.params.get_triloop(&self.seq[i..i + 5])?);
        } else if loop_size == 4 {
            thermo.add_finite(self.params.get_tetraloop(&self.seq[i..i + 6])?);
        }

        let thermo_current = self.mat.get(i, j)?;
        if thermo.dh > 0.0
            && thermo.ds > 0.0
            && (thermo_current.dh <= 0.0 || thermo_current.ds <= 0.0)
        {
            thermo = Thermo::with_inf();
        }

        if thermo_current.dg() < thermo.dg() && !traceback {
            thermo = thermo_current;
        }

        Ok(thermo)
    }

    fn init_matrix(seq1: &[u8], seq2: &[u8]) -> Result<Matrix<Thermo>> {
        let mut mat = Matrix::new(seq1.len(), seq2.len());
        for i in 1..seq1.len() {
            for j in 1..seq2.len() {
                if (j as i8 - i as i8) < (MIN_HAIRPIN_LOOP as i8 + 1)
                    || !is_base_pair(&seq1[i], &seq2[j])
                {
                    mat.set(i, j, Thermo::with_inf())?;
                }
            }
        }

        Ok(mat)
    }

    fn init_traceback(len: usize) -> Matrix<(i8, i8)> {
        Matrix::with_value(len, len, (-1, -1))
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    fn get_params() -> ThermoParams {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("thermo");
        let path = path.to_str().expect("path error").to_owned();

        let params = ThermoParams::with_file_path(&path);
        assert!(params.is_ok());
        params.unwrap()
    }

    #[test]
    fn test_hairpin() {
        let params = get_params();
        let seq = "CCCCCATCCGATCAGGGGG".as_bytes().to_vec();
        let monomer = Monomer::new(&seq, &params);
        assert!(monomer.is_ok());
        let mut monomer = monomer.unwrap();

        macro_rules! run_test {
            ($i:literal, $j:literal, $traceback:literal, $ds_i:literal, $dh_i:literal, $ds:literal, $dh:literal) => {
                assert!(monomer.mat.set($i, $j, Thermo::with_values($ds_i, $dh_i)).is_ok());
                assert_eq!(
                    monomer.hairpin($i, $j, $traceback),
                    Ok(Thermo::with_values($ds, $dh)),
                );
            };
        }

        run_test!(5, 10, false, -33324.0, 0.0, -17.18, -2600.0);
        run_test!(4, 10, false, -33324.0, 0.0, -15.74, -2100.0);
        run_test!(7, 11, false, -33324.0, 0.0, -4.379999999999999, 2200.0);
        run_test!(5, 16, false, -36.12999999998, -6000.0, -24.43, -3900.0);
        run_test!(5, 15, true, -25.2, -4300.0, -25.2, -4300.0);
    }
}
