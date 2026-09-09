use crate::{
    errors::Result,
    matrix::Matrix,
    mfe::{calc::ThermoCalc, common},
    seq::{is_base_pair_num, seq_to_num},
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
    seq_num: Vec<u8>,
    thermo_calc: ThermoCalc<'a>,
    params: &'a ThermoParams,
    max_loop: usize,
}

impl<'a> Monomer<'a> {
    pub fn new(seq: &[u8], params: &'a ThermoParams) -> Result<Self> {
        let mut seq = seq.to_vec();
        common::pad_seq(&mut seq);
        let seq_num = seq_to_num(&seq);

        let mat = Self::init_matrix(&seq_num)?;

        Ok(Self {
            mat,
            traceback: Self::init_traceback(seq.len()),
            seq,
            seq_num,
            thermo_calc: ThermoCalc::new(params),
            params,
            max_loop: 30,
        })
    }

    pub fn calculate(&mut self) -> Result<Thermo> {
        self.fill()?;
        self.find_best()
    }

    fn find_best(&self) -> Result<Thermo> {
        let oligo_len = self.seq.len() - 2;
        let mut thermo = Thermo::with_inf();

        for i in 1..=oligo_len {
            for j in 1..=oligo_len {
                if is_base_pair_num(&self.seq_num[i], &self.seq_num[j]) {
                    let current = self.mat.get(i, j)?;
                    if current.dg() < thermo.dg() {
                        thermo = current;
                    }
                }
            }
        }

        Ok(thermo)
    }

    fn fill(&mut self) -> Result<()> {
        let oligo_len = self.seq.len() - 2;

        for j in 2..=oligo_len {
            if j < MIN_HAIRPIN_LOOP + 2 {
                continue;
            }
            for i in (1..=(j - MIN_HAIRPIN_LOOP - 1)).rev() {
                if !is_base_pair_num(&self.seq_num[i], &self.seq_num[j]) {
                    continue;
                }

                self.fill_stack(i, j)?;
                self.fill_loop(i, j)?;
                self.fill_hairpin(i, j)?;
            }
        }

        Ok(())
    }

    fn fill_stack(&mut self, i: usize, j: usize) -> Result<()> {
        let mut stack = self.mat.get(i + 1, j - 1)?;
        stack.add_finite(
            self.params
                .get_stack(&self.seq[i..=i + 1], &[self.seq[j], self.seq[j - 1]])?,
        );

        if stack.dg() < self.mat.get(i, j)?.dg() {
            self.mat.set(i, j, stack)?;
            self.traceback.set(i, j, (i as i8 + 1, j as i8 - 1))?;
        }

        Ok(())
    }

    fn fill_loop(&mut self, i: usize, j: usize) -> Result<()> {
        let oligo_len = self.seq.len() - 2;
        let mut d = j as i64 - i as i64 - 3;
        let lower_bound =
            (MIN_HAIRPIN_LOOP as i64 + 1).max(j as i64 - i as i64 - 2 - self.max_loop as i64);

        while d >= lower_bound {
            let d_usize = d as usize;
            let mut ii = i + 1;
            while (ii as i64) < (j as i64 - d) && ii <= oligo_len {
                let jj = ii + d_usize;
                if is_base_pair_num(&self.seq_num[ii], &self.seq_num[jj]) {
                    let internal =
                        self.mat.get(ii, jj)? + self.optimal_internal_monomer(i, j, ii, jj)?;

                    if internal.dg() < self.mat.get(i, j)?.dg() {
                        if self.mat.get(i, j)? != internal {
                            self.traceback.set(i, j, (ii as i8, jj as i8))?;
                        }
                        self.mat.set(i, j, internal)?;
                    }
                }
                ii += 1;
            }
            d -= 1;
        }

        Ok(())
    }

    fn optimal_internal_monomer(&self, i: usize, j: usize, ii: usize, jj: usize) -> Result<Thermo> {
        let loop_size_1 = ii as i64 - i as i64 - 1;
        let loop_size_2 = j as i64 - jj as i64 - 1;
        let loop_size = loop_size_1 + loop_size_2 - 1;

        let mut thermo = Thermo::with_inf();

        if loop_size_1 == 0 || loop_size_2 == 0 {
            thermo = Thermo::new();
            thermo.add_finite(self.params.get_bulge(loop_size as usize));

            if loop_size_1 == 1 || loop_size_2 == 1 {
                thermo.add_finite(
                    self.params
                        .get_stack(&[self.seq[i], self.seq[ii]], &[self.seq[j], self.seq[jj]])?,
                );
            } else {
                thermo += ThermoParams::at_penalty(&self.seq[i], &self.seq[j])
                    + ThermoParams::at_penalty(&self.seq[ii], &self.seq[jj]);
            }
        } else if loop_size_1 == 1 && loop_size_2 == 1 {
            thermo = Thermo::new();
            thermo.add_finite(
                self.params
                    .get_stack_mm(&self.seq[i..=i + 1], &[self.seq[j], self.seq[j - 1]])?,
            );
            thermo.add_finite(self.params.get_stack_mm(
                &[self.seq[jj], self.seq[jj + 1]],
                &[self.seq[ii], self.seq[ii - 1]],
            )?);
        } else if !is_base_pair_num(&self.seq_num[ii - 1], &self.seq_num[jj + 1])
            && !is_base_pair_num(&self.seq_num[i + 1], &self.seq_num[j - 1])
        {
            thermo = Thermo::new();
            thermo.add_finite(self.params.get_internal(loop_size as usize));
            thermo.add_finite(
                self.params
                    .get_tstack(&self.seq[i..=i + 1], &[self.seq[j], self.seq[j - 1]])?,
            );
            thermo.add_finite(self.params.get_tstack(
                &[self.seq[jj], self.seq[jj + 1]],
                &[self.seq[ii], self.seq[ii - 1]],
            )?);

            let asym =
                ThermoParams::internal_loop((loop_size_1 - loop_size_2).unsigned_abs() as usize);
            if f64::is_finite(asym.dh) {
                thermo += asym;
            }
        }

        Ok(thermo)
    }

    fn fill_hairpin(&mut self, i: usize, j: usize) -> Result<()> {
        let thermo = self.hairpin(i, j, false)?;
        self.mat.set(i, j, thermo)?;
        Ok(())
    }

    fn init_matrix(seq_num: &[u8]) -> Result<Matrix<Thermo>> {
        let mut mat = Matrix::new(seq_num.len(), seq_num.len());
        for i in 1..seq_num.len() {
            for j in 1..seq_num.len() {
                if (j as i8 - i as i8) < (MIN_HAIRPIN_LOOP as i8 + 1)
                    || !is_base_pair_num(&seq_num[i], &seq_num[j])
                {
                    mat.set(i, j, Thermo::with_inf())?;
                } else {
                    mat.set(i, j, Thermo::init_base_pairs())?;
                }
            }
        }

        Ok(mat)
    }

    fn init_traceback(len: usize) -> Matrix<(i8, i8)> {
        Matrix::with_value(len, len, (-1, -1))
    }

    fn hairpin(&self, i: usize, j: usize, traceback: bool) -> Result<Thermo> {
        let mut i = i;
        let mut j = j;
        let loop_size = j - i - 1;
        let mut thermo = Thermo::with_inf();

        if loop_size < MIN_HAIRPIN_LOOP {
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
        if loop_size <= MAX_HAIRPIN_LOOP {
            thermo.add_finite(self.params.get_hairpin(loop_size - 1));
        } else {
            thermo.add_finite(self.params.get_hairpin(29));
        }

        if loop_size > MIN_HAIRPIN_LOOP {
            thermo.add_finite(
                self.params
                    .get_tstack(&self.seq[i..=i + 1], &[self.seq[j], self.seq[j - 1]])?,
            );
            if loop_size == 4 {
                thermo.add_finite(self.params.get_tetraloop(&self.seq[i..i + 6])?);
            }
        } else if loop_size == MIN_HAIRPIN_LOOP {
            thermo += ThermoParams::at_penalty(&self.seq[i], &self.seq[j]);
            thermo.add_finite(self.params.get_triloop(&self.seq[i..i + 5])?);
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hairpin() {
        let params = ThermoParams::with_defaults().unwrap();
        let seq = "CCCCCATCCGATCAGGGGG".as_bytes().to_vec();
        let monomer = Monomer::new(&seq, &params);
        assert!(monomer.is_ok());
        let mut monomer = monomer.unwrap();

        macro_rules! run_test {
            ($i:literal, $j:literal, $traceback:literal, $ds_i:literal, $dh_i:literal, $ds:literal, $dh:literal) => {
                assert!(
                    monomer
                        .mat
                        .set($i, $j, Thermo::with_values($ds_i, $dh_i))
                        .is_ok()
                );
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
