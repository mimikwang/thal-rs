use crate::{
    errors::Result,
    seq::{N_NUM, is_base_pair_num},
    thermo::{Thermo, params::ThermoParams},
};

/// ThermoCalc calculates thermodynamics
#[derive(Debug)]
pub struct ThermoCalc<'a> {
    params: &'a ThermoParams,
}

impl<'a> ThermoCalc<'a> {
    pub fn new(params: &'a ThermoParams) -> Self {
        Self { params }
    }

    /// Get the optimal terminal thermo
    ///
    /// Compare terminal stacking and dangle thermos and pick the one that is most thermodynamically
    /// favorable.
    ///
    /// b10
    ///    \
    ///    b11
    ///    b21
    ///    /
    /// b20
    pub fn optimal_terminal(
        &self,
        b10: usize,
        b11: usize,
        b20: usize,
        b21: usize,
    ) -> Result<Thermo> {
        let mut thermo = self.get_terminal_stack(b10, b11, b20, b21)?;

        if !is_base_pair_num(&b11, &b21) {
            let dangle = self.get_dangle(b10, b11, b20, b21)?;
            if dangle.dg() < thermo.dg() {
                thermo = dangle;
            }
        }

        let at_penalty = ThermoParams::at_penalty_num(b10, b20);
        if at_penalty.dg() < thermo.dg() {
            thermo = at_penalty;
        }

        Ok(thermo)
    }

    /// Get the optimal interal by comparing bulge loops and internal loops
    ///
    /// A bulge where the loop size is 2.
    ///
    ///    bXX bXX
    /// b10       b11
    /// b20       b21
    ///
    /// An internal loop where loop size 1 and loop size 2 are 1.
    ///
    ///    bXX
    /// b10    b11
    /// b20    b21
    ///    bXX
    ///
    /// An internal loop is where loop size 1 or loop size 2 are greater than 1
    ///
    ///    bXX bXX
    /// b10       b11
    /// b20       b21
    ///      bXX
    pub fn optimal_internal(
        &self,
        i: i8,
        j: i8,
        ii: i8,
        jj: i8,
        seq1_num: &[usize],
        seq2_num: &[usize],
    ) -> Result<Thermo> {
        let loop_size_1 = ii - i - 1;
        let loop_size_2 = jj - j - 1;
        let loop_size = loop_size_1 + loop_size_2 - 1;

        if loop_size_1 == 0 || loop_size_2 == 0 {
            let mut thermo = Thermo::new();
            thermo.add_finite_option(self.params.get_bulge(loop_size as usize));

            if loop_size_1 == 1 || loop_size_2 == 1 {
                thermo.add_finite(self.params.get_stack(
                    seq1_num[i as usize],
                    seq1_num[ii as usize],
                    seq2_num[j as usize],
                    seq2_num[jj as usize],
                ));
            } else {
                thermo += ThermoParams::at_penalty_num(seq1_num[i as usize], seq2_num[j as usize])
                    + ThermoParams::at_penalty_num(seq1_num[ii as usize], seq2_num[jj as usize]);
            }
            return Ok(thermo);
        }

        if loop_size_1 == 1 && loop_size_2 == 1 {
            let mut thermo = Thermo::new();
            thermo.add_finite(self.params.get_stack_mm(
                seq1_num[i as usize],
                seq1_num[i as usize + 1],
                seq2_num[j as usize],
                seq2_num[j as usize + 1],
            ));
            thermo.add_finite(self.params.get_stack_mm(
                seq2_num[jj as usize],
                seq2_num[jj as usize - 1],
                seq1_num[ii as usize],
                seq1_num[ii as usize - 1],
            ));
            return Ok(thermo);
        }

        if !is_base_pair_num(&seq1_num[ii as usize - 1], &seq2_num[jj as usize - 1])
            && !is_base_pair_num(&seq1_num[i as usize + 1], &seq2_num[j as usize + 1])
        {
            let mut thermo = Thermo::new();
            thermo.add_finite_option(self.params.get_internal(loop_size as usize));
            thermo.add_finite(self.params.get_tstack(
                seq1_num[i as usize],
                seq1_num[i as usize + 1],
                seq2_num[j as usize],
                seq2_num[j as usize + 1],
            ));
            thermo.add_finite(self.params.get_tstack(
                seq2_num[jj as usize],
                seq2_num[jj as usize - 1],
                seq1_num[ii as usize],
                seq1_num[ii as usize - 1],
            ));

            let t =
                ThermoParams::internal_loop((loop_size_1 - loop_size_2).unsigned_abs() as usize);
            if f64::is_finite(t.dh) {
                thermo += t;
            }

            return Ok(thermo);
        }

        Ok(Thermo::with_inf())
    }

    /// Get the terminal stack thermo
    ///
    /// This is calculated by adding the AT penalty to the terminal stack look up values.
    fn get_terminal_stack(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Result<Thermo> {
        let mut thermo = ThermoParams::at_penalty_num(b10, b20);
        thermo.add_finite(self.params.get_tstack(b10, b11, b20, b21));
        Ok(thermo)
    }

    /// Get the dangle thermo
    ///
    /// The bases are arranged like so:
    ///
    /// b10
    ///    \
    ///    b11
    ///    b21
    ///    /
    /// b20
    ///
    /// where b10 b11 is oriented as 5' --> 3' and b20 b21 is oriented as 3' <-- 5'
    fn get_dangle(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Result<Thermo> {
        let mut thermo = ThermoParams::at_penalty_num(b10, b20);
        let dangle_3 = self.params.get_dangle3(b10, b11, b20);
        let dangle_5 = self.params.get_dangle5(b10, b20, b21);

        if b21 == N_NUM {
            thermo.add_finite(dangle_3);
        } else if b11 == N_NUM {
            thermo.add_finite(dangle_5);
        } else {
            thermo.add_finite(dangle_3);
            thermo.add_finite(dangle_5);
        }

        Ok(thermo)
    }
}

#[cfg(test)]
mod tests {
    use crate::seq::seq_to_num;

    use super::*;

    #[test]
    fn test_optimal_internal() {
        let params = ThermoParams::with_defaults().unwrap();
        let thermo_calc = ThermoCalc::new(&params);
        let seq1 = "NCCCCCATCCGATCAGGGGGN".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();
        let seq1 = seq_to_num(&seq1);
        let seq2 = seq_to_num(&seq2);

        macro_rules! run_test {
            ($i:literal, $j:literal, $ii:literal, $jj:literal, $ds:literal, $dh:literal) => {
                assert_eq!(
                    thermo_calc.optimal_internal($i, $j, $ii, $jj, &seq1, &seq2),
                    Ok(Thermo::with_values($ds, $dh)),
                );
            };

            ($i:literal, $j:literal, $ii:literal, $jj:literal) => {
                assert_eq!(
                    thermo_calc.optimal_internal($i, $j, $ii, $jj, &seq1, &seq2),
                    Ok(Thermo::with_inf()),
                );
            };
        }

        run_test!(1, 2, 2, 4, -32.79, -8000.0);
        run_test!(1, 2, 2, 10, -11.92, 0.0);
        run_test!(1, 5, 4, 10, -35.91454779945188, -7700.0);
        run_test!(16, 18, 19, 19, -9.35, 0.0);
        run_test!(18, 15, 19, 19, -9.99, 0.0);
        run_test!(15, 16, 19, 19);
        run_test!(18, 12, 19, 19, -11.28, 0.0);
        run_test!(3, 1, 4, 10, -12.57, 0.0);
        run_test!(5, 4, 7, 6, -4.6, -1400.0);
    }

    #[test]
    fn test_optimal_terminal() {
        let params = ThermoParams::with_defaults().unwrap();
        let thermo_calc = ThermoCalc::new(&params);

        macro_rules! run_test {
            ($b10:literal, $b11:literal, $b20:literal, $b21:literal, $ds:literal, $dh:literal) => {
                assert_eq!(
                    thermo_calc.optimal_terminal($b10, $b11, $b20, $b21),
                    Ok(Thermo::with_values($ds, $dh)),
                );
            };
        }

        run_test!(2, 0, 1, 1, -27.4, -9800.0);
        run_test!(2, 2, 1, 1, -19.3, -7000.0);
        run_test!(3, 2, 0, 1, -9.200000000000001, -3800.0);
        run_test!(0, 3, 3, 0, -6.699999999999999, -2800.0);
        run_test!(1, 4, 2, 2, -12.6, -4400.0);
        run_test!(3, 2, 0, 2, 25.1, 7200.0);
    }
}
