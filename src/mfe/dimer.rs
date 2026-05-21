use crate::{
    errors::Result,
    matrix::Matrix,
    mfe::{calc::ThermoCalc, common},
    seq::is_base_pair,
    thermo::{Thermo, params::ThermoParams},
};

/// Used to calculate thermos for dimers
#[derive(Debug)]
pub struct Dimer<'a> {
    mat: Matrix<Thermo>,
    traceback: Matrix<(i8, i8)>,
    seq1: Vec<u8>,
    seq2: Vec<u8>,
    thermo_calc: ThermoCalc<'a>,
    params: &'a ThermoParams,
    max_loop: usize,
}

impl<'a> Dimer<'a> {
    pub fn new(seq1: &[u8], seq2: &[u8], params: &'a ThermoParams) -> Result<Self> {
        let mut seq1 = seq1.to_vec();
        let mut seq2 = seq2.to_vec();
        common::pad_seq(&mut seq1);
        common::pad_seq(&mut seq2);

        let mat = Self::init_matrix(&seq1, &seq2)?;

        Ok(Self {
            mat,
            traceback: Self::init_traceback(seq1.len(), seq2.len()),
            seq1,
            seq2,
            thermo_calc: ThermoCalc::new(params),
            params,
            max_loop: 30,
        })
    }

    pub fn calculate(&mut self) -> Result<Thermo> {
        self.fill()?;
        self.find_best()
    }

    fn init_matrix(seq1: &[u8], seq2: &[u8]) -> Result<Matrix<Thermo>> {
        let mut mat = Matrix::new(seq1.len(), seq2.len());

        // Fill based on base pairings
        for (i, b1) in seq1.iter().enumerate().take(seq1.len() - 1).skip(1) {
            for (j, b2) in seq2.iter().enumerate().take(seq2.len() - 1).skip(1) {
                if is_base_pair(b1, b2) {
                    mat.set(i, j, Thermo::init_base_pairs())?;
                } else {
                    mat.set(i, j, Thermo::init_not_base_pairs())?;
                }
            }
        }

        Ok(mat)
    }

    fn init_traceback(len1: usize, len2: usize) -> Matrix<(i8, i8)> {
        Matrix::with_value(len1, len2, (-1, -1))
    }

    fn find_best(&self) -> Result<Thermo> {
        let mut best = (1, 1);
        let mut thermo = Thermo::with_inf();

        for i in 1..self.seq1.len() {
            for j in 1..self.seq2.len() {
                if is_base_pair(&self.seq1[i], &self.seq2[j]) {
                    let current = self.mat.get(i, j)?
                        + Thermo::init_duplex()
                        + self.right_optimal_terminal(i, j)?;

                    if current.dg() < thermo.dg() {
                        best = (i, j);
                        thermo = current;
                    }
                }
            }
        }

        if is_base_pair(&self.seq1[best.0], &self.seq2[best.1]) {
            thermo = self.mat.get(best.0, best.1)? + self.right_optimal_terminal(best.0, best.1)?;
        }

        Ok(thermo)
    }

    fn fill(&mut self) -> Result<()> {
        for i in 1..self.seq1.len() {
            for j in 1..self.seq2.len() {
                if !is_base_pair(&self.seq1[i], &self.seq2[j]) {
                    self.mat.set(i, j, Thermo::with_inf())?;
                    continue;
                }

                self.mat.set(i, j, self.left_optimal_terminal(i, j)?)?;
                if i > 1 && j > 1 {
                    self.fill_stack(i, j)?;
                    self.fill_loop(i, j)?;
                }
            }
        }

        Ok(())
    }

    fn fill_stack(&mut self, i: usize, j: usize) -> Result<()> {
        if is_base_pair(&self.seq1[i - 1], &self.seq2[j - 1]) {
            let mut stack = self.mat.get(i - 1, j - 1)?;
            stack.add_finite(
                self.params
                    .get_stack(&self.seq1[i - 1..=i], &self.seq2[j - 1..=j])?,
            );
            self.mat.set(i, j, stack)?;
            self.traceback.set(i, j, (i as i8 - 1, j as i8 - 1))?;
        }
        Ok(())
    }

    fn fill_loop(&mut self, i: usize, j: usize) -> Result<()> {
        let mut best = self.mat.get(i, j)? + self.right_optimal_terminal(i, j)?;
        for d in 3..self.max_loop + 3 {
            let (mut ii, mut jj) = Self::get_loop_indices(i, j, d);
            while ii > 0 && jj < j as i8 {
                if is_base_pair(&self.seq1[ii as usize], &self.seq2[jj as usize]) {
                    let internal = self.mat.get(ii as usize, jj as usize)?
                        + self
                            .thermo_calc
                            .optimal_internal(ii, jj, i as i8, j as i8, &self.seq1, &self.seq2)?;

                    if internal.dg() < best.dg() {
                        best = internal;
                        if self.mat.get(i, j)? != internal {
                            self.traceback.set(i, j, (ii, jj))?;
                        }
                        self.mat.set(i, j, internal)?;
                    }
                }
                ii -= 1;
                jj += 1;
            }
        }

        Ok(())
    }

    fn left_optimal_terminal(&self, i: usize, j: usize) -> Result<Thermo> {
        self.thermo_calc.optimal_terminal(
            self.seq2[j],
            self.seq2[j - 1],
            self.seq1[i],
            self.seq1[i - 1],
        )
    }

    fn right_optimal_terminal(&self, i: usize, j: usize) -> Result<Thermo> {
        self.thermo_calc.optimal_terminal(
            self.seq1[i],
            self.seq1[i + 1],
            self.seq2[j],
            self.seq2[j + 1],
        )
    }

    fn get_loop_indices(i: usize, j: usize, d: usize) -> (i8, i8) {
        let mut ii = i as i8 - 1;
        let mut jj = j as i8 + i as i8 - ii - d as i8;
        if jj < 1 {
            ii += jj - 1;
            jj = 1;
        }
        (ii, jj)
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
    fn test_new() {
        let seq1 = &[b'A', b'T', b'C', b'G'];
        let seq2 = &[b'A', b'A'];
        let params = get_params();

        let dimer = Dimer::new(seq1, seq2, &params);
        assert!(dimer.is_ok());
        let dimer = dimer.unwrap();
        assert_eq!(dimer.mat.get(1, 1), Ok(Thermo::init_not_base_pairs()));
        assert_eq!(dimer.mat.get(1, 2), Ok(Thermo::init_not_base_pairs()));
        assert_eq!(dimer.mat.get(2, 1), Ok(Thermo::init_base_pairs()));
        assert_eq!(dimer.mat.get(2, 2), Ok(Thermo::init_base_pairs()));
        assert!(dimer.mat.get(10, 10).is_err());
    }

    #[test]
    fn test_calculate() {
        let params = get_params();
        let seq1 = "CCCCCATCCGATCAGGGGG".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();
        let dimer = Dimer::new(&seq1, &seq2, &params);
        assert!(dimer.is_ok());
        let mut dimer = dimer.unwrap();

        let res = dimer.calculate();
        assert!(res.is_ok());
        assert_eq!(res, Ok(Thermo::with_values(-288.7836433983556, -101400.0)));
    }

    #[test]
    fn test_fill() {
        let seq1 = "CCCCCATCCGATCAGGGGG".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();
        let params = get_params();

        let dimer = Dimer::new(&seq1, &seq2, &params);
        assert!(dimer.is_ok());
        let mut dimer = dimer.unwrap();
        assert!(dimer.fill().is_ok());

        macro_rules! run_test {
            ($i:literal, $j:literal, $ds:literal, $dh:literal) => {
                assert_eq!(dimer.mat.get($i, $j), Ok(Thermo::with_values($ds, $dh)),);
            };
            ($i:literal, $j:literal) => {
                assert_eq!(dimer.mat.get($i, $j), Ok(Thermo::with_inf()),);
            };
        }

        macro_rules! run_test_tb {
            ($i:literal, $j:literal, $a:literal, $b:literal) => {
                assert_eq!(dimer.traceback.get($i, $j), Ok(($a, $b)),);
            };
        }

        run_test!(0, 0, 0.0, 0.0);
        run_test!(1, 1, 0.0, 0.0);
        run_test!(1, 2, -11.2, -3900.0);
        run_test!(1, 5, -11.2, -3900.0);
        run_test!(1, 6);
        run_test!(1, 9);
        run_test!(1, 10, -3.9, -2100.0);
        run_test!(2, 1, -12.6, -4400.0);
        run_test!(2, 5, -31.099999999999998, -11900.0);
        run_test!(7, 6, -76.89999999999999, -29800.0);
        run_test!(8, 1, -10.9, -4000.0);
        run_test!(17, 17, -248.98364339835564, -85400.0);

        run_test_tb!(0, 0, -1, -1);
        run_test_tb!(1, 10, -1, -1);
        run_test_tb!(2, 2, 1, 1);
        run_test_tb!(2, 3, 1, 2);
        run_test_tb!(2, 4, 1, 3);
        run_test_tb!(2, 5, 1, 4);

        run_test_tb!(5, 2, 4, 1);
        run_test_tb!(5, 10, 4, 5);

        run_test_tb!(6, 2, -1, -1);
        run_test_tb!(6, 8, 5, 5);
        run_test_tb!(6, 13, 5, 5);

        run_test_tb!(10, 11, 9, 10);

        run_test_tb!(11, 5, -1, -1);
        run_test_tb!(11, 8, 10, 7);
        run_test_tb!(11, 13, 10, 12);

        run_test_tb!(19, 0, -1, -1);
        run_test_tb!(19, 7, 7, 6);
        run_test_tb!(19, 11, 13, 10);
        run_test_tb!(19, 12, 18, 11);
    }

    #[test]
    fn test_left_optimal_terminal() {
        let seq1 = "CCCCCATCCGATCAGGGGG".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();
        let params = get_params();

        let dimer = Dimer::new(&seq1, &seq2, &params);
        assert!(dimer.is_ok());
        let dimer = dimer.unwrap();

        macro_rules! run_test {
            ($i:literal, $j:literal, $ds:literal, $dh:literal) => {
                assert_eq!(
                    dimer.left_optimal_terminal($i, $j),
                    Ok(Thermo::with_values($ds, $dh)),
                );
            };
        }

        run_test!(19, 15, -27.4, -9800.0);
        run_test!(2, 2, -19.3, -7000.0);
        run_test!(11, 8, -9.200000000000001, -3800.0);
        run_test!(7, 14, -6.699999999999999, -2800.0);
    }

    #[test]
    fn test_right_optimal_terminal() {
        let seq1 = "CCCCCATCCGATCAGGGGG".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();
        let params = get_params();

        let dimer = Dimer::new(&seq1, &seq2, &params);
        assert!(dimer.is_ok());
        let dimer = dimer.unwrap();

        macro_rules! run_test {
            ($i:literal, $j:literal, $ds:literal, $dh:literal) => {
                assert_eq!(
                    dimer.right_optimal_terminal($i, $j),
                    Ok(Thermo::with_values($ds, $dh)),
                );
            };
        }

        run_test!(19, 15, -12.6, -4400.0);
        run_test!(2, 2, -19.3, -7000.0);
        run_test!(11, 8, -6.699999999999999, -2800.0);
        run_test!(7, 14, 25.1, 7200.0);
    }
}
