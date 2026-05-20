use crate::{
    errors::Result,
    matrix::Matrix,
    seq::is_base_pair,
    thermo::{Thermo, params::ThermoParams},
};

pub fn thal(seq1: &[u8], seq2: &[u8], params: &ThermoParams) -> Result<Thermo> {
    let mut seq1 = seq1.to_owned();
    let mut seq2 = seq2.to_owned();
    pad_seq(&mut seq1);
    pad_seq(&mut seq2);

    let mut mat = init_matrix(&seq1, &seq2)?;
    let mut traceback_mat = Matrix::with_value(seq1.len(), seq2.len(), (-1, -1));
    fill_matrix(30, &mut mat, &mut traceback_mat, &seq1, &seq2, params)?;

    let mut best_i = 1;
    let mut best_j = 1;
    let mut thermo = Thermo::with_inf();

    for i in 1..seq1.len() {
        for j in 1..seq2.len() {
            if is_base_pair(&seq1[i], &seq2[j]) {
                let mut thermo_current = mat.get(i, j)? + Thermo::with_values(-5.7, 200.0);
                if let Some(t) = rsh(&seq1, &seq2, i, j, params)? {
                    thermo_current += t;
                }

                if thermo_current.dg() < thermo.dg() {
                    thermo = thermo_current;
                    best_i = i;
                    best_j = j;
                }
            }
        }
    }

    if is_base_pair(&seq1[best_i], &seq2[best_j]) {
        let rsh_thermo = rsh(&seq1, &seq2, best_i, best_j, params)?;
        thermo = mat.get(best_i, best_j)?;
        if let Some(t) = rsh_thermo {
            thermo += t;
        }
    }

    Ok(thermo)
}

/// Add N's at the start and end of a sequence
pub fn pad_seq(seq: &mut Vec<u8>) {
    seq.push(b'N');
    seq.insert(0, b'N');
}

pub fn init_matrix(seq1: &[u8], seq2: &[u8]) -> Result<Matrix<Thermo>> {
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

pub fn fill_matrix(
    max_loop: usize,
    mat: &mut Matrix<Thermo>,
    traceback_mat: &mut Matrix<(i8, i8)>,
    seq1: &[u8],
    seq2: &[u8],
    params: &ThermoParams,
) -> Result<()> {
    for i in 1..seq1.len() {
        for j in 1..seq2.len() {
            if !is_base_pair(&seq1[i], &seq2[j]) {
                mat.set(i, j, Thermo::with_inf())?;
                continue;
            }

            let mut thermo = Thermo::new();
            if let Some(t) = lsh(seq1, seq2, i, j, params)? {
                thermo += t;
            }
            mat.set(i, j, thermo)?;

            if i > 1 && j > 1 {
                let rsh_thermo = rsh(seq1, seq2, i, j, params)?;

                if is_base_pair(&seq1[i - 1], &seq2[j - 1]) {
                    let mut thermo = mat.get(i - 1, j - 1)?;
                    let stack_thermo =
                        params.get_stack(&[seq1[i - 1], seq1[i]], &[seq2[j - 1], seq2[j]])?;
                    if let Some(t) = stack_thermo {
                        thermo += t;
                    }
                    mat.set(i, j, thermo)?;
                    traceback_mat.set(i, j, (i as i8 - 1, j as i8 - 1))?;
                }

                let mut best_thermo = mat.get(i, j)?;
                if let Some(t) = rsh_thermo {
                    best_thermo += t;
                }

                for d in 3..max_loop + 3 {
                    let mut ii = i as i8 - 1;
                    let mut jj = j as i8 + i as i8 - ii - d as i8;
                    if jj < 1 {
                        ii += jj - 1;
                        jj = 1;
                    }
                    while ii > 0 && jj < j as i8 {
                        if is_base_pair(&seq1[ii as usize], &seq2[jj as usize]) {
                            let internal_thermo_calc = get_internal_thermo(
                                ii as usize,
                                jj as usize,
                                i,
                                j,
                                seq1,
                                seq2,
                                params,
                            )?;
                            let mut internal_thermo = mat.get(ii as usize, jj as usize)?;
                            if let Some(t) = internal_thermo_calc {
                                internal_thermo += t;
                            }

                            if internal_thermo.dg() < best_thermo.dg() {
                                best_thermo = internal_thermo;
                                if mat.get(i, j)? != internal_thermo {
                                    traceback_mat.set(i, j, (ii as i8, jj as i8))?;
                                }
                                mat.set(i, j, internal_thermo)?;
                            }
                        }
                        ii -= 1;
                        jj += 1;
                    }
                }
            }
        }
    }

    Ok(())
}

pub fn lsh(
    seq1: &[u8],
    seq2: &[u8],
    i: usize,
    j: usize,
    params: &ThermoParams,
) -> Result<Option<Thermo>> {
    let b10 = seq1[i - 1];
    let b11 = seq1[i];
    let b20 = seq2[j - 1];
    let b21 = seq2[j];

    get_terminal_thermo(&[b21, b20], &[b11, b10], params)
}

pub fn rsh(
    seq1: &[u8],
    seq2: &[u8],
    i: usize,
    j: usize,
    params: &ThermoParams,
) -> Result<Option<Thermo>> {
    let b10 = seq1[i];
    let b11 = seq1[i + 1];
    let b20 = seq2[j];
    let b21 = seq2[j + 1];

    get_terminal_thermo(&[b10, b11], &[b20, b21], params)
}

/// Get the optimal terminal thermo
///
/// The seq1 is oriented as 5' to 3' and seq2 is oriented 3' to 5'
///
/// b10
///    \
///    b11
///    b21
///    /
/// b20
fn get_terminal_thermo(
    seq1: &[u8; 2],
    seq2: &[u8; 2],
    params: &ThermoParams,
) -> Result<Option<Thermo>> {
    let b10 = seq1[0];
    let b11 = seq1[1];
    let b20 = seq2[0];
    let b21 = seq2[1];

    // If b11 and b21 are base pairs
    if !is_base_pair(&b10, &b20) {
        return Ok(None);
    }

    // Get the at penalty
    let at_penalty_thermo = ThermoParams::at_penalty(&b10, &b20);

    // Initialize thermo with terminal stacking
    let mut thermo = at_penalty_thermo;
    if let Some(t) = params.get_tstack(&[b10, b11], &[b20, b21])?
        && f64::is_finite(t.dh)
    {
        thermo += t;
    }

    // If b11 and b21 are not base pairs
    if !is_base_pair(&b11, &b21) {
        let mut thermo_current = at_penalty_thermo;
        let thermo_dangle_3 = params.get_dangle(&[b10, b11], &[b20])?;
        let thermo_dangle_5 = params.get_dangle(&[b10], &[b20, b21])?;

        if b21 == b'N'
            && let Some(t) = thermo_dangle_3
            && f64::is_finite(t.dh)
        {
            thermo_current += t;
        } else if b11 == b'N'
            && let Some(t) = thermo_dangle_5
            && f64::is_finite(t.dh)
        {
            thermo_current += t;
        } else {
            if let Some(t) = thermo_dangle_3 {
                thermo_current += t;
            }
            if let Some(t) = thermo_dangle_5 {
                thermo_current += t;
            }
        }

        if thermo_current.dg() < thermo.dg() {
            thermo = thermo_current;
        }
    }

    Ok(Some(thermo))
}

fn get_internal_thermo(
    i: usize,
    j: usize,
    ii: usize,
    jj: usize,
    seq1: &[u8],
    seq2: &[u8],
    params: &ThermoParams,
) -> Result<Option<Thermo>> {
    let loop_size_1 = ii as i8 - i as i8 - 1;
    let loop_size_2 = jj as i8 - j as i8 - 1;
    let loop_size = loop_size_1 + loop_size_2 - 1;
    let mut thermo = Thermo::with_inf();
    let bulge_thermo = params.get_bulge(loop_size as usize);

    if loop_size_1 == 0 || loop_size_2 == 0 {
        thermo = Thermo::new();
        if let Some(t) = bulge_thermo
            && f64::is_finite(t.dh)
        {
            thermo += t;
        }
        if loop_size_1 == 1 || loop_size_2 == 1 {
            let stack_thermo = params.get_stack(&[seq1[i], seq1[ii]], &[seq2[j], seq2[jj]])?;
            if let Some(t) = stack_thermo
                && f64::is_finite(t.dh)
            {
                thermo += t;
            }
        } else {
            thermo += ThermoParams::at_penalty(&seq1[i], &seq2[j])
                + ThermoParams::at_penalty(&seq1[ii], &seq2[jj]);
        }
    } else if loop_size_1 == 1 && loop_size_2 == 1 {
        thermo = Thermo::new();
        let stack_thermo_1 =
            params.get_stack_mm(&[seq1[i], seq1[i + 1]], &[seq2[j], seq2[j + 1]])?;
        let stack_thermo_2 =
            params.get_stack_mm(&[seq2[jj], seq2[jj - 1]], &[seq1[ii], seq1[ii - 1]])?;

        if let Some(t) = stack_thermo_1
            && f64::is_finite(t.dh)
        {
            thermo += t;
        }
        if let Some(t) = stack_thermo_2
            && f64::is_finite(t.dh)
        {
            thermo += t;
        }
    } else {
        if !is_base_pair(&seq1[ii - 1], &seq2[jj - 1]) && !is_base_pair(&seq1[i + 1], &seq2[j + 1])
        {
            thermo = Thermo::new();
            let internal_thermo = params.get_internal(loop_size as usize);
            let tstack_thermo_1 =
                params.get_tstack(&[seq1[i], seq1[i + 1]], &[seq2[j], seq2[j + 1]])?;
            let tstack_thermo_2 =
                params.get_tstack(&[seq2[jj], seq2[jj - 1]], &[seq1[ii], seq1[ii - 1]])?;

            if let Some(t) = internal_thermo
                && f64::is_finite(t.dh)
            {
                thermo += t;
            }
            if let Some(t) = tstack_thermo_1
                && f64::is_finite(t.dh)
            {
                thermo += t;
            }
            if let Some(t) = tstack_thermo_2
                && f64::is_finite(t.dh)
            {
                thermo += t;
            }
            let internal_thermo = ThermoParams::internal_loop(
                (loop_size_1 as i8 - loop_size_2 as i8).unsigned_abs() as usize,
            );
            if f64::is_finite(internal_thermo.dh) {
                thermo += internal_thermo;
            }
        }
    }

    Ok(Some(thermo))
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn test_pad_seq() {
        let mut seq = vec![b'T', b'A', b'C'];
        pad_seq(&mut seq);
        assert_eq!(seq, vec![b'N', b'T', b'A', b'C', b'N']);
    }

    #[test]
    fn test_init_matrix() {
        let mat = init_matrix(
            &[b'N', b'A', b'T', b'C', b'G', b'N'],
            &[b'N', b'A', b'A', b'N'],
        );
        assert!(mat.is_ok());

        let mat = mat.unwrap();
        assert_eq!(mat.get(1, 1), Ok(Thermo::init_not_base_pairs()));
        assert_eq!(mat.get(1, 2), Ok(Thermo::init_not_base_pairs()));
        assert_eq!(mat.get(2, 1), Ok(Thermo::init_base_pairs()));
        assert_eq!(mat.get(2, 2), Ok(Thermo::init_base_pairs()));
        assert!(mat.get(10, 10).is_err());
    }

    fn get_thermo_params() -> ThermoParams {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("thermo");
        let path = path.to_str().expect("path error").to_owned();

        let params = ThermoParams::with_file_path(&path);
        assert!(params.is_ok());
        params.unwrap()
    }

    #[test]
    fn test_thal() {
        let params = get_thermo_params();
        let seq1 = "CCCCCATCCGATCAGGGGG".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();

        let res = thal(&seq1, &seq2, &params);
        assert!(res.is_ok());
        assert_eq!(res, Ok(Thermo::with_values(-288.7836433983556, -101400.0)));
    }

    #[test]
    fn test_fill_matrix() {
        let params = get_thermo_params();
        let seq1 = "NACN".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();
        let mut mat = Matrix::with_value(seq1.len(), seq2.len(), Thermo::new());
        let mut traceback_mat = Matrix::with_value(seq1.len(), seq2.len(), (-1, -1));
        let max_loop = 30;

        let res = fill_matrix(
            max_loop,
            &mut mat,
            &mut traceback_mat,
            &seq1,
            &seq2,
            &params,
        );
        assert!(res.is_ok());

        assert_eq!(mat.get(0, 0), Ok(Thermo::with_values(0.0, 0.0)));
        assert_eq!(mat.get(1, 1), Ok(Thermo::with_inf()));
        assert_eq!(traceback_mat, Matrix::with_value(4, 4, (-1, -1)));
    }

    #[test]
    fn test_fill_matrix_full() {
        let params = get_thermo_params();
        let seq1 = "NCCCCCATCCGATCAGGGGGN".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();
        let mut mat = Matrix::with_value(seq1.len(), seq2.len(), Thermo::new());
        let mut traceback_mat = Matrix::with_value(seq1.len(), seq2.len(), (-1, -1));
        let max_loop = 30;

        let res = fill_matrix(
            max_loop,
            &mut mat,
            &mut traceback_mat,
            &seq1,
            &seq2,
            &params,
        );
        assert!(res.is_ok());

        assert_eq!(mat.get(0, 0), Ok(Thermo::new()));
        assert_eq!(mat.get(1, 1), Ok(Thermo::new()));
        assert_eq!(mat.get(1, 2), Ok(Thermo::with_values(-11.2, -3900.0)));
        assert_eq!(mat.get(1, 5), Ok(Thermo::with_values(-11.2, -3900.0)));
        assert_eq!(mat.get(1, 6), Ok(Thermo::with_inf()));
        assert_eq!(mat.get(1, 9), Ok(Thermo::with_inf()));
        assert_eq!(mat.get(1, 10), Ok(Thermo::with_values(-3.9, -2100.0)));
        assert_eq!(mat.get(2, 1), Ok(Thermo::with_values(-12.6, -4400.0)));
        assert_eq!(
            mat.get(2, 5),
            Ok(Thermo::with_values(-31.099999999999998, -11900.0))
        );
        assert_eq!(
            mat.get(7, 6),
            Ok(Thermo::with_values(-76.89999999999999, -29800.0))
        );
        assert_eq!(mat.get(8, 1), Ok(Thermo::with_values(-10.9, -4000.0)));
        assert_eq!(
            mat.get(17, 17),
            Ok(Thermo::with_values(-248.98364339835564, -85400.0))
        );
        assert_eq!(traceback_mat.get(0, 0), Ok((-1, -1)));
        assert_eq!(traceback_mat.get(1, 10), Ok((-1, -1)));

        assert_eq!(traceback_mat.get(2, 2), Ok((1, 1)));
        assert_eq!(traceback_mat.get(2, 3), Ok((1, 2)));
        assert_eq!(traceback_mat.get(2, 4), Ok((1, 3)));
        assert_eq!(traceback_mat.get(2, 5), Ok((1, 4)));

        assert_eq!(traceback_mat.get(5, 2), Ok((4, 1)));
        assert_eq!(traceback_mat.get(5, 10), Ok((4, 5)));

        assert_eq!(traceback_mat.get(6, 2), Ok((-1, -1)));
        assert_eq!(traceback_mat.get(6, 8), Ok((5, 5)));
        assert_eq!(traceback_mat.get(6, 13), Ok((5, 5)));

        assert_eq!(traceback_mat.get(10, 11), Ok((9, 10)));

        assert_eq!(traceback_mat.get(11, 5), Ok((-1, -1)));
        assert_eq!(traceback_mat.get(11, 8), Ok((10, 7)));
        assert_eq!(traceback_mat.get(11, 13), Ok((10, 12)));

        assert_eq!(traceback_mat.get(19, 0), Ok((-1, -1)));
        assert_eq!(traceback_mat.get(19, 7), Ok((7, 6)));
        assert_eq!(traceback_mat.get(19, 11), Ok((13, 10)));
        assert_eq!(traceback_mat.get(19, 12), Ok((18, 11)));
    }

    #[test]
    fn test_lsh() {
        let params = get_thermo_params();

        let seq1 = "NCCCCCATCCGATCAGGGGGN".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();

        assert_eq!(
            lsh(&seq1, &seq2, 19, 15, &params),
            Ok(Some(Thermo::with_values(-27.4, -9800.0)))
        );
        assert_eq!(
            lsh(&seq1, &seq2, 2, 2, &params),
            Ok(Some(Thermo::with_values(-19.3, -7000.0)))
        );
        assert_eq!(
            lsh(&seq1, &seq2, 11, 8, &params),
            Ok(Some(Thermo::with_values(-9.200000000000001, -3800.0)))
        );
        assert_eq!(
            lsh(&seq1, &seq2, 7, 14, &params),
            Ok(Some(Thermo::with_values(-6.699999999999999, -2800.0)))
        );
    }

    #[test]
    fn test_rsh() {
        let params = get_thermo_params();

        let seq1 = "NCCCCCATCCGATCAGGGGGN".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();

        assert_eq!(
            rsh(&seq1, &seq2, 19, 15, &params),
            Ok(Some(Thermo::with_values(-12.6, -4400.0)))
        );
        assert_eq!(
            rsh(&seq1, &seq2, 2, 2, &params),
            Ok(Some(Thermo::with_values(-19.3, -7000.0)))
        );
        assert_eq!(
            rsh(&seq1, &seq2, 11, 8, &params),
            Ok(Some(Thermo::with_values(-6.699999999999999, -2800.0)))
        );
        assert_eq!(
            rsh(&seq1, &seq2, 7, 14, &params),
            Ok(Some(Thermo::with_values(25.1, 7200.0)))
        );
    }

    #[test]
    fn test_get_internal_thermo() {
        let params = get_thermo_params();

        let seq1 = "NCCCCCATCCGATCAGGGGGN".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();

        assert_eq!(
            get_internal_thermo(1, 2, 2, 4, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-32.79, -8000.0)))
        );
        assert_eq!(
            get_internal_thermo(1, 2, 2, 10, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-11.92, 0.0)))
        );
        assert_eq!(
            get_internal_thermo(1, 5, 4, 10, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-35.91454779945188, -7700.0)))
        );
        assert_eq!(
            get_internal_thermo(16, 18, 19, 19, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-9.35, 0.0)))
        );
        assert_eq!(
            get_internal_thermo(18, 15, 19, 19, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-9.99, 0.0)))
        );
        assert_eq!(
            get_internal_thermo(15, 16, 19, 19, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_inf()))
        );
        assert_eq!(
            get_internal_thermo(18, 12, 19, 19, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-11.28, 0.0)))
        );
        assert_eq!(
            get_internal_thermo(3, 1, 4, 10, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-12.57, 0.0)))
        );
        assert_eq!(
            get_internal_thermo(5, 4, 7, 6, &seq1, &seq2, &params),
            Ok(Some(Thermo::with_values(-4.6, -1400.0)))
        );
    }
}
