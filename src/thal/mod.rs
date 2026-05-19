use crate::{
    errors::Result,
    matrix::Matrix,
    seq::is_base_pair,
    thermo::{Thermo, params::ThermoParams},
};

/// Add N's at the start and end of a sequence
pub fn pad_seq(seq: &mut Vec<u8>) {
    seq.push(b'N');
    seq.insert(0, b'N');
}

pub fn init_matrix(seq1: &[u8], seq2: &[u8]) -> Result<Matrix> {
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

    // If not base pairs, then skip since we already initialize it to be the right values
    if !is_base_pair(&b11, &b21) {
        return Ok(None);
    }

    // Get the AT penalty
    let at_penalty_thermo = ThermoParams::at_penalty(&b11, &b21);

    // Initialize thermo with stacking
    let mut thermo = at_penalty_thermo;
    if let Some(t) = params.get_tstack(&[b21, b20], &[b11, b10])? {
        if f64::is_finite(t.dh) {
            thermo += t;
        }
    }

    // If the b10 and b20 are not base pairs
    if !is_base_pair(&b10, &b20) {
        let mut thermo_current = at_penalty_thermo;

        let thermo_dangle_3 = params.get_dangle(&[b21, b20], &[b11])?;
        let thermo_dangle_5 = params.get_dangle(&[b21], &[b11, b10])?;

        // Get the dangling 3 end
        if b10 == b'N'
            && let Some(t) = thermo_dangle_3
        {
            if f64::is_finite(t.dh) {
                thermo_current += t;
            }
        }
        // Get the dangling 5 end
        else if b20 == b'N'
            && let Some(t) = thermo_dangle_5
        {
            if f64::is_finite(t.dh) {
                thermo_current += t;
            }
        } else {
            if let Some(t) = thermo_dangle_3 {
                thermo_current += t;
            }
            if let Some(t) = thermo_dangle_5 {
                thermo_current += t;
            }
        }

        // If thermo current has lower gibbs free energy, make thermo the current
        if thermo_current.dg() < thermo.dg() {
            thermo = thermo_current;
        }
    }

    Ok(Some(thermo))
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

    // If not base pairs, then skip since we already initialize it to be the right values
    if !is_base_pair(&b10, &b20) {
        return Ok(None);
    }

    // Get the AT penalty
    let at_penalty_thermo = ThermoParams::at_penalty(&b10, &b20);

    // Initialize thermo with stacking
    let mut thermo = at_penalty_thermo;
    if let Some(t) = params.get_tstack(&[b10, b11], &[b20, b21])? {
        if f64::is_finite(t.dh) {
            thermo += t;
        }
    }

    // If the b11 and b21 are not base pairs
    if !is_base_pair(&b11, &b21) {
        let mut thermo_current = at_penalty_thermo;

        let thermo_dangle_3 = params.get_dangle(&[b10, b11], &[b20])?;
        let thermo_dangle_5 = params.get_dangle(&[b10], &[b20, b21])?;

        // Get the dangling 3 end only
        if b21 == b'N'
            && let Some(t) = thermo_dangle_3
        {
            if f64::is_finite(t.dh) {
                thermo_current += t;
            }
        }
        // Get the dangling 5 end only
        else if b11 == b'N'
            && let Some(t) = thermo_dangle_5
        {
            if f64::is_finite(t.dh) {
                thermo_current += t;
            }
        }
        // Get Both
        else {
            if let Some(t) = thermo_dangle_3 {
                thermo_current += t;
            }
            if let Some(t) = thermo_dangle_5 {
                thermo_current += t;
            }
        }

        // If thermo current has lower gibbs free energy, make thermo the current
        if thermo_current.dg() < thermo.dg() {
            thermo = thermo_current;
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

    #[test]
    fn test_lsh() {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("thermo");
        let path = path.to_str().expect("path error").to_owned();

        let params = ThermoParams::with_file_path(&path);
        assert!(params.is_ok());
        let params = params.unwrap();

        let seq1 = "NCCCCCATCCGATCAGGGGGN".as_bytes().to_vec();
        let seq2 = seq1.clone().into_iter().rev().collect::<Vec<u8>>();

        assert_eq!(
            lsh(&seq1, &seq2, 19, 15, &params),
            Ok(Some(Thermo::with_values(-27.4, -9800.0)))
        );
    }

    #[test]
    fn test_rsh() {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("thermo");
        let path = path.to_str().expect("path error").to_owned();

        let params = ThermoParams::with_file_path(&path);
        assert!(params.is_ok());
        let params = params.unwrap();

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
}
