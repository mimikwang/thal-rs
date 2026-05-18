use crate::{errors::Result, matrix::Matrix, seq::is_base_pair, thermo::Thermo};

pub fn init_matrix(seq1: &[u8], seq2: &[u8]) -> Result<Matrix> {
    // Make a matrix that has a padding on both ends
    let mut mat = Matrix::new(seq1.len() + 2, seq2.len() + 2);

    // Fill based on base pairings
    for i in 1..=seq1.len() {
        for j in 1..=seq2.len() {
            if is_base_pair(&seq1[i - 1], &seq2[j - 1]) {
                mat.set(i, j, Thermo::init_base_pairs())?;
            } else {
                mat.set(i, j, Thermo::init_not_base_pairs())?;
            }
        }
    }

    Ok(mat)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init_matrix() {
        let mat = init_matrix(&[b'A', b'T', b'C', b'G'], &[b'A', b'A']);
        assert!(mat.is_ok());

        let mat = mat.unwrap();
        assert_eq!(mat.get(1, 1), Ok(Thermo::init_not_base_pairs()));
        assert_eq!(mat.get(1, 2), Ok(Thermo::init_not_base_pairs()));
        assert_eq!(mat.get(2, 1), Ok(Thermo::init_base_pairs()));
        assert_eq!(mat.get(2, 2), Ok(Thermo::init_base_pairs()));
        assert!(mat.get(10, 10).is_err());
    }
}
