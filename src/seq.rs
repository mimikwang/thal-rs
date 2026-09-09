const A: u8 = b'A';
const G: u8 = b'G';
const C: u8 = b'C';
const T: u8 = b'T';

const A_NUM: u8 = 0;
const G_NUM: u8 = 1;
const C_NUM: u8 = 2;
const T_NUM: u8 = 3;
const N_NUM: u8 = 4;

const BASE_LUT: [u8; 256] = {
    let mut table = [N_NUM; 256];
    table[A as usize] = A_NUM;
    table[C as usize] = C_NUM;
    table[G as usize] = G_NUM;
    table[T as usize] = T_NUM;
    table
};

const BASE_PAIR_LUT: [[bool; 5]; 5] = {
    let mut table = [[false; 5]; 5];
    table[A_NUM as usize][T_NUM as usize] = true;
    table[T_NUM as usize][A_NUM as usize] = true;
    table[C_NUM as usize][G_NUM as usize] = true;
    table[G_NUM as usize][C_NUM as usize] = true;

    table
};

/// Check to see if two bases are WC pairs
pub fn is_base_pair(b1: &u8, b2: &u8) -> bool {
    if (b1 == &A && b2 == &T) || (b1 == &T && b2 == &A) {
        return true;
    }

    if (b1 == &C && b2 == &G) || (b1 == &G && b2 == &C) {
        return true;
    }

    false
}

/// Convert bases to number where
///
/// A --> 0
/// G --> 1
/// C --> 2
/// T --> 3
///
/// And anything else is 4. This helps with building fast lookups instead of using a hashmap for
/// with strings as keys. We can use an array lookup instead.
pub fn seq_to_num(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| BASE_LUT[b as usize]).collect()
}

pub fn is_base_pair_num(b1: &u8, b2: &u8) -> bool {
    BASE_PAIR_LUT[*b1 as usize][*b2 as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_base_pair() {
        assert!(is_base_pair(&b'T', &b'A'));
        assert!(is_base_pair(&b'A', &b'T'));
        assert!(is_base_pair(&b'C', &b'G'));
        assert!(is_base_pair(&b'G', &b'C'));
        assert!(!is_base_pair(&b'T', &b'T'));
        assert!(!is_base_pair(&b'G', &b'T'));
    }

    #[test]
    fn test_seq_to_num() {
        assert_eq!(
            seq_to_num(&[b'A', b'G', b'C', b'T', b'Z']),
            vec![0, 1, 2, 3, 4]
        );
    }

    #[test]
    fn test_is_base_pair_num() {
        assert!(is_base_pair_num(&3, &0));
        assert!(is_base_pair_num(&0, &3));
        assert!(is_base_pair_num(&2, &1));
        assert!(is_base_pair_num(&1, &2));
        assert!(!is_base_pair_num(&0, &0));
        assert!(!is_base_pair_num(&4, &4));
    }
}
