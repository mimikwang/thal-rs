const A: u8 = b'A';
const G: u8 = b'G';
const C: u8 = b'C';
const T: u8 = b'T';

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

/// Check to see if bases are WC pairs
pub fn is_base_pairs(b1: &[u8], b2: &[u8]) -> bool {
    if b1.len() != b2.len() {
        return false;
    }

    for (bb1, bb2) in b1.iter().zip(b2.iter()) {
        if !is_base_pair(bb1, bb2) {
            return false;
        }
    }

    true
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
    fn test_is_base_pairs() {
        assert!(is_base_pairs(
            &[b'T', b'A', b'C', b'G'],
            &[b'A', b'T', b'G', b'C']
        ));
        assert!(!is_base_pairs(
            &[b'T', b'A', b'C', b'G'],
            &[b'A', b'T', b'G', b'G']
        ));
    }
}
