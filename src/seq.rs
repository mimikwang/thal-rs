const A: u8 = b'A';
const G: u8 = b'G';
const C: u8 = b'C';
const U: u8 = b'U';

/// Check to see if two bases are WC pairs
pub fn is_base_pair(b1: &u8, b2: &u8) -> bool {
    if (b1 == &A && b2 == &U) || (b1 == &U && b2 == &A) {
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

    false
}
