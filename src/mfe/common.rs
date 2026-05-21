/// Add N's at the start and end of a sequence
pub fn pad_seq(seq: &mut Vec<u8>) {
    seq.push(b'N');
    seq.insert(0, b'N');
}
