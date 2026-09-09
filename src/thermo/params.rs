use crate::{
    errors::Result,
    seq::{A_NUM, T_NUM},
    thermo::Thermo,
};

// Generated at compile time by `build.rs` from `stack.dh` / `stack.ds` (see `write_stack`).
// Provides `STACK_ARRAY: [Thermo; 625]`, indexed the same way as `stack_array` below via
// `quad_index`. Not wired into `ThermoParams` yet — see `test_generated_stack_array_matches_runtime`.
include!(concat!(env!("OUT_DIR"), "/gen_params.rs"));

#[derive(Debug)]
pub struct ThermoParams {
    stack_array: [Thermo; 625],
    stack_mm_array: [Thermo; 625],
    tstack_array: [Thermo; 625],
    dangle3_array: [Thermo; 125],
    dangle5_array: [Thermo; 125],
    triloop_array: [Thermo; 3125],
    tetraloop_array: [Thermo; 15625],

    internal_array: [Thermo; 30],
    bulge_array: [Thermo; 30],
    hairpin_array: [Thermo; 30],
}

impl ThermoParams {
    pub fn with_defaults() -> Result<Self> {
        Ok(Self {
            stack_array: STACK_ARRAY,
            stack_mm_array: STACK_MM_ARRAY,
            tstack_array: TSTACK_ARRAY,
            dangle3_array: DANGLE3_ARRAY,
            dangle5_array: DANGLE5_ARRAY,
            triloop_array: TRILOOP_ARRAY,
            tetraloop_array: TETRALOOP_ARRAY,
            internal_array: INTERNAL_ARRAY,
            bulge_array: BULGE_ARRAY,
            hairpin_array: HAIRPIN_ARRAY,
        })
    }

    pub fn get_stack(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Thermo {
        self.stack_array[quad_index(b10, b11, b20, b21)]
    }

    pub fn get_stack_mm(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Thermo {
        self.stack_mm_array[quad_index(b10, b11, b20, b21)]
    }

    pub fn get_tstack(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Thermo {
        self.tstack_array[quad_index(b10, b11, b20, b21)]
    }

    pub fn get_dangle3(&self, b10: usize, b11: usize, b20: usize) -> Thermo {
        self.dangle3_array[triple_index(b10, b11, b20)]
    }

    pub fn get_dangle5(&self, b10: usize, b20: usize, b21: usize) -> Thermo {
        self.dangle5_array[triple_index(b10, b20, b21)]
    }

    pub fn get_triloop(&self, b0: usize, b1: usize, b2: usize, b3: usize, b4: usize) -> Thermo {
        self.triloop_array[triloop_index(b0, b1, b2, b3, b4)]
    }

    #[allow(clippy::too_many_arguments)]
    pub fn get_tetraloop(
        &self,
        b0: usize,
        b1: usize,
        b2: usize,
        b3: usize,
        b4: usize,
        b5: usize,
    ) -> Thermo {
        self.tetraloop_array[tetraloop_index(b0, b1, b2, b3, b4, b5)]
    }

    pub fn get_internal(&self, loop_length: usize) -> Option<Thermo> {
        let thermo = self.internal_array.get(loop_length)?;
        Some(thermo.to_owned())
    }

    pub fn get_bulge(&self, loop_length: usize) -> Option<Thermo> {
        let thermo = self.bulge_array.get(loop_length)?;
        Some(thermo.to_owned())
    }

    pub fn get_hairpin(&self, loop_length: usize) -> Option<Thermo> {
        let thermo = self.hairpin_array.get(loop_length)?;
        Some(thermo.to_owned())
    }

    pub fn at_penalty(b1: &u8, b2: &u8) -> Thermo {
        if (b1 == &b'A' && b2 == &b'T') || (b1 == &b'T' && b2 == &b'A') {
            return Thermo::with_values(6.9, 2200.0);
        }
        Thermo::with_values(0.0, 0.0)
    }

    pub fn at_penalty_num(b1: usize, b2: usize) -> Thermo {
        if (b1 == A_NUM && b2 == T_NUM) || (b1 == T_NUM && b2 == A_NUM) {
            return Thermo::with_values(6.9, 2200.0);
        }

        Thermo::with_values(0.0, 0.0)
    }

    pub fn internal_loop(loop_size_diff: usize) -> Thermo {
        Thermo::with_values((-300.0 / 310.15) * loop_size_diff as f64, 0.0)
    }
}

fn triple_index(b0: usize, b1: usize, b2: usize) -> usize {
    (b0 * 5 + b1) * 5 + b2
}

fn quad_index(b0: usize, b1: usize, b2: usize, b3: usize) -> usize {
    triple_index(b0, b1, b2) * 5 + b3
}

fn triloop_index(b0: usize, b1: usize, b2: usize, b3: usize, b4: usize) -> usize {
    quad_index(b0, b1, b2, b3) * 5 + b4
}

#[allow(clippy::too_many_arguments)]
fn tetraloop_index(b0: usize, b1: usize, b2: usize, b3: usize, b4: usize, b5: usize) -> usize {
    triloop_index(b0, b1, b2, b3, b4) * 5 + b5
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_at_penalty() {
        assert_eq!(
            ThermoParams::at_penalty(&b'A', &b'T'),
            Thermo::with_values(6.9, 2200.0)
        );
        assert_eq!(
            ThermoParams::at_penalty(&b'T', &b'A'),
            Thermo::with_values(6.9, 2200.0)
        );
        assert_eq!(
            ThermoParams::at_penalty(&b'A', &b'A'),
            Thermo::with_values(0.0, 0.0)
        );
    }
}
