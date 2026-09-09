use std::collections::HashMap;

use crate::{
    errors::Result,
    seq::{A_NUM, T_NUM, base_to_num},
    thermo::{
        Thermo,
        loader::{self, Loops},
    },
};

#[derive(Debug)]
pub struct ThermoParams {
    stack: HashMap<String, Thermo>,
    stack_mm: HashMap<String, Thermo>,
    tstack: HashMap<String, Thermo>,
    triloop: HashMap<String, Thermo>,
    tetraloop: HashMap<String, Thermo>,
    dangle: HashMap<String, Thermo>,
    loops: Loops,

    stack_array: Vec<Thermo>,
    stack_mm_array: Vec<Thermo>,
    tstack_array: Vec<Thermo>,
    dangle3_array: Vec<Thermo>,
    dangle5_array: Vec<Thermo>,
    triloop_array: Vec<Thermo>,
    tetraloop_array: Vec<Thermo>,
}

macro_rules! get_param {
    ($name:ident, $param:ident) => {
        pub fn $name(&self, b1: &[u8], b2: &[u8]) -> Result<Option<Thermo>> {
            let lookup = get_lookup(b1, b2)?;
            let Some(thermo) = self.$param.get(&lookup) else {
                return Ok(None);
            };
            Ok(Some(thermo.to_owned()))
        }
    };
}

macro_rules! get_loop {
    ($name:ident, $param:ident) => {
        pub fn $name(&self, loop_length: usize) -> Option<Thermo> {
            let thermo = self.loops.$param.get(loop_length)?;
            Some(thermo.to_owned())
        }
    };
}

impl ThermoParams {
    pub fn with_defaults() -> Result<Self> {
        let stack = loader::default_stack()?;
        let stack_mm = loader::default_stack_mm()?;
        let tstack = loader::default_tstack()?;
        let triloop = loader::default_triloop()?;
        let tetraloop = loader::default_tetraloop()?;
        let dangle = loader::default_dangle()?;

        let stack_array = ThermoParams::build_pair_array(&stack)?;
        let stack_mm_array = ThermoParams::build_pair_array(&stack_mm)?;
        let tstack_array = ThermoParams::build_pair_array(&tstack)?;
        let dangle3_array = ThermoParams::build_dangle3_array(&dangle)?;
        let dangle5_array = ThermoParams::build_dangle5_array(&dangle)?;
        let triloop_array = ThermoParams::build_triloop_array(&triloop)?;
        let tetraloop_array = ThermoParams::build_tetraloop_array(&tetraloop)?;

        Ok(Self {
            stack,
            stack_mm,
            tstack,
            triloop,
            tetraloop,
            dangle,
            loops: loader::default_loops()?,
            stack_array,
            stack_mm_array,
            tstack_array,
            dangle3_array,
            dangle5_array,
            triloop_array,
            tetraloop_array,
        })
    }

    pub fn with_file_path(file_path: &str) -> Result<Self> {
        let stack = loader::load_stack(file_path)?;
        let stack_mm = loader::load_stack_mm(file_path)?;
        let tstack = loader::load_tstack(file_path)?;
        let triloop = loader::load_triloop(file_path)?;
        let tetraloop = loader::load_tetraloop(file_path)?;
        let dangle = loader::load_dangle(file_path)?;

        let stack_array = ThermoParams::build_pair_array(&stack)?;
        let stack_mm_array = ThermoParams::build_pair_array(&stack_mm)?;
        let tstack_array = ThermoParams::build_pair_array(&tstack)?;
        let dangle3_array = ThermoParams::build_dangle3_array(&dangle)?;
        let dangle5_array = ThermoParams::build_dangle5_array(&dangle)?;
        let triloop_array = ThermoParams::build_triloop_array(&triloop)?;
        let tetraloop_array = ThermoParams::build_tetraloop_array(&tetraloop)?;

        Ok(Self {
            stack,
            stack_mm,
            tstack,
            triloop,
            tetraloop,
            dangle,
            loops: loader::load_loops(file_path)?,
            stack_array,
            stack_mm_array,
            tstack_array,
            dangle3_array,
            dangle5_array,
            triloop_array,
            tetraloop_array,
        })
    }

    pub fn get_stack_fast(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Thermo {
        self.stack_array[quad_index(b10, b11, b20, b21)]
    }

    pub fn get_stack_mm_fast(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Thermo {
        self.stack_mm_array[quad_index(b10, b11, b20, b21)]
    }

    pub fn get_tstack_fast(&self, b10: usize, b11: usize, b20: usize, b21: usize) -> Thermo {
        self.tstack_array[quad_index(b10, b11, b20, b21)]
    }

    pub fn get_dangle3_fast(&self, b10: usize, b11: usize, b20: usize) -> Thermo {
        self.dangle3_array[triple_index(b10, b11, b20)]
    }

    pub fn get_dangle5_fast(&self, b10: usize, b20: usize, b21: usize) -> Thermo {
        self.dangle5_array[triple_index(b10, b20, b21)]
    }

    pub fn get_triloop_fast(
        &self,
        b0: usize,
        b1: usize,
        b2: usize,
        b3: usize,
        b4: usize,
    ) -> Thermo {
        self.triloop_array[triloop_index(b0, b1, b2, b3, b4)]
    }

    #[allow(clippy::too_many_arguments)]
    pub fn get_tetraloop_fast(
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

    get_param!(get_stack, stack);
    get_param!(get_stack_mm, stack_mm);
    get_param!(get_tstack, tstack);
    get_param!(get_dangle, dangle);

    pub fn get_triloop(&self, bases: &[u8]) -> Result<Option<Thermo>> {
        let lookup = u8_to_string(bases)?;
        let Some(thermo) = self.triloop.get(&lookup) else {
            return Ok(None);
        };
        Ok(Some(thermo.to_owned()))
    }

    pub fn get_tetraloop(&self, bases: &[u8]) -> Result<Option<Thermo>> {
        let lookup = u8_to_string(bases)?;
        let Some(thermo) = self.tetraloop.get(&lookup) else {
            return Ok(None);
        };
        Ok(Some(thermo.to_owned()))
    }

    get_loop!(get_internal, internal);
    get_loop!(get_bulge, bulge);
    get_loop!(get_hairpin, hairpin);

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

    fn build_pair_array(map: &HashMap<String, Thermo>) -> Result<Vec<Thermo>> {
        let mut array = vec![Thermo::with_values(f64::INFINITY, f64::INFINITY); 5usize.pow(4)];

        for (key, thermo) in map {
            let (b1, b2) = key.split_once('_').ok_or("invalid key")?;
            let b1 = b1.as_bytes();
            let b2 = b2.as_bytes();
            if b1.len() != 2 || b2.len() != 2 {
                return Err("invalid key length");
            }

            let index = quad_index(
                base_to_num(b1[0]),
                base_to_num(b1[1]),
                base_to_num(b2[0]),
                base_to_num(b2[1]),
            );
            array[index] = *thermo;
        }

        Ok(array)
    }

    fn build_dangle3_array(dangle: &HashMap<String, Thermo>) -> Result<Vec<Thermo>> {
        let mut array = vec![Thermo::with_values(f64::INFINITY, f64::INFINITY); 5usize.pow(3)];

        for (key, thermo) in dangle {
            let (b1, b2) = key.split_once('_').ok_or("invalid dangle key")?;
            let b1 = b1.as_bytes();
            let b2 = b2.as_bytes();
            if b1.len() != 2 || b2.len() != 1 {
                continue;
            }

            let index = triple_index(base_to_num(b1[0]), base_to_num(b1[1]), base_to_num(b2[0]));
            array[index] = *thermo;
        }

        Ok(array)
    }

    fn build_dangle5_array(dangle: &HashMap<String, Thermo>) -> Result<Vec<Thermo>> {
        let mut array = vec![Thermo::with_values(f64::INFINITY, f64::INFINITY); 5usize.pow(3)];

        for (key, thermo) in dangle {
            let (b1, b2) = key.split_once('_').ok_or("invalid dangle key")?;
            let b1 = b1.as_bytes();
            let b2 = b2.as_bytes();
            if b1.len() != 1 || b2.len() != 2 {
                continue;
            }

            let index = triple_index(base_to_num(b1[0]), base_to_num(b2[0]), base_to_num(b2[1]));
            array[index] = *thermo;
        }

        Ok(array)
    }

    fn build_triloop_array(triloop: &HashMap<String, Thermo>) -> Result<Vec<Thermo>> {
        let mut array = vec![Thermo::with_values(f64::INFINITY, f64::INFINITY); 5usize.pow(5)];

        for (key, thermo) in triloop {
            let bases = key.as_bytes();
            if bases.len() != 5 {
                return Err("invalid triloop key length");
            }

            let index = triloop_index(
                base_to_num(bases[0]),
                base_to_num(bases[1]),
                base_to_num(bases[2]),
                base_to_num(bases[3]),
                base_to_num(bases[4]),
            );
            array[index] = *thermo;
        }

        Ok(array)
    }

    fn build_tetraloop_array(tetraloop: &HashMap<String, Thermo>) -> Result<Vec<Thermo>> {
        let mut array = vec![Thermo::with_values(f64::INFINITY, f64::INFINITY); 5usize.pow(6)];

        for (key, thermo) in tetraloop {
            let bases = key.as_bytes();
            if bases.len() != 6 {
                return Err("invalid tetraloop key length");
            }

            let index = tetraloop_index(
                base_to_num(bases[0]),
                base_to_num(bases[1]),
                base_to_num(bases[2]),
                base_to_num(bases[3]),
                base_to_num(bases[4]),
                base_to_num(bases[5]),
            );
            array[index] = *thermo;
        }

        Ok(array)
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

fn u8_to_string(b1: &[u8]) -> Result<String> {
    String::from_utf8(b1.to_vec()).map_err(|_| "error parsing string")
}

fn get_lookup(b1: &[u8], b2: &[u8]) -> Result<String> {
    let b1_str = u8_to_string(b1)?;
    let b2_str = u8_to_string(b2)?;
    Ok(format!("{b1_str}_{b2_str}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::base_to_num;

    #[test]
    fn test_get_stack_fast_matches_get_stack() {
        let params = ThermoParams::with_defaults().unwrap();

        macro_rules! run_test {
            ($b10:literal, $b11:literal, $b20:literal, $b21:literal) => {
                let expected = params
                    .get_stack(&[$b10, $b11], &[$b20, $b21])
                    .unwrap()
                    .unwrap_or(Thermo::with_values(f64::INFINITY, f64::INFINITY));
                let actual = params.get_stack_fast(
                    base_to_num($b10),
                    base_to_num($b11),
                    base_to_num($b20),
                    base_to_num($b21),
                );
                assert_eq!(actual, expected);
            };
        }

        run_test!(b'A', b'A', b'T', b'T');
        run_test!(b'A', b'C', b'T', b'G');
        run_test!(b'A', b'G', b'T', b'C');
        run_test!(b'A', b'A', b'A', b'A');
        run_test!(b'G', b'C', b'G', b'C');
        run_test!(b'T', b'A', b'T', b'A');
    }

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
