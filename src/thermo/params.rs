use std::collections::HashMap;

use crate::{
    errors::Result,
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
    tetraloop: HashMap<String, Thermo>,
    dangle: HashMap<String, Thermo>,
    loops: Loops,
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
    pub fn with_file_path(file_path: &str) -> Result<Self> {
        Ok(Self {
            stack: loader::load_stack(file_path)?,
            stack_mm: loader::load_stack_mm(file_path)?,
            tstack: loader::load_tstack(file_path)?,
            tetraloop: loader::load_tetraloop(file_path)?,
            dangle: loader::load_dangle(file_path)?,
            loops: loader::load_loops(file_path)?,
        })
    }

    get_param!(get_stack, stack);
    get_param!(get_stack_mm, stack_mm);
    get_param!(get_tstack, tstack);
    get_param!(get_dangle, dangle);

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

    pub fn at_penalty(&self, b1: &u8, b2: &u8) -> Thermo {
        if (b2 == &b'A' || b1 == &b'A') && b2 == &b'T' {
            return Thermo::with_values(6.9, -2200.0);
        }
        Thermo::with_values(0.0, 0.0)
    }
}

fn u8_to_string(b1: &[u8]) -> Result<String> {
    String::from_utf8(b1.to_vec()).map_err(|_| "error parsing string")
}

fn get_lookup(b1: &[u8], b2: &[u8]) -> Result<String> {
    let b1_str = u8_to_string(b1)?;
    let b2_str = u8_to_string(b2)?;
    Ok(format!("{b1_str}_{b2_str}"))
}
