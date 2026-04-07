use std::collections::HashMap;

use crate::thermo::{
    Thermo,
    loader::{self, Loops},
};

#[derive(Debug)]
pub struct ThermoParams {
    stack: HashMap<String, Thermo>,
    tstack: HashMap<String, Thermo>,
    tetraloop: HashMap<String, Thermo>,
    dangle: HashMap<String, Thermo>,
    loops: Loops,
}

impl ThermoParams {
    pub fn with_file_path(file_path: &str) -> Result<Self, &'static str> {
        Ok(Self {
            stack: loader::load_stack(file_path)?,
            tstack: loader::load_tstack(file_path)?,
            tetraloop: loader::load_tetraloop(file_path)?,
            dangle: loader::load_dangle(file_path)?,
            loops: loader::load_loops(file_path)?,
        })
    }

    pub fn get_stack(&self, b1: &[u8], b2: &[u8]) -> Result<Option<Thermo>, &'static str> {
        let lookup = get_lookup(b1, b2)?;
        let Some(thermo) = self.stack.get(&lookup) else {
            return Ok(None);
        };
        Ok(Some(thermo.to_owned()))
    }

    pub fn get_tstack(&self, b1: &[u8], b2: &[u8]) -> Result<Option<Thermo>, &'static str> {
        let lookup = get_lookup(b1, b2)?;
        let Some(thermo) = self.tstack.get(&lookup) else {
            return Ok(None);
        };
        Ok(Some(thermo.to_owned()))
    }

    pub fn get_dangle(&self, b1: &[u8], b2: &[u8]) -> Result<Option<Thermo>, &'static str> {
        let lookup = get_lookup(b1, b2)?;
        let Some(thermo) = self.dangle.get(&lookup) else {
            return Ok(None);
        };
        Ok(Some(thermo.to_owned()))
    }

    pub fn get_tetraloop(&self, bases: &[u8]) -> Result<Option<Thermo>, &'static str> {
        let lookup = u8_to_string(bases)?;
        let Some(thermo) = self.tetraloop.get(&lookup) else {
            return Ok(None);
        };
        Ok(Some(thermo.to_owned()))
    }

    pub fn get_internal(&self, loop_length: usize) -> Option<Thermo> {
        let Some(thermo) = self.loops.internal.get(loop_length) else {
            return None;
        };
        Some(thermo.to_owned())
    }

    pub fn get_bulge(&self, loop_length: usize) -> Option<Thermo> {
        let Some(thermo) = self.loops.bulge.get(loop_length) else {
            return None;
        };
        Some(thermo.to_owned())
    }

    pub fn get_hairpin(&self, loop_length: usize) -> Option<Thermo> {
        let Some(thermo) = self.loops.hairpin.get(loop_length) else {
            return None;
        };
        Some(thermo.to_owned())
    }

    pub fn at_penalty(&self, b1: &u8, b2: &u8) -> Thermo {
        if (b1 == &b'A' && b2 == &b'T') || (b2 == &b'A' && b2 == &b'T') {
            return Thermo::with_values(6.9, -2200.0);
        }
        Thermo::with_values(0.0, 0.0)
    }
}

fn u8_to_string(b1: &[u8]) -> Result<String, &'static str> {
    String::from_utf8(b1.to_vec()).map_err(|_| "error parsing string")
}

fn get_lookup(b1: &[u8], b2: &[u8]) -> Result<String, &'static str> {
    let b1_str = u8_to_string(b1)?;
    let b2_str = u8_to_string(b2)?;
    Ok(format!("{b1_str}_{b2_str}"))
}
