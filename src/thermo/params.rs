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
}
