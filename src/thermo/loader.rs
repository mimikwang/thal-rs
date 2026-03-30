use core::f64;
use std::{
    collections::HashMap,
    fs,
    path::{Path, PathBuf},
};

use crate::thermo::Thermo;

const COMMENT_PREFIX: &'static str = "#";
const STACK_DH: &'static str = "stack.dh";
const STACH_DS: &'static str = "stack.ds";
const STACKMM_DH: &'static str = "stackmm.dh";
const STACKMM_DS: &'static str = "stackmm.ds";
const TSTACK_DH: &'static str = "tstack.dh";
const TSTACK_DS: &'static str = "tstack.ds";
const DANGLE_DH: &'static str = "dangle.dh";
const DANGLE_DS: &'static str = "dangle.ds";
const LOOPS_DH: &'static str = "loops.dh";
const LOOPS_DS: &'static str = "loops.ds";
const TETRALOOP_DH: &'static str = "tetraloop.dh";
const TETRALOOP_DS: &'static str = "tetraloop.ds";

pub fn load_stack(file_path: &str) -> Result<HashMap<String, Thermo>, &'static str> {
    // Check to make sure paths exist
    let dh_path = Path::new(file_path).join(STACK_DH);
    let ds_path = Path::new(file_path).join(STACH_DS);
    let dhmm_path = Path::new(file_path).join(STACKMM_DH);
    let dsmm_path = Path::new(file_path).join(STACKMM_DS);

    if !dh_path.exists() || !ds_path.exists() || !dhmm_path.exists() || !dsmm_path.exists() {
        return Err("paths not found");
    }

    // Load files
    let mut output = HashMap::new();
    read_data(&mut output, dh_path, true)?;
    read_data(&mut output, dhmm_path, true)?;
    read_data(&mut output, ds_path, false)?;
    read_data(&mut output, dsmm_path, false)?;

    Ok(output)
}

pub fn load_tstack(file_path: &str) -> Result<HashMap<String, Thermo>, &'static str> {
    // Check to make sure paths exist
    let dh_path = Path::new(file_path).join(TSTACK_DH);
    let ds_path = Path::new(file_path).join(TSTACK_DS);

    if !dh_path.exists() || !ds_path.exists() {
        return Err("paths not found");
    }

    // Load files
    let mut output = HashMap::new();
    read_data(&mut output, dh_path, true)?;
    read_data(&mut output, ds_path, false)?;

    Ok(output)
}

pub fn load_dangle(file_path: &str) -> Result<HashMap<String, Thermo>, &'static str> {
    // Check to make sure paths exist
    let dh_path = Path::new(file_path).join(DANGLE_DH);
    let ds_path = Path::new(file_path).join(DANGLE_DS);

    if !dh_path.exists() || !ds_path.exists() {
        return Err("paths not found");
    }

    // Load files
    let mut output = HashMap::new();
    read_data(&mut output, dh_path, true)?;
    read_data(&mut output, ds_path, false)?;

    Ok(output)
}

pub fn load_tetraloop(file_path: &str) -> Result<HashMap<String, Thermo>, &'static str> {
    // Check to make sure paths exist
    let dh_path = Path::new(file_path).join(TETRALOOP_DH);
    let ds_path = Path::new(file_path).join(TETRALOOP_DS);

    if !dh_path.exists() || !ds_path.exists() {
        return Err("paths not found");
    }

    // Load files
    let mut output = HashMap::new();
    read_data(&mut output, dh_path, true)?;
    read_data(&mut output, ds_path, false)?;

    Ok(output)
}

pub fn load_loops(file_path: &str) -> Result<Loops, &'static str> {
    // Check to make sure paths exist
    let dh_path = Path::new(file_path).join(LOOPS_DH);
    let ds_path = Path::new(file_path).join(LOOPS_DS);

    if !dh_path.exists() || !ds_path.exists() {
        return Err("paths not found");
    }

    // Load files
    let mut loops = Loops::default();
    read_loop_data(&mut loops, dh_path, true)?;
    read_loop_data(&mut loops, ds_path, false)?;

    Ok(loops)
}

fn parse_line(line: &str) -> Result<(String, f64), &'static str> {
    let split = line.split(&[' ', '\t']).collect::<Vec<&str>>();
    if split.len() != 2 {
        return Err("invalid line");
    }
    Ok((split[0].to_owned(), parse_value(split[1])?))
}

fn read_data(
    output: &mut HashMap<String, Thermo>,
    path: PathBuf,
    is_enthalpy: bool,
) -> Result<(), &'static str> {
    let data = fs::read_to_string(path).map_err(|_| "error reading file")?;
    for val in data.split("\n") {
        // Skip comments and empty lines
        if val.starts_with(COMMENT_PREFIX) || val.is_empty() {
            continue;
        }

        let (key, val) = parse_line(val)?;
        if output.get(&key).is_none() {
            if is_enthalpy {
                output.insert(key, Thermo::with_values(0.0, val));
            } else {
                output.insert(key, Thermo::with_values(val, 0.0));
            }
        } else {
            if let Some(t) = output.get_mut(&key) {
                if is_enthalpy {
                    t.dh = val;
                } else {
                    t.ds = val;
                }
            }
        }
    }
    Ok(())
}

struct Loop {
    internal: f64,
    bulge: f64,
    hairpin: f64,
}

#[derive(Default, Debug)]
pub struct Loops {
    internal: [Thermo; 30],
    bulge: [Thermo; 30],
    hairpin: [Thermo; 30],
}

fn parse_loop_line(line: &str) -> Result<Loop, &'static str> {
    let split = line.split(&[' ', '\t']).collect::<Vec<&str>>();
    if split.len() != 4 {
        return Err("invalid line");
    }
    Ok(Loop {
        internal: parse_value(split[1])?,
        bulge: parse_value(split[2])?,
        hairpin: parse_value(split[3])?,
    })
}

fn read_loop_data(
    output: &mut Loops,
    path: PathBuf,
    is_enthalpy: bool,
) -> Result<(), &'static str> {
    let data = fs::read_to_string(path).map_err(|_| "error reading file")?;
    let mut loop_size = 0;
    for val in data.split("\n") {
        // Skip comments and empty lines
        if val.starts_with(COMMENT_PREFIX) || val.is_empty() {
            continue;
        }

        let loop_data = parse_loop_line(val)?;
        if is_enthalpy {
            output.internal[loop_size].dh = loop_data.internal;
            output.bulge[loop_size].dh = loop_data.bulge;
            output.hairpin[loop_size].dh = loop_data.hairpin;
        } else {
            output.internal[loop_size].ds = loop_data.internal;
            output.bulge[loop_size].ds = loop_data.bulge;
            output.hairpin[loop_size].ds = loop_data.hairpin;
        }
        loop_size += 1;
    }

    Ok(())
}

fn parse_value(val: &str) -> Result<f64, &'static str> {
    if val == "inf" {
        return Ok(f64::INFINITY);
    }

    Ok(val.parse::<f64>().map_err(|_| "value should be a number")?)
}
