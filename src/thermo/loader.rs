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

fn parse_value(line: &str) -> Result<Option<(String, f64)>, &'static str> {
    let split = line.split(" ").collect::<Vec<&str>>();
    if split.len() != 2 {
        return Err("invalid line");
    }
    if split[1] == "inf" {
        return Ok(None);
    }

    let value = split[1]
        .parse::<f64>()
        .map_err(|_| "value should be a number")?;

    Ok(Some((split[0].to_owned(), value)))
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

        if let Some((key, val)) = parse_value(val)? {
            if output.get(&key).is_none() {
                if is_enthalpy {
                    output.insert(key, Thermo::with_values(0.0, val));
                } else {
                    output.insert(key, Thermo::with_values(val, 0.0));
                }
            } else {
                if is_enthalpy {
                    if let Some(t) = output.get_mut(&key) {
                        t.dh = val;
                    }
                } else {
                    if let Some(t) = output.get_mut(&key) {
                        t.ds = val;
                    }
                }
            }
        };
    }
    Ok(())
}
