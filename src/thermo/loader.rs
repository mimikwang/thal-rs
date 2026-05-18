use core::f64;
use std::{
    collections::HashMap,
    fs,
    path::{Path, PathBuf},
};

use crate::{errors::Result, thermo::Thermo};

const COMMENT_PREFIX: &str = "#";
const STACK_DH: &str = "stack.dh";
const STACK_DS: &str = "stack.ds";
const STACKMM_DH: &str = "stackmm.dh";
const STACKMM_DS: &str = "stackmm.ds";
const TSTACK_DH: &str = "tstack.dh";
const TSTACK_DS: &str = "tstack.ds";
const DANGLE_DH: &str = "dangle.dh";
const DANGLE_DS: &str = "dangle.ds";
const LOOPS_DH: &str = "loops.dh";
const LOOPS_DS: &str = "loops.ds";
const TETRALOOP_DH: &str = "tetraloop.dh";
const TETRALOOP_DS: &str = "tetraloop.ds";

macro_rules! data_fn {
    ($name:ident, $file_name_ds:ident, $file_name_dh:ident) => {
        pub fn $name(file_path: &str) -> Result<HashMap<String, Thermo>> {
            // Check ot make sure path exists
            let dh_path = Path::new(file_path).join($file_name_dh);
            let ds_path = Path::new(file_path).join($file_name_ds);

            if !dh_path.exists() || !ds_path.exists() {
                return Err("paths not found");
            }

            let mut output = HashMap::new();
            read_data(&mut output, dh_path, true)?;
            read_data(&mut output, ds_path, false)?;

            Ok(output)
        }
    };
}

data_fn!(load_stack, STACK_DS, STACK_DH);
data_fn!(load_stack_mm, STACKMM_DS, STACKMM_DH);
data_fn!(load_tstack, TSTACK_DS, TSTACK_DH);
data_fn!(load_dangle, DANGLE_DS, DANGLE_DH);
data_fn!(load_tetraloop, TETRALOOP_DS, TETRALOOP_DH);

pub fn load_loops(file_path: &str) -> Result<Loops> {
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

fn parse_line(line: &str) -> Result<(String, f64)> {
    let split = line.split(&[' ', '\t']).collect::<Vec<&str>>();
    if split.len() != 2 {
        return Err("invalid line");
    }
    Ok((split[0].to_owned(), parse_value(split[1])?))
}

fn read_data(output: &mut HashMap<String, Thermo>, path: PathBuf, is_enthalpy: bool) -> Result<()> {
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
    pub internal: [Thermo; 30],
    pub bulge: [Thermo; 30],
    pub hairpin: [Thermo; 30],
}

fn parse_loop_line(line: &str) -> Result<Loop> {
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

fn read_loop_data(output: &mut Loops, path: PathBuf, is_enthalpy: bool) -> Result<()> {
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

fn parse_value(val: &str) -> Result<f64> {
    if val == "inf" {
        return Ok(f64::INFINITY);
    }

    val.parse::<f64>().map_err(|_| "value should be a number")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_path() -> String {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("thermo");
        path.to_str().expect("path error").to_owned()
    }

    #[test]
    fn test_load_stack() {
        let loaded = load_stack(&get_path());
        assert!(loaded.is_ok());

        let loaded = loaded.unwrap();
        assert_eq!(
            loaded.get("AA_TT"),
            Some(&Thermo {
                ds: -22.2,
                dh: -7900.0
            })
        );
        assert_eq!(
            loaded.get("AA_AA"),
            Some(&Thermo {
                ds: f64::INFINITY,
                dh: f64::INFINITY
            })
        )
    }
}
