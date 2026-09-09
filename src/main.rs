use std::env;

use crate::{errors::Result, mfe::Dimer, thermo::params::ThermoParams};

mod errors;
mod matrix;
mod mfe;
mod seq;
mod thal;
mod thermo;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        return Err("invalid args");
    }

    if &args[1] == "dimer" {
        let seq1 = &args[2];
        let seq2 = &args[3];

        let params = ThermoParams::with_file_path("thermo")?;
        let mut mfe = Dimer::new(seq1.as_bytes(), seq2.as_bytes(), &params)?;
        println!("{:?}", mfe.calculate()?);

        return Ok(());
    }

    Err("invalid args")
}
