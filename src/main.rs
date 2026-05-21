use crate::{errors::Result, mfe::Dimer, thermo::params::ThermoParams};

mod errors;
mod matrix;
mod mfe;
mod seq;
mod thal;
mod thermo;

fn main() -> Result<()> {
    let params = ThermoParams::with_file_path("thermo")?;
    let mut mfe = Dimer::new(b"AGCT", b"TCGA", &params)?;
    println!("{:?}", mfe.calculate()?);

    Ok(())
}
