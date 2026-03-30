use crate::{mfe::Mfe, thermo::params::ThermoParams};

mod matrix;
mod mfe;
mod seq;
mod thermo;

fn main() -> Result<(), &'static str> {
    let params = ThermoParams::with_file_path("thermo")?;
    let mut mfe = Mfe::new(vec![b'A'], vec![b'B'], params);
    mfe.calculate()?;

    Ok(())
}
