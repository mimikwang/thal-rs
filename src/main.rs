use crate::{mfe::Mfe, thermo::params::ThermoParams};

mod matrix;
mod mfe;
mod seq;
mod thermo;

fn main() -> Result<(), &'static str> {
    let params = ThermoParams::with_file_path("thermo")?;
    let mut mfe = Mfe::new(
        vec![b'A', b'G', b'C', b'T'],
        vec![b'T', b'C', b'G', b'A'],
        params,
    );
    mfe.calculate()?;

    Ok(())
}
