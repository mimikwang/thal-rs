use std::env;

use thal_rs::{
    errors::Result,
    mfe::{Dimer, Monomer},
    thermo::params::ThermoParams,
};

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        return Err("invalid args");
    }

    if &args[1] == "dimer" {
        if args.len() != 4 {
            return Err("dimer requires 2 seqs");
        }

        let seq1 = &args[2];
        let seq2 = &args[3];

        let params = ThermoParams::with_defaults()?;
        let mut mfe = Dimer::new(seq1.as_bytes(), seq2.as_bytes(), &params)?;
        let res = mfe.calculate()?;
        println!("dg={:?} dh={:?} ds={:?}", res.dg(), res.dh, res.ds);

        return Ok(());
    }

    if &args[1] == "monomer" {
        if args.len() != 3 {
            return Err("monomer requires 1 seq");
        }

        let seq = &args[2];
        let params = ThermoParams::with_defaults()?;
        let mut mfe = Monomer::new(seq.as_bytes(), &params)?;
        let res = mfe.calculate()?;
        println!("dg={:?} dh={:?} ds={:?}", res.dg(), res.dh, res.ds);

        return Ok(());
    }

    Err("invalid args")
}
