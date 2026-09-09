use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use ::thal_rs::mfe::{Dimer, Monomer};
use ::thal_rs::thermo::params::ThermoParams;

/// Computes the minimum free energy of a self-dimer / heterodimer.
///
/// Returns a `(dg, dh, ds)` tuple.
#[pyfunction]
fn dimer(seq1: &str, seq2: &str) -> PyResult<(f64, f64, f64)> {
    let params = ThermoParams::with_defaults().map_err(PyValueError::new_err)?;
    let mut mfe =
        Dimer::new(seq1.as_bytes(), seq2.as_bytes(), &params).map_err(PyValueError::new_err)?;
    let res = mfe.calculate().map_err(PyValueError::new_err)?;
    Ok((res.dg(), res.dh, res.ds))
}

/// Computes the minimum free energy of a monomer secondary structure.
///
/// Returns a `(dg, dh, ds)` tuple.
#[pyfunction]
fn monomer(seq: &str) -> PyResult<(f64, f64, f64)> {
    let params = ThermoParams::with_defaults().map_err(PyValueError::new_err)?;
    let mut mfe = Monomer::new(seq.as_bytes(), &params).map_err(PyValueError::new_err)?;
    let res = mfe.calculate().map_err(PyValueError::new_err)?;
    Ok((res.dg(), res.dh, res.ds))
}

#[pymodule]
fn thal_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dimer, m)?)?;
    m.add_function(wrap_pyfunction!(monomer, m)?)?;
    Ok(())
}
