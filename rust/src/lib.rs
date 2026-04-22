//! Rust-backed internals for folitools.
//!
//! The `python-extension` feature (default) exposes this crate as a Python
//! extension module named `folitools._rust_native`. The Python wrapper at
//! `src/folitools/_add_tags_entry.py` forwards `argv` into `run_add_tags`
//! below. Building with `--no-default-features` disables the pyo3 glue so
//! the standalone `foli_add_tags` binary can be built without the Python
//! dependency (used by `scripts/build_rust.sh` for local dev).

pub mod add_tags;

#[cfg(feature = "python-extension")]
mod python {
    use super::*;
    use pyo3::prelude::*;

    /// Run `foli_add_tags` with the given argv (argv[0] is the program name).
    /// Returns the process exit code the caller should propagate.
    #[pyfunction]
    #[pyo3(name = "run_add_tags")]
    fn run_add_tags(py: Python<'_>, argv: Vec<String>) -> i32 {
        // Release the GIL while the BAM processing runs; all I/O and htslib
        // calls happen entirely in Rust.
        py.allow_threads(|| add_tags::run_cli(argv))
    }

    #[pymodule]
    fn _rust_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(run_add_tags, m)?)?;
        Ok(())
    }
}
