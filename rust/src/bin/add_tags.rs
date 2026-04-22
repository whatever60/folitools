//! Thin CLI entry point. All logic lives in the `add_tags` library module so
//! the pyo3 extension can share it. For a local rebuild:
//!   cargo build --release --no-default-features --bin foli_add_tags

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    std::process::exit(folitools_rust::add_tags::run_cli(argv));
}
