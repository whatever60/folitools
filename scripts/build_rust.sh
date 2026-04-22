#!/usr/bin/env bash
# Local dev helper: build folitools' Rust binary (foli_add_tags) and drop it
# into the active environment's bin/. For release/PyPI builds this isn't
# needed — `pip install folitools` ships the binary in the wheel (built via
# maturin in CI).
#
# Usage:
#   bash scripts/build_rust.sh              # installs into $CONDA_PREFIX/bin
#   PREFIX=/somewhere/bin bash scripts/build_rust.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RUST_DIR="$REPO_ROOT/rust"

if [[ ! -d "$RUST_DIR" ]]; then
    echo "error: no rust/ directory at $RUST_DIR" >&2
    exit 1
fi

if ! command -v cargo >/dev/null 2>&1; then
    echo "error: cargo not on PATH. Install Rust (e.g. via the conda env)." >&2
    exit 1
fi

DEST_BIN="${PREFIX:-${CONDA_PREFIX:-}}"
if [[ -z "$DEST_BIN" ]]; then
    echo "error: neither PREFIX nor CONDA_PREFIX is set; pass PREFIX=..." >&2
    exit 1
fi
if [[ "$(basename "$DEST_BIN")" != "bin" ]]; then
    DEST_BIN="$DEST_BIN/bin"
fi
mkdir -p "$DEST_BIN"

echo "Building rust binary (release mode)..."
cd "$RUST_DIR"
cargo build --release --bin foli_add_tags

SRC="$RUST_DIR/target/release/foli_add_tags"
DST="$DEST_BIN/foli_add_tags"
install -m 0755 "$SRC" "$DST"
echo "Installed: $DST"

echo "Done. Verify with: foli_add_tags --help"
