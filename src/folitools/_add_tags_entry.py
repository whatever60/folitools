"""Thin CLI wrapper around the Rust-backed ``folitools._rust_native``.

Installed as the ``foli_add_tags`` console script via ``[project.scripts]``.
The Rust side owns arg parsing (via clap) and all BAM I/O; this file only
forwards ``sys.argv`` and propagates the exit code.
"""

from __future__ import annotations

import sys


def cli() -> int:
    # Imported inside the function so a broken extension import doesn't break
    # `python -c 'import folitools'` or any unrelated CLI.
    from folitools._rust_native import run_add_tags

    return int(run_add_tags(sys.argv))


if __name__ == "__main__":
    sys.exit(cli())
