"""Version-stamping helpers for foli-primer outputs.

Every foli-primer subcommand writes one or more files; each is stamped
with the folitools version that produced it so a ``grep folitools``
(text) or Excel "File Properties" inspection (xlsx) tells you which
run wrote any given file. The version is pulled from
``folitools.__version__`` at write time so release bumps flow through
automatically.

Stamp surface per format:

- TSV/CSV: leading ``# folitools <version>`` line. Internal readers in
  this package pass ``comment='#'`` so the marker is skipped on parse;
  external consumers using strict TSV/CSV parsers may need the same.
- xlsx: workbook ``properties.description`` core property. Visible
  through Excel's File Properties dialog and ``openpyxl.load_workbook``.
- FASTA: intentionally **not** stamped. cutadapt's dnaio reader rejects
  ``;`` comment lines and Biopython's default ``fasta`` reader has
  deprecated ``#`` ones, so any in-file marker risks breaking
  downstream consumers (``foli assign-probes`` consumes
  ``foli-primer recover``'s i5/i7 FASTAs through cutadapt). The sibling
  ``recover.log`` is already version-stamped, which covers
  reproducibility audits for these files.
"""

from pathlib import Path
from typing import Any

import pandas as pd

from .. import __version__


def _comment(prefix: str) -> str:
    return f"{prefix} folitools {__version__}\n"


def write_versioned_tsv(
    df: pd.DataFrame, path: str | Path, **kwargs: Any
) -> None:
    """``df.to_csv(path, sep='\\t', index=False, ...)`` with a leading version comment.

    ``mode`` is forced to ``'a'`` (append) so the version line written
    just before stays put.
    """
    kwargs.setdefault("sep", "\t")
    kwargs.setdefault("index", False)
    kwargs["mode"] = "a"
    with open(path, "w") as fh:
        fh.write(_comment("#"))
    df.to_csv(path, **kwargs)


def write_versioned_csv(
    df: pd.DataFrame, path: str | Path, **kwargs: Any
) -> None:
    """``df.to_csv(path, index=False, ...)`` with a leading version comment."""
    kwargs.setdefault("index", False)
    kwargs["mode"] = "a"
    with open(path, "w") as fh:
        fh.write(_comment("#"))
    df.to_csv(path, **kwargs)


def write_versioned_excel(
    df: pd.DataFrame, path: str | Path, **kwargs: Any
) -> None:
    """``df.to_excel(path, ...)`` followed by stamping workbook properties."""
    df.to_excel(path, **kwargs)
    stamp_excel_version(path)


def stamp_excel_version(path: str | Path) -> None:
    """Set the workbook's ``properties.description`` to the folitools version.

    Used after pandas writes an .xlsx so the file carries the version
    in a place users can inspect via Excel → File → Info → Properties.
    """
    from openpyxl import load_workbook

    wb = load_workbook(path)
    wb.properties.description = f"folitools {__version__}"
    wb.save(path)


