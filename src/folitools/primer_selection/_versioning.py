"""Version-stamping helpers for foli-primer outputs.

Each foli-primer subcommand gets at most one version stamp, placed
where it doesn't disturb the file's primary content:

- ``foli-primer summary``: stamps the workbook's
  ``properties.description`` core property on its xlsx output(s).
  Visible via Excel's "File → Info → Properties" dialog or
  ``openpyxl.load_workbook(...).properties.description``.
- ``foli-primer recover``: writes the version into ``recover.log`` via
  a single ``logger.info`` line at the top — no in-file marker on the
  CSVs/xlsx/PDF/FASTA outputs.
- ``foli-primer subset`` / ``saddle`` / ``sgad`` / ``product``: no
  natural per-command stamp surface; outputs are uncluttered TSVs
  (non-matrix, so the matrix-style ``index_label`` trick from
  ``foli get-count-mtx`` doesn't apply) or a FASTA. Skipped
  intentionally — we'd rather have no stamp than a stamp that
  breaks downstream parsers.

The version is pulled from ``folitools.__version__`` at write time so
release bumps flow through automatically.
"""

from pathlib import Path
from typing import Any

import pandas as pd

from .. import __version__


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
