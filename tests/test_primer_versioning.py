"""Roundtrip test for foli-primer's xlsx version-stamping helper."""

from pathlib import Path

import pandas as pd

from folitools import __version__
from folitools.primer_selection._versioning import write_versioned_excel


def test_write_versioned_excel_stamps_workbook_description(tmp_path: Path) -> None:
    df = pd.DataFrame({"col": [1, 2, 3]})
    out = tmp_path / "x.xlsx"

    write_versioned_excel(df, out, index=False)

    from openpyxl import load_workbook

    wb = load_workbook(out)
    assert wb.properties.description == f"folitools {__version__}"
    # Sheet content survives the re-save by the stamping step.
    parsed = pd.read_excel(out)
    pd.testing.assert_frame_equal(parsed, df)
