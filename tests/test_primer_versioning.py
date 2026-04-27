"""Roundtrip tests for foli-primer's version-stamping helpers."""

from pathlib import Path

import pandas as pd

from folitools import __version__
from folitools.primer_selection._versioning import (
    write_versioned_csv,
    write_versioned_excel,
    write_versioned_tsv,
)


def test_write_versioned_tsv_roundtrips_with_comment_skip(tmp_path: Path) -> None:
    df = pd.DataFrame({"gene": ["A", "B"], "primer_fwd": ["ACG", "TGC"]})
    out = tmp_path / "x.tsv"

    write_versioned_tsv(df, out)

    text = out.read_text()
    assert text.startswith(f"# folitools {__version__}\n"), text[:80]

    # Internal readers across foli-primer pass comment='#'; the marker is
    # transparently skipped so downstream parses the data only.
    parsed = pd.read_csv(out, sep="\t", comment="#")
    pd.testing.assert_frame_equal(parsed, df)


def test_write_versioned_csv_roundtrips_with_comment_skip(tmp_path: Path) -> None:
    df = pd.DataFrame({"a": [1, 2], "b": [3.0, 4.0]})
    out = tmp_path / "x.csv"

    write_versioned_csv(df, out)
    parsed = pd.read_csv(out, comment="#")
    pd.testing.assert_frame_equal(parsed, df)


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


