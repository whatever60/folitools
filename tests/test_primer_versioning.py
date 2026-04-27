"""Roundtrip tests for foli-primer's version-stamping helpers."""

from pathlib import Path

import pandas as pd

from folitools import __version__
from folitools.primer_selection._versioning import (
    command_log,
    default_command_log_path,
    write_versioned_excel,
)


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


def test_default_command_log_path_uses_output_dirname(tmp_path: Path) -> None:
    # Sibling of the primary output, named after the command.
    primary = tmp_path / "deep" / "primer_info.tsv"
    log = default_command_log_path(primary, "subset")
    assert log == tmp_path / "deep" / "subset.log"


def test_command_log_writes_version_first_then_user_lines(tmp_path: Path) -> None:
    log = tmp_path / "subset.log"

    with command_log("test_command_log_v1", log) as logger:
        logger.info("species: mouse")
        logger.info("target genes: 42")

    body = log.read_text()
    # First non-blank field after the timestamp is the version anchor.
    first_line = body.splitlines()[0]
    assert f"folitools version: {__version__}" in first_line
    assert "species: mouse" in body
    assert "target genes: 42" in body


def test_command_log_truncates_prior_run(tmp_path: Path) -> None:
    log = tmp_path / "saddle.log"
    log.write_text("stale\nstale\nstale\n")

    with command_log("test_command_log_v2", log) as logger:
        logger.info("fresh line")

    body = log.read_text()
    assert "stale" not in body
    assert "fresh line" in body


def test_command_log_none_path_is_a_noop(tmp_path: Path) -> None:
    # Passing None means "no log" — yielded logger is still usable but
    # writes go nowhere on disk.
    with command_log("test_command_log_v3", None) as logger:
        logger.info("disappears")

    # Sanity: no log file got created in tmp_path as a side effect.
    assert list(tmp_path.iterdir()) == []
