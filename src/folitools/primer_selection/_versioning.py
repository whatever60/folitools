"""Version-stamping helpers for foli-primer outputs.

Each foli-primer subcommand gets at most one version stamp, placed
where it doesn't disturb the file's primary content:

- ``foli-primer summary``: stamps the workbook's
  ``properties.description`` core property on its xlsx output(s).
- ``foli-primer subset`` / ``saddle`` / ``sgad`` / ``product`` /
  ``recover``: each writes a sibling ``.log`` file (path defaults to
  ``<output-dir>/<command>.log``, overridable via ``--log``). The log's
  first line is ``folitools version: <version>``, followed by the
  command's input/output paths and key parameters.

The version is pulled from ``folitools.__version__`` at write time so
release bumps flow through automatically.

Outputs that aren't stamped (intentionally):

- TSV/CSV: the non-matrix shape has no empty top-left cell to repurpose
  the way ``foli get-count-mtx`` does via ``index_label``.
- FASTA: cutadapt's dnaio reader rejects ``;`` comment lines and
  Biopython's ``fasta`` reader has deprecated ``#`` ones; either
  marker would break a downstream consumer. Reproducibility for these
  is covered by the sibling per-command log instead.
"""

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

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


def default_command_log_path(
    primary_output: str | Path, command_name: str
) -> Path:
    """Default sibling log path for a foli-primer subcommand.

    ``<dirname(primary_output)>/<command_name>.log``. Used when the user
    doesn't pass ``--log``: pipelining several foli-primer subcommands
    into the same output directory then yields one log per command,
    side-by-side with the data outputs.
    """
    return Path(primary_output).parent / f"{command_name}.log"


@contextmanager
def command_log(
    name: str, log_path: str | Path | None
) -> Iterator[logging.Logger]:
    """Open a per-command log file and yield a logger ready to write to it.

    On entry: the file is truncated (so a re-run starts clean), the
    folitools version is written as the first line, and the caller can
    record any subsequent ``logger.info(...)`` line into the same file.
    On exit: the handler is removed and closed so subsequent calls with
    a different ``log_path`` don't accidentally inherit it.

    Args:
        name: Logger name — pass ``__name__`` from the caller so each
            subcommand's logs are independently filterable in tests.
        log_path: Where to write. ``None`` short-circuits to a no-op
            (logger yielded with no file handler attached) so callers
            can opt out without branching.
    """
    logger = logging.getLogger(name)
    handler: logging.FileHandler | None = None

    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        # Truncate any prior log so cross-run noise doesn't accumulate;
        # mirrors recover()'s existing behavior.
        if log_path.exists():
            log_path.unlink()
        handler = logging.FileHandler(log_path)
        handler.setLevel(logging.INFO)
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        # Reproducibility anchor — every command log starts with this.
        logger.info(f"folitools version: {__version__}")

    try:
        yield logger
    finally:
        if handler is not None:
            logger.removeHandler(handler)
            handler.close()
