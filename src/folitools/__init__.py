"""Top-level package metadata for folitools."""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import tomllib


def _resolve_version() -> str:
    """Return the installed package version or the source-tree project version."""
    try:
        return version("folitools")
    except PackageNotFoundError:
        pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
        with pyproject.open("rb") as f:
            return tomllib.load(f)["project"]["version"]


__version__ = _resolve_version()
