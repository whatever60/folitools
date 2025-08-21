import os
import shutil
from importlib.metadata import PackageNotFoundError, version
from importlib.resources import as_file, files
from pathlib import Path

from .data_fetch import _default_cache_root, ensure_dataset


# Dataset metadata (bump when you publish a new artifact)
_DATASET_VERSION = "2025.08.18"
_ARCHIVE_TOP = f"primer_data-{_DATASET_VERSION}"
_DATASET_URL = (
    "https://github.com/whatever60/folitools/releases/download/"
    f"v{_DATASET_VERSION}/{_ARCHIVE_TOP}.tar.zst"
)
_DATASET_SHA256 = "a043e19e752ba49dd6ce2028e353650c9233c262018a46b44cd928367a8d68e9"


def _materialize_packaged_base() -> Path | None:
    """Return a filesystem path to packaged ``data/`` or copy once into cache.

    Behavior:
      - If packaged data is filesystem-backed, return it directly (no copy).
      - If packaged data is zip-backed, copy it **once** into
        ``<cache_root>/packaged/`` and reuse thereafter.
        Reuse is keyed by a simple fingerprint written to ``.source``:
        ``packaged:<distribution_version_or_unknown>``.
        If the distribution version changes, we refresh the copy.

    Returns:
        Path to the packaged data directory, or ``None`` if no packaged data exists.
    """
    trav = files(__package__).joinpath("data")
    if not trav.exists():
        return None

    # Try to use it directly if it's already on the filesystem.
    try:
        fs_path = Path(trav)
        if fs_path.is_dir():
            return fs_path
    except TypeError:
        fs_path = None  # zip-backed; must materialize

    # Zip-backed: copy once into cache/packaged and reuse via a marker.
    cache_root = _default_cache_root()
    packaged_cache = cache_root / "packaged"
    marker = packaged_cache / ".complete"
    source_meta = packaged_cache / ".source"

    try:
        dist_ver = version("folitools")
    except PackageNotFoundError:
        dist_ver = "unknown"
    fingerprint = f"packaged:{dist_ver}"

    if packaged_cache.is_dir() and marker.exists():
        try:
            if source_meta.read_text(encoding="utf-8").strip() == fingerprint:
                return packaged_cache  # up-to-date; skip copying
        except FileNotFoundError:
            pass

    packaged_cache.mkdir(parents=True, exist_ok=True)
    with as_file(trav) as tmp:
        src_dir = Path(tmp)
        shutil.copytree(src_dir, packaged_cache, dirs_exist_ok=True)
    source_meta.write_text(fingerprint + "\n", encoding="utf-8")
    marker.touch()
    return packaged_cache


def _data_dir_for_species(species: str) -> Path:
    """Resolve the data directory for a species with dev → packaged → cached fallback.

    Precedence:
      1) **Development tree**: ``.../primer_selection/data/<species>`` (if present)
      2) **Packaged data**:
         - If installed package contains ``data/`` on the real filesystem, use it directly
           (no copying).
         - If the package is **zip-backed**, copy the whole packaged ``data/`` **once**
           into ``<cache_root>/packaged/`` and read from there on subsequent runs.
      3) **Downloaded dataset** via :func:`ensure_dataset` (checksum-verified, cached).

    Platform save locations (cache root):
      - **If ``FOLITOOLS_DATA_DIR`` is set**: use that directory directly.
        Example: ``$FOLITOOLS_DATA_DIR/packaged/<species>/`` or
        ``$FOLITOOLS_DATA_DIR/<version>/<species>/``.
      - **Linux** (no env override): ``~/.cache/folitools/primer_selection``.
      - **macOS** (no env override): ``~/Library/Caches/folitools/primer_selection``.
      - **Windows** (no env override): ``%LOCALAPPDATA%\\folitools\\primer_selection``
        (fallback: ``%USERPROFILE%\\AppData\\Local\\folitools\\primer_selection``).

    Cache behavior:
      - Packaged, zip-backed data is copied once into ``<cache_root>/packaged/`` and
        a marker file ``.complete`` plus a fingerprint ``.source`` are written.
        If the installed distribution version changes (or the cache is deleted),
        the copy is refreshed; otherwise it is **not** re-copied.
      - Downloaded datasets are stored under ``<cache_root>/<version>/`` and gated by
        checksum + a ``.complete`` marker set by :func:`ensure_dataset`.

    Args:
        species: Species key (e.g., ``"mouse"``, ``"human"``); a simple name only.

    Returns:
        Path to a directory containing data for the specified species.

    Raises:
        ValueError: If ``species`` is empty, ``.``/``..``, or contains path separators.
        FileNotFoundError: If the resolved directory for ``species`` does not exist.
    """
    # Validate species name strictly (no traversal)
    name = species.strip()
    if not name or name in {".", ".."}:
        raise ValueError(f"Invalid species name: {species!r}")
    if os.sep in name or (os.altsep and os.altsep in name):
        raise ValueError(
            f"Invalid species name: {species!r}. Must not contain path separators."
        )

    # 1) Development checkout path
    dev_base = Path(__file__).parent / "data"
    if dev_base.is_dir():
        candidate = dev_base / name
        if not candidate.is_dir():
            raise FileNotFoundError(
                f"Data directory for species {name!r} not found at: {candidate}"
            )
        return candidate

    # 2) Packaged data (direct or materialized once into cache/packaged)
    pkg_base = _materialize_packaged_base()
    if pkg_base is not None:
        candidate = pkg_base / name
        if candidate.is_dir():
            return candidate
        # If packaged exists but species missing, fall through to downloader.

    # 3) Downloaded, checksum-verified dataset (cached under <cache_root>/<version>/)
    base = ensure_dataset(
        version=_DATASET_VERSION,
        url=_DATASET_URL,
        sha256=_DATASET_SHA256,
        archive_member=_ARCHIVE_TOP,
    )
    candidate = base / name
    if not candidate.is_dir():
        raise FileNotFoundError(
            f"Data directory for species {name!r} not found at: {candidate}"
        )
    return candidate
