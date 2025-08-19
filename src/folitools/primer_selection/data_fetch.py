import hashlib
import os
import sys
import tarfile
import urllib.request
from pathlib import Path


# Prefer platformdirs if you already depend on it; otherwise keep this tiny shim.
def _default_cache_root() -> Path:
    """Return a cross-platform cache directory for folitools."""
    # Respect an override first.
    env = os.getenv("FOLITOOLS_DATA_DIR")
    if env:
        p = Path(env).expanduser()
        p.mkdir(parents=True, exist_ok=True)
        return p

    # Minimal platform-aware default without extra deps.
    if sys.platform == "win32":
        base = Path(os.getenv("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        base = Path(os.getenv("XDG_CACHE_HOME", Path.home() / ".cache"))
    target = base / "folitools" / "primer_selection"
    target.mkdir(parents=True, exist_ok=True)
    return target


def _sha256sum(path: Path) -> str:
    """Compute sha256 for a file."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_dataset(
    *,
    version: str,
    url: str,
    sha256: str,
    archive_member: str | None = None,
) -> Path:
    """Ensure the primer-selection dataset is available locally.

    Args:
        version: Semantic or date version tag used for on-disk names.
        url: HTTPS URL to a release artifact (e.g., .tar.zst/.tar.gz/.zip).
        sha256: Expected sha256 of the downloaded artifact.
        archive_member: Optional single member or top-level directory inside the
            archive to extract; if None, extracts the whole archive.

    Returns:
        Path to the dataset directory on the local filesystem.

    Raises:
        RuntimeError: If checksum verification fails or extraction fails.
    """
    cache_root = _default_cache_root()
    version_dir = cache_root / version
    marker = version_dir / ".complete"

    if marker.exists():
        return version_dir

    # Download
    dl = cache_root / f"primer_data-{version}{Path(url).suffix}"
    print(f"Downloading {url} to {dl}")
    if not dl.exists():
        with urllib.request.urlopen(url) as resp, dl.open("wb") as out:
            # stream download
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)

    # Verify
    digest = _sha256sum(dl)
    if digest != sha256:
        dl.unlink(missing_ok=True)
        raise RuntimeError(f"Checksum mismatch for {dl.name}: {digest} != {sha256}")

    # Extract
    version_dir.mkdir(parents=True, exist_ok=True)
    if tarfile.is_tarfile(dl):
        with tarfile.open(dl, "r:*") as t:  # auto-detect compression
            if archive_member:
                t.extract(archive_member, path=version_dir)
                # Flatten if needed
                candidate = version_dir / archive_member
                if candidate.is_dir():
                    pass
            else:
                t.extractall(path=version_dir)
    else:
        # If the artifact is just a directory dump or single file, move it.
        # Extend here for .zip if you use zipfile.
        pass

    marker.touch()
    return version_dir
