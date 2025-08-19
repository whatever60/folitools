#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path
from typing import Final


def sha256_file(path: Path, *, chunk_mb: int = 8) -> str:
    """Compute the SHA-256 hex digest of a file by streaming.

    Args:
        path: Path to the file to hash.
        chunk_mb: Read size in MiB per chunk (tunable for very large files).

    Returns:
        The lowercase hexadecimal SHA-256 digest string.

    Raises:
        FileNotFoundError: If the file does not exist.
        PermissionError: If the file cannot be read.
    """
    h = hashlib.sha256()
    chunk_size: Final[int] = chunk_mb * 1024 * 1024
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute (and optionally verify) SHA-256 of an artifact."
    )
    parser.add_argument("file", type=Path, help="Path to the artifact (e.g., .tar.zst)")
    parser.add_argument(
        "--expect",
        type=str,
        default=None,
        help="Expected SHA-256 hex (verify mode). Exits 1 on mismatch.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit a one-line JSON object: {\"sha256\": ..., \"path\": ...}",
    )
    parser.add_argument(
        "--chunk-mb",
        type=int,
        default=8,
        help="Read size in MiB per chunk (default: 8).",
    )
    args = parser.parse_args()

    digest = sha256_file(args.file, chunk_mb=args.chunk_mb)

    if args.expect:
        if digest.lower() != args.expect.lower():
            print(
                f"SHA-256 mismatch for {args.file.name}: "
                f"{digest} != {args.expect}"
            )
            raise SystemExit(1)

    if args.json:
        print(json.dumps({"sha256": digest, "path": str(args.file)}))
    else:
        print(digest)


if __name__ == "__main__":
    main()
