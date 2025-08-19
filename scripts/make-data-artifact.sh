#!/usr/bin/env bash
# Create a reproducible .tar.zst from your data dir and emit a JSON manifest.
# Works best with GNU tar. Falls back to a portable path if BSD tar is the only option.

set -euo pipefail
IFS=$'\n\t'

# --- Config (override via env or flags) ---
SRC_DIR="${SRC_DIR:-src/folitools/primer_selection/data}"
DIST_DIR="${DIST_DIR:-dist/data}"
NAME="${NAME:-primer_data}"
# Default version tag: YYYY.MM.DD or pass --version=...
VERSION="${VERSION:-$(date -u +%Y.%m.%d)}"
# Reproducible timestamp: from SOURCE_DATE_EPOCH if set, else midnight UTC today
if [[ -z "${SOURCE_DATE_EPOCH:-}" ]]; then
  SOURCE_DATE_EPOCH="$(date -u -d "$(date -u +%Y-%m-%d) 00:00:00" +%s 2>/dev/null || date -u -j -f "%Y-%m-%d %H:%M:%S" "$(date -u +%Y-%m-%d) 00:00:00" +%s)"
fi

usage() {
  cat <<EOF
Usage: ${0##*/} [--src DIR] [--out DIR] [--name NAME] [--version TAG]
Creates: <out>/<name>-<version>.tar.zst and <out>/<name>-<version>.manifest.json

Env overrides:
  SRC_DIR   (default: $SRC_DIR)
  DIST_DIR  (default: $DIST_DIR)
  NAME      (default: $NAME)
  VERSION   (default: $VERSION)
EOF
}

for arg in "$@"; do
  case "$arg" in
    --src=*) SRC_DIR="${arg#*=}";;
    --out=*) DIST_DIR="${arg#*=}";;
    --name=*) NAME="${arg#*=}";;
    --version=*) VERSION="${arg#*=}";;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $arg" >&2; usage; exit 1;;
  esac
done

if [[ ! -d "$SRC_DIR" ]]; then
  echo "ERR: data dir not found: $SRC_DIR" >&2
  exit 1
fi

mkdir -p "$DIST_DIR"
workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

top="${NAME}-${VERSION}"
staged="$workdir/$top"
mkdir -p "$staged"

# Stage data (rsync preserves structure; adjust filters if needed)
rsync -a --delete "$SRC_DIR"/ "$staged"/

artifact="$DIST_DIR/${top}.tar.zst"
manifest="$DIST_DIR/${top}.manifest.json"

# --- Choose tar ---
tar_bin="$(command -v gtar || true)"
if [[ -z "$tar_bin" ]]; then
  tar_bin="$(command -v tar)"
fi
if [[ -z "$tar_bin" ]]; then
  echo "ERR: no tar found" >&2
  exit 1
fi

is_gnu_tar=0
if "$tar_bin" --version 2>/dev/null | grep -qi 'gnu tar'; then
  is_gnu_tar=1
fi

# --- Create artifact (deterministic when GNU tar is available) ---
if [[ "$is_gnu_tar" -eq 1 ]]; then
  # Use GNU tar with full reproducibility knobs and zstd compression
    echo $tar_bin $workdir $artifact
  ( cd "$workdir"
    TZ=UTC \
    "$tar_bin" \
      --sort=name \
      --owner=0 --group=0 --numeric-owner \
      --mtime="@$SOURCE_DATE_EPOCH" \
      -I 'zstd -T0 -19' \
      -cf "$artifact" \
      "$top"
  )
else
  echo "WARN: non-GNU tar detected; making best-effort portable artifact." >&2
  ( cd "$workdir"
    # Portable: pipe to zstd (no --sort/--owner). Still stable enough for releases.
    tar -cf - "$top" | zstd -T0 -19 -q -o "$artifact"
  )
fi

# --- Compute SHA-256 ---
if command -v sha256sum >/dev/null 2>&1; then
  sha256="$(sha256sum "$artifact" | awk '{print $1}')"
  sha_cmd="sha256sum <file>"
elif command -v shasum >/dev/null 2>&1; then
  sha256="$(shasum -a 256 "$artifact" | awk '{print $1}')"
  sha_cmd="shasum -a 256 <file>"
else
  echo "ERR: need sha256sum or shasum on PATH" >&2
  exit 1
fi

size_bytes="$(stat -c%s "$artifact" 2>/dev/null || stat -f%z "$artifact")"

# --- Quick integrity check (list archive) ---
tar -tf "$artifact" >/dev/null

# --- Write manifest ---
cat >"$manifest" <<JSON
{
  "name": "$NAME",
  "version": "$VERSION",
  "filename": "$(basename "$artifact")",
  "size_bytes": $size_bytes,
  "sha256": "$sha256",
  "sha256_cmd": "$sha_cmd",
  "source_dir_rel": "$SRC_DIR",
  "top_level_dir_in_archive": "$top",
  "created_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON

echo "OK: $artifact"
echo "OK: $manifest"
echo
echo "Use in code (example):"
cat <<EOF
ensure_dataset(
    version="$VERSION",
    url="https://<host>/<path>/$(basename "$artifact")",
    sha256="$sha256",
    archive_member="$top",  # or None to extract all
)
EOF
