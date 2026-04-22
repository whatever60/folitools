# Changelog

All notable changes to this project will be documented in this file.

Starting with version 0.3.2, releases are tracked here.

## [Unreleased]

## [0.4.0] - 2026-04-22

### Changed

- `cutadapt` is now installed from our pinned fork
  [whatever60/cutadapt@folitools-perf](https://github.com/whatever60/cutadapt/tree/folitools-perf)
  instead of the bioconda/PyPI package. The fork adds performance
  improvements that matter for `foli assign-probes` with large primer
  panels; `Match` objects are bit-identical to upstream cutadapt (verified
  50,000/50,000 on real foli-seq reads and via the cutadapt test suite).
- `environment.yaml` no longer lists `cutadapt`; pip installs it from the
  fork via the `pyproject.toml` dependency. A C compiler (`gcc`, already
  in `environment.yaml`) is required because the fork builds Cython
  extensions from source.

### Performance

- End-to-end on a representative sample (3.1 M paired reads, 542×2 i5/i7
  primers, `-j 16`), the combined effect of the cutadapt fork and the
  `xopen`/`pigz`-backed writer in `add_umi` (from 0.3.2) cuts `foli
  assign-probes` from ~7:30 to ~2:45 per sample (~2.7× faster, ~63% wall
  time reduction).

### Fork features (cutadapt side)

- `SeedMultiAdapterFilter`: seed-and-extend pre-filter for groups of 8+
  non-anchored multi-adapters. Auto-activates; transparent to callers.
- `--input-compression-threads N`: opt-in parallel decompression of
  gzipped inputs via `pigz` (default 0 = in-process via `isal`).
- `NPrefixIndexedAdapters`: fast path for `^N{k}<body>` anchored adapters.
  Present but intentionally dormant; reserved for explicit opt-in.

## [0.3.2] - 2026-04-21

### Added

- `foli assign-probes --skip-seqkit` support through the `skip_seqkit` CLI option, so output stats generation can be skipped when it is not needed.
- Source-tree package version fallback in `folitools.__version__`, so local runs report the version from `pyproject.toml`.

### Changed

- Switched compressed I/O in `folitools.add_umi` from `gzip.open` to `xopen` so it can use `pigz` when available.
