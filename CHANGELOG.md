# Changelog

All notable changes to this project will be documented in this file.

Starting with version 0.3.2, releases are tracked here.

## [Unreleased]

## [0.3.2] - 2026-04-21

### Added

- `foli assign-probes --skip-seqkit` support through the `skip_seqkit` CLI option, so output stats generation can be skipped when it is not needed.
- Source-tree package version fallback in `folitools.__version__`, so local runs report the version from `pyproject.toml`.

### Changed

- Switched compressed I/O in `folitools.add_umi` from `gzip.open` to `xopen` so it can use `pigz` when available.
