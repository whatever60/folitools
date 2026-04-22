# Changelog

All notable changes to this project will be documented in this file.

Starting with version 0.3.2, releases are tracked here.

## [Unreleased]

## [0.4.1] - 2026-04-22

### Changed

- `cutadapt` is now depended on as the PyPI distribution
  [`cutadapt-folitools`](https://pypi.org/project/cutadapt-folitools/),
  instead of a direct git URL. The installed Python module and CLI are
  unchanged (still `import cutadapt`, still the `cutadapt` command).
  Installing folitools no longer requires `git` or a C toolchain — users
  get a prebuilt wheel when one is available for their platform.
- Pinned to `cutadapt-folitools >= 5.3.1.post1`, which includes the
  q-gram lemma tightening in `SeedMultiAdapterFilter` (see below).

### Fixed (fork side)

- `SeedMultiAdapterFilter.is_acceptable` now rejects adapters whose
  pigeonhole-safe seed size falls below `MIN_SEED_SIZE`. Previously the
  computed seed could be clamped up to 3 even when the q-gram lemma
  required a smaller value, which could silently drop valid matches for
  very short non-anchored adapters with a high error rate. No observable
  effect on the `foli assign-probes` workload — the 542×2 primers yield
  seed=6, safely above the floor. Tightening ensures bit-identical
  semantics to upstream for all accepted adapter groups, not just the
  ones we measured.

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
