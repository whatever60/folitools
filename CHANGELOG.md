# Changelog

All notable changes to this project will be documented in this file.

Starting with version 0.3.2, releases are tracked here.

## [Unreleased]

### Added

- `summary_stats` / `foli summary`: `properly_mapped_depth` is renamed
  to `counted_assigned_depth` and a new `counted_depth` column is added
  in front of it. Both are per-QNAME counters emitted by
  `foli_add_tags` on its SUMMARY line so they live in the same units
  as the rest of the table:
  - `counted` = `mapped` ∩ `good_umi` — the QNAMEs that umi_tools
    emits to `group.tsv.gz` under the current
    `--chimeric-pairs use` / `--unmapped-reads discard` /
    `--unpaired-reads use` flags.
  - `counted_assigned` = `counted` ∩ `assigned` — equals the row sum
    of the pre-dedup count matrix after the `Unassigned,` filter in
    `folitools.get_matrix`.
  The DAG is now `raw → qc → long → {not_na → good_umi, mapped →
  assigned} → counted → counted_assigned → n_umi → n_genes`, with
  `mapped → counted` and `assigned → counted_assigned` as the cross
  edges. `foli summary --strict` (new flag) additionally asserts
  `counted_depth == #rows in --group-tsvs`,
  `counted_assigned_depth == #rows in --group-tsvs after the
  Unassigned, filter`, and
  `counted_assigned_depth == row sum of --count-matrix-raw`. New
  optional `--group-tsvs` arg accepts a glob/path/list of
  `<sample>.group.tsv.gz` files and is consulted only when `--strict`
  is set. `count_matrix_raw` is no longer the primary source of
  `counted_assigned_depth`; it falls back to populating that column
  only when the add_tags log is absent (preserving the previous
  behavior for callers that only have the matrix).
- `summary_stats` / `foli summary`: two new columns that capture the
  mapping/annotation funnel running parallel to the library-quality
  funnel — `mapped_depth` (QNAMEs with both primary mates aligned) and
  `assigned_depth` (QNAMEs whose `XF` doesn't start with `Unassigned`).
  Both are per-QNAME counters added to the `foli_add_tags --log`
  SUMMARY line so they share units with `not_na_adapter` and
  `good_umi` (featureCounts' `*.summary` `Assigned` row was rejected
  as a source: under foli's `-p` invocation, without
  `--countReadPairs`, it counts reads not pairs and is therefore ~2×
  the right value).
- The columns now form a DAG (`raw → qc → long → {not_na → good_umi,
  mapped → assigned} → properly_mapped → n_umi → n_genes`) and the
  per-row monotonic-non-increasing assert is replaced by an edge-by-edge
  DAG check (parent ≥ child on present-value pairs; edges with a NaN
  endpoint are skipped). When run as part of the foli pipeline, the
  QNAME intersection of `assigned_depth` and `good_umi_depth` is
  expected to equal `properly_mapped_depth` up to chimeric/unpaired
  alignments that umi_tools sends to BAM only.
- Every pipeline step now stamps its log with the folitools version that
  produced it, so a single `grep folitools_version` (JSON) or
  `grep '^# folitools'` (text) tells you which run wrote any given log.
  None of the markers are hardcoded — Python pulls from `__version__`,
  Rust pulls from `env!("CARGO_PKG_VERSION")` (kept in sync via the
  release commit). Per step:
  - `foli qc`: adds a `"folitools_version"` key to each
    `<sample>.fastp.json`.
  - `foli assign-probes`: re-enables cutadapt's `--json` report at
    `<sample>.cutadapt.json` and adds a `"folitools_version"` key.
  - `foli map`: `foli_add_tags` (Rust binary and Python reference) now
    writes a trailing `# folitools <version>` line to its `--log` file
    after the SUMMARY line.
  - `foli count`: appends `# folitools <version>` to each
    `<sample>.group.log` after `umi_tools group` finishes.
  - `foli get-read-stats`: stamps `folitools_version` into each
    parquet file's schema key-value metadata (no sidecar log).
  - `foli get-count-mtx` and `foli summary` already carry the version
    in the index label of their output tables.
  - `<dir>.stats` (seqkit `--tabular`) is intentionally **not**
    decorated — it's a strict TSV consumed by `summary_stats` and any
    extra line breaks parsing.
- `foli-primer` subcommands also stamp, but at most once per command
  and only into surfaces that don't disturb the file's primary
  content:
  - `summary`: workbook `properties.description` core property on each
    xlsx output (`primer_to_order.xlsx`, `primer_to_order_idt.xlsx`).
    Visible via Excel's "File → Info → Properties" dialog or
    `openpyxl.load_workbook(...).properties.description`.
  - `recover`: a single `logger.info` line at the top of `recover.log`.
    The CSV / xlsx / PDF / FASTA outputs are deliberately uncluttered.
  - `subset`, `saddle`, `sgad`, `product`: no per-command stamp.
    Their outputs are non-matrix TSVs or a FASTA; in either case any
    in-file marker would either pollute the table (no empty top-left
    cell to repurpose like `foli get-count-mtx` has) or break a
    downstream parser (cutadapt's dnaio rejects `;` comments;
    Biopython's `fasta` reader deprecated `#` ones). The xlsx helper
    lives in `folitools.primer_selection._versioning`.

### Changed

- `foli get-count-mtx`: the Unassigned-row filter in
  `get_matrix.process_count_file_simple` and `read_counts` now anchors
  on the actual `<gene>,<primers>` delimiter (`Unassigned,`) rather
  than the bare `Unassigned` prefix. The shape of the rows
  `foli_add_tags` writes for unassigned reads is `Unassigned,<primer>`,
  so behavior is unchanged in the common path; the anchor just keeps a
  hypothetical real gene name starting with `Unassigned` from being
  silently dropped.
- `foli map` defaults flipped to match the foli-seq protocol and the
  empirical behavior of recent runs. `--strand` now defaults to `1`
  (forward-stranded) instead of `0`; `--allow-overlap` and
  `--allow-multimapping` now default to `true` instead of `false`.
  Inspecting STAR's `ReadsPerGene.out.tab` across the
  `20260326_sparc_ibd_clinical_trial_batch2_batch3` run shows `s=2`
  classifies 50–8000× more reads as `N_noFeature` than `s=0`/`s=1`,
  while `s=0` and `s=1` are within a few percent of each other —
  reads sit on the gene's strand, so `s=1` is the right default.
  Pass `--strand 0` / `--strand 2` / `--no-allow-overlap` /
  `--no-allow-multimapping` to override.
- `foli_add_tags --log` SUMMARY semantics: `not_na_adapter` is now
  adapter-only (R1 QNAMEs whose primer pair was recognized) and
  `good_umi` requires adapter-ok **and** umi-ok — exactly the
  `criteria_ok` gate that controls `UC`-tag emission and downstream
  UMI dedup. SUMMARY field order, `summary_stats` column order, and
  the README table are flipped to match (`not_na_adapter_depth` →
  `good_umi_depth`), preserving the monotonic-non-increasing invariant
  with names that reflect the actual pipeline funnel.
- `foli count`: umi_tools group flags switched from `output` to
  `discard` for `--unmapped-reads` and to `use` for `--chimeric-pairs`
  / `--unpaired-reads`. `discard` for unmapped is TSV-equivalent to
  `output` (we don't pass `--output-bam`) but lets pysam skip the
  unmapped tail. `use` for chimeric pairs keeps cross-contig pairs in
  the bundle path so their multi-gene `XF` (R1's gene + R2's gene,
  aggregated by `foli_add_tags`) lands as a `GeneA|GeneB` co-assignment
  column — intentional, to preserve gene-fusion-like signal from
  amplicons whose primer-end and distal-end map to different genes.
  `use` for unpaired is irrelevant in foli (data is always paired) and
  set for symmetry.

## [0.6.1] - 2026-04-23

### Added

- `foli summary` CLI command wrapping
  `folitools.summary.summary_stats`. Aggregates per-sample read counts
  from every pipeline stage into a single sample × metric table and
  writes TSV/CSV chosen by the `--output` extension (`.gz` supported).

## [0.6.0] - 2026-04-23

### Added

- `folitools.summary.summary_stats`: new library function that aggregates
  per-sample read counts from every pipeline stage — seqkit stats
  (raw and fastp-trimmed), STAR `Log.final.out`, `foli_add_tags`
  SUMMARY logs, and raw/dedup count matrices — into a single
  `sample × metric` DataFrame (`Int64`). Any source left `None` yields
  an all-`NA` column. Columns are ordered so a healthy run is
  monotonically non-increasing across each row; the function asserts
  this invariant so pipeline regressions surface immediately.
- `foli get-count-mtx --output-raw`: writes a pre-UMI-dedup (read-count)
  matrix alongside or instead of `--output`. At least one of the two
  must be set; both can be written in one invocation and share the
  same input scan.
- `foli_add_tags` now writes a single-line `SUMMARY` record to the
  `--log` file at program end, capturing per-run `total_r1`,
  `good_umi`, and `not_na_adapter` counts. Consumed by `summary_stats`
  to drive the `good_umi_depth` / `not_na_adapter_depth` columns.

### Changed

- `foli_add_tags` (both the Rust binary and the Python reference
  implementation) now stamps `CB` / `US` / `PR` / `UC` / `XN` / `XT` on
  the primary R1 only. `umi_tools group`/`count` in `--paired` mode
  only inspects R1, so mirroring these onto R2 was redundant and
  caused downstream raw read counting to double the true value.
  R2 is passed through with its featureCounts-assigned `XT` preserved
  for the `XF` aggregation.
- `folitools.get_matrix.process_count_file_simple` now explicitly
  asserts that `read_id` is unique in each `group.tsv.gz`, protecting
  the new raw-count path (and the existing dedup path) from silent
  double-counting if the primary-R1-only invariant is ever broken.

## [0.5.0] - 2026-04-22

### Added

- Rust-backed `foli_add_tags` console script, a drop-in replacement for
  `python -m folitools.add_tags` used by `foli map`. The logic lives in
  `rust/src/add_tags.rs` and is exposed to Python via a pyo3 extension
  module (`folitools._rust_native`); a thin Python entry at
  `folitools._add_tags_entry` forwards `argv` into it. Matches the Python
  implementation's output bit-for-bit in the common path; same CLI
  surface (`-i`, `-o`, `--cell_tag_name`, `--cell_tag`, `--log`).
- `.github/workflows/publish.yml` now builds platform wheels with maturin
  (Linux manylinux_2_28 x86_64 + macOS universal2) plus an sdist, and
  publishes to TestPyPI then PyPI on `v*` tags via trusted publishing.

### Changed

- **Build backend**: `uv_build` → `maturin`. The wheel now bundles both
  the Python package and a compiled Rust extension
  (`folitools._rust_native`) so `pip install folitools` is sufficient —
  no separate Rust build step.
- `foli map` (specifically `foli_03_map.sh`) always calls `foli_add_tags`.
  The Python `folitools.add_tags` module is kept as the reference
  implementation and for import in tests, but is no longer wired into the
  mapping pipeline.
- `samtools collate` in `foli_03_map.sh` switched from `-l 1 --threads 1`
  to `-u --threads 4`. `add_tags` consumes the output directly, so the
  compression round-trip was pure overhead; the change saves ~45 s on a
  12.6 M-record sample.
- STAR now runs with `--chimSegmentMin 12` to emit chimeric alignments as
  supplementary records in `Aligned.out.bam`. featureCounts' `-B -C`
  already excludes them from counts, so gene-level output is unchanged;
  the extra records are useful for IGV inspection and fusion-aware
  downstream tooling.
- `add_tags.py`: per-mate SAM-compliance check (exactly one primary per
  mate per QNAME) replaced the previous "total == 2" check, with per-mate
  FLAG dumps in the log on violation. Non-compliant extras are downgraded
  from primary to secondary (`flag |= 0x100`) so downstream (umi_tools)
  sees a compliant BAM. Warnings now go to a `--log` file instead of
  stderr. The hot loop also caches `read.flag` and uses bit tests in
  place of `is_read1` / `is_secondary` / `is_supplementary` property
  lookups.
- `rust/Cargo.toml`: feature-gated the `calc_rest_stats` binary behind a
  `heavy_deps` feature so a default build doesn't compile polars/rayon.

### Performance

- On a 12.6 M-record sample (real featureCounts output, stdin → stdout
  pipe configuration matching the production pipeline): `add_tags`
  dropped from **118.5 s** (Python + pysam) to **72.4 s** (Rust
  + rust-htslib), a 1.64× wall-time and 1.85× user-CPU speedup. Peak
  memory dropped from ~34 MB to ~3 MB.

## [0.4.1] - 2026-04-22

### Changed

- Bumped the `cutadapt-folitools` pin from `==5.3.1.post0` to
  `>=5.3.1.post1`. The new upstream release includes a q-gram lemma
  tightening of `SeedMultiAdapterFilter.is_acceptable` so that adapters
  whose pigeonhole-safe seed size would fall below `MIN_SEED_SIZE` are
  correctly routed to the non-indexed `MultipleAdapters` path instead
  of being wrapped with a too-large seed. This makes the "bit-identical
  to upstream cutadapt" guarantee hold for all accepted adapter groups,
  not just ones we had measured. No observable change for the
  `foli assign-probes` workload — the 542×2 primer set already resolves
  to seed=6, safely above the floor.

## [0.4.0] - 2026-04-22

### Changed

- `cutadapt` is now installed from our fork, published on PyPI as
  [`cutadapt-folitools`](https://pypi.org/project/cutadapt-folitools/)
  (source: [whatever60/cutadapt@folitools-perf](https://github.com/whatever60/cutadapt/tree/folitools-perf)),
  instead of the bioconda/PyPI `cutadapt` package. The fork installs the
  same `cutadapt` Python module and console script as upstream, and adds
  performance improvements that matter for `foli assign-probes` with
  large primer panels; `Match` objects are bit-identical to upstream
  cutadapt (verified 50,000/50,000 on real foli-seq reads and via the
  cutadapt test suite).
- `environment.yaml` no longer lists `cutadapt`; pip installs it from
  the fork via the `pyproject.toml` dependency. A C compiler (`gcc`,
  already in `environment.yaml`) is required because `cutadapt-folitools`
  is published as an sdist and builds Cython extensions from source.

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
