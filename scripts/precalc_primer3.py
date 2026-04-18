import csv
import gzip
import os
from pathlib import Path

import primer3
import pyarrow as pa
import pyarrow.parquet as pq
from joblib import Parallel, delayed
from tqdm.auto import tqdm, trange


def _normalize_primer_seq(seq: str) -> str:
    return seq.strip().upper().replace(" ", "")


def _primer_index_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("primer_id", pa.uint32()),
            pa.field("primer_name", pa.string()),  # <primerIndex>_[fwd|rev]
            pa.field("primerIndex", pa.string()),  # original pair name
            pa.field("transcriptID", pa.string()),
            pa.field("orientation", pa.string()),  # "fwd" / "rev"
            pa.field("sequence", pa.string()),
        ]
    )


def _p3_schema() -> pa.Schema:
    # Filtered output only (structure_found=True and dg<=cutoff), so we don't store structure_found.
    return pa.schema(
        [
            pa.field("primer1_id", pa.uint32()),
            pa.field("primer2_id", pa.uint32()),
            pa.field("tm", pa.float64()),
            pa.field("dg", pa.float64()),
            pa.field("dh", pa.float64()),
            pa.field("ds", pa.float64()),
        ]
    )


def _open_parquet_writer(
    file_path: Path, schema: pa.Schema, *, compression: str
) -> pq.ParquetWriter:
    file_path.parent.mkdir(parents=True, exist_ok=True)
    return pq.ParquetWriter(
        where=str(file_path),
        schema=schema,
        compression=compression,
        use_dictionary=True,
        write_statistics=True,
    )


def _write_table(
    writer: pq.ParquetWriter,
    schema: pa.Schema,
    columns: dict[str, list[object]],
    *,
    row_group_size: int,
) -> None:
    arrays = [
        pa.array(columns[name], type=schema.field(name).type) for name in schema.names
    ]
    table = pa.Table.from_arrays(arrays, names=schema.names)
    writer.write_table(table, row_group_size=row_group_size)


def _write_csv_gz_atomic(
    out_path: Path,
    header: list[str],
    rows: list[list[object]],
) -> None:
    tmp_path = out_path.with_suffix(out_path.suffix + ".inprogress")
    tmp_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(tmp_path, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

    os.replace(tmp_path, out_path)


def _compute_block_primer3_to_parquet(
    i_start: int,
    i_end: int,
    primer_seqs: list[str],
    out_calc_dir: str,
    out_counts_parts_dir: str,
    *,
    dg_store_cutoff_cal: float,
    batch_rows: int,
    row_group_size: int,
    compression: str,
    skip_if_exists: bool,
) -> tuple[str, str]:
    """
    Compute primer3 heterodimer for primer1_id in [i_start, i_end) vs all primer2_id,
    filter by (structure_found AND dg <= cutoff), and write only passing rows.

    Also writes per-primer (primer1_id) pass counts as a gzipped CSV for this block.

    Args:
        i_start: Block start primer1_id (inclusive).
        i_end: Block end primer1_id (exclusive).
        primer_seqs: All primer sequences (indexed by primer_id).
        out_calc_dir: Directory for Parquet output.
        out_counts_parts_dir: Directory for per-block count CSV parts.
        dg_store_cutoff_cal: Filter cutoff (cal/mol). Keep if dg <= cutoff.
        batch_rows: Flush buffer around this many kept rows.
        row_group_size: Parquet row group size.
        compression: Parquet compression codec.
        skip_if_exists: If True, skip work if both outputs for this block exist.

    Returns:
        (calc_parquet_path, counts_csv_gz_path)
    """
    n = len(primer_seqs)
    schema = _p3_schema()

    calc_path = Path(out_calc_dir) / f"part-{i_start:09d}-{i_end - 1:09d}.parquet"
    counts_path = (
        Path(out_counts_parts_dir)
        / f"primer1_pass_counts-{i_start:09d}-{i_end - 1:09d}.csv.gz"
    )

    if skip_if_exists and calc_path.exists() and counts_path.exists():
        return str(calc_path), str(counts_path)

    calc_tmp_path = calc_path.with_suffix(calc_path.suffix + ".inprogress")

    # Use primer3.calc_heterodimer if present; else fall back.
    calc_heterodimer = getattr(
        primer3, "calc_heterodimer", primer3.bindings.calc_heterodimer
    )

    writer = _open_parquet_writer(calc_tmp_path, schema, compression=compression)

    try:
        buf: dict[str, list[object]] = {name: [] for name in schema.names}

        def flush(force: bool = False) -> None:
            if force or len(buf["primer1_id"]) >= batch_rows:
                if buf["primer1_id"]:
                    _write_table(writer, schema, buf, row_group_size=row_group_size)
                    for k in buf:
                        buf[k].clear()

        # Per-primer counts for this block
        count_rows: list[list[object]] = []
        for i in trange(i_start, i_end):
            s1 = primer_seqs[i]
            passed_for_i = 0

            for j in range(n):
                s2 = primer_seqs[j]

                # No ASCII output stored; only need structure_found + thermodynamics.
                r = calc_heterodimer(s1, s2, output_structure=False)
                if not bool(getattr(r, "structure_found", False)):
                    continue

                dg = float(getattr(r, "dg"))
                if dg > dg_store_cutoff_cal:
                    continue

                # Keep row
                buf["primer1_id"].append(i)
                buf["primer2_id"].append(j)
                buf["tm"].append(float(getattr(r, "tm")))
                buf["dg"].append(dg)
                buf["dh"].append(float(getattr(r, "dh")))
                buf["ds"].append(float(getattr(r, "ds")))

                passed_for_i += 1
                flush(force=False)

            count_rows.append(
                [
                    i,
                    passed_for_i,
                    n,
                    (passed_for_i / n) if n else 0.0,
                ]
            )

        flush(force=True)

    finally:
        writer.close()

    os.replace(calc_tmp_path, calc_path)

    _write_csv_gz_atomic(
        out_path=counts_path,
        header=["primer1_id", "passed_pairs", "n_total_primers", "passed_fraction"],
        rows=count_rows,
    )

    return str(calc_path), str(counts_path)


if __name__ == "__main__":
    # --- user-specified config ---
    input_pairs_parquet = (
        "/home/ubuntu/dev/folitools/src/folitools/primer_selection/data/human/"
        "02.SADDLE_input.320_380.primer_0.parquet"
    )
    out_dir = "/home/ubuntu/data/folitools/precalc_primer3/human_320_380"

    cores = 12
    dg_store_cutoff_cal = -9000.0  # cal/mol (=-9 kcal/mol)

    # --- tuning knobs (adjust as needed) ---
    # Target ~this many rows per part file (kept rows are fewer; this controls i_block_size for scheduling/files).
    rows_per_part_target = 50_000_000
    i_block_size = 0  # 0 => auto from rows_per_part_target and n

    # Buffering / parquet layout
    batch_rows = 200_000
    row_group_size = 1_000_000
    compression = "zstd"

    # Resume behavior
    skip_existing_parts = True

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    primers_index_path = out_path / "primers.parquet"
    calc_dir = out_path / "calc_heterodimer"
    counts_parts_dir = out_path / "counts_parts"
    final_counts_path = out_path / "primer1_pass_counts.csv.gz"

    calc_dir.mkdir(parents=True, exist_ok=True)
    counts_parts_dir.mkdir(parents=True, exist_ok=True)

    # --- build primer index (2 rows per input row) ---
    pairs_table = pq.read_table(
        input_pairs_parquet,
        columns=["transcriptID", "primerIndex", "sequenceLeft", "sequenceRight"],
    )

    transcript_ids = pairs_table.column("transcriptID").to_pylist()
    pair_names = pairs_table.column("primerIndex").to_pylist()
    left_seqs = pairs_table.column("sequenceLeft").to_pylist()
    right_seqs = pairs_table.column("sequenceRight").to_pylist()

    m = len(pair_names)
    primer_ids = list(range(2 * m))

    primer_names = [f"{p}_fwd" for p in pair_names] + [f"{p}_rev" for p in pair_names]
    primer_seqs = [_normalize_primer_seq(s) for s in left_seqs] + [
        _normalize_primer_seq(s) for s in right_seqs
    ]
    orientations = ["fwd"] * m + ["rev"] * m
    transcript_ids_out = transcript_ids + transcript_ids
    pair_names_out = pair_names + pair_names

    idx_schema = _primer_index_schema()
    idx_table = pa.Table.from_arrays(
        [
            pa.array(primer_ids, type=pa.uint32()),
            pa.array(primer_names, type=pa.string()),
            pa.array(pair_names_out, type=pa.string()),
            pa.array(transcript_ids_out, type=pa.string()),
            pa.array(orientations, type=pa.string()),
            pa.array(primer_seqs, type=pa.string()),
        ],
        names=idx_schema.names,
    )
    pq.write_table(idx_table, str(primers_index_path), compression="zstd")

    n = len(primer_seqs)
    if i_block_size <= 0:
        i_block_size = max(1, rows_per_part_target // max(1, n))

    blocks: list[tuple[int, int]] = []
    for i_start in range(0, n, i_block_size):
        i_end = min(n, i_start + i_block_size)
        blocks.append((i_start, i_end))

    # --- compute primer3 filtered pairs to parquet + per-block counts ---
    results = Parallel(n_jobs=cores, backend="loky")(
        delayed(_compute_block_primer3_to_parquet)(
            i_start=i_start,
            i_end=i_end,
            primer_seqs=primer_seqs,
            out_calc_dir=str(calc_dir),
            out_counts_parts_dir=str(counts_parts_dir),
            dg_store_cutoff_cal=dg_store_cutoff_cal,
            batch_rows=batch_rows,
            row_group_size=row_group_size,
            compression=compression,
            skip_if_exists=skip_existing_parts,
        )
        for (i_start, i_end) in tqdm(
            blocks, desc="Computing blocks (primer3, multiprocessing)"
        )
    )

    # --- merge per-block count CSVs into one csv.gz (streaming; no pandas) ---
    # We write header once, then append all block rows (skipping each part header).
    part_count_files = sorted({counts_path for _, counts_path in results})
    tmp_final = final_counts_path.with_suffix(final_counts_path.suffix + ".inprogress")

    with gzip.open(tmp_final, "wt", newline="") as out_f:
        wrote_header = False
        for part_path in tqdm(part_count_files, desc="Merging count parts"):
            with gzip.open(part_path, "rt", newline="") as in_f:
                for line_idx, line in enumerate(in_f):
                    if line_idx == 0:
                        if not wrote_header:
                            out_f.write(line)
                            wrote_header = True
                        continue
                    out_f.write(line)

    os.replace(tmp_final, final_counts_path)

    print(
        {
            "primers_index": str(primers_index_path),
            "calc_heterodimer_dir": str(calc_dir),
            "primer1_pass_counts": str(final_counts_path),
            "n_primers": n,
            "i_block_size": i_block_size,
            "n_blocks": len(blocks),
        }
    )
