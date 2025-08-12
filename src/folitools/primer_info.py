import os
import gzip

from typing import Iterable
from tqdm.auto import tqdm
from Bio import SeqIO
from joblib import Parallel, delayed
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

from .utils import expand_path_to_list


class IncrementalParquetWriter:
    """
    Buffer rows and write them incrementally as row groups to a single Parquet file
    via a kept-open pyarrow.ParquetWriter.

    Args:
        path: Output parquet file path.
        schema: Polars schema or pyarrow.Schema to enforce on each batch.
        batch_size: Number of rows to buffer before flushing as a row group.
        **parquet_kwargs: Additional kwargs forwarded to pyarrow.parquet.ParquetWriter
            (e.g., compression="snappy").
    """

    def __init__(
        self,
        path: str,
        schema: pl.Schema | pa.Schema,
        batch_size: int,
        **parquet_kwargs,
    ) -> None:
        if isinstance(schema, pl.Schema):
            df = pl.DataFrame({k: pl.Series(k, [], dtype=v) for k, v in schema.items()})
            arrow_schema = df.to_arrow().schema
        elif isinstance(schema, pa.Schema):
            arrow_schema = schema
        else:
            raise TypeError("schema must be pl.Schema or pyarrow.Schema")

        self._batch_size = batch_size
        self._buffer: list[dict] = []
        self._writer = pq.ParquetWriter(path, arrow_schema, **parquet_kwargs)

    def write_row(self, row: dict) -> None:
        """
        Add a single row (mapping column name to value) and flush if batch_size reached.
        """
        self._buffer.append(row)
        if len(self._buffer) >= self._batch_size:
            self.flush()

    def write_rows(self, rows: Iterable[dict]) -> None:
        """
        Add multiple rows at once.
        """
        for row in rows:
            self.write_row(row)

    def flush(self) -> None:
        """
        Flush buffered rows as one row group to the Parquet file.
        """
        if not self._buffer:
            return
        df = pl.DataFrame(self._buffer)
        table = df.to_arrow()
        self._writer.write_table(table)
        self._buffer.clear()

    def close(self) -> None:
        """
        Flush remaining rows and close writer.
        """
        self.flush()
        self._writer.close()

    def __enter__(self) -> "IncrementalParquetWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def process_sample(
    read1_path: str, read2_path: str, output_path: str, batch_size: int = 50_000
) -> None:
    """
    Parse paired FASTQ (gzipped) and incrementally write read metadata to a Parquet file.

    Args:
        read1_path: Path to the R1 FASTQ file (.fq.gz).
        read2_path: Path to the R2 FASTQ file (.fq.gz).
        output_path: Path where the parquet output will be saved. Always written/overwritten.
        batch_size: Number of rows to buffer before flushing to parquet.
    """
    schema = {
        "read_id": pl.Utf8,
        "primer_fwd": pl.Utf8,
        "primer_rev": pl.Utf8,
        "r1_length": pl.Int64,
        "r2_length": pl.Int64,
        "umi5": pl.Utf8,
        "umi3": pl.Utf8,
    }

    # pyarrow compression default is snappy, but polars default is zstd
    with (
        IncrementalParquetWriter(
            output_path, pl.Schema(schema), batch_size=batch_size, compression="zstd"
        ) as writer,
        gzip.open(read1_path, "rt") as r1_handle,
        gzip.open(read2_path, "rt") as r2_handle,
    ):
        for r1, r2 in zip(
            SeqIO.parse(r1_handle, "fastq"),
            SeqIO.parse(r2_handle, "fastq"),
            strict=True,
        ):
            if r1.id != r2.id:
                raise ValueError(f"Read ID mismatch: {r1.id} vs {r2.id}")

            # parse primers from description, expecting last token like "FWD+REV"
            try:
                primer_fwd, primer_rev = r1.description.split(" ")[-1].split("+")
            except ValueError:
                raise ValueError(
                    f"Unexpected primer format in description: {r1.description!r}"
                )

            # parse read_id and UMIs: expecting format readid_umi5_umi3
            parts = r1.id.split("_")
            if len(parts) != 3:
                raise ValueError(f"Unexpected read ID format: {r1.id}")
            read_id, umi5, umi3 = parts

            row = {
                "read_id": read_id,
                "primer_fwd": primer_fwd,
                "primer_rev": primer_rev,
                "r1_length": len(r1.seq),
                "r2_length": len(r2.seq),
                "umi5": umi5,
                "umi3": umi3,
            }
            writer.write_row(row)


def get_read_stats(
    read1_pattern: str | list[str],
    output_dir: str,
    cores: int = 1,
    overwrite: bool = False,
    skip: int = 0,
) -> None:
    """
    Process a single pair of FASTQ files (R1 and R2), extracting UMI and primer info,
    and write to a parquet file if it does not already exist.
    """
    os.makedirs(output_dir, exist_ok=True)
    read1_files = expand_path_to_list(read1_pattern)
    args = []
    for i, read1_path in enumerate(tqdm(read1_files)):
        if i < skip:
            continue
        sample = os.path.basename(read1_path).split("_")[0]
        out_path = os.path.join(output_dir, f"{sample}.parquet")
        read2_path = read1_path.replace("_1.fq.gz", "_2.fq.gz")
        output_path = os.path.join(output_dir, f"{sample}.parquet")
        if os.path.exists(out_path) and not overwrite:
            print(f"Output file {out_path} already exists. Skipping sample {sample}.")
            continue
        args.append((read1_path, read2_path, output_path))
    Parallel(n_jobs=cores)(
        delayed(process_sample)(read1_path, read2_path, output_path)
        for read1_path, read2_path, output_path in args
    )
