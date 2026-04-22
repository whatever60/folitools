import sys
from typing import Annotated, TextIO

import pysam
from cyclopts import run, Parameter


def add_tags_wo_fastq(
    bam_input: str,
    bam_output: str,
    umi_tag_s_name: str = "US",  # single-read UMI
    umi_tag_s_correct_name: str = "UC",  # single-read UMI (correct ones only)
    cell_tag_name: str = "CB",
    primer_tag_name: str = "PR",
    gene_with_primer_tag_name: str = "XF",
    cell_tag: str | None = None,
    umi_length: int = 6,
    log_fh: TextIO | None = None,
) -> None:
    """
    Reads SAM from stdin, parses embedded UMI sequences from the read names,
    adds UMI and optional cell tags, and writes a BAM file.

    VERY IMPORTANT NOTE:
        The read IDs must be formatted as:
            <original_read_id>_{umi_r1}_{umi_r2}
        where the two UMI sequences are appended to the original read name,
        separated by underscores. These UMIs will be split out and tagged
        appropriately.

    Args:
        bam_input: Standard input stream for SAM data.
        bam_output: Path to the output BAM file ('-' or None for stdout).
        umi_tag_s_name: Name of the single-read UMI tag (default: 'US').
        umi_tag_p_name: Name of the paired-read UMI tag (default: 'UP').
        cell_tag_name: Name of the cell barcode tag (default: 'CB').
        adapter_tag_name: Name of the adapter tag (unused here) (default: 'XA').
        cell_tag: Optional cell barcode value to add to every read.
        umi_length: Expected length of each individual UMI (default: 6).
    """
    if bam_input == "-":
        bam_input = "/dev/stdin"
    mode_out = "w"
    if bam_output == "-":
        mode_out = "wb0"  # same as wbu, i.e., uncompressed BAM
    elif bam_output.endswith(".bam"):
        mode_out = "wb"
    elif bam_output.endswith(".sam"):
        mode_out = "w"
    else:
        raise ValueError("Output file must end with .bam or .sam")

    def process_and_write_pair():
        """Write the current primary read pair with a gene-with-primer tag."""
        # SAM spec: a record is primary iff FLAG & 0x900 == 0. Because mates
        # share a QNAME, the invariant is per-mate: at most one primary R1 and
        # one primary R2. STAR can, rarely, emit more than one such record per
        # mate (non-compliant output). Treat that as a warning, not a fatal
        # error: log it, and downgrade extras to secondary (flag |= 0x100) so
        # the output BAM stays SAM-compliant and downstream (umi_tools) sees
        # exactly one primary per mate.
        r1_entries = [r for k, r in primary_alignments if k == "R1"]
        r2_entries = [r for k, r in primary_alignments if k == "R2"]
        if len(r1_entries) != 1 or len(r2_entries) != 1:
            if log_fh is not None:
                flag_dump = ", ".join(
                    f"{k} flag=0x{r.flag:x} ref={r.reference_name}:{r.reference_start}"
                    for k, r in primary_alignments
                )
                log_fh.write(
                    f"WARNING non-compliant primary count for {read_id}: "
                    f"R1={len(r1_entries)} R2={len(r2_entries)} [{flag_dump}] "
                    f"(next read: {query_name_current})\n"
                )
            keep = set()
            if r1_entries:
                keep.add(id(r1_entries[0]))
            if r2_entries:
                keep.add(id(r2_entries[0]))
            for _, r in primary_alignments:
                if id(r) not in keep:
                    r.flag |= 0x100

        # Join unique XS tags, filter out any "Unassigned" tags if others exist
        genes = {
            gene for tag in xt_tags if tag != "Unassigned" for gene in tag.split(",")
        }
        if genes:
            gene_with_primer = ",".join(sorted(genes) + [primers])
        else:
            gene_with_primer = ",".join(["Unassigned", primers])
        for _, read in primary_alignments:
            if not (read.is_secondary or read.is_supplementary):
                read.set_tag(
                    gene_with_primer_tag_name, gene_with_primer, value_type="Z"
                )
            bam_out.write(read)

    with (
        # NOTE: pysam read can handle both BAM and SAM format automatically if opened with "rb".
        pysam.AlignmentFile(bam_input, "rb") as sam_in,
        pysam.AlignmentFile(bam_output, mode_out, header=sam_in.header) as bam_out,
    ):
        primary_alignments: list[tuple[str, pysam.AlignedSegment]] = []
        xt_tags: list[str] = []
        read_id: str | None = None
        primers: str = ""
        # Pull per-call knobs into locals — trivially faster than repeated
        # closure/attr lookup in the hot loop.
        cell_tag_set = cell_tag is not None
        for read in sam_in.fetch(until_eof=True):
            query_name_current = read.query_name
            if read_id is not None and read_id != query_name_current:
                # New QNAME: flush the previous group.
                process_and_write_pair()
                primary_alignments = []
                xt_tags = []
                read_id = None

            if read_id is None:
                read_id = query_name_current

            # Parse QNAME once. All four splits are cheap; the goal below is
            # to avoid the per-record pysam attribute roundtrips, not this.
            query_name_simple, umi1, umi2, primers = query_name_current.split("_", 3)
            primer_fwd, primer_rev = primers.split("+")
            read.query_name = query_name_simple

            # One flag read; derive mate + primary-ness from bits. Avoids
            # read.is_read1 / is_secondary / is_supplementary (each is a
            # property that re-reads .flag).
            flag = read.flag
            is_read1 = (flag & 0x40) != 0
            is_not_primary = (flag & 0x900) != 0

            if cell_tag_set:
                read.set_tag(cell_tag_name, cell_tag, value_type="Z")

            read.set_tag(
                umi_tag_s_name, umi1 if is_read1 else umi2, value_type="Z"
            )
            read.set_tag(primer_tag_name, primers, value_type="Z")

            # Short-circuit chain beats building a list + all().
            if (
                len(umi1) == umi_length
                and len(umi2) == umi_length
                and "N" not in umi1
                and "N" not in umi2
                and primer_fwd != "no_adapter"
                and primer_rev != "no_adapter"
            ):
                read.set_tag(umi_tag_s_correct_name, umi1 + umi2, value_type="Z")

            # XT fallback. We also cache the value we'll need for XF so we
            # don't call get_tag a second time at the bottom of the loop.
            if read.has_tag("XT"):
                xt_value = read.get_tag("XT")
            else:
                # featureCounts did not assign — should be an "Unassigned_*"
                # case. Debug assertion was dropped; trust featureCounts.
                read.set_tag("XN", -1, value_type="i")
                read.set_tag("XT", "Unassigned", value_type="Z")
                xt_value = "Unassigned"

            if is_not_primary:
                bam_out.write(read)
            else:
                primary_alignments.append(("R1" if is_read1 else "R2", read))
            xt_tags.append(xt_value)

        # Wrap up the last pair
        if read_id is not None:
            process_and_write_pair()


def main(
    input_: Annotated[
        str,
        Parameter(
            name=["--input", "-i"],
            help="Path to the input BAM file ('-' or omit for stdin)",
        ),
    ] = "-",
    output: Annotated[
        str,
        Parameter(
            name=["--output", "-o"],
            help="Path to the output BAM file ('-' or omit for stdout)",
        ),
    ] = "-",
    cell_tag_name: Annotated[
        str,
        Parameter(
            name="--cell_tag_name",
            help="Name of the cell tag to add (default: 'CB')",
        ),
    ] = "CB",
    cell_tag: Annotated[
        str | None,
        Parameter(
            name="--cell_tag",
            help="Value of the cell tag to add (default: None)",
        ),
    ] = None,
    log: Annotated[
        str | None,
        Parameter(
            name="--log",
            help=(
                "Path to a log file for non-fatal warnings (e.g. SAM-compliance "
                "anomalies). If omitted, such warnings are silently dropped. "
                "Fatal errors still go to stderr."
            ),
        ),
    ] = None,
) -> int:
    """
    Dispatch to add_tags or add_tags_wo_fastq based on presence of --index_fastq.

    All existing comments and logic are preserved.
    """
    log_fh = open(log, "a") if log else None
    try:
        try:
            add_tags_wo_fastq(
                bam_input=input_,
                bam_output=output,
                cell_tag_name=cell_tag_name,
                cell_tag=cell_tag,
                log_fh=log_fh,
            )
        except Exception as e:
            print(e, file=sys.stderr)
            return 1
    finally:
        if log_fh is not None:
            log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(run(main))
