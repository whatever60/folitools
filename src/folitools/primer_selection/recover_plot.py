"""Plotting functionality for primer selection reports."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl


def make_report(
    out_pdf: Path,
    amplicon_all: pd.DataFrame,
    amplicon_sub: pd.DataFrame,
    locate_final: pd.DataFrame,
    length_range: tuple[int, int],
) -> None:
    """Create a multi-page PDF report with sanity-check figures.

    Pages:
      1. Amplicon length histogram (log-log) with target range lines.
      2. Venn diagram: primers involved in valid amplicons vs all primers detected.
      3. Histogram: number of transcripts per primer pair (within target length).

    Args:
        out_pdf: Destination PDF path.
        amplicon_all: DataFrame of all enumerated amplicons.
        amplicon_sub: DataFrame filtered by `length_range`.
        locate_final: Locate results joined with primer metadata.
        length_range: Inclusive range used for filtering.
    """
    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    with (
        mpl.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none"}),
        PdfPages(out_pdf) as pdf,
    ):
        # Page 1: amplicon length histogram
        fig1, ax1 = plt.subplots(figsize=(6, 4.5))
        if not amplicon_all.empty:
            x = amplicon_all["amplicon_length"].astype(int).to_numpy()
            hist, bins = np.histogram(x, bins=100)
            logbins = np.logspace(
                np.log10(max(1, bins[0])), np.log10(max(bins[-1], 10)), len(bins)
            )
            ax1.hist(x, bins=logbins)
            ax1.set_xscale("log")
            ax1.set_yscale("log")
        ax1.set_xlabel("Amplicon length (bp)")
        ax1.set_ylabel("Frequency")
        ax1.set_title("Amplicon length distribution")
        # Only draw vertical lines for finite bounds (not -1)
        for v in length_range:
            if v != -1:
                ax1.axvline(v, linestyle="--", alpha=0.5)
        pdf.savefig(fig1)
        plt.close(fig1)

        # Page 2: venn of primers
        # Sanity check: The venn diagram is to check whether all primers are involved in at least one amplicon.
        # fig2, ax2 = plt.subplots(figsize=(5, 5))
        # primers_amp = set(amplicon_sub["primer_seq_fwd"].tolist()) | set(
        #     amplicon_sub["primer_seq_rev"].tolist()
        # )
        # primers_all = set(locate_final["primer_seq"].tolist())
        # venn2((primers_amp, primers_all), set_labels=("Amplicon", "All"), ax=ax2)
        # ax2.set_title("Primers participating in valid amplicons")
        # pdf.savefig(fig2)
        # plt.close(fig2)

        # # Page 3: histogram of number of transcripts per pair
        # fig3, ax3 = plt.subplots(figsize=(6, 4.5))
        # if not amplicon_sub.empty:
        #     counts = amplicon_sub.groupby(["primer_seq_fwd", "primer_seq_rev"])[
        #         "transcript_id"
        #     ].nunique()
        #     ax3.hist(counts.tolist(), bins=np.arange(0.5, counts.max() + 1.5).tolist())
        # ax3.set_xlabel("# transcripts per primer pair (within range)")
        # ax3.set_ylabel("Count")
        # ax3.set_title("Mapping multiplicity")
        # pdf.savefig(fig3)
        # plt.close(fig3)
