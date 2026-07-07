"""Tests for UMI annotation of cutadapt-renamed FASTQ records."""

from io import StringIO

from folitools.add_umi import add_umi


def test_add_umi_keeps_0_7_header_shape(tmp_path) -> None:
    """UMI annotation should keep the old underscore ID plus space comment."""
    stream = StringIO(
        "\n".join(
            [
                "@READ_with_underscore FWD ACGTAA",
                "GGGACGTAA",
                "+",
                "IIIIIIIII",
                "@READ_with_underscore REV TTGGCC",
                "CCCTTGGCC",
                "+",
                "JJJJJJJJJ",
                "",
            ]
        )
    )
    out1 = tmp_path / "out_1.fq"
    out2 = tmp_path / "out_2.fq"

    add_umi(
        stream,
        str(out1),
        str(out2),
        compression_threads=0,
    )

    assert out1.read_text().splitlines() == [
        "@READ_with_underscore_GGG_CCC FWD+REV",
        "ACGTAA",
        "+",
        "IIIIII",
    ]
    assert out2.read_text().splitlines() == [
        "@READ_with_underscore_GGG_CCC FWD+REV",
        "TTGGCC",
        "+",
        "JJJJJJ",
    ]


def test_add_umi_normalizes_legacy_mate_suffixes(tmp_path) -> None:
    """R1/R2 IDs ending in `/1` and `/2` should still pair cleanly."""
    stream = StringIO(
        "\n".join(
            [
                "@READ_with_underscore/1 FWD ACGTAA",
                "GGGACGTAA",
                "+",
                "IIIIIIIII",
                "@READ_with_underscore/2 REV TTGGCC",
                "CCCTTGGCC",
                "+",
                "JJJJJJJJJ",
                "",
            ]
        )
    )
    out1 = tmp_path / "out_1.fq"
    out2 = tmp_path / "out_2.fq"

    add_umi(
        stream,
        str(out1),
        str(out2),
        compression_threads=0,
    )

    assert (
        out1.read_text().splitlines()[0]
        == "@READ_with_underscore_GGG_CCC FWD+REV"
    )
    assert (
        out2.read_text().splitlines()[0]
        == "@READ_with_underscore_GGG_CCC FWD+REV"
    )
