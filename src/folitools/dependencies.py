"""External dependency reporting and capability checks."""

import re
import shutil
import subprocess

from . import __version__


PROGRAMS = (
    {
        "name": "fastp",
        "command": "fastp",
        "version_args": ("--version",),
        "checks": (
            (
                ("--help",),
                ("--cut_tail", "--correction", "--stdout"),
            ),
        ),
        "required": True,
    },
    {
        "name": "FastQC",
        "command": "fastqc",
        "version_args": ("--version",),
        "checks": ((("--help",), ("--threads", "--outdir")),),
        "required": True,
    },
    {
        "name": "seqkit",
        "command": "seqkit",
        "version_args": ("version",),
        "checks": (
            (("stats", "--help"), ("--all", "--tabular", "--threads")),
            (
                ("locate", "--help"),
                ("--ignore-case", "--max-mismatch", "--use-fmi", "--pattern"),
            ),
        ),
        "required": True,
    },
    {
        "name": "STAR",
        "command": "STAR",
        "version_args": ("--version",),
        "checks": (
            (
                ("--help",),
                ("genomeLoad", "readFilesCommand", "outSAMorder", "chimOutType"),
            ),
        ),
        "required": True,
    },
    {
        "name": "featureCounts",
        "command": "featureCounts",
        "version_args": ("-v",),
        "checks": ((("-h",), ("--fraction", "--Rpath", "--donotsort")),),
        "required": True,
    },
    {
        "name": "samtools",
        "command": "samtools",
        "version_args": ("--version",),
        "checks": (
            (("collate", "--help"), ("-O", "-u", "--threads")),
            (("fastq", "--help"), ("-1", "-2")),
        ),
        "required": True,
    },
    {
        "name": "sambamba",
        "command": "sambamba",
        "version_args": ("--version",),
        "checks": ((("sort",), ("--nthreads", "--memory-limit", "--out")),),
        "required": True,
    },
    {
        "name": "cutadapt",
        "command": "cutadapt",
        "version_args": ("--version",),
        "checks": (
            (
                ("--help",),
                (
                    "--interleaved",
                    "--json",
                    "--rename",
                    "--action",
                    "--minimum-length",
                ),
            ),
        ),
        "required": True,
    },
    {
        "name": "UMI-tools",
        "command": "umi_tools",
        "version_args": ("--version",),
        "checks": (
            (
                ("group", "--help"),
                (
                    "--paired",
                    "--group-out",
                    "--gene-tag",
                    "--cell-tag-split",
                    "--skip-tags-regex",
                ),
            ),
        ),
        "required": True,
    },
    {
        "name": "pigz",
        "command": "pigz",
        "version_args": ("--version",),
        "checks": (),
        "required": False,
    },
    {
        "name": "bwa-mem2",
        "command": "bwa-mem2",
        "version_args": ("version",),
        "checks": (),
        "required": False,
    },
)


def _command_output(path: str, args: tuple[str, ...]) -> str:
    """Run one dependency probe and combine its standard output and error."""
    result = subprocess.run(
        [path, *args],
        check=False,
        capture_output=True,
        text=True,
    )
    return f"{result.stdout}\n{result.stderr}"


def check_dependencies() -> bool:
    """Print resolved dependency paths and capabilities, returning overall health."""
    healthy = True
    print(f"folitools {__version__}")

    for program in PROGRAMS:
        path = shutil.which(program["command"])
        if path is None:
            if program["required"]:
                healthy = False
                print(f"MISSING {program['name']}")
            else:
                print(f"OPTIONAL-MISSING {program['name']}")
            continue

        version_output = _command_output(path, program["version_args"])
        version_match = re.search(r"\d+(?:\.\d+)+(?:[A-Za-z0-9.+-]*)?", version_output)
        version = version_match.group(0) if version_match else "unknown"
        missing_capabilities = []
        for args, capabilities in program["checks"]:
            help_output = _command_output(path, args)
            missing_capabilities.extend(
                capability
                for capability in capabilities
                if capability not in help_output
            )

        if missing_capabilities:
            healthy = False
            missing = ", ".join(missing_capabilities)
            print(f"INCOMPATIBLE {program['name']} {version} {path} missing: {missing}")
        else:
            label = "OK" if program["required"] else "OPTIONAL"
            print(f"{label} {program['name']} {version} {path}")

    return healthy
