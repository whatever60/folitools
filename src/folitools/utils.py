import os

from wcmatch import glob
from natsort import natsorted, ns


def expand_path_to_list(
    arg: str | list[str], suffix: str | list[str] = ""
) -> list[str]:
    """
    Normalize `arg` into a list of strings. If `arg` is a string, treat it as a glob
    pattern, expand it, and return the matching paths. If `arg` is a directory,
    add all files in that directory. If no matches are found, raise.

    Args:
        arg: Either a list of strings (returned as-is after validation) or a glob-style
            pattern string to expand, or a directory path.
        suffix: If provided, filter results to only include files ending with this suffix.

    Returns:
        A list of strings.

    Raises:
        ValueError: If `arg` is a string and the glob pattern matches nothing.
        TypeError: If `arg` is neither `str` nor `list[str]`, or if a list contains non-str.
    """
    if isinstance(arg, str):
        arg = [arg]
    if not all(isinstance(s, str) for s in arg):
        raise TypeError("All elements of the list must be strings.")

    matches = []
    for a in arg:
        # Check if the path is an existing directory
        if os.path.isdir(a):
            a = os.path.join(a, "*")  # Expand to all files in the directory
        glob_matches = sorted(glob.glob(a, flags=glob.EXTGLOB))
        matches.extend(glob_matches)

    # Remove duplicates and sort
    matches = natsorted(set(matches), alg=ns.INT | ns.PATH)

    # Apply suffix filter if provided
    if suffix:
        if isinstance(suffix, str):
            suffix = [suffix]
        matches = [
            m
            for m in matches
            if any(m.endswith(s) for s in suffix) and os.path.isfile(m)
        ]

    # Check if we found any matches
    if not matches:
        raise ValueError(
            f"No files found for pattern(s) {arg!r}"
            + (f" with suffix '{suffix}'" if suffix else "")
        )

    return matches
