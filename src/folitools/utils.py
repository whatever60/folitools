from wcmatch import glob


def expand_path_to_list(arg: str | list[str]) -> list[str]:
    """
    Normalize `arg` into a list of strings. If `arg` is a string, treat it as a glob
    pattern, expand it, and return the matching paths. If no matches are found, raise.

    Args:
        arg: Either a list of strings (returned as-is after validation) or a glob-style
            pattern string to expand.

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
        matches.extend(sorted(glob.glob(a, flags=glob.EXTGLOB)))
        if not matches:
            raise ValueError(f"Glob pattern {a!r} did not match any files.")
    return sorted(set(matches))