import math

from cyclopts import run


def allocate_pipeline_cores(total_cores: int) -> tuple[int, int]:
    """
    Calculates an optimal core allocation for a bioinformatics pipeline.

    This function distributes a total number of cores among two processes:
    1. STAR (a multi-threaded aligner)
    2. featureCounts (a multi-threaded feature counter)

    Args:
        total_cores: The total number of CPU cores available for the pipeline.

    Returns:
        A tuple of two integers:
        (star_cores, featurecounts_cores).

    Raises:
        ValueError: If total_cores is less than 2, as this is insufficient
                    to run the pipeline with this allocation strategy.
    """
    if total_cores < 2:
        raise ValueError(
            "A minimum of 2 cores is required to run the pipeline "
            "(1 for STAR, 1 for featureCounts)."
        )

    # Proportional split
    star_proportion = 0.7
    featurecounts_proportion = 0.3

    star_cores = math.floor(total_cores * star_proportion)
    featurecounts_cores = math.ceil(total_cores * featurecounts_proportion)

    # Sanity checks to ensure at least one core each
    if star_cores < 1:
        star_cores = 1
        featurecounts_cores = total_cores - star_cores

    if featurecounts_cores < 1:
        featurecounts_cores = 1
        star_cores = total_cores - featurecounts_cores

    # output to stdout
    print(star_cores, featurecounts_cores)
    return star_cores, featurecounts_cores


if __name__ == "__main__":
    run(allocate_pipeline_cores)
