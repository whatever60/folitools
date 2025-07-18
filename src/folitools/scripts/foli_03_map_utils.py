import math

from cyclopts import run

def allocate_pipeline_cores(total_cores: int) -> tuple[int, int, int]:
    """
    Calculates an optimal core allocation for a bioinformatics pipeline.

    This function distributes a total number of cores among three processes:
    1. STAR (a multi-threaded aligner)
    2. A custom Python script (assumed to be single-threaded)
    3. featureCounts (a multi-threaded feature counter)

    Args:
        total_cores: The total number of CPU cores available for the pipeline.

    Returns:
        A tuple of three integers:
        (star_cores, python_script_cores, featurecounts_cores).

    Raises:
        ValueError: If total_cores is less than 3, as this is insufficient
                    to run the pipeline with this allocation strategy.
    """
    if total_cores < 3:
        raise ValueError(
            "A minimum of 3 cores is required to run the pipeline "
            "(1 for STAR, 1 for Python, 1 for featureCounts)."
        )

    # Reserve one core for the Python script
    python_script_cores = 1
    available_for_apps = total_cores - python_script_cores

    # Proportional split
    star_proportion = 0.7
    featurecounts_proportion = 0.3

    star_cores = math.floor(available_for_apps * star_proportion)
    featurecounts_cores = math.ceil(available_for_apps * featurecounts_proportion)

    # Sanity checks to ensure at least one core each
    if star_cores < 1:
        star_cores = 1
        featurecounts_cores = available_for_apps - star_cores

    if featurecounts_cores < 1:
        featurecounts_cores = 1
        star_cores = available_for_apps - featurecounts_cores

    # output to stdout
    print(star_cores, python_script_cores, featurecounts_cores)
    return star_cores, python_script_cores, featurecounts_cores


if __name__ == "__main__":
    run(allocate_pipeline_cores)
