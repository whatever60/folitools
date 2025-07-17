import math

def allocate_pipeline_cores(total_cores: int) -> dict:
    """
    Calculates an optimal core allocation for a bioinformatics pipeline.

    This function distributes a total number of cores among three processes:
    1. STAR (a multi-threaded aligner)
    2. A custom Python script (assumed to be single-threaded)
    3. featureCounts (a multi-threaded feature counter)

    The user is responsible for ensuring the OS has sufficient resources.

    The strategy is as follows:
    - A minimum of 3 cores is required to run the pipeline.
    - 1 core is always assigned to the single-threaded Python script.
    - The remaining cores are distributed proportionally between STAR (≈70%)
      and featureCounts (≈30%), ensuring each gets at least one core.

    Args:
        total_cores: The total number of CPU cores available for the pipeline.

    Returns:
        A dictionary containing the number of cores allocated to 'star',
        'python_script', and 'featurecounts'.
        
    Raises:
        ValueError: If total_cores is less than 3, as this is insufficient
                    to run the pipeline with this allocation strategy.
    """
    if total_cores < 3:
        raise ValueError("A minimum of 3 cores is required to run the pipeline (1 for STAR, 1 for Python, 1 for featureCounts).")

    # --- Core Allocation Logic ---

    # 1. Reserve one core for the single-threaded Python script.
    python_script_cores = 1
    
    available_for_apps = total_cores - python_script_cores

    # 2. Define the proportional split for the multi-threaded applications.
    # STAR is typically the most demanding part of the pipeline.
    star_proportion = 0.7  # 70% for STAR
    featurecounts_proportion = 0.3  # 30% for featureCounts

    # 3. Calculate the core counts based on the proportions.
    # Use math.floor for one and math.ceil for the other to ensure the sum
    # of the parts equals the whole available pool, preventing off-by-one errors.
    star_cores = math.floor(available_for_apps * star_proportion)
    featurecounts_cores = math.ceil(available_for_apps * featurecounts_proportion)

    # 4. Sanity check: Ensure both STAR and featureCounts get at least one core.
    # This is important for low total_cores values (e.g., 3 or 4).
    if star_cores < 1:
        star_cores = 1
        # Recalculate the other to not exceed the total
        featurecounts_cores = available_for_apps - star_cores

    if featurecounts_cores < 1:
        featurecounts_cores = 1
        # Recalculate the other
        star_cores = available_for_apps - featurecounts_cores

    # --- Construct the final allocation dictionary ---
    allocation = {
        'star': star_cores,
        'python_script': python_script_cores,
        'featurecounts': featurecounts_cores,
    }

    return allocation

# --- Example Usage ---
if __name__ == "__main__":
    # Test with different total core counts
    core_scenarios = [3, 4, 8, 16, 32, 64]

    for cores in core_scenarios:
        try:
            core_allocation = allocate_pipeline_cores(cores)
            print(f"--- Allocation for {cores} total cores ---")
            print(f"  STAR:          {core_allocation['star']} cores")
            print(f"  Python Script: {core_allocation['python_script']} core")
            print(f"  featureCounts: {core_allocation['featurecounts']} cores")
            total_assigned = sum(core_allocation.values())
            print("  ---------------------------------")
            print(f"  Total Assigned:  {total_assigned} cores\n")
        except ValueError as e:
            print(f"Could not allocate for {cores} cores: {e}\n")

    # Example of handling the exception for an invalid number of cores
    try:
        allocate_pipeline_cores(2)
    except ValueError as e:
        print("--- Testing error handling ---")
        print(f"Correctly caught error for 2 cores: {e}")
