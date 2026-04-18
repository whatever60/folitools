import numpy as np
import pandas as pd
from joblib import Parallel, delayed


def _rarefy_array(arr: np.ndarray, n: int, k: int, seed: int = 42) -> np.ndarray:
    """Rarefy one sample-level count vector.

    Args:
        arr: One-dimensional count vector for a single sample.
        n: Target sequencing depth after rarefaction.
        k: Number of independent rarefaction repeats to generate.
        seed: Seed used to initialize the NumPy random number generator.

    Returns:
        A 2D array with one rarefied count vector per row. The output shape is
        ``(k, arr.size)`` when rarefaction is performed. If ``n`` is greater
        than or equal to the observed depth, the original counts are returned
        once with shape ``(1, arr.size)``.
    """
    depth = arr.sum()
    rng = np.random.default_rng(seed=seed)
    if n >= depth:  # don't rarefy if expected depth is larger than the actual depth.
        return arr.reshape(1, -1)
    # all_elements = np.repeat(np.arange(arr.size), arr)
    # return np.stack(
    #     [
    #         np.bincount(
    #             rng.choice(all_elements, size=n, replace=False), minlength=arr.size
    #         )
    #         for _ in range(k)
    #     ],
    #     axis=0,
    # )
    return rng.multivariate_hypergeometric(arr, n, size=k)


def _rarefy(
    df: pd.DataFrame, ref: list[int], repeat_num: int = 20, cores: int = 4
) -> tuple[pd.DataFrame, list]:
    """Rarefy each sample in a count table to its corresponding target depth.

    Args:
        df: Count table with samples in rows and features in columns.
        ref: Target rarefaction depth for each row in ``df``. Must have the
            same length as ``df``.
        repeat_num: Number of rarefaction repeats for samples that need to be
            downsampled.
        cores: Number of worker threads used to process rows in parallel.

    Returns:
        A tuple ``(df_rarefied, original_names)``. ``df_rarefied`` contains the
        rarefied count table with row names in the form
        ``<sample_name>_rarefied_<repeat_index>``. Samples whose target depth is
        greater than or equal to their observed depth appear only once.
        ``original_names`` records the original sample name for each returned
        row, in the same order as ``df_rarefied``.

    Raises:
        ValueError: If ``ref`` does not have the same length as ``df``.
    """
    # make sure all columns in df are positive integers and are in ref.
    if not len(ref) == len(df):
        raise ValueError(
            "The length of ref should be the same as the number of rows in df."
        )

    values = df.to_numpy()

    # get a list of 2d numpy array using joblib parallisim
    res = Parallel(n_jobs=cores, prefer="threads")(
        delayed(_rarefy_array)(values[idx], n, repeat_num) for idx, n in enumerate(ref)
    )
    assert isinstance(res, list)
    sample_names_new_orig = [
        (f"{sample_name}_rarefied_{j}", sample_name)
        for idx, sample_name in enumerate(df.index)
        for j in range(res[idx].shape[0])
    ]
    res = np.concatenate(res, axis=0)
    idx_new, idx_orig = zip(*sample_names_new_orig)
    return pd.DataFrame(res, index=idx_new, columns=df.columns), idx_orig


def rarefy(
    df_otu_count: pd.DataFrame,
    df_meta: pd.DataFrame,
    rarefying_repeat: int = 0,
    rarefying_value: None | int = None,
    rarefying_key: None | str = None,
    cores: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rarefy an OTU count table and expand metadata to match the new rows.

    Args:
        df_otu_count: OTU count table with samples in rows and features in
            columns.
        df_meta: Metadata table indexed by sample name.
        rarefying_repeat: Number of rarefaction repeats to generate for each
            sample when rarefaction is enabled. Rarefaction is skipped when this
            value is less than or equal to 1.
        rarefying_value: Default target depth used when ``rarefying_key`` does
            not provide a per-sample reference. If not provided, the default is
            one less than the minimum observed sample depth.
        rarefying_key: Column in ``df_meta`` that defines the rarefaction
            target for each sample. String values refer to another sample whose
            depth should be used. Integer values are treated as explicit target
            depths. Missing values fall back to ``rarefying_value``.
        cores: Number of worker threads used during rarefaction.

    Returns:
        A tuple ``(df_otu_count_rarefied, df_meta_rarefied)``. When rarefaction
        is enabled, both outputs are expanded so that repeated rarefied samples
        and their metadata stay aligned. When rarefaction is disabled, the
        original inputs are returned unchanged.

    Raises:
        ValueError: If ``rarefying_key`` refers to a sample name that is not
            present in ``df_otu_count``.
        ValueError: If a value in ``rarefying_key`` is neither a string, an
            integer, nor a missing value.
    """
    if rarefying_repeat > 1:
        # Rarefying takes place by the following order:
        # - When rarefying_key is specified:
        #    - When this key is a string, rarefy to the depth of the sample corresponding to the key.
        #    - When this key is a int, rarefy to this value.
        #    - When this key is None, rarefy to `rarefying_value` if it is
        #        specified, otherwise rarefy to the minimum depth of all samples.
        # - When rarefying_key isn't specified, treat as if rarefying_key is None
        #     for all samples.
        depth = df_otu_count.sum(axis=1)
        if rarefying_value is None:
            rarefying_value = int(depth.min()) - 1

        ref = []
        if rarefying_key is None:
            rarefying_keys = [np.nan] * len(df_otu_count)
        else:
            rarefying_keys = df_meta[rarefying_key].tolist()

        for i in rarefying_keys:
            if isinstance(i, str):
                if i not in df_otu_count.index:
                    raise ValueError(f"Reference sample {i} is not available.")
                ref.append(depth[i])
            elif isinstance(i, int):
                ref.append(i)
            elif np.isnan(i):
                ref.append(rarefying_value)
            else:
                raise ValueError(
                    f"Unknown data type for rarefying_key: {i} ({type(i)})"
                )

        # rarefy
        df_otu_count, names_orig = _rarefy(
            df_otu_count, ref, rarefying_repeat, cores=cores
        )
        # NOTE:
        # Must create the dummy dataframe as a column, cannot be empty dataframe
        # with index, otherwise order of merged dataframe index will be slightly
        # wrong.
        df_meta = pd.merge(
            pd.DataFrame({"original_sample_name": names_orig}),
            df_meta,
            left_on="original_sample_name",
            right_index=True,
            validate="many_to_one",
        ).reset_index(names=df_meta.index.name)
        df_meta.index = df_otu_count.index
    return df_otu_count, df_meta
