import numpy as np
import pandas as pd
from joblib import Parallel, delayed


def _rarefy_array(arr: np.ndarray, n: int, k: int, seed: int = 42) -> np.ndarray:
    """Rarefy one row k times so that each row sum up to n."""
    depth = arr.sum()
    rng = np.random.default_rng(seed=seed)
    if n >= depth:  # don't rarefy if expected depth is larger than the actual depth.
        return arr.reshape(1, -1)
    all_elements = np.repeat(np.arange(arr.size), arr)
    return np.stack(
        [
            np.bincount(
                rng.choice(all_elements, size=n, replace=False), minlength=arr.size
            )
            for _ in range(k)
        ],
        axis=0,
    )


def _rarefy(
    df: pd.DataFrame, ref: list[int], repeat_num: int = 20
) -> tuple[pd.DataFrame, list]:
    """Rarefy the dataframe to the reference list or integer.

    The input dataframe should be sample x features, with sample name as index.

    In the returned output, each sample in the input dataframe will be rarefied
    `repeat_num` times, resulting in a dataframe with shape (sample x features x
    repeat_num), except for samples whose reference is themselves, which will be
    repeated only once. Sample name in the returned output will be like
    <original_sample_name>_rarefied_<repeat_num>.
    """
    # make sure all columns in df are positive integers and are in ref.
    if not len(ref) == len(df):
        raise ValueError(
            "The length of ref should be the same as the number of rows in df."
        )

    # get a list of 2d numpy array using joblib parallisim
    res = Parallel(n_jobs=4)(
        delayed(_rarefy_array)(df.iloc[idx].to_numpy(), n, repeat_num)
        for idx, n in enumerate(ref)
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
) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        df_otu_count, names_orig = _rarefy(df_otu_count, ref, rarefying_repeat)
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
