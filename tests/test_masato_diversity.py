import numpy as np
import pandas as pd

from src.folitools.masato.diversity import _rarefy
from src.folitools.masato.diversity import _rarefy_array
from src.folitools.masato.diversity import rarefy


def test_rarefy_array_preserves_requested_depth():
    """Rarefied rows should keep the requested total count."""
    arr = np.array([5, 3, 2, 0], dtype=np.int64)
    res = _rarefy_array(arr, n=4, k=6, seed=42)

    assert res.shape == (6, 4)
    assert np.all(res.sum(axis=1) == 4)
    assert np.all(res <= arr)


def test_rarefy_repeats_only_samples_that_need_rarefaction():
    """Samples at or below the target depth should only appear once."""
    df = pd.DataFrame(
        [[4, 3, 1], [2, 1, 0]],
        index=["sample_a", "sample_b"],
        columns=["otu_1", "otu_2", "otu_3"],
    )

    res, names_orig = _rarefy(df, ref=[4, 10], repeat_num=3)

    assert res.shape == (4, 3)
    assert list(res.index) == [
        "sample_a_rarefied_0",
        "sample_a_rarefied_1",
        "sample_a_rarefied_2",
        "sample_b_rarefied_0",
    ]
    assert list(names_orig) == ["sample_a", "sample_a", "sample_a", "sample_b"]
    assert (res.sum(axis=1).tolist()) == [4, 4, 4, 3]


def test_rarefy_accepts_cores_argument():
    """Public rarefy should let callers choose how many workers to use."""
    df_otu_count = pd.DataFrame(
        [[4, 3, 1], [2, 1, 0]],
        index=["sample_a", "sample_b"],
        columns=["otu_1", "otu_2", "otu_3"],
    )
    df_meta = pd.DataFrame({"group": ["x", "y"]}, index=df_otu_count.index)

    res_count, res_meta = rarefy(df_otu_count, df_meta, rarefying_repeat=3, cores=1)

    assert res_count.shape == (6, 3)
    assert res_count.sum(axis=1).tolist() == [2, 2, 2, 2, 2, 2]
    assert res_meta["group"].tolist() == ["x", "x", "x", "y", "y", "y"]


def test_rarefy_seed_is_reproducible():
    """Public rarefy should give stable output for the same seed."""
    df_otu_count = pd.DataFrame(
        [[20, 15, 10, 5], [18, 9, 6, 3]],
        index=["sample_a", "sample_b"],
        columns=["otu_1", "otu_2", "otu_3", "otu_4"],
    )
    df_meta = pd.DataFrame({"group": ["x", "y"]}, index=df_otu_count.index)

    res_count_1, res_meta_1 = rarefy(
        df_otu_count, df_meta, rarefying_repeat=4, cores=1, seed=123
    )
    res_count_2, res_meta_2 = rarefy(
        df_otu_count, df_meta, rarefying_repeat=4, cores=1, seed=123
    )
    res_count_3, _ = rarefy(df_otu_count, df_meta, rarefying_repeat=4, cores=1, seed=456)

    pd.testing.assert_frame_equal(res_count_1, res_count_2)
    pd.testing.assert_frame_equal(res_meta_1, res_meta_2)
    assert not res_count_1.equals(res_count_3)
