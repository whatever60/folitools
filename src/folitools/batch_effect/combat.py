from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray


def map_counts_with_combat(
    ref_counts: NDArray[np.integer],
    query_counts: NDArray[np.integer],
    *,
    method: Literal["combat_ref", "combat_seq"] = "combat_ref",
    ref_batch_label: str = "ref",
    query_batch_label: str = "query",
    # Optional biological group labels (preserved signal). If None, treated as one group.
    ref_group: NDArray | None = None,
    query_group: NDArray | None = None,
    # ComBat-ref option
    genewise_disp: bool = False,
    # ComBat-seq options (sva::ComBat_seq)
    full_mod: bool = False,
    shrink: bool = False,
    shrink_disp: bool = False,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """
    Batch-correct two bulk RNA-seq count matrices using ComBat-ref or ComBat-seq, returning
    corrected matrices with the same shapes as the inputs (num_samples x num_genes).

    Args:
        ref_counts: Reference raw counts, shape (n_ref_samples, n_genes).
        query_counts: Query raw counts, shape (n_query_samples, n_genes).
        method: "combat_ref" (source local helper_seq.R + Combat_ref.R) or
            "combat_seq" (call sva::ComBat_seq).
        ref_batch_label: Batch label assigned to reference samples in the combined run.
        query_batch_label: Batch label assigned to query samples in the combined run.
        ref_group: Optional group labels for reference samples (length n_ref_samples).
        query_group: Optional group labels for query samples (length n_query_samples).
        genewise_disp: Passed to ComBat_ref(..., genewise.disp=...). Repo examples use FALSE.
        full_mod: Passed to sva::ComBat_seq (only relevant if covariates are provided).
        shrink: Passed to sva::ComBat_seq.
        shrink_disp: Passed to sva::ComBat_seq.

    Returns:
        (ref_corrected, query_corrected): int64 arrays with shapes matching the inputs.

    Notes:
        - R expects a matrix with genes as rows and samples as columns; this wrapper transposes
          for you and then transposes results back.
        - ComBat_seq does NOT have a "reference batch" mode; it adjusts all batches jointly.
        - ComBat_ref (as implemented in that repo) chooses a reference batch internally; this
          wrapper cannot force it unless the R implementation exposes such an argument.
    """
    if method not in {"combat_ref", "combat_seq"}:
        raise ValueError("method must be 'combat_ref' or 'combat_seq'.")

    ref = np.asarray(ref_counts)
    query = np.asarray(query_counts)

    if ref.ndim != 2 or query.ndim != 2:
        raise ValueError("ref_counts and query_counts must both be 2D arrays (samples x genes).")
    if ref.shape[1] != query.shape[1]:
        raise ValueError("ref_counts and query_counts must have the same number of genes (axis=1).")
    if np.any(ref < 0) or np.any(query < 0):
        raise ValueError("Counts must be non-negative.")

    ref_i64 = ref.astype(np.int64, copy=False)
    query_i64 = query.astype(np.int64, copy=False)

    n_ref, n_genes = ref_i64.shape
    n_query = query_i64.shape[0]

    # Combine samples (samples x genes) -> (genes x samples) for R
    combined_sxg = np.concatenate([ref_i64, query_i64], axis=0)
    combined_gxs = combined_sxg.T

    batch_labels = [ref_batch_label] * n_ref + [query_batch_label] * n_query

    if ref_group is None or query_group is None:
        group_labels = [1] * (n_ref + n_query)
    else:
        ref_group_arr = np.asarray(ref_group)
        query_group_arr = np.asarray(query_group)
        if ref_group_arr.shape[0] != n_ref or query_group_arr.shape[0] != n_query:
            raise ValueError("ref_group/query_group must match the number of samples in ref/query.")
        group_labels = np.concatenate([ref_group_arr, query_group_arr]).tolist()

    # Import rpy2 *inside* the function (as requested)
    import rpy2.robjects as ro
    from rpy2.robjects import default_converter
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import FactorVector, IntVector, StrVector
    from rpy2.robjects.numpy2ri import converter as numpy_converter

    def _as_r_int_matrix(x: NDArray[np.int64]) -> ro.Matrix:
        flat = x.flatten(order="F")  # R fills matrices column-major
        return ro.r["matrix"](IntVector(flat), nrow=x.shape[0], ncol=x.shape[1])

    counts_r = _as_r_int_matrix(combined_gxs)
    batch_r = FactorVector(StrVector(batch_labels))
    group_r = FactorVector(StrVector([str(g) for g in group_labels]))

    if method == "combat_seq":
        sva = importr("sva")
        adjusted_r = sva.ComBat_seq(
            counts_r,
            batch=batch_r,
            group=group_r,
            covar_mod=ro.NULL,
            full_mod=full_mod,
            shrink=shrink,
            shrink_disp=shrink_disp,
            gene_subset_n=ro.NULL,
        )

        with localconverter(default_converter + numpy_converter):
            adjusted_gxs = np.asarray(adjusted_r)

        adjusted_sxg = np.maximum(adjusted_gxs, 0).astype(np.int64, copy=False).T
        return adjusted_sxg[:n_ref], adjusted_sxg[n_ref:]

    # method == "combat_ref"
    # Assume helper_seq.R and Combat_ref.R live next to the Python module.
    try:
        module_dir = Path(__file__).resolve().parent
    except NameError:
        module_dir = Path.cwd()

    helper_path = module_dir / "helper_seq.R"
    combat_ref_path = module_dir / "Combat_ref.R"
    if not helper_path.exists():
        raise FileNotFoundError(f"Missing R script: {helper_path}")
    if not combat_ref_path.exists():
        raise FileNotFoundError(f"Missing R script: {combat_ref_path}")

    # Use POSIX paths for R even on Windows
    helper_r = str(helper_path).replace("\\", "/")
    combat_ref_r = str(combat_ref_path).replace("\\", "/")

    ro.r(f'source("{helper_r}")')
    ro.r(f'source("{combat_ref_r}")')

    if "ComBat_ref" not in ro.globalenv:
        raise RuntimeError("R globalenv does not contain ComBat_ref after sourcing Combat_ref.R")

    combat_ref_fun = ro.globalenv["ComBat_ref"]
    adjusted_r = combat_ref_fun(
        counts=counts_r,
        batch=batch_r,
        group=group_r,
        genewise_disp=bool(genewise_disp),
    )

    with localconverter(default_converter + numpy_converter):
        adjusted_gxs = np.asarray(adjusted_r)

    adjusted_sxg = np.maximum(adjusted_gxs, 0).astype(np.int64, copy=False).T
    return adjusted_sxg[:n_ref], adjusted_sxg[n_ref:]
