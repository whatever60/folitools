from typing import Literal

import numpy as np
import pandas as pd


def deseq2_vst_transform(
    ref_counts: np.ndarray,
    query_counts: np.ndarray,
    *,
    method: Literal["A", "B"] = "B",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply DESeq2's variance stabilizing transformation (VST) to map a query bulk RNA-seq
    count matrix onto a reference, returning transformed matrices with the same shapes
    as the inputs (num_samples x num_genes).

    Method A ("A"): Fit size factors + dispersions on the combined (ref + query) dataset
    with design ~ dataset, then VST and split back.
    Method B ("B"): Fit on reference only (design ~ 1), freeze the dispersion function,
    then apply VST to query using the frozen dispersion trend.

    Args:
        ref_counts: Reference raw counts with shape (n_ref_samples, n_genes).
        query_counts: Query raw counts with shape (n_query_samples, n_genes).
        method: Either "A" or "B" (see above).

    Returns:
        A tuple (ref_vst, query_vst), each a float numpy array with the same shape as the
        corresponding input matrix (num_samples x num_genes).

    Raises:
        ValueError: If inputs are not 2D, have different gene counts, contain negatives,
            or method is not "A"/"B".
        ImportError/RuntimeError: If rpy2/R/DESeq2 are not available or R execution fails.

    Notes:
        - Inputs must be non-negative raw counts. Values are cast to int64 for DESeq2.
        - This is intended for visualization / clustering / downstream ML feature spaces,
          not differential expression testing.
    """
    if method not in {"A", "B"}:
        raise ValueError("method must be 'A' or 'B'.")

    ref_counts = np.asarray(ref_counts)
    query_counts = np.asarray(query_counts)

    if ref_counts.ndim != 2 or query_counts.ndim != 2:
        raise ValueError("ref_counts and query_counts must both be 2D arrays.")
    if ref_counts.shape[1] != query_counts.shape[1]:
        raise ValueError(
            "ref_counts and query_counts must have the same num_genes (axis=1)."
        )
    if np.any(ref_counts < 0) or np.any(query_counts < 0):
        raise ValueError("Counts must be non-negative.")

    # Cast to integers (DESeq2 expects counts)
    ref_counts_int = ref_counts.astype(np.int64, copy=False)
    query_counts_int = query_counts.astype(np.int64, copy=False)

    n_ref_samples, n_genes = ref_counts_int.shape
    n_query_samples, _ = query_counts_int.shape

    gene_names = [f"g{i}" for i in range(n_genes)]
    ref_sample_names = [f"ref{i}" for i in range(n_ref_samples)]
    query_sample_names = [f"query{i}" for i in range(n_query_samples)]

    # Build gene x sample matrices for R (DESeq2 uses genes as rows, samples as cols).
    # Flatten in column-major order so R's matrix(fill-by-column) matches.
    ref_vec = ref_counts_int.T.reshape(-1, order="F")
    query_vec = query_counts_int.T.reshape(-1, order="F")

    # Import rpy2 *inside* the function as requested
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects import numpy2ri, pandas2ri

    r_map_vst = ro.r(
        """
        function(ref_vec, query_vec, n_genes, n_ref, n_query,
                 gene_names, ref_names, query_names, method) {

          suppressPackageStartupMessages(library(DESeq2))

          ref <- matrix(ref_vec, nrow=n_genes, ncol=n_ref)
          query <- matrix(query_vec, nrow=n_genes, ncol=n_query)

          rownames(ref) <- gene_names
          colnames(ref) <- ref_names
          rownames(query) <- gene_names
          colnames(query) <- query_names

          if (method == "A") {
            counts <- cbind(ref, query)
            col_data <- data.frame(
              dataset = factor(c(rep("ref", n_ref), rep("query", n_query))),
              row.names = colnames(counts)
            )

            dds <- DESeqDataSetFromMatrix(countData=counts, colData=col_data, design=~dataset)
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersions(dds, quiet=TRUE)

            vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
            mat <- t(assay(vsd))  # samples x genes

            ref_out <- mat[ref_names, , drop=FALSE]
            query_out <- mat[query_names, , drop=FALSE]
            return(list(ref=ref_out, query=query_out))
          }

          if (method == "B") {
            # Fit on reference only (frozen dispersion trend)
            col_data_ref <- data.frame(row.names = ref_names)
            dds_ref <- DESeqDataSetFromMatrix(countData=ref, colData=col_data_ref, design=~1)
            dds_ref <- estimateSizeFactors(dds_ref)
            dds_ref <- estimateDispersions(dds_ref, quiet=TRUE)
            disp_fun <- dispersionFunction(dds_ref)

            vsd_ref <- varianceStabilizingTransformation(dds_ref, blind=FALSE)
            ref_out <- t(assay(vsd_ref))  # samples x genes

            # Apply frozen trend to query
            col_data_q <- data.frame(row.names = query_names)
            dds_q <- DESeqDataSetFromMatrix(countData=query, colData=col_data_q, design=~1)
            dds_q <- estimateSizeFactors(dds_q)
            dispersionFunction(dds_q) <- disp_fun

            vsd_q <- varianceStabilizingTransformation(dds_q, blind=FALSE)
            query_out <- t(assay(vsd_q))  # samples x genes

            return(list(ref=ref_out, query=query_out))
          }

          stop("Unknown method")
        }
        """
    )

    with localconverter(ro.default_converter + numpy2ri.converter):
        r_ref_vec = ro.conversion.py2rpy(ref_vec)
        r_query_vec = ro.conversion.py2rpy(query_vec)

    res = r_map_vst(
        r_ref_vec,
        r_query_vec,
        int(n_genes),
        int(n_ref_samples),
        int(n_query_samples),
        ro.StrVector(gene_names),
        ro.StrVector(ref_sample_names),
        ro.StrVector(query_sample_names),
        method,
    )

    with localconverter(ro.default_converter + pandas2ri.converter):
        ref_vst: np.ndarray = ro.conversion.rpy2py(res.rx2("ref"))
        query_vst: np.ndarray = ro.conversion.rpy2py(res.rx2("query"))

    if ref_vst.shape != ref_counts.shape or query_vst.shape != query_counts.shape:
        raise RuntimeError(
            "Unexpected output shape from DESeq2 VST (this should not happen)."
        )

    return ref_vst, query_vst


def deseq2_vst_transform_baseline(
    ref_counts: np.ndarray,
    query_counts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Baseline DESeq2 VST that is *not* batch-aware.

    This fits a single VST on the combined (ref + query) counts using design ~ 1 and
    blind=TRUE, then returns ref/query transformed matrices with the same shapes as
    the inputs (num_samples x num_genes).

    Args:
        ref_counts: Reference raw counts with shape (n_ref_samples, n_genes).
        query_counts: Query raw counts with shape (n_query_samples, n_genes).

    Returns:
        (ref_vst, query_vst): float64 arrays with shapes matching the inputs.

    Raises:
        ValueError: If shapes are incompatible or counts are negative.
    """
    ref_counts = np.asarray(ref_counts)
    query_counts = np.asarray(query_counts)

    if ref_counts.ndim != 2 or query_counts.ndim != 2:
        raise ValueError("ref_counts and query_counts must both be 2D arrays.")
    if ref_counts.shape[1] != query_counts.shape[1]:
        raise ValueError(
            "ref_counts and query_counts must have the same num_genes (axis=1)."
        )
    if np.any(ref_counts < 0) or np.any(query_counts < 0):
        raise ValueError("Counts must be non-negative.")

    ref_counts_int = ref_counts.astype(np.int64, copy=False)
    query_counts_int = query_counts.astype(np.int64, copy=False)

    n_ref_samples, n_genes = ref_counts_int.shape
    n_query_samples, _ = query_counts_int.shape

    gene_names = [f"g{i}" for i in range(n_genes)]
    ref_sample_names = [f"ref{i}" for i in range(n_ref_samples)]
    query_sample_names = [f"query{i}" for i in range(n_query_samples)]

    # R expects genes x samples; flatten in column-major order.
    ref_vec = ref_counts_int.T.reshape(-1, order="F")
    query_vec = query_counts_int.T.reshape(-1, order="F")

    # Import rpy2 inside the function
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    r_vst_baseline = ro.r(
        """
        function(ref_vec, query_vec, n_genes, n_ref, n_query,
                 gene_names, ref_names, query_names) {

          suppressPackageStartupMessages(library(DESeq2))

          ref <- matrix(ref_vec, nrow=n_genes, ncol=n_ref)
          query <- matrix(query_vec, nrow=n_genes, ncol=n_query)

          rownames(ref) <- gene_names
          colnames(ref) <- ref_names
          rownames(query) <- gene_names
          colnames(query) <- query_names

          counts <- cbind(ref, query)
          col_data <- data.frame(row.names = colnames(counts))

          # Baseline: no batch/condition in design, and blind=TRUE
          dds <- DESeqDataSetFromMatrix(countData=counts, colData=col_data, design=~1)
          dds <- estimateSizeFactors(dds)
          dds <- estimateDispersions(dds, quiet=TRUE)

          vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
          mat <- t(assay(vsd))  # samples x genes

          ref_out <- mat[ref_names, , drop=FALSE]
          query_out <- mat[query_names, , drop=FALSE]
          return(list(ref=ref_out, query=query_out))
        }
        """
    )

    with localconverter(ro.default_converter + numpy2ri.converter):
        r_ref_vec = ro.conversion.py2rpy(ref_vec)
        r_query_vec = ro.conversion.py2rpy(query_vec)

    res = r_vst_baseline(
        r_ref_vec,
        r_query_vec,
        int(n_genes),
        int(n_ref_samples),
        int(n_query_samples),
        ro.StrVector(gene_names),
        ro.StrVector(ref_sample_names),
        ro.StrVector(query_sample_names),
    )

    with localconverter(ro.default_converter + pandas2ri.converter):
        ref_vst: np.ndarray = ro.conversion.rpy2py(res.rx2("ref"))
        query_vst: np.ndarray = ro.conversion.rpy2py(res.rx2("query"))

    if ref_vst.shape != ref_counts.shape or query_vst.shape != query_counts.shape:
        raise RuntimeError(
            "Unexpected output shape from DESeq2 VST (this should not happen)."
        )
    return ref_vst, query_vst
