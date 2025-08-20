import pandas as pd


def order(
    input_: str,
    primer_selection: str,
    primer_info: str,
    output: str,
    has_linker: bool = False,
):
    """
    Generate an Excel file for primer ordering.
    Args:
        input_ (str): Path to gene info tsv file.
        primer_selection (str): Path to selected primers tsv file.
        primer_info (str): Path to candidate primer info tsv file.
        output (str): Path to output Excel file (should end with .xlsx).
    """
    selection = pd.read_table(primer_selection)
    primer_info_df = pd.read_table(primer_info)
    df_gene = pd.read_table(input_).rename(
        {"group": "Group", "gene": "geneSymbol"}, axis=1
    )

    # primer_info["amplicon_index"].is_unique, selection["primerIndex"].is_unique

    assert df_gene["geneSymbol"].is_unique, "geneSymbol should be unique"

    res_df = pd.merge(
        selection,
        primer_info_df,
        left_on="primerIndex",
        right_on="amplicon_index",
        how="left",
    ).merge(df_gene, on="geneSymbol", how="left")
    res_df["Chosen Index"] = res_df.primerIndex.map(lambda x: int(x.split("_")[1]))
    # reorder rows
    # res_df = (
    #     res_df.set_index("geneSymbol")
    #     .loc[[i for i in gene_names if i in res_df.geneSymbol.tolist()]]
    #     .reset_index()
    # )
    # take columns
    res_df = res_df[
        [
            "Group",
            "geneSymbol",
            "geneID",
            "Chosen Index",
            "amplicon_index",
            "L_seq",
            "R_seq",
            "pool",
        ]
    ]

    # format to excel like:
    # amplicon_index  L_seq  R_seq  primer_sequence_to_order_forward  primer_sequence_to_order_reverse
    # ENSMUST00000001184_0  TTTTCACCTACCCCTGCTGTGTTCG  GCAGGGAGACAAAACACTGTAGCCA  CCTACACGACGCTCTTCCGATCTNNNNNNTTTTCACCTACCCCTGCTGTGTTCG  TCAGACGTGTGCTCTTCCGATCTNNNNNNGCAGGGAGACAAAACACTGTAGCCA
    # ENSMUST00000001184_1  CATCAGGTCTGTTCACTGCACCCAA  ACAGTGACTCAGTCCAGCTTCCAGA  CCTACACGACGCTCTTCCGATCTNNNNNNCATCAGGTCTGTTCACTGCACCCAA  TCAGACGTGTGCTCTTCCGATCTNNNNNNACAGTGACTCAGTCCAGCTTCCAGA
    if has_linker:
        FWD_PREFIX = "CCTACACGACGCTCTTCCGATCTNNNNNNACATCA"
        REV_PREFIX = "TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTT"
    else:
        FWD_PREFIX = "CCTACACGACGCTCTTCCGATCTNNNNNN"
        REV_PREFIX = "TCAGACGTGTGCTCTTCCGATCTNNNNNN"

    res_df["primer_sequence_to_order_forward"] = FWD_PREFIX + res_df["L_seq"]
    res_df["primer_sequence_to_order_reverse"] = REV_PREFIX + res_df["R_seq"]
    res_df = res_df.sort_values(by=["pool", "Group", "geneSymbol", "Chosen Index"])
    res_df.to_excel(output, index=False)
    print(res_df["Group"].value_counts())
    return res_df
