import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import logging

logger = logging.getLogger(__name__)


def signif(x, digits=2):
    return float(f"%.{digits}g" % x)


def record_top_sps(df: pd.DataFrame, columns: list[str], output_file: str) -> None:
    """
    Record the top SPs using multiple rank columns
    """

    # Round only float columns to 2 significant figures
    float_columns = df.select_dtypes(include=["float64", "float32"]).columns
    df.loc[:, float_columns] = df.loc[:, float_columns].map(signif)

    # Get SPs in top 20 of any rank column
    top_20_any_idx = []
    for col in columns:
        top_20_any_idx.extend(df[col].nsmallest(20).index)
    # Get unique indices
    top_20_any_idx = list(set(top_20_any_idx))
    top_20_any = df.loc[top_20_any_idx]

    # Top 100 in any rank column
    top_100_any_idx = []
    for col in columns:
        top_100_any_idx.extend(df[col].nsmallest(100).index)
    top_100_any_idx = list(set(top_100_any_idx))
    top_100_any = df.loc[top_100_any_idx]

    # Get SPs in top 100 of chosen rank columns
    top_100_ef_cohen_d_idx = df.loc[:, ["EF rank", "Cohen d rank"]].apply(
        lambda row: all(row <= 100), axis=1
    )
    top_100_ef_cohen_d = df.loc[top_100_ef_cohen_d_idx]

    top_100_ef_p_value_idx = df.loc[:, ["EF rank", "p_value rank"]].apply(
        lambda row: all(row <= 100), axis=1
    )
    top_100_ef_p_value = df.loc[top_100_ef_p_value_idx]

    top_100_hf_percent_cohen_d_idx = df.loc[
        :, ["HF mean percent rank", "Cohen d rank"]
    ].apply(lambda row: all(row <= 100), axis=1)
    top_100_hf_percent_cohen_d = df.loc[top_100_hf_percent_cohen_d_idx]

    top_100_hf_percent_p_value_idx = df.loc[
        :, ["HF mean percent rank", "p_value rank"]
    ].apply(lambda row: all(row <= 100), axis=1)
    top_100_hf_percent_p_value = df.loc[top_100_hf_percent_p_value_idx]

    top_100_all_idx = df.loc[:, columns].apply(lambda row: all(row <= 100), axis=1)
    top_100_all = df.loc[top_100_all_idx]

    # Save to Excel
    with pd.ExcelWriter(output_file) as writer:
        df.to_excel(writer, sheet_name="All", index=False)
        top_20_any.to_excel(writer, sheet_name="Top 20 Any", index=False)
        top_100_any.to_excel(writer, sheet_name="Top 100 Any", index=False)
        top_100_ef_cohen_d.to_excel(
            writer, sheet_name="Top 100 EF and Cohen d", index=False
        )
        top_100_ef_p_value.to_excel(
            writer, sheet_name="Top 100 EF and p-value", index=False
        )
        top_100_hf_percent_cohen_d.to_excel(
            writer, sheet_name="Top 100 HF percent and Cohen d", index=False
        )
        top_100_hf_percent_p_value.to_excel(
            writer, sheet_name="Top 100 HF percent and p", index=False
        )
        top_100_all.to_excel(writer, sheet_name="Top 100 All", index=False)


def add_cohen_d(
    df: pd.DataFrame, hf_columns: list[str], lf_columns: list[str]
) -> pd.DataFrame:
    """
    Calculate Cohen's d for two groups.
    """
    logger.info("Calculating Cohen's d")

    def cohen_d(group1, group2):
        mean1, mean2 = np.mean(group1), np.mean(group2)
        var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
        n1, n2 = len(group1), len(group2)

        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
        if pooled_std == 0:
            return np.nan
        d = (mean1 - mean2) / pooled_std

        return d

    df["Cohen d"] = df.apply(
        lambda row: cohen_d(row[hf_columns], row[lf_columns]),
        axis=1,
    )
    df["Cohen d rank"] = df["Cohen d"].rank(ascending=False, method="min")

    return df


def add_hf_percentages(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the percentage of reads that map to each HF
    """
    logger.info("Calculating mean percentagesof HF reads")
    df["HF mean percent"] = df["HF mean"] / (df["HF mean"] + df["LF mean"])

    df["HF mean percent rank"] = df["HF mean percent"].rank(
        ascending=False, method="min"
    )

    return df


def add_t_test_p_value(
    df: pd.DataFrame, hf_columns: list[str], lf_columns: list[str]
) -> pd.DataFrame:
    """
    Compute the p-value for each SP using a t-test
    """
    logger.info("Calculating p-values of t-test between HFs and LFs")
    try:
        df["p_value"] = df.apply(
            lambda row: ttest_ind(
                row[hf_columns].astype(float),  # Convert to float
                row[lf_columns].astype(float),  # Convert to float
                equal_var=False,
            ).pvalue,
            axis=1,
        )
        df["p_value rank"] = (
            -np.log10(df["p_value"]) * np.sign(df["HF mean"] - df["LF mean"])
        ).rank(ascending=False, method="min")
    except Exception as e:
        print("Could not compute p-values")
        print(e)

    return df
