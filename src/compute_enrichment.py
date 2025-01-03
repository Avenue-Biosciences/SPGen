import pandas as pd
import argparse
import re
import itertools
import os
import json
from sqlalchemy import text, Engine
import logging
from visualisations import *
from multiple_replicates import *
from top_sp_stats import *
from db_utils import *

logger = logging.getLogger(__name__)


def read_input(input_path: str) -> pd.DataFrame:
    """
    Read the input JSON file and validate required fields

    Args:
        input_path: Path to the JSON input file

    Returns:
        dict: Validated input configuration

    Raises:
        ValueError: If required fields are missing
    """
    with open(input_path, "r") as f:
        input = json.load(f)

    # Check required fields
    required_fields = ["protein", "target_directory", "library"]
    missing_fields = [field for field in required_fields if field not in input]

    if missing_fields:
        raise ValueError(
            f"Missing required fields in input file: {', '.join(missing_fields)}"
        )

    return input


def normalise_read_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise read counts by:
    1. Computing sum of reads in LF and HF columns
    2. Dividing each column by its sum
    3. Multiplying by the average of original column sums

    Args:
        df: DataFrame containing HF and LF columns

    Returns:
        DataFrame with normalised read counts
    """
    logger.info("Normalising read counts")
    # Create a copy to avoid modifying the original
    df = df.copy()

    # Get HF and LF columns
    hf_columns = [col for col in df.columns if col.startswith("HF")]
    lf_columns = [col for col in df.columns if col.startswith("LF")]

    # Calculate column sums and their average
    column_sums = df[hf_columns + lf_columns].sum()
    average_sum = column_sums.mean()

    # Normalize each column and multiply by average sum
    for col in hf_columns + lf_columns:
        df[col + "_norm"] = df[col] / column_sums[col] * average_sum

    return df


def get_common_sps(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """
    For each pair of EF columns, get the number of common SPs in top 20 and top 100
    """
    logger.info("Recording common SPs in top ranks between replicates")
    # Generate pairs of columns using itertools
    column_pairs = list(itertools.combinations(columns, 2))

    # Create a DataFrame with the pairs
    results = pd.DataFrame(column_pairs, columns=["Column1", "Column2"])
    results["Comparison"] = results.apply(
        lambda row: f"{row['Column1']} vs {row['Column2']}", axis=1
    )

    # Function to count common elements
    def count_common_elements(row, threshold):
        ef1 = df[row["Column1"]]
        ef2 = df[row["Column2"]]
        return ((ef1 <= threshold) & (ef2 <= threshold)).sum()

    # Apply the function to get n_common_top_20 and n_common_top_100
    results["n_common_top_20"] = results.apply(
        lambda row: count_common_elements(row, 20), axis=1
    )
    results["n_common_top_100"] = results.apply(
        lambda row: count_common_elements(row, 100), axis=1
    )

    # Drop the temporary columns
    results = results.drop(["Column1", "Column2"], axis=1)

    return results


def add_enrichment_factors(
    df: pd.DataFrame, hf_columns: list[str], lf_columns: list[str]
) -> pd.DataFrame:
    """
    Compute the enrichment factors and their ranks using normalized read counts
    """
    logger.info("Computing enrichment factors")
    if not hf_columns or not lf_columns:
        raise ValueError(
            "No normalized HF or LF columns found. Run normalise_read_counts first."
        )
    if len(hf_columns) != len(lf_columns):
        raise ValueError("Number of HF columns and LF columns must be the same")

    for hf_column, lf_column in zip(hf_columns, lf_columns):
        df[f"EF{hf_column[2:3]}"] = (df[hf_column] + 1) / (df[lf_column] + 1)

    ef_columns = [col for col in df.columns if col.startswith("EF")]
    for ef_column in ef_columns:
        df[f"{ef_column} rank"] = df[ef_column].rank(ascending=False, method="min")

    df["HF mean"] = df[hf_columns].mean(axis=1)
    df["LF mean"] = df[lf_columns].mean(axis=1)
    df["EF mean"] = (df["HF mean"] + 1) / (df["LF mean"] + 1)
    df["EF rank"] = df["EF mean"].rank(ascending=False, method="min")

    return df


def get_sp_origins(sp_names: list[str], engine: Engine) -> pd.DataFrame:
    logger.info("Getting SP origins from database")
    names_joined = ",".join(f"'{name}'" for name in sp_names)
    query = text(f"SELECT name, origin FROM sp WHERE name IN ({names_joined})")
    df = pd.read_sql(query, engine)
    return df


def get_read_counts(target_dir: str) -> pd.DataFrame:
    logger.info("Getting read counts")
    filepath = os.path.join(target_dir, "6_read_counts", "read_counts.csv")
    if not os.path.exists(filepath):
        raise ValueError(
            "Read counts file not found: should be saved under target directory as 6_read_counts/read_counts.csv"
        )
    df = pd.read_csv(filepath)

    # Check column names
    if "name" not in df.columns:
        raise ValueError("Read counts file must contain a 'name' column")

    # Check for HF and LF columns using regex pattern
    valid_pattern = re.compile(r"^(HF|LF)\d+$")
    count_columns = [col for col in df.columns if col != "name"]
    invalid_columns = [col for col in count_columns if not valid_pattern.match(col)]

    if invalid_columns:
        raise ValueError(
            f"Invalid column names found: {', '.join(invalid_columns)}. "
            "Column names should be 'name' and then only 'HF' or 'LF' followed by a number (e.g., HF1, LF2)"
        )
    if len(count_columns) < 2:
        raise ValueError("Read counts file must contain at least two columns")

    # Add check for equal number of HF and LF columns
    hf_columns = [col for col in count_columns if col.startswith("HF")]
    lf_columns = [col for col in count_columns if col.startswith("LF")]
    if len(hf_columns) != len(lf_columns):
        raise ValueError(
            f"Unequal number of HF ({len(hf_columns)}) and LF ({len(lf_columns)}) columns"
        )

    if "unmapped" in df["name"].values:
        # Remove unmapped SPs
        logger.info("Removing unmapped counts")
        df = df.loc[df["name"] != "unmapped", :]

    return df


def compute_enrichment(
    target_dir: str,
    db_config_file: str,
    protein: str,
    library: str,
    top_n_similarity: int,
    top_n_signalp: int,
):

    output_dir = os.path.join(target_dir, "7_enrichment")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Connect to db
    db_config = get_db_config(db_config_file)
    engine = get_db_engine(db_config)

    # Read read counts
    df = get_read_counts(target_dir)
    # Get signal peptide origins
    sp_origins = get_sp_origins(df["name"].tolist(), engine)
    # Merge
    df = pd.merge(sp_origins, df, on="name", how="right")

    # Normalise read counts
    df = normalise_read_counts(df)

    # Compute enrichment factors and their ranks
    hf_columns = [
        col for col in df.columns if col.startswith("HF") and col.endswith("_norm")
    ]
    lf_columns = [
        col for col in df.columns if col.startswith("LF") and col.endswith("_norm")
    ]
    df = add_enrichment_factors(df, hf_columns, lf_columns)
    df = df.sort_values(by="EF rank")

    # Add HF percentages
    df = add_hf_percentages(df)

    rank_columns = [
        "EF rank",
        "HF mean percent rank",
    ]

    if len(hf_columns) > 1:
        # Add Cohen's d
        df = add_cohen_d(df, hf_columns, lf_columns)

        # Compute p-values
        df = add_t_test_p_value(df, hf_columns, lf_columns)

        rank_columns += ["p_value rank", "Cohen d rank"]

    common_sps_ranks = get_common_sps(df, rank_columns)

    plot_rank_pairplot(df, rank_columns, os.path.join(output_dir, "rank_pairplot.pdf"))

    # Compute SignalP predictions
    signalp_results = run_signalp(
        df.head(top_n_signalp),
        library,
        protein,
        engine,
        os.path.join(output_dir, "signalp"),
    )
    df = pd.merge(df, signalp_results, on="name", how="left")

    save_df = df.copy()
    float_columns = save_df.select_dtypes(include=["float64", "float32"]).columns
    save_df.loc[:, float_columns] = save_df.loc[:, float_columns].map(
        lambda x: signif(x, 3)
    )

    save_df.to_excel(os.path.join(output_dir, "enrichment_factors.xlsx"), index=False)

    # Create Excel writer object
    with pd.ExcelWriter(os.path.join(output_dir, "common_sps.xlsx")) as writer:
        # Get common SPs
        ef_rank_columns = [col for col in df.columns if re.match(r"^EF\d+ rank$", col)]
        if len(ef_rank_columns) > 1:
            common_sps_ef = get_common_sps(df, ef_rank_columns)
            common_sps_ef.to_excel(
                writer, sheet_name="Replicate Comparisons", index=False
            )
        common_sps_ranks.to_excel(writer, sheet_name="Ranking Comparisons", index=False)

    plot_rank_vs_ef(
        df,
        os.path.join(output_dir, "rank_vs_ef.pdf"),
        "EF mean",
        "EF rank",
        protein,
    )

    if len(hf_columns) > 1:
        try:
            volcano_plot(
                df, os.path.join(output_dir, "volcano.pdf"), "EF mean", protein
            )
        except Exception as e:
            print(e)
            print("Could not plot volcano plot")

    plot_correlation_heatmap(
        df,
        hf_columns + ["HF mean"],
        lf_columns + ["LF mean"],
        os.path.join(output_dir, "correlation_heatmap.pdf"),
    )

    # Compute similarities and draw phylogenetic tree
    compute_similarities_with_tree(
        df.head(top_n_similarity), library, engine, output_dir
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to the input JSON file")
    parser.add_argument(
        "db_config", help="Path to the database configuration JSON file"
    )
    args = parser.parse_args()

    input = read_input(args.input)

    compute_enrichment(
        input["target_directory"],
        args.db_config,
        input["protein"],
        input["library"],
        input["top_n_similarity"],
        input["top_n_signalp"],
    )


if __name__ == "__main__":
    main()
