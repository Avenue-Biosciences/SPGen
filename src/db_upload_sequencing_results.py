import pandas as pd
import argparse
import json
from run_pipeline import get_db_config, get_db_engine
from sqlalchemy import text, Engine
from typing import Iterable


def get_screening_id(screening_name: str, engine: Engine) -> int:
    """
    Get the id of a screening from the database
    If the screening does not exist, create it
    """
    with engine.connect() as conn:
        print(f"Getting screening id for {screening_name}")
        fetch_query = text(f"SELECT id FROM screening WHERE name = '{screening_name}'")
        id = conn.execute(fetch_query).scalar()
        if id is None:
            print(f"Screening {screening_name} not found in database, creating...")
            insert_query = text(
                f"INSERT INTO screening (name) VALUES ('{screening_name}')"
            )
            conn.execute(insert_query)
            conn.commit()
            id = conn.execute(fetch_query).scalar()
        else:
            print(f"Screening {screening_name} found in database")
    return id


def get_protein_id(protein_name: str, engine: Engine) -> int:
    """
    Get the id of a protein from the database
    If the protein does not exist, create it
    """
    with engine.connect() as conn:
        print(f"Getting protein id for {protein_name}")
        fetch_query = text(f"SELECT id FROM protein WHERE name = '{protein_name}'")
        id = conn.execute(fetch_query).scalar()
        if id is None:
            print(f"Protein {protein_name} not found in database, creating...")
            insert_query = text(f"INSERT INTO protein (name) VALUES ('{protein_name}')")
            conn.execute(insert_query)
            conn.commit()
            id = conn.execute(fetch_query).scalar()
        else:
            print(f"Protein {protein_name} found in database")
    return id


def get_sp_ids(sp_names: Iterable[str], engine: Engine) -> pd.DataFrame:
    """
    Get the ids of a list of sps from the database
    """
    print(f"Getting SP ids for {len(sp_names)} SPs")
    fetch_query = text(
        f"SELECT name, id as sp_id FROM sp WHERE name IN ({','.join(f"'{name}'" for name in sp_names)})"
    )
    ids = pd.read_sql(fetch_query, engine)
    print(f"Found {ids.shape[0]} SPs in the database")
    if ids.shape[0] != len(sp_names):
        raise ValueError(
            f"{len(sp_names) - ids.shape[0]} SPs do not exist in the database"
        )

    return ids


def get_existing_results(
    screening_id: int, protein_id: int, engine: Engine
) -> pd.DataFrame:
    screening_result_id_query = text(
        "SELECT sp_id, replicate_number "
        "FROM screening_result "
        "JOIN screening_result_replicate ON screening_result.id = screening_result_replicate.screening_result_id "
        f"WHERE screening_result.screening_id = {screening_id} AND screening_result.protein_id = {protein_id}"
    )
    screening_result_ids = pd.read_sql(screening_result_id_query, engine)
    return screening_result_ids


def insert_results(
    read_counts: pd.DataFrame, screening_id: int, protein_id: int, engine: Engine
):

    # Pivot to longer
    hf_columns = [col for col in read_counts.columns if col.startswith("HF")]
    lf_columns = [col for col in read_counts.columns if col.startswith("LF")]

    read_counts_melted = read_counts.melt(
        id_vars=["sp_id"],
        value_vars=hf_columns + lf_columns,
        var_name="variable",
        value_name="count",
    )
    read_counts_melted["replicate_number"] = (
        read_counts_melted["variable"].str.extract(r"(\d)").astype(int)
    )
    read_counts_melted["variable"] = read_counts_melted["variable"].str.extract(
        r"([A-Z]+)"
    )
    read_counts_melted = read_counts_melted.pivot(
        index=["sp_id", "replicate_number"], columns="variable", values="count"
    ).reset_index()
    read_counts_melted.columns = ["sp_id", "replicate_number", "hf_count", "lf_count"]

    # Check for existing results
    print("Checking for existing results")
    existing_results = get_existing_results(screening_id, protein_id, engine)

    if existing_results.empty:
        print("No existing results found")
        # Create screening results
        print("Creating screening results")
        screening_result = pd.DataFrame(
            {
                "sp_id": read_counts_melted["sp_id"].unique(),
                "screening_id": screening_id,
                "protein_id": protein_id,
            }
        )
        screening_result.to_sql(
            "screening_result", engine, if_exists="append", index=False
        )

        read_counts_insert = read_counts_melted
        print(f"Inserting {read_counts_insert.shape[0]} results into the database")
    else:
        # Remove existing results
        read_counts_insert = pd.merge(
            read_counts_melted,
            existing_results,
            on=["sp_id", "replicate_number"],
            how="outer",
            indicator=True,
        )
        matches_existing = read_counts_insert.loc[
            read_counts_insert["_merge"] == "both", :
        ]
        print(
            f"Found {matches_existing.shape[0]} existing results from replicate numbers {matches_existing['replicate_number'].unique()}"
        )
        read_counts_insert = read_counts_insert.loc[
            read_counts_insert["_merge"] == "left_only", :
        ].drop(columns=["_merge"])

        if read_counts_insert.empty:
            print("No new results to insert")
            return

    # Get IDs of screening results
    screening_result_id_query = text(
        f"SELECT id as screening_result_id, sp_id FROM screening_result WHERE screening_id = {screening_id} AND protein_id = {protein_id}"
    )
    screening_result_ids = pd.read_sql(screening_result_id_query, engine)

    if screening_result_ids.shape[0] != read_counts_insert["sp_id"].nunique():
        raise ValueError(
            f"Number of screening results ({screening_result_ids.shape[0]}) does not match number of SPs to insert results for ({read_counts_insert['sp_id'].nunique()})"
        )

    read_counts_insert = pd.merge(
        read_counts_insert, screening_result_ids, left_on="sp_id", right_on="sp_id"
    )

    read_counts_insert[
        ["screening_result_id", "replicate_number", "hf_count", "lf_count"]
    ].to_sql("screening_result_replicate", engine, if_exists="append", index=False)


def read_input(input_path: str) -> dict:
    with open(input_path, "r") as f:
        input = json.load(f)

    # Check required fields
    required_fields = ["input_file", "screening_name", "protein_name"]
    missing_fields = [field for field in required_fields if field not in input]

    if missing_fields:
        raise ValueError(
            f"Missing required fields in input file: {', '.join(missing_fields)}"
        )

    return input


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to the input JSON file")
    parser.add_argument(
        "db_config", help="Path to the database configuration JSON file"
    )
    parser.add_argument(
        "--empty", action="store_true", help="Empty the database before inserting"
    )
    args = parser.parse_args()

    db_config = get_db_config(args.db_config)
    engine = get_db_engine(db_config)

    input = read_input(args.input)

    screening_id = get_screening_id(input["screening_name"], engine)
    protein_id = get_protein_id(input["protein_name"], engine)

    read_counts_df = pd.read_csv(input["input_file"])
    if "unmapped" in read_counts_df["name"].values:
        # Remove unmapped SPs
        print("Removing unmapped counts")
        read_counts_df = read_counts_df.loc[read_counts_df["name"] != "unmapped", :]

    # Remove unmapped SPs
    sp_ids = get_sp_ids(read_counts_df["name"], engine)
    read_counts_df = pd.merge(read_counts_df, sp_ids, on="name")

    # Development:
    if args.empty:
        print("Emptying the database before inserting")
        with engine.connect() as conn:
            conn.execute(text("DELETE FROM screening_result_replicate"))
            conn.execute(text("DELETE FROM screening_result"))
            conn.commit()

    insert_results(read_counts_df, screening_id, protein_id, engine)


if __name__ == "__main__":
    main()
