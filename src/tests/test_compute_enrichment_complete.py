import pytest
import pandas as pd
import os
from unittest.mock import patch
import sys
import tempfile
import shutil
import json
from sqlalchemy import create_engine
from compute_enrichment import main


def signif(x, digits=2):
    return float(f"%.{digits}g" % x)


@pytest.fixture
def test_engine():
    """Initialize a test database using SQLite in-memory database"""
    engine = create_engine("sqlite:///:memory:")

    # Insert test data
    insert_data = pd.read_csv("tests/complete_enrichment_files/db_insert_sps.csv")

    # Insert sp data
    sp_data = insert_data[["id", "name", "origin", "amino_acid_sequence"]]
    sp_data.to_sql("sp", engine, if_exists="replace", index=False)

    # Create a new library
    pd.DataFrame({"id": [1], "name": ["integration_test"], "type": ["test"]}).to_sql(
        "sp_library", engine, if_exists="replace", index=False
    )

    # Link sps to the new library
    sp_to_sp_library_data = insert_data[
        ["dna_sequence_sp_only", "dna_sequence_full"]
    ].copy()
    sp_to_sp_library_data["sp_id"] = insert_data["id"]
    sp_to_sp_library_data["sp_library_id"] = 1
    sp_to_sp_library_data.to_sql(
        "sp_to_sp_library", engine, if_exists="replace", index=False
    )

    return engine


@pytest.fixture
def temp_test_dir():
    """Create a temporary directory for testing and clean it up afterward."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy all files and directories from complete_run_files directory
        source_dir = "tests/complete_enrichment_files"
        for item in os.listdir(source_dir):
            source_path = os.path.join(source_dir, item)
            dest_path = os.path.join(temp_dir, item)
            if os.path.isdir(source_path):
                shutil.copytree(source_path, dest_path)
            else:
                shutil.copy(source_path, dest_path)

        # Update input.json with temp_dir as output
        with open(os.path.join(temp_dir, "input.json"), "r+") as f:
            data = json.loads(f.read())
            data["target_directory"] = temp_dir
            f.seek(0)
            json.dump(data, f, indent=4)
            f.truncate()

        # Store original working directory
        original_dir = os.getcwd()

        yield temp_dir

        # Change back to original directory after test
        os.chdir(original_dir)


def test_complete_run(temp_test_dir, test_engine, monkeypatch):
    """Test the complete run execution."""
    # Modify sys.argv to pass the temp directory as argument
    monkeypatch.setattr(
        "sys.argv",
        [
            "compute_enrichment.py",
            temp_test_dir,
            os.path.join(temp_test_dir, "db_config.json"),
        ],
    )
    monkeypatch.setattr("compute_enrichment.get_db_engine", lambda config: test_engine)

    main()

    # Assert output files exist
    output_dir = os.path.join(temp_test_dir, "7_enrichment")
    assert os.path.exists(os.path.join(output_dir, "enrichment_factors.xlsx"))
    assert os.path.exists(os.path.join(output_dir, "common_sps.xlsx"))
    assert os.path.exists(os.path.join(output_dir, "rank_vs_ef.pdf"))
    assert os.path.exists(os.path.join(output_dir, "volcano.pdf"))
    assert os.path.exists(os.path.join(output_dir, "correlation_heatmap.pdf"))

    # Read and verify output data
    results = pd.read_excel(os.path.join(output_dir, "enrichment_factors.xlsx"))
    # TODO: testing p_value
    results = results.drop(
        columns=[
            "p_value",
            "p_value rank",
            "Cohen d",
            "Cohen d rank",
            "HF mean percent",
            "HF mean percent rank",
            "SignalP_prediction",
            "SignalP_P_Other",
            "SignalP_P_SP",
            "SP_length",
            "CS_start",
            "CS_end",
            "CS_prob",
        ]
    )
    expected_results = pd.read_excel(
        "tests/complete_enrichment_files/expected_output.xlsx"
    )
    float_columns = expected_results.select_dtypes(
        include=["float64", "float32"]
    ).columns
    expected_results.loc[:, float_columns] = expected_results.loc[:, float_columns].map(
        lambda x: signif(x, 3)
    )

    pd.testing.assert_frame_equal(results, expected_results)
