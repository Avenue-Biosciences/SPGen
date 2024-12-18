import os
import shutil
import tempfile
import pytest
import pandas as pd
from run_pipeline import main
from sqlalchemy import create_engine, text


@pytest.fixture
def test_engine():
    """Initialize a test database using SQLite in-memory database"""
    engine = create_engine("sqlite:///:memory:")

    # Insert test data
    insert_data = pd.read_csv("tests/complete_pipeline_files/db_insert.csv")

    # Insert sp data
    sp_data = insert_data[
        ["id", "name", "uniprot_entry", "origin", "amino_acid_sequence"]
    ]
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
        # Copy all files from complete_pipeline_files directory
        source_dir = "tests/complete_pipeline_files"
        for file in os.listdir(source_dir):
            shutil.copy(os.path.join(source_dir, file), temp_dir)

        # Store original working directory
        original_dir = os.getcwd()

        yield temp_dir

        # Change back to original directory after test
        os.chdir(original_dir)


def test_complete_pipeline(temp_test_dir, test_engine, monkeypatch):
    """Test the complete pipeline execution."""
    # Modify sys.argv to pass the temp directory as argument
    monkeypatch.setattr("sys.argv", ["run_pipeline.py", temp_test_dir])
    monkeypatch.setattr("run_pipeline.get_db_engine", lambda config: test_engine)

    expected = (
        pd.read_csv("tests/complete_pipeline_files/expected_read_counts.csv")
        .sort_values("name")
        .reset_index(drop=True)
    )
    # Run the pipeline
    main()

    # Assert expected output directories exist
    expected_dirs = [
        "1_fastqc",
        "2_adapter_removed",
        "3_forward_reverse_reads_merged",
        "4_filler_restriction_removed",
        "5_alignments",
        "6_read_counts",
    ]

    for dir_name in expected_dirs:
        assert os.path.exists(os.path.join(temp_test_dir, dir_name))

    # Assert final output file exists
    assert os.path.exists(
        os.path.join(temp_test_dir, "6_read_counts", "read_counts.csv")
    )

    # Compare output with expected results
    output = (
        pd.read_csv(os.path.join(temp_test_dir, "6_read_counts", "read_counts.csv"))
        .sort_values("name")
        .reset_index(drop=True)
    )

    pd.testing.assert_frame_equal(output, expected)
