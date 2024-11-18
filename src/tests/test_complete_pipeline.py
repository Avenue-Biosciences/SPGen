import os
import shutil
import tempfile
import pytest
from run_pipeline import main


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


def test_complete_pipeline(temp_test_dir, monkeypatch):
    """Test the complete pipeline execution."""
    # Modify sys.argv to pass the temp directory as argument
    monkeypatch.setattr("sys.argv", ["run_pipeline.py", temp_test_dir])

    # Run the pipeline
    main()

    # Assert expected output directories exist
    expected_dirs = [
        "1_fastqc_output",
        "2_adapter_removed",
        "3_forward_reverse_reads_merged",
        "4_filler_removed",
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
    output_file = os.path.join(temp_test_dir, "6_read_counts", "read_counts.csv")
    expected_file = "expected_read_counts.csv"

    with open(output_file, "r") as f1, open(expected_file, "r") as f2:
        assert f1.read() == f2.read(), "Output does not match expected results"
