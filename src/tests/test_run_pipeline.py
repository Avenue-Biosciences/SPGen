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

    # Create test tables
    with engine.connect() as conn:
        conn.execute(
            text(
                """
            CREATE TABLE screening (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL UNIQUE
            )
        """
            )
        )

        conn.execute(
            text(
                """
            CREATE TABLE protein (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL UNIQUE
            )
        """
            )
        )

        conn.execute(
            text(
                """
            CREATE TABLE screening_result (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                screening_id INTEGER NOT NULL,
                protein_id INTEGER NOT NULL,
                sp_id INTEGER NOT NULL,
                FOREIGN KEY (screening_id) REFERENCES screening (id),
                FOREIGN KEY (protein_id) REFERENCES protein (id),
                FOREIGN KEY (sp_id) REFERENCES sp (id)
            )
        """
            )
        )

        conn.execute(
            text(
                """
            CREATE TABLE screening_result_replicate (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                screening_result_id INTEGER NOT NULL,
                replicate_number INTEGER NOT NULL,
                hf_count INTEGER NOT NULL,
                lf_count INTEGER NOT NULL,
                FOREIGN KEY (screening_result_id) REFERENCES screening_result (id)
            )
        """
            )
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


def verify_db_contents(expected_read_counts, n_replicates, engine):

    with engine.connect() as conn:
        n_rows_screening_result = conn.execute(
            text("SELECT COUNT(*) FROM screening_result")
        ).scalar()
        assert n_rows_screening_result == expected_read_counts.shape[0]

        n_rows_screening_result_replicate = conn.execute(
            text("SELECT COUNT(*) FROM screening_result_replicate")
        ).scalar()
        assert (
            n_rows_screening_result_replicate
            == expected_read_counts.shape[0] * n_replicates
        )

    query = text(
        """
        SELECT 
            sp.name,
            MAX(CASE WHEN replicate_number = 1 THEN hf_count END) as HF1,
            MAX(CASE WHEN replicate_number = 1 THEN lf_count END) as LF1,
            MAX(CASE WHEN replicate_number = 2 THEN hf_count END) as HF2,
            MAX(CASE WHEN replicate_number = 2 THEN lf_count END) as LF2
        FROM screening_result_replicate
        JOIN screening_result ON screening_result_replicate.screening_result_id = screening_result.id
        JOIN sp ON screening_result.sp_id = sp.id
        WHERE screening_result.screening_id = 1
        AND screening_result.protein_id = 1
        GROUP BY sp.name
        """
    )

    inserted = (
        pd.read_sql(query, engine)[expected_read_counts.columns]
        .sort_values(by="name")
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(
        inserted,
        expected_read_counts.sort_values(by="name").reset_index(drop=True),
    )


def test_complete_pipeline(temp_test_dir, test_engine, monkeypatch):
    """Test the complete pipeline execution."""
    # Modify sys.argv to pass the temp directory as argument
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_pipeline.py",
            temp_test_dir,
            os.path.join(temp_test_dir, "db_config.json"),
        ],
    )
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
        "7_enrichment",
    ]

    for dir_name in expected_dirs:
        assert os.path.exists(os.path.join(temp_test_dir, dir_name))

    # Assert read counts file exists
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

    # Verify db contents
    verify_db_contents(expected.loc[expected["name"] != "unmapped", :], 2, test_engine)

    # Assert enrichment files exits
    expected_enrichment_files = set(
        [
            "common_sps.xlsx",
            "enrichment_factors.xlsx",
            "correlation_heatmap.pdf",
            "rank_pairplot.pdf",
            "rank_vs_ef.pdf",
            "SP similarities heatmap.pdf",
            "SP similarities wide.xlsx",
            "SPs.fasta",
            "signalp",
        ]
    )

    found_files = set(os.listdir(os.path.join(temp_test_dir, "7_enrichment")))

    assert found_files == expected_enrichment_files
