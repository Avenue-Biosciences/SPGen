import pytest
import pandas as pd
from sqlalchemy import create_engine, text
from db_upload_sequencing_results import (
    get_screening_id,
    get_protein_id,
    get_sp_ids,
    insert_results,
    get_existing_results,
)


@pytest.fixture
def test_engine():
    """Create a test database connection"""
    engine = create_engine("sqlite:///:memory:")

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
            CREATE TABLE sp (
                id INTEGER PRIMARY KEY,
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

    # Insert some test data
    db_insert = pd.read_csv("tests/complete_pipeline_files/db_insert.csv")[
        ["id", "name"]
    ]
    db_insert.to_sql("sp", engine, if_exists="replace", index=False)

    return engine


def test_get_screening_id(test_engine):
    """Test getting/creating screening ID"""
    # Test creating new screening
    screening_id = get_screening_id("Test Screening", test_engine)
    assert screening_id == 1

    # Test getting existing screening
    screening_id_2 = get_screening_id("Test Screening", test_engine)
    assert screening_id_2 == screening_id


def test_get_protein_id(test_engine):
    """Test getting/creating protein ID"""
    # Test creating new protein
    protein_id = get_protein_id("Test Protein", test_engine)
    assert protein_id == 1

    # Test getting existing protein
    protein_id_2 = get_protein_id("Test Protein", test_engine)
    assert protein_id_2 == protein_id


def test_get_sp_ids(test_engine):
    """Test getting SP IDs"""
    sp_names = ["Q9BZM4", "P36888"]
    sp_ids = get_sp_ids(sp_names, test_engine)
    expected_ids = pd.DataFrame({"name": sp_names, "sp_id": [1, 3]})
    pd.testing.assert_frame_equal(sp_ids, expected_ids)


def verify_results(
    expected_read_counts, n_replicates, screening_id, protein_id, engine
):

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
        f"""
        SELECT 
            sp.name,
            MAX(CASE WHEN replicate_number = 1 THEN hf_count END) as HF1,
            MAX(CASE WHEN replicate_number = 1 THEN lf_count END) as LF1,
            MAX(CASE WHEN replicate_number = 2 THEN hf_count END) as HF2,
            MAX(CASE WHEN replicate_number = 2 THEN lf_count END) as LF2
        FROM screening_result_replicate
        JOIN screening_result ON screening_result_replicate.screening_result_id = screening_result.id
        JOIN sp ON screening_result.sp_id = sp.id
        WHERE screening_result.screening_id = {screening_id}
        AND screening_result.protein_id = {protein_id}
        GROUP BY sp.name
        """
    )

    inserted = pd.read_sql(query, engine).sort_values(by="name").reset_index(drop=True)
    pd.testing.assert_frame_equal(
        inserted, expected_read_counts.sort_values(by="name").reset_index(drop=True)
    )


def verify_results(
    expected_read_counts, n_replicates, screening_id, protein_id, engine
):

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
        f"""
        SELECT 
            sp.name,
            MAX(CASE WHEN replicate_number = 1 THEN hf_count END) as HF1,
            MAX(CASE WHEN replicate_number = 1 THEN lf_count END) as LF1,
            MAX(CASE WHEN replicate_number = 2 THEN hf_count END) as HF2,
            MAX(CASE WHEN replicate_number = 2 THEN lf_count END) as LF2
        FROM screening_result_replicate
        JOIN screening_result ON screening_result_replicate.screening_result_id = screening_result.id
        JOIN sp ON screening_result.sp_id = sp.id
        WHERE screening_result.screening_id = {screening_id}
        AND screening_result.protein_id = {protein_id}
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


def test_insert_results(test_engine):
    """Test inserting sequencing results"""
    # Create test data
    screening_id = get_screening_id("Test Screening", test_engine)
    protein_id = get_protein_id("Test Protein", test_engine)

    read_counts = pd.read_csv(
        "tests/complete_pipeline_files/expected_read_counts.csv"
    ).loc[:3, :]

    # Get SP IDs
    sp_ids = get_sp_ids(read_counts["name"], test_engine)
    read_counts_merged = pd.merge(read_counts, sp_ids, on="name")

    # Insert results
    insert_results(read_counts_merged, screening_id, protein_id, test_engine)

    # Verify results
    verify_results(read_counts, 2, screening_id, protein_id, test_engine)

    # Test inserting duplicate results (should not raise error)
    insert_results(read_counts_merged, screening_id, protein_id, test_engine)

    # Verify no duplicates were added
    verify_results(read_counts, 2, screening_id, protein_id, test_engine)


def test_insert_results_with_new_replicates(test_engine):
    """Test inserting additional replicates"""
    screening_id = get_screening_id("Test Screening", test_engine)
    protein_id = get_protein_id("Test Protein", test_engine)

    # First insertion
    read_counts_1 = pd.read_csv(
        "tests/complete_pipeline_files/expected_read_counts.csv"
    ).loc[:3, ["name", "HF1", "LF1"]]
    sp_ids = get_sp_ids(read_counts_1["name"], test_engine)
    read_counts_1_merged = pd.merge(read_counts_1, sp_ids, on="name")
    insert_results(read_counts_1_merged, screening_id, protein_id, test_engine)

    # Verify results
    verify_results(read_counts_1, 1, screening_id, protein_id, test_engine)

    # Second insertion with new replicates
    read_counts_2 = pd.read_csv(
        "tests/complete_pipeline_files/expected_read_counts.csv"
    ).loc[:3, ["name", "HF2", "LF2"]]
    read_counts_2_merged = pd.merge(read_counts_2, sp_ids, on="name")
    insert_results(read_counts_2_merged, screening_id, protein_id, test_engine)
    # Verify results
    verify_results(
        pd.merge(read_counts_1, read_counts_2, on="name"),
        2,
        screening_id,
        protein_id,
        test_engine,
    )


def test_insert_results_with_new_sps(test_engine):
    """Test inserting results with new SPs"""
    screening_id = get_screening_id("Test Screening", test_engine)
    protein_id = get_protein_id("Test Protein", test_engine)

    # First insertion
    read_counts = pd.read_csv(
        "tests/complete_pipeline_files/expected_read_counts.csv"
    ).loc[:3, :]
    sp_ids = get_sp_ids(read_counts["name"], test_engine)
    read_counts_merged = pd.merge(read_counts, sp_ids, on="name")
    insert_results(read_counts_merged, screening_id, protein_id, test_engine)

    # Verify results
    verify_results(read_counts, 2, screening_id, protein_id, test_engine)

    # Second insertion with new screening id
    screening_id_2 = get_screening_id("Test Screening 2", test_engine)
    insert_results(read_counts_merged, screening_id_2, protein_id, test_engine)
    # Verify results
    with test_engine.connect() as conn:
        n_rows_screening_result = conn.execute(
            text("SELECT COUNT(*) FROM screening_result")
        ).scalar()
        assert n_rows_screening_result == 8

        n_rows_screening_result_replicate = conn.execute(
            text("SELECT COUNT(*) FROM screening_result_replicate")
        ).scalar()
        assert n_rows_screening_result_replicate == 16

    query = text(
        f"""
        SELECT 
            sp.name,
            screening_result.screening_id,
            MAX(CASE WHEN replicate_number = 1 THEN hf_count END) as HF1,
            MAX(CASE WHEN replicate_number = 1 THEN lf_count END) as LF1,
            MAX(CASE WHEN replicate_number = 2 THEN hf_count END) as HF2,
            MAX(CASE WHEN replicate_number = 2 THEN lf_count END) as LF2
        FROM screening_result_replicate
        JOIN screening_result ON screening_result_replicate.screening_result_id = screening_result.id
        JOIN sp ON screening_result.sp_id = sp.id
        WHERE screening_result.screening_id IN ({screening_id}, {screening_id_2})
        AND screening_result.protein_id = {protein_id}
        GROUP BY sp.name, screening_result.screening_id
        """
    )

    inserted = (
        pd.read_sql(query, test_engine).sort_values(by="name").reset_index(drop=True)
    )
    expected = pd.concat([read_counts, read_counts])
    expected["screening_id"] = [screening_id] * 4 + [screening_id_2] * 4
    pd.testing.assert_frame_equal(
        inserted,
        expected[inserted.columns]
        .sort_values(by=["name", "screening_id"])
        .reset_index(drop=True),
    )
