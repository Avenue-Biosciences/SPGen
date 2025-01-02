import pandas as pd
import numpy as np
import pytest
from compute_enrichment import add_enrichment_factors, get_common_sps


@pytest.fixture
def sample_df():
    # Create a sample DataFrame with rank columns
    data = {
        "EF1 rank": [1, 2, 3, 50, 150],
        "EF2 rank": [2, 1, 4, 60, 160],
        "EF3 rank": [1, 3, 2, 70, 170],
        "EF rank": [1, 2, 3, 50, 150],
        "other_column": [10, 20, 30, 40, 50],  # Should be ignored
    }
    return pd.DataFrame(data)


def test_get_common_sps_basic(sample_df):
    result = get_common_sps(sample_df, columns=["EF1 rank", "EF2 rank", "EF3 rank"])

    # Check the shape of the result
    assert result.shape == (3, 3)

    # Check column names
    assert list(result.columns) == ["Comparison", "n_common_top_20", "n_common_top_100"]

    # Check the pairs are correctly formed
    expected_pairs = [
        "EF1 rank vs EF2 rank",
        "EF1 rank vs EF3 rank",
        "EF2 rank vs EF3 rank",
    ]
    assert list(result["Comparison"]) == expected_pairs


def test_get_common_sps_counts(sample_df):
    result = get_common_sps(sample_df, columns=["EF1 rank", "EF2 rank", "EF3 rank"])

    # For top 20, all three samples have ranks 1-3, so they should share 3 elements
    assert all(result["n_common_top_20"] == 3)

    # For top 100, all five samples have ranks below 100
    assert all(result["n_common_top_100"] == 4)


def test_add_enrichment_factors():
    # Create sample input data
    data = {
        "HF1": [99, 50],
        "LF1": [9, 10],
        "HF2": [39, 40],
        "LF2": [7, 8],
    }
    df = pd.DataFrame(data)

    # Process the data
    result = add_enrichment_factors(
        df, hf_columns=["HF1", "HF2"], lf_columns=["LF1", "LF2"]
    )

    # Test that all expected columns are present
    expected_columns = [
        "HF1",
        "LF1",
        "HF2",
        "LF2",
        "EF1",
        "EF2",
        "EF1 rank",
        "EF2 rank",
        "HF mean",
        "LF mean",
        "EF mean",
        "EF rank",
    ]
    assert all(col in result.columns for col in expected_columns)

    # Test enrichment factor calculations
    np.testing.assert_almost_equal(result["EF1"][0], 10)
    np.testing.assert_almost_equal(result["EF2"][0], 5)

    # Test mean calculations
    np.testing.assert_almost_equal(result["HF mean"][0], 69)
    np.testing.assert_almost_equal(result["LF mean"][0], 8)
    np.testing.assert_almost_equal(result["EF mean"][0], (69 + 1) / (8 + 1))

    # Test ranking
    assert result["EF1 rank"][0] == 1  # Highest EF should be rank 1
    assert result["EF rank"][0] == 1  # Highest mean EF should be rank 1


def test_add_enrichment_factors_zero_values():
    # Test with zero values in LF columns
    data = {
        "HF1": [10, 5, 0],
        "LF1": [0, 0, 0],
    }
    df = pd.DataFrame(data)
    result = add_enrichment_factors(df, hf_columns=["HF1"], lf_columns=["LF1"])

    # Check that division by zero is handled (we add 1 to denominator)
    np.testing.assert_almost_equal(result["EF1"][0], 11)  # (10 + 1) / (0 + 1)
