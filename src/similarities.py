import itertools
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage


def compute_similarity(seq1: str, seq2: str) -> int:
    # Compute pairwise similarity
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    similarity = aligner.score(seq1, seq2)
    return similarity


def compute_similarities(sps: pd.DataFrame) -> pd.DataFrame:
    # Create all possible combinations of sequences
    comb_df = pd.DataFrame(itertools.combinations(sps.index, 2), columns=["i", "j"])

    # Compute similarities
    comb_df["Similarity"] = comb_df.apply(
        lambda x: compute_similarity(
            sps.loc[x["i"], "amino_acid_sequence"],
            sps.loc[x["j"], "amino_acid_sequence"],
        ),
        axis=1,
    )

    return comb_df


def plot_similarity_heatmap(similarities_df: pd.DataFrame, output_path: str):
    # Compute distances from similarities
    similarities_df = similarities_df[["SP1", "SP2", "Similarity"]].copy()
    similarities_df["Distance"] = 1 - (
        similarities_df["Similarity"] - similarities_df["Similarity"].min()
    ) / (similarities_df["Similarity"].max() - similarities_df["Similarity"].min())

    # Fill to symmetric matrix
    similarities_df2 = similarities_df.copy()
    similarities_df2["SP1"] = similarities_df["SP2"]
    similarities_df2["SP2"] = similarities_df["SP1"]

    similarities_df = pd.concat(
        [similarities_df, similarities_df2], ignore_index=True
    ).reset_index(drop=True)

    # Convert to wide format
    similarities = similarities_df.drop(columns=["Distance"]).pivot(
        index="SP1", columns="SP2", values="Similarity"
    )
    similarities_wide = similarities.copy()

    distances = similarities_df.drop(columns=["Similarity"]).pivot(
        index="SP1", columns="SP2", values="Distance"
    )
    # Add 1 to the diagonal
    np.fill_diagonal(distances.values, 1)

    # Perform hierarchical clustering using the row distances
    linkage_matrix = linkage(distances, method="average")
    row_order = dendrogram(linkage_matrix, no_plot=True)["leaves"]

    # Reorder rows and columns based on clustering
    similarities = similarities.iloc[row_order, row_order]

    plt.figure(figsize=(10, 6))
    sns.heatmap(similarities, annot=True, cmap="viridis")
    plt.suptitle("Similarity heatmap")
    plt.title("Higher values indicate more similar sequences")
    plt.savefig(output_path, bbox_inches="tight")

    return similarities_wide
