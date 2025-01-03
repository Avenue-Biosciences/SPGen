import itertools
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Iterable
from sqlalchemy import Engine
import logging
from scipy.cluster.hierarchy import dendrogram, linkage
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from db_utils import get_library_sequences
from utils import run_command

logger = logging.getLogger(__name__)


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
    logger.info("Computing pairwise similarities between SP amino acid sequences")
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


def draw_phylogenetic_tree(
    fasta_file,
    output_dir,
    matrix="blosum62",
):
    """
    Creates a phylogenetic tree from a FASTA file using the following steps:
    1. Multiple sequence alignment using MUSCLE
    2. Distance matrix calculation using amino acid substitution matrix
    3. Neighbor-joining tree construction
    4. Tree visualization using graphviz

    Args:
        fasta_file (str): Path to input FASTA file
        fig: matplotlib figure object
        ax: matplotlib axis object
        matrix (str): Substitution matrix to use ('blosum62', 'pam250', etc.)
    """
    logger.info("Performing multiple sequence alignment")
    output = os.path.join(output_dir, "alignment_" + os.path.basename(fasta_file))
    # Run MUSCLE alignment using subprocess
    muscle_cmd = [
        "muscle",
        "-align",
        fasta_file,
        "-output",
        output,
    ]
    try:
        run_command(muscle_cmd)
    except Exception as e:
        logger.error(f"Error running MUSCLE, could not draw phylogenetic tree")
        return

    # Read the alignment
    logger.info("Reading alignment and drawing phylogenetic tree")
    alignment = AlignIO.read(output, "fasta")

    # Calculate distance matrix using specified substitution matrix
    calculator = DistanceCalculator(model=matrix)
    dm = calculator.get_distance(alignment)

    tree_methods = ["nj", "upgma"]

    # Create a single PDF with all plots

    fig, ax = plt.subplots(1, 2, figsize=(15, 10))
    for i, tree_method in enumerate(tree_methods):

        # Construct tree using neighbor-joining
        constructor = DistanceTreeConstructor(calculator, tree_method)
        tree = constructor.build_tree(alignment)
        tree.ladderize()

        # Save the tree to a file
        Phylo.draw(
            tree,
            axes=ax[i],
            do_show=False,
            label_func=lambda x: x.name if "Inner" not in x.name else "",
        )
        ax[i].set_title(
            os.path.basename(fasta_file).replace(".fasta", "")
            + " Distance Matrix: "
            + matrix
            + " Tree Method: "
            + tree_method
        )

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Top_SPS_phylogenetic_trees.pdf"))
    plt.close()


def write_fasta(ids: Iterable[str], sequences: Iterable[str], output_path: str):
    records = [
        SeqRecord(Seq(seq), id=id, description="") for id, seq in zip(ids, sequences)
    ]
    SeqIO.write(records, output_path, "fasta")


def compute_similarities_with_tree(
    df: pd.DataFrame, library: str, engine: Engine, output_dir: str
):
    # Get library sequences
    library_sequences = get_library_sequences(df["name"].tolist(), library, engine)
    df = pd.merge(df, library_sequences, on="name", how="left")

    # Compute similarities
    top_sps_similarities = compute_similarities(df)
    top_sps_similarities["SP1"] = df.loc[top_sps_similarities["i"], "name"].reset_index(
        drop=True
    )
    top_sps_similarities["SP2"] = df.loc[top_sps_similarities["j"], "name"].reset_index(
        drop=True
    )
    similarities_wide = plot_similarity_heatmap(
        top_sps_similarities,
        os.path.join(output_dir, f"SP similarities heatmap.pdf"),
    )
    similarities_wide.to_excel(
        os.path.join(output_dir, f"SP similarities wide.xlsx"), index=True
    )

    # Write to FASTA and draw phylogenetic tree
    fasta_path = os.path.join(output_dir, f"SPs.fasta")
    write_fasta(df["name"], df["amino_acid_sequence"], fasta_path)
    draw_phylogenetic_tree(fasta_path, output_dir)


def get_mature_sequence(protein: str):
    logger.info(f"Getting mature sequence of {protein} from database")
    return "LLLLL"


def run_signalp(
    df: pd.DataFrame, library: str, protein: str, engine: Engine, output_dir: str
):

    os.makedirs(output_dir, exist_ok=True)
    # Get library sequences
    library_sequences = get_library_sequences(df["name"].tolist(), library, engine)
    df = pd.merge(df, library_sequences, on="name", how="left")

    # Add mature sequence to the amino acid sequence
    mature_seq = get_mature_sequence(protein)
    df["amino_acid_sequence_with_mature"] = df["amino_acid_sequence"] + mature_seq

    # Write to FASTA
    fasta_path = os.path.join(output_dir, f"SPs_with_mature.fasta")

    logger.info("Running SignalP")
    try:
        write_fasta(df["name"], df["amino_acid_sequence_with_mature"], fasta_path)

        signalp_cmd = [
            "signalp6",
            "--fastafile",
            fasta_path,
            "--organism",
            "eukarya",
            "--output_dir",
            output_dir,
            "--format",
            "none",
            "--mode",
            "slow-sequential",
        ]
        run_command(signalp_cmd)
    except Exception as e:
        logger.error(f"Error running SignalP: {e}")
        raise e
    finally:
        if os.path.exists(fasta_path):
            os.remove(fasta_path)

    logger.info("Processing SignalP results")
    results = pd.read_csv(
        os.path.join(output_dir, "prediction_results.txt"), sep="\t", skiprows=1
    )
    results.columns = [
        "name",
        "SignalP_prediction",
        "SignalP_P_Other",
        "SignalP_P_SP",
        "CS_position",
    ]
    results[["CS_start", "CS_end", "CS_prob"]] = results["CS_position"].str.extract(
        r"CS pos\: (\d+)-(\d+)\. Pr: (\d+\.\d+)"
    )
    results["CS_prob"] = results["CS_prob"].astype(float)
    results["CS_start"] = results["CS_start"].astype(float)
    results["CS_end"] = results["CS_end"].astype(float)

    df["SP_length"] = df["amino_acid_sequence"].str.len()
    results = pd.merge(df[["name", "SP_length"]], results, on="name", how="left")

    return results[
        [
            "name",
            "SignalP_prediction",
            "SignalP_P_Other",
            "SignalP_P_SP",
            "SP_length",
            "CS_start",
            "CS_end",
            "CS_prob",
        ]
    ]
