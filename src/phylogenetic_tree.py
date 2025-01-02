from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os
import subprocess
import matplotlib.pyplot as plt


def draw_phylogenetic_tree(
    fasta_file,
    output_dir,
    matrix="blosum62",
    muscle_path="/Users/anton/Projects/screening_platform/src/muscle-osx-arm64.v5.3",
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

    output = os.path.join(output_dir, "alignment_" + os.path.basename(fasta_file))
    # Run MUSCLE alignment using subprocess
    muscle_cmd = [
        muscle_path,
        "-align",
        fasta_file,
        "-output",
        output,
    ]
    subprocess.run(muscle_cmd, check=True, capture_output=True)

    # Read the alignment
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
