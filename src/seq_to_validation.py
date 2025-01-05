import argparse
import json
import logging
import os
import pandas as pd
from db_utils import get_db_config, get_db_engine, get_library_sequences
from top_sp_stats import write_fasta

logger = logging.getLogger(__name__)


def read_input(input_file):
    with open(input_file, "r") as f:
        input = json.load(f)

    required_fields = ["SP library"]
    for field in required_fields:
        if field not in input:
            raise ValueError(
                f"Error: {field} was not found. Please check the input file."
            )
    return input


def main():
    parser = argparse.ArgumentParser(
        description="""
        Write FASTA files of chosen SPs for validation.
        Expects a validation_sps.csv file with SP names in the 'name' column in the 8_validation subdirectory.
        """
    )
    parser.add_argument("directory", help="Directory containing input files")
    parser.add_argument(
        "db_config", help="Path to the database configuration JSON file"
    )
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    input = read_input(os.path.join(args.directory, "input.json"))
    # Read chosen SPs
    chosen_sps = pd.read_csv(
        os.path.join(args.directory, "8_validation", "validation_sps.csv")
    )

    db_config = get_db_config(args.db_config)
    engine = get_db_engine(db_config)

    # Get DNA sequences for chosen SPs
    seqs = get_library_sequences(chosen_sps["name"], input["SP library"], engine)
    chosen_sps = pd.merge(chosen_sps, seqs, on="name")

    # Write FASTA files
    logger.info(f"Writing FASTA files to {args.directory}/8_validation")
    write_fasta(
        chosen_sps["name"],
        chosen_sps["amino_acid_sequence"],
        os.path.join(args.directory, "8_validation", "validation_sps_amino_acid.fasta"),
    )
    write_fasta(
        chosen_sps["name"],
        chosen_sps["dna_sequence_sp_only"],
        os.path.join(args.directory, "8_validation", "validation_sps_dna.fasta"),
    )


if __name__ == "__main__":
    main()
