import argparse
import logging
import os
from db_utils import get_db_config, get_db_engine
from process_reads import process_reads
from db_upload_sequencing_results import upload_sequencing_results
from compute_enrichment import compute_enrichment

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run the sequencing pipeline.")
    parser.add_argument("directory", help="Directory containing input files")
    parser.add_argument(
        "db_config", help="Path to the database configuration JSON file"
    )
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        filename=os.path.join(args.directory, "run_pipeline.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    db_config = get_db_config(args.db_config)
    engine = get_db_engine(db_config)

    process_reads(args.directory, engine)
    upload_sequencing_results(args.directory, engine)
    compute_enrichment(args.directory, engine)


if __name__ == "__main__":
    main()
