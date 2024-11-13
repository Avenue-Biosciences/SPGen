import os
import sys
import json
import subprocess
import argparse
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def check_and_create_directories(output_dirs):
    """Create output directories after checking they don't exist."""
    for dir_name in output_dirs:
        if os.path.exists(dir_name):
            raise ValueError(
                f"Error: '{dir_name}' directory already exists. Please clean the directory before running the pipeline."
            )
        os.makedirs(dir_name)


def read_config(config_file):
    """Read and validate the configuration file."""
    with open(config_file) as f:
        config = json.load(f)

    # Validate values
    required_fields = ["Sample ID", "R1 filename", "R2 filename"]
    for sample_type in ["High fluorescence", "Low fluorescence"]:
        for field in required_fields:
            value = config[sample_type][field]
            if not value:
                raise ValueError(
                    f"Error: {sample_type} {field} was not found. Please check the config file."
                )

    return config


def run_command(command):
    """
    Wrapper function to run subprocess commands and log their output.

    Args:
        command (list): Command to run as a list of strings

    Returns:
        subprocess.CompletedProcess: Result of the subprocess run
    """
    logger.info(f"Running command: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        logger.info(result.stdout)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with error: {e}")
        logger.error(f"Error output:\n{e.stderr}")
        raise


def run_fastqc(r1_file, r2_file, output_dir):
    """Run FastQC on the sequencing files."""
    logger.info("Running FastQC")
    return run_command(["fastqc", r1_file, r2_file, "-o", output_dir])


def cut_adapters(sample_id, r1_file, r2_file, output_dir):
    """Remove adapters from the sequencing files."""
    logger.info("Removing adapters")
    run_command(
        [
            "cutadapt",
            "-q 25",
            "-m 50",
            "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGTGTAGATCTCGGTGGTCGCCGTATCATT",
            r1_file,
            r2_file,
            "-o",
            f"{output_dir}/{sample_id}_adapter_removed_R1.fastq",
            "-p",
            f"{output_dir}/{sample_id}_adapter_removed_R2.fastq",
        ]
    )
    return (
        f"{output_dir}/{sample_id}_adapter_removed_R1.fastq",
        f"{output_dir}/{sample_id}_adapter_removed_R2.fastq",
    )


def merge_reads(sample_id, r1_file, r2_file, output_dir):
    """Merge forward and reverse reads."""
    logger.info("Merging forward and reverse reads")

    os.chdir(output_dir)
    run_command(
        [
            "flash",
            "-M",
            "350",
            r1_file,
            r2_file,
        ]
    )
    # Rename flash output files
    flash_files = {
        "out.extendedFrags.fastq": f"{sample_id}_merged_reads.fastq",
        "out.hist": f"{sample_id}.hist",
        "out.histogram": f"{sample_id}.histogram",
        "out.notCombined_1.fastq": f"{sample_id}_notCombined_R1.fastq",
        "out.notCombined_2.fastq": f"{sample_id}_notCombined_R2.fastq",
    }
    for old, new in flash_files.items():
        os.rename(old, new)
    os.chdir("..")
    return f"{output_dir}/{sample_id}_merged_reads.fastq"


def remove_filler_sequence(sample_id, input_file, output_dir, kozak_seq):
    """Remove filler sequence from the merged reads."""
    logger.info("Removing filler sequence")

    run_command(
        [
            "cutadapt",
            "-e",
            "0",
            "-g",
            kozak_seq,
            input_file,
            "-o",
            f"{output_dir}/{sample_id}_reads_no_trimmer.fastq",
        ],
    )
    return f"{output_dir}/{sample_id}_reads_no_trimmer.fastq"


def remove_restriction_site(sample_id, input_file, output_dir):
    """Remove MLuI restriction site from the reads."""
    logger.info("Removing MLuI restriction site and mature protein chain")
    run_command(
        [
            "cutadapt",
            "-e",
            "0",
            "-a",
            "ACGCGT",
            input_file,
            "-o",
            f"{output_dir}/{sample_id}_reads_final.fastq",
        ],
    )
    return f"{output_dir}/{sample_id}_reads_final.fastq"


def align_reads(sample_id, input_file, output_dir, splib_idx):
    """Align the reads to the SP library."""
    logger.info("Aligning reads to SP library")

    os.chdir(output_dir)
    run_command(
        [
            "bowtie2",
            "-p",
            "1",
            "--score-min",
            "C,0,-1",
            "-x",
            splib_idx,
            input_file,
            "-S",
            f"{sample_id}_reads_bowtie_exact.sam",
        ],
    )

    # Sort the reads and convert to binary alignment format
    run_command(
        [
            "samtools",
            "view",
            "-bS",
            f"{sample_id}_reads_bowtie_exact.sam",
            "-o",
            f"{sample_id}_reads_bowtie_exact.bam",
        ],
    )
    run_command(
        [
            "samtools",
            "sort",
            f"{sample_id}_reads_bowtie_exact.bam",
            "-o",
            f"{sample_id}_reads_bowtie_exact.sorted.bam",
        ],
    )

    # Index the sorted BAM file
    run_command(
        ["samtools", "index", f"{sample_id}_reads_bowtie_exact.sorted.bam"],
    )

    # Get the number of reads that aligned to each SP in the reference index
    subprocess.run(
        [
            "samtools",
            "idxstats",
            f"{sample_id}_reads_bowtie_exact.sorted.bam",
        ],
        stdout=open(f"{output_dir}/{sample_id}_matches_bowtie_exact.txt", "w"),
    )

    os.chdir("..")
    return f"{output_dir}/{sample_id}_matches_bowtie_exact.txt"


def process_sample(sample_id, r1_file, r2_file, output_dirs, sp_library):
    """Process a single sample through the pipeline."""
    if sp_library in ["R2B", "integration_test"]:
        kozak_seq = "GCTAGCCCACC"
    elif sp_library in ["v2"]:
        kozak_seq = "GCCGCCACC"
    else:
        raise ValueError(f"Unknown SP library: {sp_library}")

    # Run FastQC
    run_fastqc(r1_file, r2_file, output_dirs["fastqcout"])

    # Remove adapters
    r1_file, r2_file = cut_adapters(
        sample_id, r1_file, r2_file, output_dirs["adaptout"]
    )

    # Merge forward and reverse reads
    merged_reads_file = merge_reads(
        sample_id, r1_file, r2_file, output_dirs["flashout"]
    )

    # Remove filler sequence
    filler_removed_file = remove_filler_sequence(
        sample_id, merged_reads_file, output_dirs["fillerout"], kozak_seq
    )

    # Remove MLuI restriction site
    restriction_removed_file = remove_restriction_site(
        sample_id, filler_removed_file, output_dirs["fillerout"]
    )

    # Alignment steps
    sp_lib_idx = f"/sp_library_refs/{sp_library}"
    alignment_counts_file = align_reads(
        sample_id, restriction_removed_file, output_dirs["alignout"], sp_lib_idx
    )
    return alignment_counts_file


def read_alignment_counts(filename):
    """Read the alignment counts for a single sample."""
    df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["SP", "length", "count", "unmapped"],
    )
    df.loc[df["SP"] == "*", "count"] = df.loc[df["SP"] == "*", "unmapped"]
    df.loc[df["SP"] == "*", "SP"] = "unmapped"
    return df.drop(columns=["length", "unmapped"])


def process_read_counts(hf_file, lf_file):
    """Process the read counts for the high and low fluorescence samples."""
    logger.info("Processing read counts")

    hf_read_counts = read_alignment_counts(hf_file)
    lf_read_counts = read_alignment_counts(lf_file)

    # Merge the two dataframes
    merged_read_counts = pd.merge(
        hf_read_counts, lf_read_counts, on="SP", suffixes=["_hf", "_lf"]
    )
    merged_read_counts.columns = ["SP", "HF_count", "LF_count"]
    return merged_read_counts


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run the sequencing pipeline.")
    parser.add_argument("directory", help="Directory containing input files")
    args = parser.parse_args()

    # Change to specified directory
    abs_directory = os.path.abspath(args.directory)
    try:
        os.chdir(abs_directory)
    except OSError as e:
        print(f"Error: Unable to change to directory '{abs_directory}'")
        sys.exit(1)

    # Set up logging
    logging.basicConfig(filename="run_pipeline.log", level=logging.INFO)
    logger.info("Starting pipeline")

    # Define output directories
    output_dirs = {
        "fastqcout": "1_fastqc_output",
        "adaptout": "2_adapter_removed",
        "flashout": "3_forward_reverse_reads_merged",
        "fillerout": "4_filler_removed",
        "alignout": "5_alignments",
        "countout": "6_read_counts",
    }
    output_dirs = {
        key: os.path.join(abs_directory, value) for key, value in output_dirs.items()
    }

    # Create output directories
    check_and_create_directories(output_dirs.values())

    # Read and validate config
    config = read_config("config.json")

    alignment_counts = {}
    for sample in ["High fluorescence", "Low fluorescence"]:
        logger.info(f"Processing {sample}")
        alignment_counts_file = process_sample(
            config[sample]["Sample ID"],
            config[sample]["R1 filename"],
            config[sample]["R2 filename"],
            output_dirs,
            config["SP library"],
        )
        alignment_counts[sample] = alignment_counts_file
    # Process read counts
    merged_read_counts = process_read_counts(
        alignment_counts["High fluorescence"],
        alignment_counts["Low fluorescence"],
    )
    # Write to file
    merged_read_counts.to_csv(f"{output_dirs['countout']}/read_counts.csv", index=False)


if __name__ == "__main__":
    main()
