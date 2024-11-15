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


def run_command(command, output_file=None):
    """
    Wrapper function to run subprocess commands and log their output.

    Args:
        command (list): Command to run as a list of strings
        output_file (str): Path to output file, if any
    """
    logger.info(f"Running command: {' '.join(command)}")
    try:
        if output_file:
            with open(output_file, "w") as f:
                subprocess.run(
                    command,
                    check=True,
                    stdout=f,
                    stderr=f,
                )
        else:
            result = subprocess.run(command, capture_output=True, check=True, text=True)
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

    r1_output = f"{output_dir}/{sample_id}_adapter_removed_R1.fastq"
    r2_output = f"{output_dir}/{sample_id}_adapter_removed_R2.fastq"
    log_output = f"{output_dir}/{sample_id}_adapter_removed.log"

    run_command(
        [
            "cutadapt",
            "-q 25",
            "-m 50",
            "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGTAGATCTCGGTGGTCGCCGTATCATT",
            r1_file,
            r2_file,
            "-o",
            r1_output,
            "-p",
            r2_output,
        ],
        log_output,
    )

    # Parse the cutadapt output file
    stats = {}
    with open(log_output) as f:
        for line in f:
            if "Total read pairs processed:" in line:
                stats["Total reads"] = int(line.split()[-1].replace(",", ""))
            elif "Read 1 with adapter:" in line:
                stats["Read 1 with adapter"] = int(line.split()[-2].replace(",", ""))
                stats["Read 1 with adapter %"] = float(line.split()[-1].strip("()%"))
            elif "Read 2 with adapter:" in line:
                stats["Read 2 with adapter"] = int(line.split()[-2].replace(",", ""))
                stats["Read 2 with adapter %"] = float(line.split()[-1].strip("()%"))
            elif "Pairs that were too short:" in line:
                stats["Pairs too short"] = int(line.split()[-2].replace(",", ""))
                stats["Pairs too short %"] = float(line.split()[-1].strip("()%"))
            elif "Pairs written (passing filters):" in line:
                stats["Pairs retained after cutting adapters"] = int(
                    line.split()[-2].replace(",", "")
                )
                stats["Pairs retained after cutting adapters %"] = float(
                    line.split()[-1].strip("()%")
                )

    return (
        r1_output,
        r2_output,
        pd.DataFrame(stats, index=[sample_id]),
    )


def merge_reads(sample_id, r1_file, r2_file, output_dir):
    """Merge forward and reverse reads."""
    logger.info("Merging forward and reverse reads")

    output = f"{output_dir}/{sample_id}_merged_reads.fastq"
    log_output = f"{output_dir}/{sample_id}_flash.log"

    os.chdir(output_dir)
    run_command(
        [
            "flash",
            "-M",
            "350",
            r1_file,
            r2_file,
        ],
        log_output,
    )

    # Parse the flash output file
    stats = {}
    with open(log_output) as f:
        for line in f:
            if "Combined pairs:" in line:
                stats["Merged pairs"] = int(line.split()[-1].replace(",", ""))
            elif "Uncombined pairs:" in line:
                stats["Unmerged pairs"] = int(line.split()[-1].replace(",", ""))
            elif "Percent combined:" in line:
                stats["Merged %"] = float(line.split()[-1].strip("%"))

    # Rename flash output files
    flash_files = {
        "out.extendedFrags.fastq": output,
        "out.hist": f"{sample_id}.hist",
        "out.histogram": f"{sample_id}.histogram",
        "out.notCombined_1.fastq": f"{sample_id}_notCombined_R1.fastq",
        "out.notCombined_2.fastq": f"{sample_id}_notCombined_R2.fastq",
    }
    for old, new in flash_files.items():
        os.rename(old, new)
    os.chdir("..")

    return (
        output,
        pd.DataFrame(stats, index=[sample_id]),
    )


def remove_filler_sequence(sample_id, input_file, output_dir, kozak_seq):
    """Remove filler sequence from the merged reads."""
    logger.info("Removing filler sequence")

    output = f"{output_dir}/{sample_id}_reads_no_trimmer.fastq"
    log_output = f"{output_dir}/{sample_id}_cutadapt_filler.txt"

    run_command(
        [
            "cutadapt",
            "-e",
            "0",
            "-g",
            kozak_seq,
            input_file,
            "-o",
            output,
        ],
        log_output,
    )

    # Parse the cutadapt output file
    stats = {}
    with open(log_output) as f:
        for line in f:
            if "Reads with adapter:" in line:
                stats["Reads with filler sequence"] = int(
                    line.split()[-2].replace(",", "")
                )
                stats["Reads with filler sequence %"] = float(
                    line.split()[-1].strip("()%")
                )

    return output, pd.DataFrame(stats, index=[sample_id])


def remove_restriction_site(sample_id, input_file, output_dir):
    """Remove MLuI restriction site from the reads."""
    logger.info("Removing MLuI restriction site and mature protein chain")

    output = f"{output_dir}/{sample_id}_reads_final.fastq"
    log_output = f"{output_dir}/{sample_id}_cutadapt_restriction.txt"

    run_command(
        [
            "cutadapt",
            "-e",
            "0",
            "-a",
            "ACGCGT",
            input_file,
            "-o",
            output,
        ],
        log_output,
    )

    # Parse the cutadapt output file
    stats = {}
    with open(log_output) as f:
        for line in f:
            if "Reads with adapter:" in line:
                stats["Reads with restriction site"] = int(
                    line.split()[-2].replace(",", "")
                )
                stats["Reads with restriction site %"] = float(
                    line.split()[-1].strip("()%")
                )

    return output, pd.DataFrame(stats, index=[sample_id])


def align_reads(sample_id, input_file, output_dir, splib_idx):
    """Align the reads to the SP library."""
    logger.info("Aligning reads to SP library")

    log_output = f"{output_dir}/{sample_id}_bowtie.log"

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
        log_output,
    )

    # Parse the bowtie2 output file
    stats = {}
    with open(log_output) as f:
        for line in f:
            if "reads; of these:" in line:
                stats["Total reads to alignment"] = int(line.split()[0])
            elif "aligned 0 times" in line:
                stats["Unaligned reads"] = int(line.split()[0])
                stats["Unaligned reads %"] = float(line.split("(")[1].split("%")[0])
            elif "aligned exactly 1 time" in line:
                stats["Uniquely aligned reads"] = int(line.split()[0])
                stats["Uniquely aligned reads %"] = float(
                    line.split("(")[1].split("%")[0]
                )
            elif "aligned >1 times" in line:
                stats["Multi-mapped reads"] = int(line.split()[0])
                stats["Multi-mapped reads %"] = float(line.split("(")[1].split("%")[0])
            elif "overall alignment rate" in line:
                stats["Overall alignment rate %"] = float(line.split("%")[0])

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
    return (
        f"{output_dir}/{sample_id}_matches_bowtie_exact.txt",
        pd.DataFrame(stats, index=[sample_id]),
    )


def process_sample(sample_id, r1_file, r2_file, output_dirs, sp_library):
    """Process a single sample through the pipeline."""
    if sp_library in ["R2B", "R2B_full", "integration_test"]:
        kozak_seq = "GCTAGCCCACC"
    elif sp_library in ["v2"]:
        kozak_seq = "GCCGCCACC"
    else:
        raise ValueError(f"Unknown SP library: {sp_library}")

    # Run FastQC
    run_fastqc(r1_file, r2_file, output_dirs["fastqcout"])

    # Remove adapters
    r1_file, r2_file, cutadapt_stats = cut_adapters(
        sample_id, r1_file, r2_file, output_dirs["adaptout"]
    )

    # Merge forward and reverse reads
    merged_reads_file, flash_stats = merge_reads(
        sample_id, r1_file, r2_file, output_dirs["flashout"]
    )

    # Remove filler sequence
    filler_removed_file, filler_stats = remove_filler_sequence(
        sample_id, merged_reads_file, output_dirs["fillerout"], kozak_seq
    )

    # Remove MLuI restriction site
    restriction_removed_file, restriction_stats = remove_restriction_site(
        sample_id, filler_removed_file, output_dirs["fillerout"]
    )

    # Alignment steps
    sp_lib_idx = f"/sp_library_refs/{sp_library}"
    alignment_counts_file, alignment_stats = align_reads(
        sample_id, restriction_removed_file, output_dirs["alignout"], sp_lib_idx
    )
    return alignment_counts_file, pd.concat(
        [cutadapt_stats, flash_stats, filler_stats, restriction_stats, alignment_stats],
        axis=1,
    )


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
    stats = {}
    for sample in ["High fluorescence", "Low fluorescence"]:
        logger.info(f"Processing {sample}")
        alignment_counts_file, stats_sample = process_sample(
            config[sample]["Sample ID"],
            config[sample]["R1 filename"],
            config[sample]["R2 filename"],
            output_dirs,
            config["SP library"],
        )
        alignment_counts[sample] = alignment_counts_file
        stats[sample] = stats_sample

    full_stats = pd.concat(stats.values(), axis=0)
    full_stats.to_csv(f"{output_dirs['countout']}/stats.csv", index=True)

    # Process read counts
    logger.info("Processing read counts")
    merged_read_counts = process_read_counts(
        alignment_counts["High fluorescence"],
        alignment_counts["Low fluorescence"],
    )
    # Write to file
    merged_read_counts.to_csv(f"{output_dirs['countout']}/read_counts.csv", index=False)

    logger.info("Pipeline completed")


if __name__ == "__main__":
    main()
