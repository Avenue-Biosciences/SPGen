import os
import sys
import json
import subprocess
import argparse
import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sqlalchemy import create_engine, text, Engine

logger = logging.getLogger(__name__)


def check_and_create_directories(output_dirs):
    """Create output directories after checking they don't exist."""
    for dir_name in output_dirs:
        if os.path.exists(dir_name):
            logger.warning(
                f"Warning: '{dir_name}' directory already exists. Existing files may be overwritten."
            )
        else:
            os.makedirs(dir_name)


def read_input(input_file):
    """Read and validate the configuration file."""
    with open(input_file) as f:
        input = json.load(f)

    required_fields = ["Replicates", "SP library"]
    for field in required_fields:
        if field not in input:
            raise ValueError(
                f"Error: {field} was not found. Please check the input file."
            )

    # Validate values
    required_rep_fields = ["High fluorescence", "Low fluorescence"]
    required_sample_fields = ["Sample ID", "R1 filename", "R2 filename"]
    for i, replicate in enumerate(input["Replicates"]):
        for field in required_rep_fields:
            value = replicate[field]
            if not value:
                raise ValueError(
                    f"Error: replicate {i} {field} was not found. Please check the input file."
                )
            for s_field in required_sample_fields:
                value = replicate[field][s_field]
                if not value:
                    raise ValueError(
                        f"Error: replicate {i} {field} {s_field} was not found. Please check the input file."
                    )
    # Check that all sample IDs are unique
    sample_ids = [
        replicate[sample]["Sample ID"]
        for replicate in input["Replicates"]
        for sample in ["High fluorescence", "Low fluorescence"]
    ]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(
            "Error: Duplicate sample IDs found. Please check the input file."
        )

    return input


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

    output = f"{output_dir}/{sample_id}_reads_no_filler.fastq"
    log_output = f"{output_dir}/{sample_id}_cutadapt_filler.log"

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
            if "Reads with adapters:" in line:
                stats["Reads with filler sequence"] = int(
                    line.split()[-2].replace(",", "")
                )
                stats["Reads with filler sequence %"] = float(
                    line.split()[-1].strip("()%")
                )

    return output, pd.DataFrame(stats, index=[sample_id])


def remove_restriction_site(sample_id, input_file, min_length, output_dir):
    """Remove MLuI restriction site from the reads."""
    logger.info("Removing MLuI restriction site and mature protein chain")

    output = f"{output_dir}/{sample_id}_reads_final.fastq"
    log_output = f"{output_dir}/{sample_id}_cutadapt_restriction.log"

    run_command(
        [
            "cutadapt",
            "-e",
            "0",
            "-m",
            str(min_length),
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
            if "Reads with adapters:" in line:
                stats["Reads with restriction site"] = int(
                    line.split()[-2].replace(",", "")
                )
                stats["Reads with restriction site %"] = float(
                    line.split()[-1].strip("()%")
                )
            elif "Reads written (passing filters):" in line:
                stats["Reads retained after restriction site removal"] = int(
                    line.split()[-2].replace(",", "")
                )
                stats["Reads retained after restriction site removal %"] = float(
                    line.split()[-1].strip("()%")
                )

    return output, pd.DataFrame(stats, index=[sample_id])


def get_db_config(config_file: str):
    with open(config_file) as f:
        config = json.load(f)

    # Check that all required fields are present
    required_fields = ["db_name", "db_user", "db_password", "db_host", "db_port"]
    for field in required_fields:
        if field not in config:
            raise ValueError(
                f"Error: {field} was not found. Please check the db_config.json file."
            )

    return config


def get_db_engine(db_config: dict) -> Engine:
    return create_engine(
        f"postgresql://{db_config['db_user']}:{db_config['db_password']}@{db_config['db_host']}:{db_config['db_port']}/{db_config['db_name']}"
    )


def build_library_reference(sp_library: str, db_config: dict):
    """Build the SP library reference."""
    logger.info("Building SP library reference")

    # Download the SP library from the database
    engine = get_db_engine(db_config)

    sp_library_df = pd.read_sql(
        text(
            "SELECT sp.name, sp_to_sp_library.dna_sequence_sp_only "
            "FROM sp_library "
            "JOIN sp_to_sp_library ON sp_library.id = sp_to_sp_library.sp_library_id "
            "JOIN sp ON sp_to_sp_library.sp_id = sp.id "
            f"WHERE sp_library.name = '{sp_library}'"
        ),
        engine,
    )

    records = [
        SeqRecord(Seq(row["dna_sequence_sp_only"]), id=row["name"], description="")
        for _, row in sp_library_df.iterrows()
    ]

    os.makedirs("/sp_library_fastas", exist_ok=True)
    os.makedirs("/sp_library_refs", exist_ok=True)
    fasta_file = os.path.join("/sp_library_fastas", f"{sp_library}.fasta")
    SeqIO.write(records, fasta_file, "fasta")
    run_command(
        [
            "bowtie2-build",
            fasta_file,
            os.path.join("/sp_library_refs", f"{sp_library}"),
        ],
        f"{sp_library}_bowtie_build.log",
    )


def align_reads(sample_id, input_file, output_dir, splib_idx):
    """Align the reads to the SP library."""
    logger.info("Aligning reads to SP library")

    log_output = f"{output_dir}/{sample_id}_bowtie.log"

    os.chdir(output_dir)
    run_command(
        [
            "bowtie2",
            "-p",
            "4",
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
        min_length = 21
    elif sp_library in ["v2"]:
        kozak_seq = "GCCGCCACC"
        min_length = 21
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
        sample_id, filler_removed_file, min_length, output_dirs["fillerout"]
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


def read_alignment_counts(filename, label):
    """Read the alignment counts for a single sample."""
    df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["SP", "length", f"{label}", "unmapped"],
    )
    df.loc[df["SP"] == "*", f"{label}"] = df.loc[df["SP"] == "*", "unmapped"]
    df.loc[df["SP"] == "*", "SP"] = "unmapped"
    return df.drop(columns=["length", "unmapped"])


def process_read_counts(files):
    """Process the read counts for the high and low fluorescence samples."""
    logger.info("Processing read counts")

    read_counts = [read_alignment_counts(file, label) for label, file in files.items()]

    # Start with the first dataframe
    merged_read_counts = read_counts[0]

    # Merge with remaining dataframes one by one
    for df in read_counts[1:]:
        merged_read_counts = pd.merge(merged_read_counts, df, on="SP")

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

    # Read and validate input
    input = read_input("input.json")

    # Set up logging
    logging.basicConfig(
        filename="run_pipeline.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger.info("Starting pipeline")

    # Define output directories
    output_dirs = {
        "fastqcout": "1_fastqc",
        "adaptout": "2_adapter_removed",
        "flashout": "3_forward_reverse_reads_merged",
        "fillerout": "4_filler_restriction_removed",
        "alignout": "5_alignments",
        "countout": "6_read_counts",
    }
    output_dirs = {
        key: os.path.join(abs_directory, value) for key, value in output_dirs.items()
    }

    # Create output directories
    check_and_create_directories(output_dirs.values())

    db_config = get_db_config("db_config.json")

    # Build the SP library reference
    build_library_reference(input["SP library"], db_config)

    alignment_counts = {}
    stats = []
    for i, replicate in enumerate(input["Replicates"]):
        for sample in ["High fluorescence", "Low fluorescence"]:
            logger.info(f"Processing replicate {i+1} {sample}")
            alignment_counts_file, stats_sample = process_sample(
                replicate[sample]["Sample ID"],
                replicate[sample]["R1 filename"],
                replicate[sample]["R2 filename"],
                output_dirs,
                input["SP library"],
            )
            label = f"{'HF' if sample == 'High fluorescence' else 'LF'}{i+1}"
            alignment_counts[label] = alignment_counts_file
            stats.append(stats_sample)

    full_stats = pd.concat(stats, axis=0)
    full_stats.to_csv(f"{output_dirs['countout']}/stats.csv", index=True)

    # Process read counts
    logger.info("Processing read counts")
    merged_read_counts = process_read_counts(alignment_counts)
    # Write to file
    merged_read_counts.to_csv(f"{output_dirs['countout']}/read_counts.csv", index=False)

    logger.info("Pipeline completed")


if __name__ == "__main__":
    main()
