import subprocess
import logging

logger = logging.getLogger(__name__)


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
