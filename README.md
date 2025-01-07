# NGS data analysis pipeline

### Overview

This repo contains a complete NGS data analysis pipeline for our amplicon sequencing. For detailed documentation on how to use the pipeline, please see the [Notion page](https://www.notion.so/NGS-data-analysis-pipeline-13c8d1326f288056a743d2204cc3dde8)


The complete logic for the pipeline is in the `src/run_pipeline.py` file.

For now, the Python script calls command line tools directly, instead of using their Python API. This is done to maintain easy compatibility with running commands on the command line. All the tools except Flash seem to have a Python API.


### Docker image

The pipeline is dockerized. To build the image, simply run

```
docker build -t ngs-pipeline:tag .
```

### Testing the pipeline

To test the pipeline, run

```
docker run --rm ngs-pipeline:tag bash -c "PYTHONPATH=. pytest tests"
```

### Installation of tools locally

If for some reason you want to install the tools locally, here are the instructions:

`bowtie2`, `fastqc` and `samtools` can be installed with Homebrew: 

```
brew install bowtie2 fastqc samtools
```

To install `flash`:
- Download the installable from https://github.com/ebiggers/flash/releases/download/v1.2.11/FLASH-1.2.11.tar.gz
```
tar xzf FLASH-1.2.11.tar.gz
cd FLASH-1.2.11
make
sudo cp flash /usr/local/bin
```

To install `cutadapt`:
```
python3 -m pip install --user --upgrade cutadapt
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc
```

To install MUSCLE:

```
wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-osx-arm64.v5.3
# Add to path
```
