FROM python:3.12-slim

SHELL ["/bin/bash", "-c"]

RUN apt update && apt install -y wget unzip
RUN mkdir /tools

# FastQC
RUN apt install -y default-jre perl
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
RUN unzip fastqc_v0.12.1.zip -d fastqc_v0.12.1 && mv fastqc_v0.12.1/FastQC /tools && rm fastqc_v0.12.1.zip

# cutadapt
RUN apt install -y cutadapt=4.2-1

# FLASH
RUN apt install -y build-essential zlib1g-dev
RUN wget https://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz
RUN tar -xzf FLASH-1.2.11.tar.gz
RUN cd FLASH-1.2.11 &&  make && mv flash /tools

# bowtie2
RUN apt install -y bowtie2=2.5.0-3+b2

# samtools
RUN apt install -y build-essential libncurses-dev libbz2-dev liblzma-dev zlib1g-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
RUN tar -xvjf samtools-1.21.tar.bz2
RUN cd samtools-1.21 && ./configure && make && make install

# MUSCLE
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-aarch64.v5.3
RUN mv muscle-aarch64.v5.3 /tools/muscle && chmod +x /tools/muscle

# Add to path
RUN echo "export PATH=/tools:/tools/FastQC:\$PATH"  >> ~/.bashrc
ENV PATH=/tools:/tools/FastQC:$PATH

# Cleanup
RUN rm -r FLASH* samtools*


# Python packages
COPY requirements.txt /
RUN pip install --upgrade pip && \
    pip install -r requirements.txt

# SignalP
COPY SignalP-6.0/signalp6_slow_sequential/signalp-6-package /signalp-6-package
RUN pip install /signalp-6-package/
RUN SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" ) && \
    cp -r /signalp-6-package/models/* $SIGNALP_DIR/model_weights/
RUN rm -r /signalp-6-package


# Copy pipeline code and tests
COPY src/ /
