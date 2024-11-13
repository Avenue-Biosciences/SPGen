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

# Add to path
RUN echo "export PATH=/tools:/tools/FastQC:\$PATH"  >> ~/.bashrc

# Cleanup
RUN rm -r FLASH* samtools*

# Python packages
RUN pip install --upgrade pip && \
    pip install pandas==2.2.3 pytest==8.3.3

# Copy pipeline code and tests
COPY src /src
COPY tests /tests

# Create reference indices
COPY create_reference_indices.sh /
COPY sp_library_fasta /sp_library_fasta
RUN chmod +x /create_reference_indices.sh && \
    /create_reference_indices.sh
