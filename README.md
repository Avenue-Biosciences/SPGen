# NGS data analysis pipeline

### Installation of tools locally

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

