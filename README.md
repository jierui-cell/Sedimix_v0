# Sedimix_v0

**Sedimix**: An Automated Pipeline for Genomic Analysis  
This project was submitted in partial fulfillment of the MCB honors program.

## Overview
Sedimix is currently under development and will undergo significant improvements over the next few months. This version, v0, primarily serves as a demo test.

**Key Features:**
- Utilizes Snakemake, a workflow management system for Python
- Identifies ancient hominin reads in the samples through mapping and taxonomic classification
- Generates a report file with summary statistics (number of reads, deamination percentage, lineage site information, etc.)

## Requirements
Ensure the following tools are installed and configured before running the pipeline:

### Centrifuge
```bash
git clone https://github.com/DaehwanKimLab/centrifuge
make -C centrifuge
export PATH=$PATH:$(pwd)/centrifuge
```

### Kraken2
```bash
git clone https://github.com/DerrickWood/kraken2.git
./kraken2/install_kraken2.sh kraken2
export PATH=$PATH:$(pwd)/kraken2
```

### Seqtk
```bash
git clone https://github.com/lh3/seqtk.git
make -C seqtk
export PATH=$PATH:$(pwd)/seqtk
```

### BWA
```bash
git clone https://github.com/lh3/bwa.git
make -C bwa
export PATH=$PATH:$(pwd)/bwa
```

### Samtools
```bash
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -xvjf samtools-1.20.tar.bz2
rm samtools-1.20.tar.bz2
./samtools-1.20/configure
make -C samtools-1.20
export PATH=$PATH:$(pwd)/samtools-1.20
```

Add these tools to your `PATH` environment variable permanently by adding the `export PATH` commands to your shell profile file (e.g., `.bashrc` or `.zshrc`).

**Index Files**  
Download index files for Centrifuge and Kraken2 from the following:
- [AWS Indexes for Centrifuge](https://benlangmead.github.io/aws-indexes/centrifuge)
- [AWS Indexes for Kraken2](https://benlangmead.github.io/aws-indexes/k2)

**Human Reference Genome**  
Download the human reference genome hg19.fq.gz from the following:  
- [UCSC Genome Browser hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/)

Alternatively, you can build these indexes yourself by following the instructions provided on their respective GitHub repositories.

## Python and Other Dependencies
To ensure all necessary dependencies are installed, create a conda environment using the provided `environment.yaml` file.

### Conda Environment

Create a conda environment with the following command:

```bash
conda env create -f environment.yaml
```

Alternatively, you can use mamba for faster environment creation:

```bash
mamba env create -f environment.yaml
```

Activate the environment with:

```bash
conda activate sedimix
```

## Pipeline Functionality
1. Takes raw FASTQ files as input.
2. Processes these through BWA to identify reads of interest.
3. Classifies Homo sapiens reads using Centrifuge.
4. Outputs a comprehensive taxonomic report and a folder containing classified hominin reads.

## Usage Instructions
1. Place your input FASTQ files in the `data` folder.
2. Run the pipeline with the following command:

   ```bash
   snakemake -s ../rules/centrifuge.smk --cores {n_cores}
   ```

## Retrieve Your Results
- **Classified Reads**: Located in the `final_reads` folder
- **Data Summary Report**: Located in the `final_report` folder
- **Example Folder**: An example folder can be found in `SIM_Set3_centrifuge`.