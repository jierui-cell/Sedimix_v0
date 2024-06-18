# Sedimix_v0

**Sedimix**: An Automated Pipeline for Genomic Analysis  
This project was submitted in partial fulfillment of the MCB honors program.

## Overview
Sedimix is currently under development and will undergo significant improvements over the next few months. This version, v0, primarily serves as a demo test.

**Key Features:**
- Utilizes Snakemake, a workflow management system for Python
- Identifies ancient hominin reads in the samples through mapping and taxonomic classification  
- Generate report file with summary statistics (number of reads, deamination percentage, lineage site information etc.)

## Requirements
Ensure the following tools are installed and configured before running the pipeline:
- **Centrifuge**: [DaehwanKimLab/centrifuge](https://github.com/DaehwanKimLab/centrifuge)
- **Kraken2**: [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2)
- **Seqtk**: [lh3/seqtk](https://github.com/lh3/seqtk)

**Index Files**  
Download index files for Centrifuge and Kraken2 from the following:
- [AWS Indexes for Centrifuge](https://benlangmead.github.io/aws-indexes/centrifuge)
- [AWS Indexes for Kraken2](https://benlangmead.github.io/aws-indexes/k2)

Alternatively, you can build these indexes yourself by following the instructions provided on their respective GitHub repositories.

## Python and Other Dependencies
To ensure all necessary dependencies are installed, create a conda environment using the provided `environment.yaml` file.

### Conda Environment

Create a conda environment with the following command:

```bash
conda env create -f environment.yaml
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
   snakemake -s rules/centrifuge.smk --cores {n_cores} 
   ```
## Retrieve Your Results
- **Classified Reads**: Located in the `final_reads` folder
- **Data Summary Report**: Located in the `final_report` folder
- **Example Folde**r: An example folder can be found in `SIM_Set3_centrifuge`. 