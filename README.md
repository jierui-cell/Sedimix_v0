# Sedimix_v0

**Sedimix**: An Automated Pipeline for Genomic Analysis  
This project was submitted in partial fulfillment of the MCB honors program.

## Overview
Sedimix is currently under development and will undergo significant improvements over the next few months. This version, v0, primarily serves as a demo test.

**Key Features:**
- Written in Bash and Python
- Utilizes Snakemake, a workflow management system for Python
- Outputs a taxonomic report of the genomic samples
- Identifies and classifies hominin reads in the samples

## Requirements
Ensure the following tools are installed and configured before running the pipeline:
- **Centrifuge**: [DaehwanKimLab/centrifuge](https://github.com/DaehwanKimLab/centrifuge)
- **Kraken2**: [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2)
- **Snakemake**: A Python workflow management system
- **Seqtk**: [lh3/seqtk](https://github.com/lh3/seqtk)
- **BWA**: A fast, accurate read mapping software
- **Samtools**: Utilities for manipulating alignments in the SAM/BAM format

**Index Files**  
Download index files for Centrifuge and Kraken2 from the following:
- [AWS Indexes for Centrifuge](https://benlangmead.github.io/aws-indexes/centrifuge)
- [AWS Indexes for Kraken2](https://benlangmead.github.io/aws-indexes/k2)

Alternatively, you can build these indexes yourself by following the instructions provided on their respective GitHub repositories.

## Pipeline Functionality
1. Takes raw FASTQ files as input.
2. Processes these through BWA to identify reads of interest.
3. Classifies Homo sapiens reads using Centrifuge.
4. Outputs a comprehensive taxonomic report and a folder containing classified hominin reads.

## Usage Instructions
1. Place your input FASTQ files in the `data` folder.
2. Run the pipeline with the following command:

   ```bash
   snakemake -s rules/centrifuge.smk -c1
   ```
## Retrieve Your Results
- **Classified Reads**: Located in the `final_reads` folder
- **Data Summary Report**: Located in the `final_report` folder
- **Example Folde**r: An example folder can be found in SIM_Set3_centrifuge. 