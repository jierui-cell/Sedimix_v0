# write the default snakemake workflow 
# mapping first strategy 

# put all the files in data   

## ------ 
## Small functions to use 
## ------ 
import glob
import os 

# List all .fq and .fq.gz files in the data folder
input_files = glob.glob("data/*.fq") + glob.glob("data/*.fq.gz")        

## ------
## Snakemake workflow 
## ------

# define the final output files
# define the final output files
rule all:
    input:
        expand([
            "final_reads/{sample}_classified_homo_sapiens.fq", 
            "final_reads/{sample}_reads.lst",
            "mapdamage_result/{sample}",  
            "mapdamage_result_after_classification/{sample}"
        ], sample=[os.path.basename(f).split('.fq')[0] for f in input_files])

# Compress all fastq files if not done before
rule gzip_fq_files:
    input:
        lambda wildcards: "data/" + wildcards.filename if not wildcards.filename.endswith('.gz') else "data/" + wildcards.filename[:-3]
    output:
        "data/{filename}.fq.gz"
    benchmark:
        "benchmarks/gzip/{filename}.txt"
    shell:
        """
        gzip -c {input} > {output}
        """

# Map reads using bwa to the hg19 genome, convert to BAM, sort, and index
rule map_sort_index:
    input:
        "data/{sample}.fq.gz"
    output:
        bam="mapping/{sample}.bam",
        bam_bai="mapping/{sample}.bam.bai",
        sai="temp/{sample}.sai",
        sam="temp/{sample}.sam"
    params:
        ref_genome="/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa",
        threads=8
    benchmark:
        "benchmarks/map_sort_index/{sample}.txt"
    shell:
        """
        bwa aln {params.ref_genome} -t {params.threads} {input} > {output.sai}
        bwa samse {params.ref_genome} {output.sai} {input} > {output.sam}
        samtools view -h -b {output.sam} | samtools sort -T {wildcards.sample} -o {output.bam}
        samtools index {output.bam}
        """

# Filter mapped reads with sequence length of at least 30 and mapping quality of at least 25
rule mapping_quality_filter:
    input:
        "mapping/{sample}.bam"
    output:
        "mapping/{sample}_L30MQ25.bam"
    benchmark:
        "benchmarks/mapping_quality_filter/{sample}.txt"
    shell:
        """
        samtools view -h -q 25 {input} |
        awk 'length($10) > 29 || $1 ~ /^@/' |
        samtools view -bS - > {output}
        """

# Detour: add statistical analysis for ancient DNA damage
rule run_mapdamage:
    input:
        sam_file="mapping/{sample}_L30MQ25.bam",
        ref_file="/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa"
    output:
        directory("mapdamage_result/{sample}")
    benchmark:
        "benchmarks/mapdamage/{sample}.txt"
    shell:
        """
        mapDamage --merge-libraries -i {input.sam_file} -r {input.ref_file} -d {output}
        """

# Convert BAM to fastq
rule bam_to_fastq:
    input:
        "mapping/{sample}_L30MQ25.bam"
    output:
        "mapping/{sample}_L30MQ25.fq"
    benchmark:
        "benchmarks/bam_to_fastq/{sample}.txt"
    shell:
        """
        samtools fastq {input} > {output}
        """

# Classification
rule classification:
    input:
        "mapping/{sample}_L30MQ25.fq"
    output:
        kuniq_results="classification/{sample}_L30MQ25.kuniq", 
        kuniq_report="classification/{sample}_L30MQ25.kuniqreport" 
    # change index to nt_db (710GB) if needed
    params:
        kuniq_path="/global/home/users/jieruixu/jieruixu/sediment_dna/peerj_replication/krakenuniq-1.0.4/", 
        kuniq_db="standard_db" # database size: 387GB  
    benchmark:
        "benchmarks/classification/{sample}.txt"
    shell:
        # TODO: add support for multithread (--threads)
        # just try memory mapping, --preload-size 64GB 
        """ 
        {params.kuniq_path}krakenuniq --threads 24 --db {params.kuniq_path}{params.kuniq_db} --output {output.kuniq_results} --report-file {output.kuniq_report} {input}
        """ 

# Generate classified reads
rule generate_classified_reads:
    input:
        kuniq_result="classification/{sample}_L30MQ25.kuniq",
        fq="mapping/{sample}_L30MQ25.fq"
    output:
        fq="final_reads/{sample}_classified_homo_sapiens.fq",
        lst="final_reads/{sample}_reads.lst"
    params:
        taxID=9606  # Homo sapiens
    benchmark:
        "benchmarks/generate_classified_reads/{sample}.txt"
    shell:
        """
        python ../scripts/select_reads_from_kraken2.py {params.taxID} {input.kuniq_result} {wildcards.sample}
        mv {wildcards.sample}_reads.lst final_reads/
        /global/home/users/jieruixu/jieruixu/sediment_dna/seqtk/seqtk subseq {input.fq} {output.lst} > {output.fq}
        """

# Map again, and run mapDamage
rule map_after_classification:
    input:
        "final_reads/{sample}_classified_homo_sapiens.fq"
    output:
        bam="mapping/{sample}_classified_homo_sapiens.bam",
        bam_bai="mapping/{sample}_classified_homo_sapiens.bam.bai",
        sai="temp/{sample}_classified_homo_sapiens.sai",
        sam="temp/{sample}_classified_homo_sapiens.sam"
    params:
        ref_genome="/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa",
        threads=8
    benchmark:
        "benchmarks/map_after_classification/{sample}.txt"
    shell:
        """
        bwa aln {params.ref_genome} -t {params.threads} {input} > {output.sai}
        bwa samse {params.ref_genome} {output.sai} {input} > {output.sam}
        samtools view -h -b {output.sam} | samtools sort -T {wildcards.sample} -o {output.bam}
        samtools index {output.bam}
        """

# Run mapDamage after classification
rule run_mapdamage_after_classification:
    input:
        sam_file="mapping/{sample}_classified_homo_sapiens.bam",
        ref_file="/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa"
    output:
        directory("mapdamage_result_after_classification/{sample}")
    benchmark:
        "benchmarks/run_mapdamage_after_classification/{sample}.txt"
    shell:
        """
        mapDamage --merge-libraries -i {input.sam_file} -r {input.ref_file} -d {output}
        """