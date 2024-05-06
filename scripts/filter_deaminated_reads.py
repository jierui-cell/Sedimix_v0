import pysam
from pyfaidx import Fasta

def is_deaminated(read, reference):
    # Extract the sequence of the read and its reference positions
    query_sequence = read.query_sequence
    positions = read.get_reference_positions(full_length=True)
    
    # Check for mutations in the first 3 and last 3 positions of the read
    check_positions = positions[:3] + positions[-3:] if len(positions) >= 3 else positions
    check_sequence = query_sequence[:3] + query_sequence[-3:] if len(query_sequence) >= 3 else query_sequence
    
    for pos, base in zip(check_positions, check_sequence):
        if pos is None:
            continue  # Skip positions that do not align to the reference
        ref_base = reference[read.reference_name][pos].seq
        if (base == 'A' and ref_base == 'G') or (base == 'T' and ref_base == 'C'):
            return True
    return False

def filter_reads(input_bam, output_deaminated, output_non_deaminated, reference_path):
    # Open the reference genome
    reference = Fasta(reference_path)
    
    # Open the input BAM file and create output BAM files
    with pysam.AlignmentFile(input_bam, "rb") as bamfile, \
         pysam.AlignmentFile(output_deaminated, "wb", template=bamfile) as deaminated_bam, \
         pysam.AlignmentFile(output_non_deaminated, "wb", template=bamfile) as non_deaminated_bam:
        
        # Iterate over each read in the BAM file
        for read in bamfile:
            if read.is_unmapped:
                continue  # Skip unmapped reads

            # Determine if the read is deaminated
            if is_deaminated(read, reference):
                deaminated_bam.write(read)
            else:
                non_deaminated_bam.write(read)


# Usage
input_bam_path = '../SIM_Set3_centrifuge/mapping/stripped_simulation_s_classified_homo_sapiens.bam'
ref_genome_path = '/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa'
deaminated_bam_path = 'test_deaminated.bam'
non_deaminated_bam_path = 'test_non_deaminated.bam'

filter_reads(input_bam_path, deaminated_bam_path, non_deaminated_bam_path, ref_genome_path)

# mapDamage --merge-libraries -l 200 -y 1 -i test_non_deaminated.bam -r ../../reference_data/human/hg19.fa -d non_deaminated
# mapDamage --merge-libraries -l 200 -y 1 -i test_deaminated.bam -r ../../reference_data/human/hg19.fa -d deaminated