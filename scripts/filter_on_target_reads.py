import argparse
import pysam

def load_snp_positions(snp_file):
    snp_positions = {}
    with open(snp_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split()
            chrom = 'chr' + parts[0]
            pos = int(parts[2])  # 1-based position
            snp_positions.setdefault(chrom, set()).add(pos)
    return snp_positions

def filter_reads(bam_file, snp_positions, output_bam):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    out_bam = pysam.AlignmentFile(output_bam, 'wb', header=bam.header)
    
    for read in bam.fetch():
        chrom = bam.get_reference_name(read.reference_id)
        if chrom in snp_positions:
            for pos in read.get_reference_positions():
                if (pos + 1) in snp_positions[chrom]:  # pos + 1 to convert to 1-based
                    out_bam.write(read)
                    break
    
    bam.close()
    out_bam.close()

def main():
    parser = argparse.ArgumentParser(description="Filter reads that overlap with SNP positions")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("snp_file", help="SNP positions file")
    parser.add_argument("output_bam", help="Output BAM file")
    args = parser.parse_args()

    snp_positions = load_snp_positions(args.snp_file)
    filter_reads(args.bam_file, snp_positions, args.output_bam)

if __name__ == "__main__":
    main()
