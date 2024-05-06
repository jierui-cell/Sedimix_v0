import pandas as pd 
import sys 

_, taxID, centrifuge_result, sample = sys.argv 
taxID = int(taxID) 
output_file = f"{sample}_reads.lst"

def extract_homo_sapiens_reads(taxID, centrifuge_file, output_file):
    # Read centrifuge output into pandas DataFrame
    centrifuge_df = pd.read_csv(centrifuge_file, sep='\t', header=0)
    # Filter reads classified as Homo sapiens
    # TODO: not sure why same read map to same taxID multiple times 
    homo_sapiens_reads = set(centrifuge_df[centrifuge_df['taxID'] == taxID]['readID'].tolist())
    
    print (f'Number of total reads: {len(centrifuge_df)}')
    print (f'Number of homo sapien reads classified: {len(homo_sapiens_reads)}')
    
    # prepend \n to output 
    homo_sapiens_reads_output = [f'{read_id}\n' for read_id in homo_sapiens_reads]

    # Write homo_sapiens_reads to homo_sapiens.lst file
    with open(output_file, 'w') as file:
        file.writelines(homo_sapiens_reads_output) 

extract_homo_sapiens_reads(taxID, centrifuge_result, output_file) 
