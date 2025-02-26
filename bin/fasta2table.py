#!/usr/bin/env python
import argparse
import pandas as pd
import os.path
from Bio import SeqIO


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--fasta", type=str, required=True, help='provide fasta file')
    parser.add_argument("--sample", type=str, required=True, help='provide sample name')
    parser.add_argument("--tophits", type=str, required=True, help='provide blast top hits')
    args = parser.parse_args()
    fasta_file = args.fasta
    sample_name = args.sample
    blast = args.tophits

    fasta_df = fasta_to_dataframe(fasta_file)

    blastn_results = pd.read_csv(blast, sep="\t", header=0)
    blastn_results.drop(['sample_name'], axis=1, inplace=True)
    merged_df = pd.merge(fasta_df, blastn_results, on = ['qseqid'], how = 'outer')
    merged_df.insert(0, "sample_name", sample_name)




   # merged_df.to_csv(str(sample_name) + "_blastn_top_hits.txt", index=None, sep="\t")
    merged_df.to_csv(os.path.basename(blast).replace("_top_hits_tmp.txt", "_top_hits.txt"), index=None, sep="\t")

# Function to convert FASTA file to DataFrame
def fasta_to_dataframe(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    
    # List to hold the sequence data
    data = []
    
    for record in records:
        # Append ID and sequence to the list
        data.append([record.id, str(record.seq)])
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=["qseqid", "consensus_seq"])
    return df


if __name__ == "__main__":
    main()
