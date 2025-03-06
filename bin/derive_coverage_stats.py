#!/usr/bin/env python
import argparse
import pandas as pd
from functools import reduce
from glob import glob
from subprocess import PIPE
import numpy as np

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--sample", type=str, required=True, help='provide sample name')
    parser.add_argument("--blastn_results", type=str)
    parser.add_argument("--nanostat", type=str)
    parser.add_argument("--bed", type=str)
    parser.add_argument("--coverage", type=str)
    parser.add_argument("--target_size", type=str)

    args = parser.parse_args()
    sample_name = args.sample
    blast = args.blastn_results
    nanostat = args.nanostat
    bed = args.bed
    coverage = args.coverage
    target_size = args.target_size

    
    samtools_cov = pd.DataFrame()
    blast_df = pd.DataFrame()

    with open(nanostat) as f:
        a = " "
        while(a):
            a = f.readline()
            l = a.find("number_of_reads") #Gives a non-negative value when there is a match
            if ( l >= 0 ):
                elements = a.split("\t")
                filtered_read_counts = int(float(elements[1].strip()))
                print(filtered_read_counts)
    f.close()
    
    blastn_results = pd.read_csv(blast, sep="\t", header=0)
    blast_df = blastn_results.copy()
    blast_df.rename(columns={"length": "alignment_length"}, inplace=True)
#        blast_df = blastn_results[["sample_name", "qseqid", "consensus_seq", "qlen", "qseq", "stitle", "sacc", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "sstrand", "species", "sskingdoms", "FullLineage", "target_organism_match"]]
#        blast_df.columns = ["sample_name", "qseqid", "full_consensus_seq", "full_consensus_length", "qseq", "reference_title", "reference_accession", "reference_length", "pc_ident", "mismatch", "gapopen", "evalue", "bitscore", "query_coverage", "orientation", "species", "kingdom", "full_lineage", "target_organism_match"]
#        print(blast_df)
#        blast_df.loc[:,"qseqid2"]  = '>' + blast_df["qseqid"].astype(str).copy()
#        blast_df.loc[:,"qseq"] = blast_df['qseq'].str.replace('-','').copy()
#        blast_df["query_match_seq"] = blast_df[["qseqid2", "qseq"]].apply("\n".join, axis=1)
    #samtools_cov = pd.read_csv(sample_name + "_coverage.txt", sep="\t", usecols=["#rname", "endpos", "numreads", "meandepth"], header=0)
    samtools_cov = pd.read_csv(coverage, sep="\t", usecols=["#rname", "endpos", "numreads", "meandepth"], header=0)
    #samtools_cov.drop(["startpos", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"] , axis=1, inplace=True)
    samtools_cov2 = samtools_cov.copy()
    samtools_cov2.rename(columns={"#rname": "qseqid", 
                                  "endpos": "query_match_length", 
                                  "numreads": "qseq_mapping_read_count", 
                                  "meandepth": "qseq_mean_depth"}, inplace=True)
    samtools_cov2['qseq_pc_mapping_read'] = samtools_cov2['qseq_mapping_read_count'] / filtered_read_counts * 100
    print(samtools_cov2)
    #sum_aligned_read_counts = samtools_cov2['read_count'].sum()
    #sum_unaligned= filtered_read_counts - sum_aligned_read_counts
    #sum_unaligned_pc = sum_unaligned / filtered_read_counts * 100

    mosdepth = pd.read_csv(bed, sep="\t", header=0)
    
    mosdepth.columns = ["qseqid", "start", "end", "region", "base_counts_at_depth_30X"]
    #address why the base_counts_at_depth_30X can be higher than length of the sequence, soft clipping?
    print(mosdepth)
    mosdepth['qseq_pc_depth_30X'] = np.where(mosdepth['base_counts_at_depth_30X'] > mosdepth['end'], 100, mosdepth['base_counts_at_depth_30X'] / mosdepth['end'] * 100)
    mosdepth['qseq_pc_depth_30X'] = mosdepth['qseq_pc_depth_30X'].apply(lambda x: float("{:.1f}".format(x)))
    
    mosdepth_final = mosdepth[["qseqid", "qseq_pc_depth_30X"]]

    summary_dfs = (blast_df, samtools_cov2, mosdepth_final)
    blast_df_final = reduce(lambda left,right: pd.merge(left,right,on=["qseqid"],how='outer').fillna(0), summary_dfs)    

    blast_df_final = blast_df_final.sort_values(["qseq_pc_mapping_read", "target_organism_match"], ascending = (False, False))
    #######30X_DEPTH_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the qseq_pc_depth_30X is >=90
    #ORANGE: If sgi != 0 and the qseq_pc_depth_30X is between 75 and 90.
    #RED: If sgi != 0 and the qseq_pc_depth_30X is <75.
    #GREY: If sgi == 0.

    blast_df_final['30X_DEPTH_FLAG'] = np.where(
        (blast_df_final['sgi'] != 0) & 
        (blast_df_final['qseq_pc_depth_30X'] >= 90), 
        "GREEN",
        np.where((blast_df_final['sgi'] != 0) & 
            (blast_df_final['qseq_pc_depth_30X'] >= 75) & 
            (blast_df_final['qseq_pc_depth_30X'] < 90), 
            "ORANGE",
            np.where((blast_df_final['sgi'] != 0) & 
                (blast_df_final['qseq_pc_depth_30X'] < 75), 
                "RED",
                np.where((blast_df_final['sgi'] == 0) & 
                    (blast_df_final['qseq_pc_depth_30X'] == 0), 
                    "GREY", 
                    ""
                )
            )
        )
    )
    #######MAPPED_READ_COUNT_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the qseq_mapping_read_count is >=200
    #ORANGE: If sgi != 0 and the qseq_mapping_read_count is between 100 and 200.
    #RED: If sgi != 0 and the qseq_mapping_read_count is <100.
    #GREY: If sgi == 0.
    
    blast_df_final['MAPPED_READ_COUNT_FLAG'] = np.where(
        (blast_df_final['sgi'] != 0) & 
        (blast_df_final['qseq_mapping_read_count'] >= 200), 
        "GREEN",
        np.where((blast_df_final['sgi'] != 0) & 
            (blast_df_final['qseq_mapping_read_count'] >= 100) & 
            (blast_df_final['qseq_mapping_read_count'] < 200), 
            "ORANGE",
            np.where((blast_df_final['sgi'] != 0) & 
                (blast_df_final['qseq_mapping_read_count'] < 100), 
                "RED",
                np.where((blast_df_final['sgi'] == 0) & 
                    (blast_df_final['qseq_mapping_read_count'] == 0), 
                    "GREY", 
                    ""
                )
            )
        )
    )
    #######TARGET_ORGANISM_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the target_organism_match is Y
    #RED: If sgi != 0 and the target_organism_match is N
    #GREY: If sgi == 0.

    blast_df_final["TARGET_ORGANISM_FLAG"] = np.where(
        (blast_df_final['sgi'] != 0) & 
        (blast_df_final.target_organism_match.str.match("Y")), 
        "GREEN",
        np.where((blast_df_final['sgi'] != 0) & 
            (blast_df_final.target_organism_match.str.match("N")), 
            "RED",
            np.where((blast_df_final['sgi'] == 0), 
                "GREY",
                ""
            )
        )
    )
    #######TARGET_SIZE_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the query_match_length is within ±20% of the target_size.
    #ORANGE: If sgi != 0 and the query_match_length is between:
    #    Target size + 20% to 40%, OR
    #    Target size - 20% to -40%.
    #RED: If sgi != 0 and the query_match_length is outside the range of ±40% of the target_size.
    #GREY: If sgi == 0.
    
    blast_df_final["TARGET_SIZE_FLAG"] = np.where(
        (blast_df_final['sgi'] != 0) &
        (blast_df_final['query_match_length'] <= float(target_size) + (0.2 * float(target_size))) &
        (blast_df_final['query_match_length'] >= float(target_size) - (0.2 * float(target_size))),
        "GREEN",
         np.where(
            ((blast_df_final['sgi'] != 0) &
            (blast_df_final['query_match_length'] > float(target_size) + (0.2 * float(target_size))) &
            (blast_df_final['query_match_length'] <= float(target_size) + (0.4 * float(target_size)))) |
            ((blast_df_final['sgi'] != 0) &
            (blast_df_final['query_match_length'] < float(target_size) - (0.2 * float(target_size))) &
            (blast_df_final['query_match_length'] >= float(target_size) - (0.4 * float(target_size)))),
            "ORANGE",
            np.where(
                (blast_df_final['sgi'] != 0) &
                ((blast_df_final['query_match_length'] < float(target_size) - (0.4 * float(target_size))) |
                (blast_df_final['query_match_length'] > float(target_size) + (0.4 * float(target_size)))),
                "RED",
                np.where(
                    blast_df_final['sgi'] == 0,
                    "GREY",
                    ""
                )
            )
        )
    )

    blast_df_final['qseq_pc_mapping_read'] = blast_df_final['qseq_pc_mapping_read'].apply(lambda x: float("{:.1f}".format(x)))
    blast_df_final.to_csv(str(sample_name) + "_top_blast_with_cov_stats.txt", index=None, sep="\t")
    print(blast_df_final)
if __name__ == "__main__":
    main()