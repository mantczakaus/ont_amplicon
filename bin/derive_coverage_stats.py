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
    args = parser.parse_args()
    sample_name = args.sample
    
    samtools_cov = pd.DataFrame()
    blast_df = pd.DataFrame()

    nanostatfile = (sample_name + "_filtered_NanoStats.txt")
    with open(nanostatfile) as f:
        a = " "
        while(a):
            a = f.readline()
            l = a.find("number_of_reads") #Gives a non-negative value when there is a match
            if ( l >= 0 ):
                elements = a.split("\t")
                filtered_read_counts = int(float(elements[1].strip()))
                print(filtered_read_counts)
    f.close()

    for blast in glob("*_megablast_top_hits.txt"):
        blastn_results = pd.read_csv(blast, sep="\t", header=0)
        blast_df = blastn_results.copy()
        blast_df.rename(columns={"length": "alignment_length"}, inplace=True)
#        blast_df = blastn_results[["sample_name", "qseqid", "consensus_seq", "qlen", "qseq", "stitle", "sacc", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "sstrand", "species", "sskingdoms", "FullLineage", "target_organism_match"]]
#        blast_df.columns = ["sample_name", "qseqid", "full_consensus_seq", "full_consensus_length", "qseq", "reference_title", "reference_accession", "reference_length", "pc_ident", "mismatch", "gapopen", "evalue", "bitscore", "query_coverage", "orientation", "species", "kingdom", "full_lineage", "target_organism_match"]
#        print(blast_df)
#        blast_df.loc[:,"qseqid2"]  = '>' + blast_df["qseqid"].astype(str).copy()
#        blast_df.loc[:,"qseq"] = blast_df['qseq'].str.replace('-','').copy()
#        blast_df["query_match_seq"] = blast_df[["qseqid2", "qseq"]].apply("\n".join, axis=1)
    samtools_cov = pd.read_csv(sample_name + "_coverage.txt", sep="\t", usecols=["#rname", "endpos", "numreads", "meandepth"], header=0)
    #samtools_cov.drop(["startpos", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"] , axis=1, inplace=True)
    samtools_cov2 = samtools_cov.copy()
    samtools_cov2.rename(columns={"#rname": "qseqid", "endpos": "query_match_length", "numreads": "qseq_mapping_read_count", "meandepth": "qseq_mean_depth"}, inplace=True)
    samtools_cov2['qseq_pc_mapping_read'] = samtools_cov2['qseq_mapping_read_count'] / filtered_read_counts * 100
    print(samtools_cov2)
    #sum_aligned_read_counts = samtools_cov2['read_count'].sum()
    #sum_unaligned= filtered_read_counts - sum_aligned_read_counts
    #sum_unaligned_pc = sum_unaligned / filtered_read_counts * 100

    mosdepth = pd.read_csv(sample_name + ".thresholds.bed", sep="\t", header=0)
    
    mosdepth.columns = ["qseqid", "start", "end", "region", "base_counts_at_depth_30X"]
    #address why the base_counts_at_depth_30X can be higher than length of the sequence, soft clipping?
    print(mosdepth)
    mosdepth['qseq_pc_depth_30X'] = np.where(mosdepth['base_counts_at_depth_30X'] > mosdepth['end'], 100, mosdepth['base_counts_at_depth_30X'] / mosdepth['end'] * 100)
    mosdepth['qseq_pc_depth_30X'] = mosdepth['qseq_pc_depth_30X'].apply(lambda x: float("{:.1f}".format(x)))
    
    mosdepth_final = mosdepth[["qseqid", "qseq_pc_depth_30X"]]

    summary_dfs = (blast_df, samtools_cov2, mosdepth_final)
    blast_df_final = reduce(lambda left,right: pd.merge(left,right,on=["qseqid"],how='outer').fillna(0), summary_dfs)    

    blast_df_final = blast_df_final.sort_values(["qseq_pc_mapping_read", "target_organism_match"], ascending = (False, False))
    #blast_df_final['status'] = np.where((blast_df_final['pc_depth_30X'] < 90), "FAIL", "PASS")
    
    #Unmapped = {"sample_name" : sample_name, "qseqid" : "Unmapped", "consensus_seq" : "NA", "read_counts" : "NA", "sgi" : "NA", "sacc" : "NA", "length": "NA", "nident": "NA", "pident": "NA", "mismatch" : "NA", "gapopen" : "NA", "qstart" : "NA", "qend": "NA", "qlen" : "NA", "sstart": "NA", "send": "NA", "slen": "NA", "sstrand": "NA", "evalue": "NA", "bitscore": "NA", "query_coverage": "NA", "orientation": "NA", "species" : "NA", "kingdom": "NA", "full_lineage": "NA", "query_match_seq" : "NA", "query_match_length" : "NA", "read_count": sum_unaligned, "pc_read" : sum_unaligned_pc, "mean_depth" : "NA", "pc_depth_30X": "NA", "target_organism_match": "NA"}
    #blast_df_final.loc[len(blast_df_final)] = Unmapped
    #blast_df_final = blast_df_final.reset_index(drop=True)
    blast_df_final['qseq_pc_mapping_read'] = blast_df_final['qseq_pc_mapping_read'].apply(lambda x: float("{:.1f}".format(x)))
    blast_df_final.to_csv(str(sample_name) + "_top_blast_with_cov_stats.txt", index=None, sep="\t")

    print(blast_df_final)


if __name__ == "__main__":
    main()