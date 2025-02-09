#!/usr/bin/env python
import argparse
import pandas as pd
from functools import reduce
from glob import glob
from subprocess import PIPE

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--sample", type=str, required=True, help='provide sample name')
    args = parser.parse_args()
    sample_name = args.sample
    
    samtools_cov = pd.DataFrame()
    blast_df = pd.DataFrame()
    PCTs_all = pd.DataFrame()

    nanostatfile = (sample_name + "_raw_NanoStats.txt")
    with open(nanostatfile) as f:
        a = " "
        while(a):
            a = f.readline()
            l = a.find("number_of_reads") #Gives a non-negative value when there is a match
            if ( l >= 0 ):
                elements = a.split("\t")
                raw_read_counts = int(float(elements[1].strip()))
                print(raw_read_counts)
    f.close()

    #for blast_results in glob("*_blastn_top_hits.txt"):
    blastn_results = pd.read_csv(sample_name + "_final_polished_consensus_megablast_blastn_top_hits.txt", sep="\t", header=0)
    blast_df = blastn_results[["Sample_name", "qseqid", "stitle", "sacc", "length", "pident", "evalue", "bitscore", "qcovs", "sstrand", "species", "sskingdoms", "FullLineage", "target_organism_match"]]
    blast_df.columns = ["Sample_name", "qseqid", "reference_title", "reference_accession", "reference_length", "pc_ident", "evalue", "bitscore", "query_coverage", "orientation", "species", "kingdom", "full_lineage", "target_organism_match"]
    print(blast_df)
        
    samtools_cov = pd.read_csv(sample_name + "_coverage.txt", sep="\t", usecols=["#rname", "endpos", "numreads"], header=0)
    #samtools_cov.drop(["startpos", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"] , axis=1, inplace=True)
    samtools_cov2 = samtools_cov.copy()
    samtools_cov2.rename(columns={"#rname": "qseqid", "endpos": "length", "numreads": "read_aligning"}, inplace=True)
    samtools_cov2['PCR'] = samtools_cov2['read_aligning'] / raw_read_counts * 100
    print(samtools_cov2)

    mosdepth = pd.read_csv(sample_name + ".regions.bed", sep="\t", header=None)
    mosdepth.columns = ["qseqid", "start", "end", "pc_depth_30X"]
    mosdepth_final = mosdepth[["qseqid", "pc_depth_30X"]]

    
    
    summary_dfs = (blast_df, samtools_cov2, mosdepth_final)
    blast_df_final = reduce(lambda left,right: pd.merge(left,right,on=["qseqid"],how='outer').fillna(0), summary_dfs)
    blast_df_final.to_csv(str(sample_name) + "_top_blast_with_cov_stats.txt", index=None, sep="\t")
    print(blast_df_final)

    



            
if __name__ == "__main__":
    main()