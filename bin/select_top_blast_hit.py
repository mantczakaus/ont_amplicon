#!/usr/bin/env python
import pandas as pd
import argparse
import os.path
from functools import reduce
import pytaxonkit
import numpy as np


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--blastn_results", type=str)
    parser.add_argument("--sample_name", type=str)
    parser.add_argument("--mode", type=str)
    parser.add_argument("--spp_targets", type=str)
    args = parser.parse_args()
    
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    spp_targets = args.spp_targets
    mode = args.mode

    if mode == "ncbi":
        blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species", "sskingdoms"], dtype={"stitle": 'str', "staxids": 'str', "species": 'str'})
        #remove synthetic construct hits
        blastn_results = blastn_results[~blastn_results["species"].str.contains("synthetic construct", na=False)]

    #elif mode == "localdb":
        #retrieve spp name and accession from local db fasta header
        #rearrange column so it matches the one for NCBI
    #    blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "seq_desc", "length", "pident", "mismatch", "gaps", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"], dtype={"stitle": 'str', "staxids": 'int'})
    #    blastn_results['sacc'] = blastn_results['seq_desc'].str.split('|').str[0]
    #    blastn_results['species'] = blastn_results['seq_desc'].str.split('|').str[1]
    #    blastn_results['species'] = blastn_results['species'].str.replace("Species:","")
    #    blastn_results = blastn_results[["qseqid", "sacc", "length", "pident", "mismatch", "gaps", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"]]
 
    blastn_top_hit = blastn_results.drop_duplicates(subset=["qseqid"], keep="first").copy()
    blastn_top_hit["staxids"] = blastn_top_hit["staxids"].str.split(";").str[0].astype(int)

    staxids_list = blastn_top_hit['staxids'].unique().tolist()
    result = pytaxonkit.lineage(staxids_list)
    FullLineage_df = result[['TaxID', 'FullLineage']]
    
    FullLineage_df2 = FullLineage_df.copy()
    FullLineage_df2.rename(columns={"TaxID": "staxids"}, inplace=True)
    FullLineage_df2['staxids']=FullLineage_df2['staxids'].astype(int)
    blast_with_full_phylo_desc_df = blastn_top_hit.merge(FullLineage_df2,  how='right', on='staxids')
    blast_with_full_phylo_desc_df.insert(0, "sample_name", sample_name)
    blast_with_full_phylo_desc_df = blast_with_full_phylo_desc_df[["sample_name", "qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species", "sskingdoms", "FullLineage"]]
    
    blast_with_full_phylo_desc_df['FullLineage'] = blast_with_full_phylo_desc_df['FullLineage'].str.lower().replace(" ","_")
    blast_with_full_phylo_desc_df['FullLineage'] = blast_with_full_phylo_desc_df['FullLineage'].str.replace(" ","_")

    blast_with_full_phylo_desc_df['sskingdoms'] = blast_with_full_phylo_desc_df['sskingdoms'].str.lower()
    organism_target_lower = spp_targets.lower().replace(" ","_")
    blast_with_full_phylo_desc_df["target_organism_match"] = np.where((blast_with_full_phylo_desc_df.sskingdoms.str.contains(organism_target_lower) | (blast_with_full_phylo_desc_df.FullLineage.str.contains(organism_target_lower))), "Y", "N")
    blast_with_full_phylo_desc_df.to_csv(os.path.basename(blastn_results_path).replace("_top_10_hits.txt", "_top_hits_tmp.txt"), index=False, sep="\t")

if __name__ == "__main__":
    main()