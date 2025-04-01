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
#        blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species", "sskingdoms"], dtype={"stitle": 'str', "staxids": 'str', "species": 'str'})
        columns = ["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart",
                  "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle",
                  "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"]
        blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, header=0, usecols=columns,
                                      dtype={"qseqid": 'str',"sgi": 'str',"sacc": 'str',"length": 'int64',"nident": 'int64',
                                            "pident": 'float64',"mismatch": 'int64',"gaps": 'int64',"gapopen": 'int64',"qstart": 'int64',
                                            "qend": 'int64',"qlen": 'int64',"sstart": 'int64',"send": 'int64',"slen": 'int64',"sstrand": 'str',
                                            "evalue": 'float64',"bitscore": 'float64',"qcovhsp": 'int64',"stitle": 'str',"staxids": 'str',
                                            "qseq": 'str',"sseq": 'str',"sseqid": 'str',"qcovs": 'int64',"qframe": 'int64',"sframe": 'int64'})


    #elif mode == "localdb":
        #retrieve spp name and accession from local db fasta header
        #rearrange column so it matches the one for NCBI
    #    blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "seq_desc", "length", "pident", "mismatch", "gaps", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"], dtype={"stitle": 'str', "staxids": 'int'})
    #    blastn_results['sacc'] = blastn_results['seq_desc'].str.split('|').str[0]
    #    blastn_results['species'] = blastn_results['seq_desc'].str.split('|').str[1]
    #    blastn_results['species'] = blastn_results['species'].str.replace("Species:","")
    #    blastn_results = blastn_results[["qseqid", "sacc", "length", "pident", "mismatch", "gaps", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"]]


    blastn_results["staxids"] = blastn_results["staxids"].str.split(";").str[0].astype(int)
    print(blastn_results)
    staxids_list = blastn_results['staxids'].unique().tolist()
    result = pytaxonkit.lineage(staxids_list)
    FullLineage_df = result[['TaxID', 'FullLineage']]
    FullLineage_df2 = FullLineage_df.copy()
    FullLineage_df2.rename(columns={"TaxID": "staxids"}, inplace=True)
    FullLineage_df2['staxids']=FullLineage_df2['staxids'].astype(int)

    spp_name_results = pytaxonkit.filter(staxids_list, equal_to='species', higher_than='species')
    spp_name = pytaxonkit.name(spp_name_results)
    spp_df = spp_name[['TaxID', 'Name']]
    spp_df2 = spp_df.copy()
    spp_df2.rename(columns={"TaxID": "staxids", "Name": "species"}, inplace=True)
    spp_df2['staxids']=spp_df2['staxids'].astype(int)

    result2 = pytaxonkit.lineage((staxids_list),formatstr="{k}")
    # Check if the expected columns 'TaxID' and 'Lineage' are in the DataFrame
    if 'TaxID' in result2.columns and 'Lineage' in result2.columns:
        # Extract the relevant columns and rename them
        super_kingdom_df2 = result2[['TaxID', 'Lineage']].copy()
        super_kingdom_df2.rename(columns={"TaxID": "staxids", "Lineage": "sskingdoms"}, inplace=True)
    
        # Convert 'staxids' to integers while handling missing values (if any)
        super_kingdom_df2['staxids'] = pd.to_numeric(super_kingdom_df2['staxids'], errors='coerce').astype('Int64')
    else:
        raise KeyError("Columns 'TaxID' and 'Lineage' not found in the result2 DataFrame")

    # Merge the DataFrames
    dfs = [blastn_results, spp_df2, super_kingdom_df2, FullLineage_df2]
    #blast_with_full_phylo_desc_df = blastn_top_hit.merge(dfs,  how='right', on='staxids')
    blast_with_full_phylo_desc_df = reduce(lambda left,right: pd.merge(left,right,on=["staxids"],how='outer'), dfs)
    
    
    blast_with_full_phylo_desc_df.insert(0, "sample_name", sample_name)
    blast_with_full_phylo_desc_df = blast_with_full_phylo_desc_df[~blast_with_full_phylo_desc_df["species"].str.contains("synthetic construct", na=False)]
    blast_with_full_phylo_desc_df = blast_with_full_phylo_desc_df[["sample_name", "qseqid", "sgi", "sacc", "length", "nident",
                                                                   "pident", "mismatch", "gaps", "gapopen", "qstart", "qend",
                                                                   "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore",
                                                                   "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs",
                                                                   "qframe", "sframe", "species", "sskingdoms", "FullLineage"]]
    blast_with_full_phylo_desc_df['FullLineage'] = blast_with_full_phylo_desc_df['FullLineage'].str.lower().replace(" ","_")
    blast_with_full_phylo_desc_df['FullLineage'] = blast_with_full_phylo_desc_df['FullLineage'].str.replace(" ","_")
    blast_with_full_phylo_desc_df['sskingdoms'] = blast_with_full_phylo_desc_df['sskingdoms'].str.lower()
    organism_target_lower = spp_targets.lower().replace(" ","_")

    blastn_top_hit = blast_with_full_phylo_desc_df.drop_duplicates(subset=["qseqid"], keep="first").copy()

    blastn_top_hit["target_organism_match"] = np.where((blastn_top_hit.sskingdoms.str.contains(organism_target_lower) | (blastn_top_hit.FullLineage.str.contains(organism_target_lower))), "Y", "N")
    blastn_top_hit.to_csv(os.path.basename(blastn_results_path).replace("_top_10_hits.txt", "_top_hits_tmp.txt"), index=False, sep="\t")

if __name__ == "__main__":
    main()