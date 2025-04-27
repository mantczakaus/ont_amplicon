#!/usr/bin/env python

import pandas as pd
import argparse
import os
from functools import reduce
import pytaxonkit
import numpy as np

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Load and enrich BLASTn results.")
    parser.add_argument("--blastn_results", required=True, type=str)
    parser.add_argument("--sample_name", required=True, type=str)
    parser.add_argument("--mode", required=True, choices=["ncbi", "localdb"], type=str)
    parser.add_argument("--target_organism", required=True, type=str)
    parser.add_argument("--taxonkit_database_dir", required=True, type=str)
    return parser.parse_args()

def load_blast_results(path, mode):
    """Load BLASTn results based on mode."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"BLASTn results file not found: {path}")

    if mode == "ncbi":
        columns = ["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart",
                   "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle",
                   "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"]

        dtype = {
            "qseqid": 'str', "sgi": 'str', "sacc": 'str', "length": 'int64', "nident": 'int64',
            "pident": 'float64', "mismatch": 'int64', "gaps": 'int64', "gapopen": 'int64', "qstart": 'int64',
            "qend": 'int64', "qlen": 'int64', "sstart": 'int64', "send": 'int64', "slen": 'int64', "sstrand": 'str',
            "evalue": 'float64', "bitscore": 'float64', "qcovhsp": 'int64', "stitle": 'str', "staxids": 'str',
            "qseq": 'str', "sseq": 'str', "sseqid": 'str', "qcovs": 'int64', "qframe": 'int64', "sframe": 'int64'
        }

        df = pd.read_csv(path, sep="\t", header=0, usecols=columns, dtype=dtype)
        df["staxids"] = pd.to_numeric(df["staxids"].str.split(";").str[0], errors='coerce').fillna(0).astype(int)
        top_hit = df.drop_duplicates(subset=["qseqid"], keep="first").copy()

        return top_hit
    else:
        raise NotImplementedError("Mode 'localdb' not implemented yet.")

def enrich_with_taxonomy(df, taxonkit_dir):
    """Add taxonomy information to the dataFrame."""
    #retain unique staxids
    staxids_l = df["staxids"].unique().tolist()

    lineage_df = pytaxonkit.lineage(staxids_l, data_dir=taxonkit_dir)[['TaxID', 'FullLineage']]
    lineage_df.columns = ["staxids", "FullLineage"]
    lineage_df["staxids"] = lineage_df["staxids"].astype(int)
    lineage_df["FullLineage"] = lineage_df["FullLineage"].str.lower().str.replace(" ", "_", regex=False)
    lineage_df["broad_taxonomic_category"] = np.where(
        lineage_df["FullLineage"].str.contains(";virus;"),
        "virus",
        np.where(
            lineage_df["FullLineage"].str.contains(";candidatus_phytoplasma;"),
            "bacteria;phytoplasma",
            np.where(
            (lineage_df["FullLineage"].str.contains(";bacteria;")) &
            (~lineage_df["FullLineage"].str.contains(";candidatus_phytoplasma;")),
            "bacteria;other",
                np.where(
                    lineage_df["FullLineage"].str.contains(";archaea;"),
                    "archaea",
                    np.where(
                        lineage_df["FullLineage"].str.contains(";erysiphaceae;"),
                        "eukaryota;fungi;powdery_mildew",
                        np.where(
                            (lineage_df["FullLineage"].str.contains(";fungi;")) &
                            (~lineage_df["FullLineage"].str.contains(";erysiphaceae;")),
                            "eukaryota;fungi",
                            np.where(
                                (lineage_df["FullLineage"].str.contains(";deuterostomia;")),
                                "eukaryota;deuterostomia",
                                np.where(
                                    (lineage_df["FullLineage"].str.contains(";protostomia;")),
                                    "eukaryota;protostomia",
                                    np.where(
                                        (lineage_df["FullLineage"].str.contains(";eukaryota;")) &
                                        (~lineage_df["FullLineage"].str.contains(";fungi;")) &
                                        (~lineage_df["FullLineage"].str.contains(";deuterostomia;")) &
                                        (~lineage_df["FullLineage"].str.contains(";protostomia;")),
                                        "eukaryota;other",
                                        "unknown"
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
    print(lineage_df)
    #kingdom_df = pytaxonkit.lineage(staxids_l, data_dir=taxonkit_dir, formatstr="{r}")
    #if "TaxID" not in kingdom_df or "Lineage" not in kingdom_df:
    #    raise KeyError("Missing expected columns in lineage data.")
    #kingdom_df = kingdom_df[["TaxID", "Lineage"]]
    #kingdom_df.columns = ["staxids", "sskingdoms"]
    #kingdom_df["staxids"] = pd.to_numeric(kingdom_df["staxids"], errors='coerce').astype('Int64')

    #names_df = pytaxonkit.filter(staxids_l, equal_to="species", rank_file=rank, debug=True)
    names_df = pytaxonkit.name(staxids_l, data_dir=taxonkit_dir)[['TaxID', 'Name']]
    names_df.columns = ["staxids", "species"]
    names_df["staxids"] = names_df["staxids"].astype(int)

    return [df, names_df, lineage_df]

def merge_taxonomy(dfs):
    """Merge taxonomy-enriched data."""
    return reduce(lambda left, right: pd.merge(left, right, on="staxids", how="outer"), dfs)

def filter_and_format(df, sample_name, target_organism):
    """Final formatting, filtering, and matching."""
    df.insert(0, "sample_name", sample_name)
    df = df[~df["species"].str.contains("synthetic construct", na=False)]

    
    target_organism_clean = target_organism.lower().replace(" ", "_")

    #top_hit = df.drop_duplicates(subset=["qseqid"], keep="first").copy()
    df["target_organism_match"] = np.where(
        df["broad_taxonomic_category"].str.contains(target_organism_clean) |
        df["FullLineage"].str.contains(target_organism_clean),
        "Y", "N"
    )

    final_columns = ["sample_name", "qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore",
                     "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe",
                     "species", "broad_taxonomic_category", "FullLineage", "target_organism_match"]

    return df[final_columns]

def main():
    args = parse_arguments()
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    target_organism = args.target_organism
    mode = args.mode
    tk_db_dir = args.taxonkit_database_dir
    #rank = os.path.join(tk_db_dir, "ranks.txt")

    if not os.path.isfile(blastn_results_path):
        raise FileNotFoundError(f"{blastn_results_path} does not exist.")

    blastn_results = load_blast_results(blastn_results_path, mode)
    enriched_dfs = enrich_with_taxonomy(blastn_results, tk_db_dir)
    merged_df = merge_taxonomy(enriched_dfs)
    final_df = filter_and_format(merged_df, sample_name, target_organism)
    out_file = os.path.basename(args.blastn_results).replace("_top_10_hits.txt", "_top_hits_tmp.txt")
    final_df.to_csv(out_file, sep="\t", index=False)
    print(f"Results saved to {out_file}")

if __name__ == "__main__":
    main()