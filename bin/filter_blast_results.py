#!/usr/bin/env python
import pandas as pd
import argparse
import numpy as np
from functools import reduce
from Bio.Seq import Seq





def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--blastn_results", type=str)
    parser.add_argument("--sample_name", type=str)
    parser.add_argument("--acc_phylo_info", type=str)
    parser.add_argument("--spp_targets", type=str)
    parser.add_argument("--gene_targets", type=str)
    parser.add_argument("--mode", type=str)
    parser.add_argument("--target_size", type=str)
    args = parser.parse_args()
    
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    mode = args.mode
    acc_phylo_info = args.acc_phylo_info
    spp_targets = args.spp_targets
    gene_targets = args.gene_targets
    target_size = args.target_size

    blast = pd.read_csv(blastn_results_path, sep="\t", index_col=False, header=0, names=["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species", "sskingdoms"], dtype={"stitle": 'str', "staxids": 'str', "species": 'str'})
    #, index_col=False,  
    #deal with header duplication, fix later
    blast = blast.drop_duplicates()

    phylo = pd.read_csv(acc_phylo_info, sep="\t", index_col=False, names=["phylogenetic_description", "sacc"])
    #blast_phylo = pd.merge(blast, phylo, how='left', on='sacc')
    #blast_phylo = blast.reset_index().merge(phylo,  how='left', on='sacc').set_index('qseqid')
    blast_phylo = blast.reset_index().merge(phylo,  how='left', on='sacc')
    blast_phylo = blast.merge(phylo,  how='left', on='sacc')
    blast_phylo['stitle'] = blast_phylo['stitle'].str.replace(" ","_")
    #blast_phylo['qlen'] = blast_phylo['qlen'].astype(int)
    blast_phylo['qlen'] = blast_phylo['qlen'].apply(pd.to_numeric, errors='coerce')
    #df['Age'] = df['Age'].apply(pd.to_numeric, errors='coerce')
    blast_phylo['phylogenetic_description'] = blast_phylo['phylogenetic_description'].str.replace(", ","_")
    print(blast_phylo)
    #blast_phylo_f = blast_phylo[~blast_phylo["phylogenetic_description"].str.contains(spp_targets, na=False)]
    if gene_targets == "COI":
        gene_terms_to_search = ['mitochondrion', 'cytochrome_oxidase_subunit_1', 'cytochrome_c_oxidase_subunit_I', 'COX1', 'COI']
        #gene_terms_to_search[gene_terms_to_search.str.contains('|'.join(gene_terms_to_search))]
    elif gene_targets == "16s":
        gene_terms_to_search = ['chromosome']
    if spp_targets == 'insect':
        spp_terms_to_search = ('Insecta')
    elif spp_targets == 'bacteria':
        spp_terms_to_search = ('Bacteria')
        #blast_phylo["target_consensus"] = np.where((blast_phylo.stitle.str.contains("mitochondrion|cytochrome_oxidase_subunit_1|cytochrome_c_oxidase_subunit_I|(COX1)|(COI)")&(blast_phylo.phylogenetic_description.str.contains("Insecta"))), "Y", "N")
    blast_phylo["target_type_match"] = np.where((blast_phylo.stitle.str.contains('|'.join(gene_terms_to_search))&(blast_phylo.phylogenetic_description.str.contains(spp_terms_to_search))), "Y", "N")
    blast_phylo["target_size_match"] = np.where(abs(blast_phylo['qlen'] - (int(target_size))) < 50, "Y", "N")
    
    #blastn_top_hit.to_csv(sample_name + "_blastn_top_hits.txt", index=False, sep="\t")
    #blast_phylo["real_orientation"] = np.where(blast_phylo.sstrand.str.contains("minus"), reverse_complement(blast_phylo.sstrand) , blast_phylo.sstrand)
    blast_phylo["real_orientation"] = blast_phylo.apply(lambda x: reverse_complement(x['sstrand'], x['qseq']), axis=1)

    
    print(blast_phylo)
    #derive read/contig count per spps
    blast_phylo.to_csv(sample_name + "_blastn_top_hits_filtered.txt", index=False, sep="\t")
    #reorder columns before saving

#def reverse_complement(x):
#    seq = Seq(x)
#    y = seq.reverse_complement(x)
#    return y

def reverse_complement(x,y):
    seq = Seq(y).replace("-","")
    if "minus" in str(x):
        print(seq)
        yrc = seq.reverse_complement()
        return str(yrc)
    elif "plus" in str(x):
        return y


#cutadapt -g "ATAAAGATATTGG"  -a "ATTTTTTGGTCAC" --times 2  -o  trial.fa barcode01_VE24-1279_COI_medaka.consensus.fasta

if __name__ == "__main__":
    main()