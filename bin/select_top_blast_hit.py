#!/usr/bin/env python
import pandas as pd
import argparse
from functools import reduce


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--blastn_results", type=str)
    parser.add_argument("--sample_name", type=str)
    parser.add_argument("--mode", type=str)
    parser.add_argument("--ids", type=str)
    args = parser.parse_args()
    
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    mode = args.mode
    ids = args.ids

    abundant_ids = pd.read_csv(ids, index_col=False, names=["qseqid"])
    if mode == "ncbi":
        blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"], dtype={"stitle": 'str', "staxids": 'str', "species": 'str'})
        #remove synthetic construct hits
        blastn_results = blastn_results[~blastn_results["species"].str.contains("synthetic construct", na=False)]

    elif mode == "localdb":
        #retrieve spp name and accession from local db fasta header
        #rearrange column so it matches the one for NCBI
        blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "seq_desc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"], dtype={"stitle": 'str', "staxids": 'str'})
        blastn_results['sacc'] = blastn_results['seq_desc'].str.split('|').str[0]
        blastn_results['species'] = blastn_results['seq_desc'].str.split('|').str[1]
        blastn_results['species'] = blastn_results['species'].str.replace("Species:","")
        
        blastn_results = blastn_results[["qseqid", "sgi", "sacc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"]]
    #print(blastn_results.dtypes)
    #Fields: query acc., subject acc., subject title, evalue, q. start, q. end, s. start, s. end, bit score, alignment length, mismatches, % identity, identical, % query coverage per subject, subject seq, query seq

    blastn_results_f = blastn_results[["qseqid", "sacc", "species", "stitle", "evalue", "qstart", "qend", "sstart", "send", "bitscore", "length", "mismatch",  "pident", "nident", "qcovs", "qseq", "sseq"]]
    #dfs = (blastn_results_f, abundant_ids)
    #blastn_results_f_red = reduce(lambda left,right: pd.merge(left,right,on=["qseqid"],how='right'), dfs)
    blastn_results_f_red = pd.merge(blastn_results_f, abundant_ids, how='right', on='qseqid')
    blastn_results_f_red.to_csv(sample_name + "_blastn_report.txt", index=False, sep="\t")
    blastn_top_hit = blastn_results.drop_duplicates(subset=["qseqid"], keep="first").copy()
    blastn_top_hit.to_csv(sample_name + "_blastn_top_hits.txt", index=False, sep="\t")
    
    #derive read/contig count per spps
    summary_per_spp = blastn_top_hit['species'].value_counts().rename_axis('species').reset_index(name='count')
    summary_per_spp.to_csv(sample_name + "_spp_abundance.txt", index=False, sep="\t")

    spp = blastn_top_hit[['sacc','species','qseqid']].copy()
    spp['count'] = spp.groupby(['species', 'sacc'])['qseqid'].transform('size')
    #collapse all contigs to given accession number and species
    f = lambda x: x.tolist() if len(x) > 1 else x
    spp = spp.groupby(['species','sacc', 'count'])['qseqid'].agg(f).reset_index().reindex(spp.columns, axis=1)
    #reorder columns before saving
    spp = spp[["species", "sacc", "count", "qseqid"]].sort_values(["count"], ascending=[False])
    spp.to_csv(sample_name + "_queryid_list_with_spp_match.txt", index=False, sep="\t")
    
    blastn_top_hit_f = blastn_top_hit[["qseqid", "qlen", "species", "sacc", "stitle", "slen", "length", "pident", "sstrand", "evalue", "bitscore", "qcovs"]]
    blastn_top_hit_spp = blastn_top_hit_f.sort_values(["evalue", "bitscore"], ascending=[True, False]).groupby("species", as_index=False).first().copy()
    blastn_top_hit_spp.to_csv(sample_name + "_blastn_top_spp_hits.txt", index=False, sep="\t")
    
    summary_per_spp = summary_per_spp.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    blastn_top_hit_spp = blastn_top_hit_spp.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    spp = spp.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    blastn_top_hit_f = blastn_top_hit_f.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }
            .collapsible {
            background-color: #777;
            color: white;
            cursor: pointer;
            padding: 18px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 15px;
            }
            .active, .collapsible:hover {
            background-color: #555;
            }
            .content {
            padding: 0 18px;
            display: none;
            overflow: hidden;
            background-color: #f1f1f1;
            }
            scroll-box {
                        width: 1250px;
                        height: 1250px;
                        overflow: auto;
                        border: 1px solid #ccc;
            }
            </style>
        </head>

        <body>
            <!-- Body Container -->
            <div class="scroll-box">
                <div class="container">
                    <h1>Homology blast results</h1>
                    <button type="button" class="collapsible"> All blast results</button>
                    <div class="content">
                        ''' + blastn_top_hit_f + '''
                    </div>

                    <button type="button" class="collapsible"> Total number of matches to a given species</button>
                    <div class="content">
                        ''' + summary_per_spp + '''
                    </div>

                    <button type="button" class="collapsible"> Total number of matches to specific accession number </h2></button>
                    <div class="content">
                        ''' + spp + '''
                    </div>

                    <button type="button" class="collapsible"> Top match per species based on evalue</h2></button>
                    <div class="content">
                        <p>If for a given species, more than one match share the same top evalue, then bitscore is considered next</p>
                        ''' + blastn_top_hit_spp + '''
                    </div>
                </div>
             </div>

        <script>
        var coll = document.getElementsByClassName("collapsible");
        var i;

        for (i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.display === "block") {
            content.style.display = "none";
            } else {
            content.style.display = "block";
            }
        });
        }
        </script>

        </body>
    </html>'''

    report = open(sample_name + "_blast_report.html", "w")
    report.write(html_string)
    report.close()

if __name__ == "__main__":
    main()