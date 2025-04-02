# ont_amplicon
## Introduction
## Pipeline overview
## Installation
## Running the pipeline  
## Output files
## Authors



## Introduction

ont_amplicon is a Nextflow-based bioinformatics pipeline designed to derive consensus sequences from:
 amplicon sequencing data that was generated using rapid library preparation kit from Oxford nanopore Technologies. 

It takes compressed fastq files as input.


## Pipeline overview
- Data quality check (QC) and preprocessing
  - Merge fastq files (Fascat, optional)
  - Raw fastq file QC (Nanoplot)
  - Trim adaptors (PoreChop ABI - optional)
  - Filter reads based on length and/or quality (Chopper - optional)
  - Reformat fastq files so read names are trimmed after the first whitespace (bbmap)
  - Processed fastq file QC (if PoreChop and/or Chopper is run) (Nanoplot)
- QC report
  - Derive read counts recovered pre and post data processing and post host filtering
- Clustering mode
  - Read clustering (Rattle)
  - Convert fastq to fasta format (seqtk)
  - Polishing (Minimap, Racon, Medaka2, samtools - optional)
  - Remove adapters if provided (Cutadapt)
  - Megablast homology search against COI database (if COI is targetted) and reverse complement where required
  - Megablast homology search against NCBI database
  - Derive top candidate hits, assign preliminary taxonomy and target organism flag (pytaxonkit)
  - Map reads back to segment of consensus sequence that align to reference and derive BAM file and alignment statistics (minimap, samtools and mosdepth)
  - Map reads to segment of NCBI reference sequence that align to consensus and derive BAM file and consensus (minimap, samtools)


## Installation
### Requirements  
1. Install Java if not already on your system. Follow the java download instructions provided on this page [`page`](https://www.nextflow.io/docs/latest/getstarted.html#installation).

2. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

3. Install [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) to suit your environment. The pipeline has been validated using singularity version 3.10.2-1 and apptainer version 1.3.6-1.el9 but has not yet been tested with version 4.

3. Install taxonkit using the script install_taxonkit.sh or follow the steps described on this page [`page`](https://bioinf.shenwei.me/taxonkit/download/).


4. Install NCBI NT or coreNT.  
Download a local copy of the NCBI database of interest, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:
  ```
  mkdir blastDB/20230930
  ```
  You will need to use a current update_blastdb.pl script from the blast+ version used with the pipeline (ie 2.16.0).
  For example:
  ```
  singularity exec -B /scratch https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4  update_blastdb.pl --decompress nt
  singularity exec -B /scratch https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4  update_blastdb.pl --decompress taxdb
  tar -xzf taxdb.tar.gz
  ```
  
  Specify the path of your local NCBI blast directories in your nextflow command using ```--blastn_db = '/full/path/to/blastDB/20230930/nt'``` or specify the following lines in a user config file.
  For instance:
  ```
  params {
    --blastn_db = '/full/path/to/blastDB/20230930/nt'
  }
  ```

5. Download the Cytochrome oxydase 1 (COI1) database if you are planning to analyse COI samples.
  ```
  git clone https://github.com/bachob5/MetaCOXI.git
  #extract MetaCOXI_Seqs.fasta from the MetaCOXI_Seqs.tar.gz file
  tar -xvf MetaCOXI_Seqs.tar.gz
  #make a blast database from the fasta file
  singularity exec -B /scratch https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4 makeblastdb -in MetaCOXI_Seqs.fasta -parse_seqids -dbtype prot
  ```
Specify the path of your COI database in your nextflow command using ```--blastn_COI = '/full/path/to/MetaCOXI_Seqs.fasta'``` or specify the following lines in a user config file.
  For instance:
  ```
  params {
    --blastn_COI = '/full/path/to/MetaCOXI_Seqs.fasta'
  }
  ```

## Running the pipeline  

### Run test data
- Run the command:
  ```
  nextflow run maelyg/ont_amplicon -profile singularity --samplesheet index.csv
  ```
  The first time the command runs, it will download the pipeline into your assets.  

  The source code can also be downloaded directly from GitHub using the git command:
  ```
  git clone https://github.com/maelyg/ont_amplicon
  ```

- Provide an index.csv file.  
  Create a **comma separated file (csv)** that will be the input for the workflow. 
  
  **Please note**: it is best to edit the csv file with an editor that does not add special characters/symbols (e.g. VSCode or Atom). If using other editors, check your files and if necessary, run dos2unix[`dos2unix`](https://www.linuxfromscratch.org/blfs/view/git/general/dos2unix.html) on the file to get rid of these unwanted characters/symbols as they will cause the pipeline to fail.  
  
  By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the ```--samplesheet [filename]``` in the nextflow run command, as long as it has a **.csv** suffix. This text file requires the following columns (which need to be included as a header): ```sampleid,sample_files,spp_targets,gene_targets,target_size,fwd_primer,rev_primer``` 

   - **sampleid** will be the sample name that will be given to the files created by the pipeline (required).  
   - **sample_path** is the full path to the fastq files that the pipeline requires as starting input (required).  
   - **spp_targets** is the organism targetted by the PCR (required).  
   - **gene_targets** is the gene targetted by the PCR (optional).  
   - **target_size** is the expected size of the amplicon (required).  
   - **fwd_primer** is the nucleotide sequence of the FWD primer (optional).  
   - **rev_primer** is the nucleotide sequence of the REV primer (optional).  



  For the fastq files path, the pipeline is currently expecting either 1) multiple fastq.gz files per sample located within one folder or 2) a single fastq.gz file per sample.  
  If there are **multiple fastq.gz files per sample**, their full path can be specified on one line using **an asterisk (*fastq.gz)** and you will need to specify the parameter ```--merge``` either on the command line or a config file.  
  See an example of an index.csv file for 2 MTDT samples:  
  ```
  sampleid,sample_files,spp_targets,gene_targets,target_size,fwd_primer,rev_primer
  VE24-1279_COI,tests/mtdt_data/barcode01_VE24-1279_COI/*fastq.gz,drosophilidae,COI,711,GGTCAACAAATCATAAAGATATTGG,ATTTTTTGGTCACCCTGAAGTTTA
  MP24-1051A_16S,tests/mtdt_data/barcode06_MP24-1051A_16S/*fastq.gz,bacteria,16s,1509,AGAGTTTGATCATGGCTCAG,AAGTCGTAACAAGGTAACCGT
  MP24-1096B_gyrB,tests/mtdt_data/barcode19_MP24-1096B_gyrB/*fastq.gz,bacteria,gyrB,1258,GAAGTCATCATGACCGTTCTGCAYGCNGGNGGNAARTTYGA,ATGACNGAYGCNGAYGTNGAYGGCTCGCACATCCGTACCCTGCT
  ```
  For samples with a single fastq.gz file, specify **the full path to the fastq.gz file.**


- Specify a profile:
  ```
  nextflow run maelyg/ont_amplicon -profile singularity --samplesheet index_example.csv
  ```
  setting the profile parameter to one of ```docker``` or ```singularity``` to suit your environment.
  
- Specify the analysis mode: ```--analysis_mode clustering``` (this is set to clustering by default; haivng this parameter in place will enable us to add other analysis modes like mapping to ref down the track if required).  

- Specify the ``--analyst_name`` and the ``--facility`` either on the nextflow command or in a config file.  The analysis cannot proceed without these being set.  

- To set additional parameters, you can either include these in your nextflow run command:
  ```
  nextflow run maelyg/ont_amplicon -profile singularity --samplesheet index_example.csv --adapter_trimming
  ```
  or set them to true in the nextflow.config file.
  ```
  params {
    adapter_trimming = true
  }
  ```

- Two tests are currently provided to check if the pipeline was successfully installed and demonstrate the current outputs generated by the pipeline prototype. The mtdt_test runs three samples provided by  MTDT (barcode01_VE24-1279_COI, barcode06_MP24-1051A_16S and barcode19_MP24-1096B_gyrB). The peq_test runs two samples provided by PEQ (ONT141 and ONT142).  

To use the tests, change directory to ont_amplicon and run the following command for the MTDT test:
  ```
  nextflow run main.nf -profile mtdt_test,singularity
  ```
  and this command for the PEQ test:  
  ```
  nextflow run main.nf -profile peq_test,singularity
  ```
The tests should take less than 5 minutes to run to completion.

If the installation is successful, it will generate a results/test folder with the following structure:
```
results/
├── barcode01_VE24-1279_COI
│   ├── clustering
│   │   └── barcode01_VE24-1279_COI_rattle.fasta
│   ├── html_report
│   │   ├── bam-alignment.html
│   │   ├── example_report_context.json
│   │   └── report.html
│   ├── mapping_to_consensus
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam.bai
│   │   ├── barcode01_VE24-1279_COI_coverage.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_match.fasta
│   │   ├── barcode01_VE24-1279_COI_histogram.txt
│   │   ├── barcode01_VE24-1279_COI.per-base.bed
│   │   └── barcode01_VE24-1279_COI_top_blast_with_cov_stats.txt
│   ├── mapping_to_ref
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam.bai
│   │   ├── barcode01_VE24-1279_COI_coverage.txt
│   │   ├── barcode01_VE24-1279_COI_histogram.txt
│   │   ├── barcode01_VE24-1279_COI_reference_match.fasta
│   │   └── barcode01_VE24-1279_COI_samtools_consensus_from_ref.fasta
│   ├── megablast
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_match.fasta
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_megablast_COI_top_hit.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc.fasta
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_10_hits_temp.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_10_hits.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_hits.txt
│   │   └── barcode01_VE24-1279_COI_reference_match.fasta
│   ├── polishing
│   │   ├── barcode01_VE24-1279_COI_cutadapt.log
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_medaka_consensus.bam
│   │   ├── barcode01_VE24-1279_COI_medaka_consensus.bam.bai
│   │   ├── barcode01_VE24-1279_COI_medaka_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_preprocessed.fastq.gz
│   │   ├── barcode01_VE24-1279_COI_racon_polished.fasta
│   │   └── barcode01_VE24-1279_COI_samtools_consensus.fasta
│   ├── preprocessing
│   │   ├── barcode01_VE24-1279_COI_basecalling_model_inference.txt
│   │   ├── barcode01_VE24-1279_COI_preprocessed.fastq.gz
│   │   ├── chopper
│   │   │   └── barcode01_VE24-1279_COI_chopper.log
│   │   └── porechop
│   │       └── barcode01_VE24-1279_COI_porechop.log
│   └── qc
│       ├── fastcat
│       │   ├── barcode01_VE24-1279_COI.fastq.gz
│       │   ├── barcode01_VE24-1279_COI_stats.tsv
│       │   └── histograms
│       │       ├── length.hist
│       │       └── quality.hist
│       └── nanoplot
│           ├── barcode01_VE24-1279_COI_filtered_LengthvsQualityScatterPlot_dot.html
│           ├── barcode01_VE24-1279_COI_filtered_NanoPlot-report.html
│           ├── barcode01_VE24-1279_COI_filtered_NanoStats.txt
│           ├── barcode01_VE24-1279_COI_raw_LengthvsQualityScatterPlot_dot.html
│           ├── barcode01_VE24-1279_COI_raw_NanoPlot-report.html
│           └── barcode01_VE24-1279_COI_raw_NanoStats.txt
└── qc_report
    ├── run_qc_report_20250401-210340.html
    └── run_qc_report_20250401-210340.txt
```

### QC step
By default the pipeline will run a quality control check of the raw reads using NanoPlot.

- It is recommended to first run only the quality control step to have a preliminary look at the data before proceeding with downstream analyses by specifying the ```--qc_only``` parameter.

The command you would run would look like this:
```
nextflow run maelyg/ont_amplicon -profile singularity \
                            --merge \
                            --qc_only
```

### Preprocessing reads
If multiple fastq files exist for a single sample, they will first need to be merged using the `--merge` option using [`Fascat`](https://github.com/epi2me-labs/fastcat).
Then the read names of the fastq file created will be trimmed after the first whitespace, for compatiblity purposes with all downstream tools.  

Reads can also be optionally trimmed of adapters and/or quality filtered:  
- Search for presence of sequencing adapters in sequences reads using [`Porechop ABI`](https://github.com/rrwick/Porechop) by specifying the ``--adapter_trimming`` parameter. Porechop ABI parameters can be specified using ```--porechop_options '{options} '```, making sure you leave a space at the end before the closing quote. Please refer to the Porechop manual.  

  **Special usage:**  
  To limit the search to known adapters listed in [`adapter.py`](https://github.com/bonsai-team/Porechop_ABI/blob/master/porechop_abi/adapters.py), just specify the ```--adapter_trimming``` option.  
  To search ab initio for adapters on top of known adapters, specify ```--adapter_trimming --porechop_options '-abi '```.  
T  o limit the search to custom adapters, specify ```--adapter_trimming --porechop_custom_primers --porechop_options '-ddb '``` and list the custom adapters in the text file located under bin/adapters.txt following the format:  
    ```
     line 1: Adapter name
     line 2: Start adapter sequence
     line 3: End adapter sequence
     --- repeat for each adapter pair---
    ```

- Perform a quality filtering step using [`Chopper`](https://github.com/wdecoster/chopper) by specifying the ```--qual_filt``` parameter. The following parameters can be specified using the ```--chopper_options '{options}'```. Please refer to the Chopper manual.  
For instance to filter reads shorter than 1000 bp and longer than 20000 bp, and reads with a minimum Phred average quality score of 10, you would specify: ```--qual_filt --chopper_options '-q 10 -l 1000 --maxlength 20000'```.  **Based on our benchmarking, we recommend using the following parameters ```--chopper_options '-q 8 -l 100'``` as a first pass**.  

  If you are analysing samples that are of poor quality (i.e. failed the QC_FLAG) or amplifying a very short amplicon (e.g. <150 bp), then we recommend using the following setting ```--chopper_options '-q 8 -l 25'``` to retain reads of all lengths.  

A zipped copy of the resulting **preprocessed** and/or **quality filtered fastq file** will be saved in the preprocessing folder.  

If you trim raw read of adapters and/or quality filter the raw reads, an additional quality control step will be performed and a qc report will be generated summarising the read counts recovered before and after preprocessing for all samples listed in the index.csv file.

A qc report will be generated in text and html formats summarising the read counts recovered after the pre-processing step.  

If the user wants to check the data after preprocessing before performing downstream analysis, they can apply the parameter ``--preprocessing_only``.

### Clustering step (RATTLE)

In the clustering mode, the tool [`RATTLE`](https://github.com/comprna/RATTLE#Description-of-clustering-parameters) will be run. 

- The ont_amplicon pipeline will automatically set a **lower read length** of **100** bp during the RATTLE clustering step if the amplicon target_size specified in the csv file is **<=300 bp**.  
- If the amplicon target_size specified in the csv file is **>300 bp**, the default lower read length of **150 bp** will be applied at the RATTLE clustering step instead.  
- For poor quality samples (i.e. failed the QC_FLAG) or if your amplicon is known to be shorter than 150 bp, use the parameter ```--rattle_raw``` to use all the reads without any length filtering during the RATTLE clustering step.  
- Finally, the ``rattle_clustering_max_variance`` is set by default to 10000. It is recommended to drop it to 10 if analysing fastq files that were generated using a **fast** basecalling model.  

  **Special usage:**
  The parameters ``--rattle_clustering_min_length [number]``` (by default: 150) and ```--rattle_clustering_max_length [number]``` (by default: 100,000) can also be specified on the command line to restrict more strictly read size.  
  Additional parameters (other than raw, lower-length, upper-length and max-variance) can be set using the parameter ```--rattle_clustering_options '[additional paramater]'```.  

Example in which all reads will be retained during the clustering step:  
```
nextflow run maelyg/ont_amplicon -resume -profile singularity \
                            --analysis_mode clustering \
                            --adapter_trimming \
                            --qual_filt \
                            --chopper_options '-q 8 -l 25' \
                            --rattle_raw \
                            --blast_threads 2 \
                            --blastn_db /path/to/ncbi_blast_db/nt
```

Example in which reads are first quality filtered using the tool chopper (only reads with a Phread average quality score above 10 are retained). Then for the clustering step, only reads ranging between 500 and 2000 bp will be retained:  
```
nextflow run maelyg/ont_amplicon -resume -profile singularity \
                            --qual_filt \
                            --chopper_options chopper_options = '-q 8 -l 100' \
                            --analysis_mode clustering \
                            --rattle_clustering_min_length 200 \
                            --rattle_clustering_max_length 2000 \
                            --blast_threads 2 \
                            --blastn_db /path/to/ncbi_blast_db/nt
```

### Polishing step (optional)
The clusters derived using RATTLE can be polished. The reads are first mapped back to the clusters using Minimap2 and then the clusters are polished using Racon, Medaka2 and Samtools consensus. 
This step is performed by default by the pipeline but can be skipped by specifying the paramater ``--polishing false``.  

### Primer search
If the fwd_primer and the rev_primer have been provided in the csv file, clusters are then searched for primers using Cutadapt.  

### Blast homology search against NCBI
A preliminary megablast homology search against a Cytochrome oxidase I (COI) database will be performed if the gene targetted is COI and based on the strandedness of the consensus in the blast results, they will be reverse complemented where required.  
Blast homology search of the consensuses are then performed.  
A first step 
All the top hits derived for each contig are listed under the file **SampleName_assembly_blastn_top_hits.txt**. This file contains the following 26 columns:
```
- qseqid
- sgi
- sacc
- length
- pident
- mismatch
- gapopen
- qstart
- qend
- qlen
- sstart
- send
- slen
- sstrand
- evalue
- bitscore
- qcovhsp
- stitle
- staxids
- qseq
- sseq
- sseqid
- qcovs
- qframe
- sframe
- species


