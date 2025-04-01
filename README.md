# ont_amplicon
## Introduction
## Pipeline overview
## Installation
## Running the pipeline  
## Output files
## Authors



Introduction

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
  - Polishing (Minimap, Racon, Medaka2 - optional)
  - Remove adapters if provided (Cutadapt)
  - Megablast homology search against COI database (if COI is targetted) and reverse complement where required
  - Megablast homology search against NCBI database
  - Derive top candidate hits, assign preliminary taxonomy and target organism flag(pytaxonkit)
  - Map reads back to segment of consensus sequence that align to reference and derive BAM file and alignment statistics (minimap, samtools and mosdepth)
  - Map reads to segment of NCBI reference sequence that align to consensus and derive BAM file and consensus (minimap, samtools)


## Installation
### Requirements  
1. Install Java if not already on your system. Follow these instructions on this page [`page`](https://www.nextflow.io/docs/latest/getstarted.html#installation).

2. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

3. Install [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) to suit your environment. The pipeline has been validated using singularity version 3.10.2-1 and apptainer version 1.3.6-1.el9 but has not yet been tested with version 4.

3. Install taxonkit using the script install_taxonkit.sh or follow the steps described on this page [`page`](https://bioinf.shenwei.me/taxonkit/download/).


Download a local copy of the NCBI NT database, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:
  ```
  mkdir blastDB/20230930
  ```
  You will need to use a current update_blastdb.pl script from the blast+ version used with the pipeline (ie 2.16.0).
  For example:
  ```
  perl update_blastdb.pl --decompress nt
  perl update_blastdb.pl taxdb
  tar -xzf taxdb.tar.gz
  ```

  Specify the path of your local NCBI blast nt directories in the nextflow.config file or your enxtflow command.
  For instance:
  ```
  params {
    --blastn_db = '/work/hia_mt18005_db/blastDB/20230930/nt'
  }
  ```