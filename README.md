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
  Create a comma separated file that will be the input for the workflow. By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the ```--samplesheet [filename]``` in the nextflow run command. This text file requires the following columns (which need to be included as a header): ```sampleid,sample_files,spp_targets,gene_targets,target_size,fwd_primer,rev_primer``` 

  **sampleid** will be the sample name that will be given to the files created by the pipeline (required).  
  **sample_path** is the full path to the fastq files that the pipeline requires as starting input (required).  
  **spp_targets** is the organism targetted by the PCR (required).
  **gene_targets** is the gene targetted by the PCR (optional).
  **target_size** is the expected size of the amplicon (required).
  **fwd_primer** is the nucleotide sequence of the FWD primer (optional).
  **rev_primer** is the nucleotide sequence of the REV primer (optional).

  This is an example of an index.csv file which specifies the name and path of fastq.gz files for 2 samples. Specify the full path length for samples with a single fastq.gz file. If there are multiple fastq.gz files per sample, place them all in a single folder and the path can be specified on one line using an asterisk:
  ```
  sampleid,sample_files
  MT212,/path_to_fastq_file_folder/*fastq.gz
  MT213,/path_to_fastq_file_folder/*fastq.gz
  ```

- Specify a profile:
  ```
  nextflow run eresearchqut/ontvisc -profile {singularity, docker} --samplesheet index_example.csv
  ```
  setting the profile parameter to one of ```docker``` or ```singularity``` to suit your environment.
  
- Specify one analysis mode: ```--analysis_mode {read_classification, clustering, denovo_assembly, map2ref}``` (see below for more details)

- To set additional parameters, you can either include these in your nextflow run command:
  ```
  nextflow run eresearchqut/ontvisc -profile {singularity, docker} --samplesheet index_example.csv --adapter_trimming
  ```
  or set them to true in the nextflow.config file.
  ```
  params {
    adapter_trimming = true
  }
  ```

- A test is provided to check if the pipeline was successfully installed. The test.fastq.gz file is derived from of a plant infected with Miscanthus sinensis mosaic virus. To use the test, run the following command, selecting the adequate profile (singularity/docker):
  ```
  nextflow run eresearchqut/ontvisc -profile test,{singularity, docker}
  ```
  If this command does not run properly, you might need to also specify thte test.config file in the command line:  
  ```
  nextflow -c conf/test.config run eresearchqut/ontvisc -profile test,{singularity, docker}
  ```
  The test requires 2 cpus at least 16Gb of memory to run and can be executed locally.  
