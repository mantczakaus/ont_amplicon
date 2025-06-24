# ont_amplicon

## Introduction

ont_amplicon is a Nextflow-based bioinformatics pipeline designed to derive consensus sequences from **amplicon sequencing data** that were generated using **rapid library preparation kit** from **Oxford nanopore Technologies**. The pipeline expects the fastq files to have been generated using a **high accuracy (HAC)** basecalling model.  **The pipeline will fail to run if the fastq files provided are generated with a Fast basecalling model.**

It takes compressed fastq files (i.e. fastq.gz) as input.

## Table of contents
1. [Pipeline overview](#pipeline-overview)  
2. [Installation](#installation)  
a. [Requirements](#requirements)  
3. [Running the pipeline](#running-the-pipeline)  
a. [Quick start](#quick-start)  
b. [Run the pipeline for the first time](#run-the-pipeline-for-the-first-time)  
c. [Run test data](#run-test-data)  
d. [QC and preprocessing steps](#qc-and-preprocessing-steps)  
e. [Subsample reads](#subsample-reads)  
e. [Clustering step (RATTLE)](#clustering-step-rattle)  
f. [Polishing step (optional)](#polishing-step-optional)  
g. [Primer search](#primer-search)  
h. [Blast homology search against NCBI](#blast-homology-search-against-NCBI)  
i. [Mapping back to consensus](#mapping-back-to-consensus)  
j. [Mapping back to reference (optional)](#mapping-back-to-reference-optional)  
k. [HTML report](#HTML-report)  
5. [Output files](#output-files)  
a. [Nextflow reports](#nextflow-reports)  
b. [Preprocessing and quality check outputs](#preprocessing-and-quality-check-outputs)  
c. [Clustering step outputs](#clustering-step-outputs)  
d. [Polishing step outputs](#polishing-step-outputs)  
e. [Blast search outputs](#blast-search-outputs)  
f. [Outputs from mapping reads back to consensus matches step](#outputs-from-mapping-reads-back-to-consensus-matches-step)  
g. [Outputs from mapping reads back to reference matches step](#outputs-from-mapping-reads-back-to-reference-matches-step)  
h. [HTML report output](#html-report-output)  
7. [Authors](#authors) 

## Pipeline overview

<p><img src="docs/images/ont_amplicon_workflow.png" width="625"></p>

- Data quality check (QC) and preprocessing
  - Merge fastq.gz files ([Fascat](https://github.com/epi2me-labs/fastcat)) - optional
  - Quality check of raw fastq file ([NanoPlot](https://github.com/wdecoster/NanoPlot))
  - Trim adaptors ([PoreChop ABI](https://github.com/bonsai-team/Porechop_ABI))
  - Filter reads based on length and/or mean quality ([Chopper](https://github.com/wdecoster/chopper)) - optional
  - Reformat fastq files so read names are trimmed after the first whitespace ([bbmap](https://github.com/BioInfoTools/BBMap))
  - Quality check of processed fastq file ([NanoPlot](https://github.com/wdecoster/NanoPlot))
  - Subsample reads ([Seqkit](https://bioinf.shenwei.me/seqkit/usage/)) - optional
- QC report
  - Derive read counts recovered pre and post data processing
- Clustering mode
  - Read clustering ([Rattle](https://github.com/comprna/RATTLE))
  - Convert fastq to fasta format ([seqtk](https://github.com/lh3/seqtk))
  - Polishing ([Minimap2](https://lh3.github.io/minimap2/minimap2.html), [Racon](https://github.com/lbcb-sci/racon), [Medaka2](https://github.com/nanoporetech/medaka), [Samtools](http://www.htslib.org/doc/samtools.html)) - optional
  - Remove adapters, if provided (Cutadapt](https://cutadapt.readthedocs.io/en/stable/reference.html))
  - Megablast homology search against COI database (if COI is targetted) and reverse complement where required ([Blast+](https://www.ncbi.nlm.nih.gov/books/NBK279690/))
  - Megablast homology search against NCBI database ([Blast+](https://www.ncbi.nlm.nih.gov/books/NBK279690/))
  - Derive top candidate hits, assign preliminary taxonomy and set target organism flag ([pytaxonkit](https://github.com/bioforensics/pytaxonkit))  
  - Map reads back to segment of consensus sequence that aligns to reference and derive BAM file and alignment statistics ([Minimap2](https://lh3.github.io/minimap2/minimap2.html), [Samtools](http://www.htslib.org/doc/samtools.html) and [Mosdepth)](https://github.com/brentp/mosdepth))  
  - Map reads to segment of NCBI reference sequence that aligns to consensus and derive BAM file and consensus ([Minimap2](https://lh3.github.io/minimap2/minimap2.html), [Samtools](http://www.htslib.org/doc/samtools.html)) - optional


## Installation
### Requirements  

Make sure you are using Windows subsystem for Linux if you are on a Windows machine. Please follow these [`steps`](https://www.ssl.com/how-to/enable-linux-subsystem-install-ubuntu-windows-10/).

If the pipeline is run on a local machine, it will require between 300-800Gb of space for the installation of required containers and databases alone. This includes:  
- 80 Mb ont_amplicon pipeline  
- ~ 3.8Gb for containers  
- 600Mb for taxonkit databases  
- 280Gb/760Gb for the blast NCBI database coreNT/NT
- 3.4Gb for the MetaCOXI database  

To run, the pipeline will also require at least 2 cores and ~40Gb of memory per sample.  
The pipeline will generate ~5-100Mb of files per sample, depending on the number of consensuses recovered per sample and if mapping back to reference is required. Make sure you have enough space available on your local machine before running several samples at the same time.  

**1. Install Java** if not already on your system. Follow the java download instructions provided on this [`page`](https://www.nextflow.io/docs/latest/getstarted.html#installation).

**2. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)**.

  Nextflow memory requirements:

  In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in ~/.bashrc or ~./bash_profile): 
  ```
  NXF_OPTS='-Xms1g -Xmx4g'
  ```

**3. Install [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps)** to suit your environment. The pipeline has been validated using singularity version 3.10.2-1 and apptainer version 1.3.6-1.el9 but has not yet been tested with singularity version 4.

**4. Important: To avoid nextflow installing the singularity containers each time you are running the pipeline, set your singularity container cache directory.** For example, the command below will set a cache directory called $HOME/.nextflow/NXF_SINGULARITY_CACHEDIR and this directs Nextflow to always save the containers into this centralised location:  
```
[[ -d $HOME/.nextflow ]] || mkdir -p $HOME/.nextflow

mkdir $HOME/.nextflow/NXF_SINGULARITY_CACHEDIR

cat <<EOF > $HOME/.nextflow/config
singularity {
    cacheDir = '$HOME/.nextflow/NXF_SINGULARITY_CACHEDIR'
    autoMounts = true
}
```
Specify a different `cacheDir` location if space is limited in your home directory.  


**5. Install the taxonkit databases** using the script install_taxonkit.sh located in the bin folder or follow the steps described on this [`page`](https://bioinf.shenwei.me/taxonkit/download/).


**6. Install NCBI NT or coreNT.**  
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
  
  Specify the full path of your local NCBI blast name in your parameter file; **the path should include nt/core_nt**.
  For instance:
  ```
  blastn_db: /full/path/to/blastDB/20230930/nt
  ```

**7. Download the Cytochrome oxydase 1 (COI1) database** if you are planning to analyse COI samples.
  ```
  wget https://zenodo.org/record/6246634/files/MetaCOXI_Seqs_1.tar.gz
  #extract MetaCOXI_Seqs.fasta from the MetaCOXI_Seqs.tar.gz file
  tar -xvzf MetaCOXI_Seqs.tar.gz
  #make a blast database from the fasta file
  singularity exec https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4 makeblastdb -in MetaCOXI_Seqs.fasta -parse_seqids -dbtype nucl
  ```
Specify the full path to your COI database name in your in your parameter file.
  For instance:
  ```
  blastn_COI: /full/path/to/MetaCOXI_Seqs.fasta
  ```

## Running the pipeline  

### Quick start
A typical command for running the pipeline is as follows:
```
nextflow run main.nf -profile singularity -params-file params/params_example.yml
```
With the following parmaters specified in the params_sample.yml. **Please update paths to databases to match your local set up**:
```
{
samplesheet: /full/path/to/index.csv
merge: true
qual_filt: true
chopper_options: -q 8 -l 100
polishing: true
blastn_db: /full/path/to/NCBI/core_nt
blastn_COI: /full/path/to/MetaCOXI_Seqs.fasta
blast_threads: 2
taxdump: ~/.taxonkit
analyst_name: Maely Gauthier
facility: QUT
mapping_back_to_ref: true
}
```

And below is an example of an **index.file**:
```
sampleid,fastq_path,target_organism,target_gene,target_size,fwd_primer,rev_primer
VE24-1279_COI,/work/tests/mtdt_data/barcode01_VE24-1279_COI/*fastq.gz,drosophilidae,COI,711,GGTCAACAAATCATAAAGATATTGG,ATTTTTTGGTCACCCTGAAGTTTA
MP24-1051A_16S,/work/tests/mtdt_data/barcode06_MP24-1051A_16S/*fastq.gz,bacteria,16s,1509,AGAGTTTGATCATGGCTCAG,AAGTCGTAACAAGGTAACCGT
MP24-1096B_gyrB,/work/tests/mtdt_data/barcode19_MP24-1096B_gyrB/*fastq.gz,bacteria,gyrB,1258,GAAGTCATCATGACCGTTCTGCAYGCNGGNGGNAARTTYGA,ATGACNGAYGCNGAYGTNGAYGGCTCGCACATCCGTACCCTGCT
```

### Run the pipeline for the first time
-  To download the pipeline, you will need your github username and password.  
  
-  The source code can also be downloaded directly from GitHub using the git command:
  ```
  git clone https://github.com/maelyg/ont_amplicon.git
  ```
  Once you have downloaded the pipeline, you can either move into the pipeline **ont_amplicon** directory or direct Nextflow to the pipeline code by providing the full path to the **main.nf** file within the **ont_amplicon**. **Please note that you can only run one Nextflow analysis at a time from a given folder.**  

- Provide an index.csv file.  
  Create a **comma separated file (csv)** that will be the input for the workflow. 
  
  **Please note**: it is best to edit the csv file with an editor that does not add special characters/symbols (e.g. VSCode or Atom). If using other editors, check your files and if necessary, run [`dos2unix`](https://www.linuxfromscratch.org/blfs/view/git/general/dos2unix.html) on the file to get rid of these unwanted characters/symbols as they will cause the pipeline to fail.  
  
  By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the ```--samplesheet [filename]``` in the parameter file, as long as it has a **.csv** suffix. This text file requires the following columns (which need to be included as a header): ```sampleid,fastq_path,target_organism,target_gene,target_size,fwd_primer,rev_primer``` 

   - **sampleid** will be the sample name that will be given to the files created by the pipeline **(required)**.  **Please do not use space or special characters when naming your samples (e.g. &()*|\/) as this will break the code.**  
   - **fastq_path** is the full path to the fastq.gz files that the pipeline requires as starting input **(required)**.  
   - **target_organism** is the organism targetted by the PCR **(required)**. If more than one organism is targetted, use a pipe '|' as separator between the targets with no space.
   - **target_gene** is the gene targetted by the PCR **(required)**.  
   - **target_size** is the expected size of the amplicon **(required)**.  
   - **fwd_primer** is the nucleotide sequence of the FWD primer (optional).  
   - **rev_primer** is the nucleotide sequence of the REV primer (optional). Please note that the reverse primer has to be reverse complemented so it reads in the 5-3 direction.  
   - **test** is a field specific to NAQS LIMS (optional)
   - **method** is a field specific to NAQS LIMS (optional)

  For the **fastq files path**, the pipeline can process 1) multiple fastq.gz files per sample located within one folder (default) or 2) a single fastq.gz file per sample.  
  If there are **multiple fastq.gz files per sample**, their full path can be specified on one line using **an asterisk (i.e. *fastq.gz)** and you will need to specify the parameter ```--merge``` in the parameter file (default setting).  
  See an example of an index.csv file for 3 MTDT samples:  
  ```
  sampleid,fastq_path,target_organism,target_gene,target_size
  VE24-1279_COI,/work/tests/mtdt_data/barcode01_VE24-1279_COI/*fastq.gz,drosophilidae,COI,711,GGTCAACAAATCATAAAGATATTGG,ATTTTTTGGTCACCCTGAAGTTTA
  MP24-1051A_16S,/work/tests/mtdt_data/barcode06_MP24-1051A_16S/*fastq.gz,bacteria,16s,1509,AGAGTTTGATCATGGCTCAG,AAGTCGTAACAAGGTAACCGT
  MP24-1096B_gyrB,/work/tests/mtdt_data/barcode19_MP24-1096B_gyrB/*fastq.gz,bacteria,gyrB,1258,GAAGTCATCATGACCGTTCTGCAYGCNGGNGGNAARTTYGA,ATGACNGAYGCNGAYGTNGAYGGCTCGCACATCCGTACCCTGCT
  ```
  For samples with **a single fastq.gz** file, specify **the full path to the fastq.gz file**  and set the merge parameter as **false**.  

  A python script is provided in the **bin** folder called derive_sample_sheet.py that will create an index.csv file with header and automatically populate the sample name and the fastq.gz file path.
  It requires python to be installed on your local system.  
  It expects fastq files to be in gzip format (i.e. fastq.gz).
  
  To run the script, specify the directory where the fastq.gz files are located with **-d [full/path/to/fastqgz/directory** and the output file name with  **-o index.csv**:  
  ```
  python derive_sample_sheet.py -d /full/path/to/fastqgz/directory -o index_example.csv
  ```
  If there is only a single fastq.gz file per sample, specify the single option with **-s**:  
  ```
  python derive_sample_sheet.py -d /full/path/to/fastqgz/directory -o index_example.csv -s
  ```
  **Please note that the script will derive sample names based on the fastq.gz file prefix. Adjust the sample name accordingly.**
  
- Specify your container engine ```singularity``` as the profile:
  ```
  nextflow run main.nf -profile singularity
  ```
  
- To specify parameters, we recommend that you provide a **params.yml** file with all your parameters. Please refer to this document for additional information on how to pass parameters and how parameter passing priority works (https://software.pixelgen.com/nf-core-pixelator/1.3.x/usage/passing-parameters/).  

```
nextflow run main.nf  -profile singularity -params-file params/params_example.yml
```
The default parameters are provided in the default_params.yml under the params folder:
```
{
samplesheet: index.csv
merge: true
qual_filt: true
chopper_options: -q 8
polishing: true
blastn_db: null
taxdump: null
blast_threads: 2
analyst_name: null
facility: null
mapping_back_to_ref: true
outdir: results
help: false
qc_only: false
preprocessing_only: false
porechop_options: null
porechop_custom_primers: false
porechop_custom_primers_path: ~/ont_amplicon/bin/adapters.txt
chopper_options: null
analysis_mode: clustering
rattle_clustering_options: null
rattle_clustering_min_length: null 
rattle_clustering_max_length: null
rattle_raw: false
rattle_clustering_max_variance: 1000000
rattle_polishing_options: null
blastn_COI: null
subsample: false
reads_downsampling_size: 10000
}
```
- Specify the full path to your blast and taxonkit databases in your parameter file.  The analysis cannot proceed without these being set.

- Specify the ``--analyst_name`` and the ``--facility`` in your parameter file.  The analysis cannot proceed without these being set.

- Specify the full path to your COI database name in your in your parameter file if some of your gene targets are COI. The analysis will not process if one of your target gene is COI and the blastn_COI parameter is not set.   

- If you are invoking the pipeline from another folder, create a config file (for example local.config) in which you specify the full path of your local ont_amplicon repository:  
```
singularity {
  runOptions = '-B /full/path/to/ont_amplicon:/mnt/ont_amplicon'
  }
```
And update your command to:  
```
nextflow run /full/path/to/main.nf  -profile singularity -resume  -params-file  params.yml -c local.config
```

- By default, nextflow saves the run log to a file called **.nextflow.log** in the folder from which the analysis is run. Add the **`-log` option** to your nextflow command to specify a different log file name and location.  

### Run test data
Two tests are currently provided to check if the pipeline was successfully installed. 
- The mtdt_test runs three samples provided by  MTDT (barcode01_VE24-1279_COI, barcode06_MP24-1051A_16S and barcode19_MP24-1096B_gyrB). 
- The peq_test runs two samples provided by PEQ (ONT141 and ONT142). 
A small NCBI blast database and COI database have been derived to speed up the analysis run in test mode.  

To use the tests, change directory to the ont_amplicon github repository and run the following command for the MTDT test:
  ```
  nextflow run main.nf -profile singularity,mtdt_test -resume -params-file params/params_mtdt_test.yml
  ```
  and this command for the PEQ test:  
  ```
  nextflow run main.nf -profile singularity,peq_test -resume -params-file params/params_peq_test.yml
  ```

Please note that you will need to specify the location of the taxonkit databases folder in your command if it is not located where expected, at `~/.taxonkit` using the parameter `--taxdump path/to/taxonkit_db_folder`. For example:
  ```
  nextflow run main.nf -profile singularity,peq_test -resume -params-file params/params_peq_test.yml --taxdump /full/path/to/taxonkit_db/folder 
  ```

The tests should take less than 5 minutes to run to completion.

If the installation is successful, your screen should output something similar to this, with a completion and success status at the bottom:
```
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [hungry_hilbert] DSL2 - revision: cea782866c

executor >  pbspro (76)
[8a/42602e] TIMESTAMP_START                                     [100%] 1 of 1 ✔
[f9/ad32cb] FASTCAT (barcode19_MP24-1096B_gyrB)                 [100%] 3 of 3 ✔
[76/8c97a7] QC_PRE_DATA_PROCESSING (barcode19_MP24-1096B_gyrB)  [100%] 3 of 3 ✔
[e7/cfb6bf] PORECHOP_ABI (barcode19_MP24-1096B_gyrB)            [100%] 3 of 3 ✔
[e7/f281b9] CHOPPER (barcode19_MP24-1096B_gyrB)                 [100%] 3 of 3 ✔
[f8/987357] REFORMAT (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[d4/4333ad] QC_POST_DATA_PROCESSING (barcode19_MP24-1096B_gyrB) [100%] 3 of 3 ✔
[1b/de3961] QCREPORT                                            [100%] 1 of 1 ✔
[ae/563754] RATTLE (barcode19_MP24-1096B_gyrB)                  [100%] 3 of 3 ✔
[d2/482ec5] CLUSTER2FASTA (barcode19_MP24-1096B_gyrB)           [100%] 3 of 3 ✔
[61/6cc7fe] MINIMAP2_RACON (barcode19_MP24-1096B_gyrB)          [100%] 3 of 3 ✔
[d1/5025de] RACON (barcode19_MP24-1096B_gyrB)                   [100%] 3 of 3 ✔
[fa/4d7987] MEDAKA2 (barcode19_MP24-1096B_gyrB)                 [100%] 3 of 3 ✔
[80/d7647f] CUTADAPT (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[74/f9c567] BLASTN_COI (barcode01_VE24-1279_COI)                [100%] 1 of 1 ✔
[ed/aa205d] REVCOMP (barcode01_VE24-1279_COI)                   [100%] 1 of 1 ✔
[e0/d44563] BLASTN (barcode01_VE24-1279_COI)                    [100%] 1 of 1 ✔
[73/541c3f] BLASTN2 (barcode19_MP24-1096B_gyrB)                 [100%] 2 of 2 ✔
[13/fc6720] EXTRACT_BLAST_HITS (barcode06_MP24-1051A_16S)       [100%] 3 of 3 ✔
[b3/13c5ea] FASTA2TABLE (barcode19_MP24-1096B_gyrB)             [100%] 3 of 3 ✔
[e0/2e5e64] MINIMAP2_CONSENSUS (barcode06_MP24-1051A_16S)       [100%] 3 of 3 ✔
[78/2822af] SAMTOOLS_CONSENSUS (barcode06_MP24-1051A_16S)       [100%] 3 of 3 ✔
[5c/9e082b] PYFAIDX (barcode06_MP24-1051A_16S)                  [100%] 3 of 3 ✔
[34/ee513d] MOSDEPTH (barcode06_MP24-1051A_16S)                 [100%] 3 of 3 ✔
[e7/b555f5] SEQTK (barcode06_MP24-1051A_16S)                    [100%] 3 of 3 ✔
[7d/cb9d44] COVSTATS (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
[70/c85da3] HTML_REPORT (2)                                     [100%] 3 of 3 ✔
[96/3ff7d4] MINIMAP2_REF (barcode19_MP24-1096B_gyrB)            [100%] 3 of 3 ✔
[29/f193f6] SAMTOOLS (barcode19_MP24-1096B_gyrB)                [100%] 3 of 3 ✔
Completed at: 23-Jun-2025 11:14:02
Duration    : 6m 18s
CPU hours   : 0.3
Succeeded   : 76
```

By default, the output files will be saved under the **results** folder (this can be changed by setting the **`--outdir`** parameter to something else).  

The results folder has the following structure:  
```
├── 00_QC_report
│   ├── run_qc_report_20250428-091847.html
│   └── run_qc_report_20250428-091847.txt
├── 01_pipeline_info
│   ├── 20250428091444_nextflow_start_timestamp.txt
│   ├── 20250428092534_nextflow_start_timestamp.txt
│   ├── execution_report_2025-04-28_09-14-26.html
│   ├── execution_timeline_2025-04-28_09-14-26.html
│   ├── execution_trace_2025-04-28_09-14-26.txt
│   ├── execution_trace_2025-04-28_09-25-22.txt
│   └── pipeline_dag_2025-04-28_09-14-26.html
├── barcode01_VE24-1279_COI
│   ├── 00_preprocessing
│   │   ├── barcode01_VE24-1279_COI_preprocessed.fastq.gz
│   │   ├── chopper
│   │   │   └── barcode01_VE24-1279_COI_chopper.log
│   │   └── porechop
│   │       └── barcode01_VE24-1279_COI_porechop.log
│   ├── 01_QC
│   │   ├── fastcat
│   │   │   ├── barcode01_VE24-1279_COI_basecalling_model_inference.txt
│   │   │   └── barcode01_VE24-1279_COI.fastq.gz
│   │   └── nanoplot
│   │       ├── barcode01_VE24-1279_COI_filtered_LengthvsQualityScatterPlot_dot.html
│   │       ├── barcode01_VE24-1279_COI_filtered_NanoPlot-report.html
│   │       ├── barcode01_VE24-1279_COI_filtered_NanoStats.txt
│   │       ├── barcode01_VE24-1279_COI_raw_LengthvsQualityScatterPlot_dot.html
│   │       ├── barcode01_VE24-1279_COI_raw_NanoPlot-report.html
│   │       └── barcode01_VE24-1279_COI_raw_NanoStats.txt
│   ├── 02_clustering
│   │   ├── barcode01_VE24-1279_COI_rattle.fasta
│   │   └── barcode01_VE24-1279_COI_rattle.log
│   │   └── barcode01_VE24-1279_COI_rattle_status.txt
│   ├── 03_polishing
│   │   ├── barcode01_VE24-1279_COI_cutadapt.log
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_medaka_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_medaka.log
│   │   ├── barcode01_VE24-1279_COI_racon_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_racon.log
│   │   ├── barcode01_VE24-1279_COI_samtools_consensus.fasta
│   │   ├── barcode01_VE24-1279_COI_samtools_consensus.fastq
│   │   └── barcode01_VE24-1279_COI_samtools_consensus.log
│   ├── 04_megablast
│   │   ├── barcode01_VE24-1279_COI_blast_status.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc.fasta
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_10_hits.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_rc_megablast_top_hits.txt
│   ├── 05_mapping_to_consensus
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam
│   │   ├── barcode01_VE24-1279_COI_aln.sorted.bam.bai
│   │   ├── barcode01_VE24-1279_COI_coverage.txt
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_match.fasta
│   │   ├── barcode01_VE24-1279_COI_final_polished_consensus_match.fastq
│   │   └── barcode01_VE24-1279_COI_top_blast_with_cov_stats.txt
│   ├── 06_mapping_to_ref
│   │   ├── barcode01_VE24-1279_COI_coverage.txt
│   │   ├── barcode01_VE24-1279_COI_ref_aln.sorted.bam
│   │   ├── barcode01_VE24-1279_COI_ref_aln.sorted.bam.bai
│   │   ├── barcode01_VE24-1279_COI_reference_match.fasta
│   │   └── barcode01_VE24-1279_COI_samtools_consensus_from_ref.fasta
│   └── 07_html_report
│       ├── barcode01_VE24-1279_COI_bam-alignment.html
│       ├── barcode01_VE24-1279_COI_report.html
│       ├── default_params.yml
│       ├── example_report_context.json
│       ├── run_qc_report.html
│       └── versions.yml
```

### QC and preprocessing steps
By defaut, the pipeline expects that multiple fastq files exist for a single sample and the `merge: true` option is set. These fastq files will be merged using [`Fascat`](https://github.com/epi2me-labs/fastcat).  
Then the read names of the fastq file created will be trimmed after the first whitespace in the read header, for compatiblity purposes with all downstream tools.  

The pipeline will run a quality control check of the raw reads using [NanoPlot](https://github.com/wdecoster/NanoPlot).  

**It is recommended to first run only the quality control step to have a preliminary look at the data before proceeding with downstream analyses by specifying the `qc_only: true` parameter.**

Reads are then trimmed of adapters and optionally quality filtered:  
- Reads are searched for the presence of sequencing adapters using [`Porechop ABI`](https://github.com/rrwick/Porechop). Other porechop ABI parameters can be specified using ```porechop_options: '{options} '```, making sure you leave a space at the end before the closing quote. Please refer to the Porechop manual.  

  **Special usage:**  
  By default, Porechop limits the search to known adapters listed in [`adapter.py`](https://github.com/bonsai-team/Porechop_ABI/blob/master/porechop_abi/adapters.py).  
  To search ab initio for adapters on top of known adapters, specify:
  ```
  porechop_options: '-abi '
  ```  
  If you want instead to limit the search to custom adapters, specify: 
  ```
  porechop_custom_primers: true
  porechop_options '-ddb '
  ```
  and list the custom adapters in the text file located under bin/adapters.txt following the format:  
  ```
  line 1: Adapter name
  line 2: Start adapter sequence
  line 3: End adapter sequence
  --- repeat for each adapter pair---
  ```

- The user can perform a quality filtering step using [`Chopper`](https://github.com/wdecoster/chopper) by specifying  the ```qual_filt: true``` parameter. Chopper parameters to apply will need to be specified separately using the ```chopper_options: {options}```. Please refer to the Chopper manual.  
  For instance to filter reads shorter than 1000 bp and longer than 20000 bp, and reads with a minimum Phred average quality score of 10, you would specify in your parameter file: 
  ```
  qual_filt: true
  chopper_options: -q 10 -l 1000 --maxlength 20000
  ```
**Based on the benchmarking performed, we recommend using the following parameters ```chopper_options: -q 8 -l 100``` as a first pass**.  

  If you are analysing samples that are of poor quality (i.e. failed the QC_FLAG) or amplifying a very short amplicon (e.g. <150 bp), then we recommend using the following setting `qual_filt: true` and `chopper_options: -q 8 -l 25`, or skip the quality trimming step altogther (i.e. `qual_filt: false`) if you want to retain all reads.  

A zipped copy of the resulting **preprocessed** and/or **quality filtered fastq file** will be saved in the preprocessing folder.  

After processing raw reads, an additional quality control step will be performed.  

A **qc report** will be generated in text and html formats summarising the read counts recovered before and after the pre-processing step for all samples listed in the index.csv file.

If the user wants to check the data after preprocessing before performing downstream analysis, they can apply the parameter `preprocessing_only: true`.

### Subsample reads 
Pre-processed reads can be subsampled by specifying `subsample: true` and the number of reads to subsample are specified using the parameter `reads_downsampling_size: [number of reads]` (set to 10,000 by default).  

### Clustering step (RATTLE)

In the clustering mode, the tool [`RATTLE`](https://github.com/comprna/RATTLE#Description-of-clustering-parameters) will be run. 

- The ont_amplicon pipeline will automatically set a **lower read length** of **100** bp during the RATTLE clustering step if the amplicon target_size specified in the csv file is **<=300 bp**.  
- If the amplicon target_size specified in the csv file is **>300 bp**, the lower read length of **150 bp** (Rattle default) will be applied at the RATTLE clustering step instead.  
- For poor quality samples (i.e. failed the QC_FLAG) or if your amplicon is known to be shorter than 150 bp, use the parameter `rattle_raw: true` to use all the reads without any length filtering during the RATTLE clustering step.  
- Finally, the `rattle_clustering_max_variance` is set by default to 1,000,000.  

  **Special usage:**
  The parameters `rattle_clustering_min_length: [number]` (by default: 150) and `rattle_clustering_max_length: [number]` (by default: 100,000) can also be specified in the parameter file to restrict lower and upper read size length.  
  Additional parameters (other than raw, lower-length, upper-length and max-variance) can be set using the parameter `rattle_clustering_options: [additional paramater]`.  

Example of parameter file in which all reads will be retained during the quality filtering and clustering steps (polishing is also set to false):  
```
{
samplesheet: tests/index_mtdt.csv
merge: true
qual_filt: false
rattle_raw: true
polishing: false
blastn_db: $HOME/ont_amplicon/tests/blastdb/reference.fasta
blastn_COI: $HOME/ont_amplicon/tests/COIdb/MetaCOXI_Seqs.fasta
taxdump: ~/.taxonkit
blast_threads: 2
analyst_name: John Smith
facility: MTDT
mapping_back_to_ref: true
}
```

Example in which reads are first quality filtered using the tool Chopper (only reads with a Phread average quality score above 8 and length >100 bp are retained). Then for the clustering step, only reads ranging between 500 and 2000 bp are retained:  
```
{
samplesheet: tests/index_mtdt.csv
merge: true
qual_filt: true
chopper_options: -q 8 -l 100
rattle_clustering_min_length: 500
rattle_clustering_max_length: 2000
polishing: true
blastn_db: $HOME/ont_amplicon/tests/blastdb/reference.fasta
blastn_COI: $HOME/ont_amplicon/tests/COIdb/MetaCOXI_Seqs.fasta
taxdump: ~/.taxonkit
blast_threads: 2
analyst_name: John Smith
facility: MTDT
mapping_back_to_ref: true
}
```

### Polishing step (optional)
The clusters derived using RATTLE can be polished. The reads are first mapped back to the clusters using [Minimap2](https://lh3.github.io/minimap2/minimap2.html) and then the clusters are polished using [Racon](https://github.com/lbcb-sci/racon) and [Medaka2](https://github.com/nanoporetech/medaka). Samtools consensus is then used to identify positions showing poor base and mapping qualities, using the predefined sets of configuration parameters that have been optimised for ONT reads (i.e. r10.4_sup) (please see the configuration section at the bottom of the [`Samtools consensus documentation`](https://www.htslib.org/doc/samtools-consensus.html)). Any stretches of Ns at the 5' and 3' end of the consensuses are then removed with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/reference.html). If a polishing step fails, it will be skipped.  
This polishing step is performed by default by the pipeline but can be skipped by specifying the paramater ``--polishing false``.  

### Primer search
If the fwd_primer and the rev_primer have been provided in the samplesheet, clusters are then searched for primers using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/reference.html).  

### Blast homology searches
If the gene targetted is Cytochrome oxidase I (COI), a preliminary megablast homology search against a COI database will be performed; then based on the strandedness of the blast results for the consensuses , some will be reverse complemented where required.  

Blast homology search of the consensuses against NCBI is then performed and the top 10 hits are returned.
A separate blast output is then derived using [pytaxonkit](https://github.com/bioforensics/pytaxonkit), to output preliminary taxonomic assignment to the top blast hit for each consensus. The nucleotide sequence of qseq **(i.e. consensus match)** and sseq **(i.e. reference match)** are extracted to use when mapping reads back to consensus and reference respectively (see steps below).  

### Mapping back to consensus
The quality filtered reads derived during the pre-processing step are mapped back to the consensus matches using Mimimap2. Samtools and Mosdepth are then used to derive bam files and coverage statistics. A summary of the blast results, preliminary taxonomic assignment, coverage statistics and associated **flags** and **confidence scores** are then derived for each consensus using python.  

### Mapping back to reference (optional)
By default the quality filtered reads derived during the pre-processing step are also mapped back to the reference blast match and [Samtools consensus](http://www.htslib.org/doc/samtools-consensus.html) is used to derive independent guided-reference consensuses. Their nucleotide sequences can be compared to that of the original consensuses to resolve ambiguities (ie low complexity and repetitive regions).  

### HTML report
An html summary report is generated for each sample, incorporating sample metadata, QC before and after 
preprocessing, blast results and coverage statistics. It also provides a link to the bam files generated when mapping back to consensus.  

## Output files
The output files will be saved by default under the **results** folder. This can be changed by setting the **`--outdir` parameter**.  

### Nextflow reports
Nextflow generates several outputs which are stored under the **01_pipeline_info** folder. Please find detailed information about these on this [page](https://www.nextflow.io/docs/latest/reports.html). All of these have the date and time appended as a suffix.

#### - HTML execution report
Nextflow outputs an **HTML execution report** which includes general metrics about the run. The report is organised into 3 main sections:  
- The **Summary** section reports the execution status, the launch command, overall execution time and some other workflow metadata.  
- The **Resources** section plots the distribution of resource usage for each workflow process. Plots are shown for CPU, memory, job duration and disk I/O. They have two (or three) tabs with the raw values and a percentage representation showing what proportion of the requested resources were used. These plots are very helpful to check that task resources are used efficiently.  
- The **Tasks** section lists all executed tasks, reporting for each of them the status, the actual command script, and several other metrics.  

#### - Trace file
Nextflow creates an execution tracing text file that contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used.  

#### - Execution timeline
Nextflow can render an HTML timeline for all processes executed in your pipeline. Each bar represents a process run in the pipeline execution. The bar length represents the task duration time (wall-time). The colored area in each bar represents the real execution time. The grey area to the left of the colored area represents the task scheduling wait time. The grey area to the right of the colored area represents the task termination time (clean-up and file un-staging). The numbers on the x-axis represent the time in absolute units e.g. minutes, hours, etc.  

Each bar displays two numbers: the task duration time and the virtual memory size peak.  

As each process can spawn many tasks, colors are used to identify those tasks belonging to the same process.  

#### - Workflow diagram
The pipeline executed is represented as an HTML diagram in direct acyclic graph format. The vertices in the graph represent the pipeline’s processes and operators, while the edges represent the data dependencies (i.e. channels) between them.

### Preprocessing and quality check outputs  
The files are located under the **Sample_name/00_preprocessing**, the **Sample_name/01_QC** & the **00_QC_report** folders.  

If the fastq.gz files need to be merged for a sample, the resulting fastq.g file will be stored under the **SampleName/01_QC/fastcat** folder. The basecallign model used will also be captured in this folder in a file called **basecalling_model_inference.txt**.  
A quality check will be performed on the raw fastq file using [NanoPlot](https://github.com/wdecoster/NanoPlot) which is a tool that can be used to produce general quality metrics e.g. quality score distribution, read lengths and other general stats. A NanoPlot-report.html file will be saved under the **SampleName/01_QC/nanoplot** folder with the prefix **raw**. This report displays 6 plots as well as a table of summary statistics.  

<p align="center"><img src="docs/images/Example_Statistics.png" width="1000"></p>

Example of output plots:
<p align="center"><img src="docs/images/Example_raw_WeightedHistogramReadlength.png" width="750"></p>
<p align="center"><img src="docs/images/Example_LengthvsQualityScatterPlot.png" width="750"></p>

A preprocessed fastq file will be saved in the **SampleName/00_preprocessing** output directory which will have its read names trimmed after the first whitespace, for compatiblity purposes with all downstream tools. This fastq file will also be trimmed of adapters (PoreChopABI) and optionally filtered based on quality and/or length (if Chopper was run).  

After adapter trimming, a PoreChopABI log will be saved under the **SampleName/00_preprocessing/porechop** folder.  

After quality/length trimming, a Chopper log file will be saved under the **SampleName/preprocessing/chopper** folder.  

A second quality check will be performed on the processsed fastq file and a NanoPlot-report.html file will be saved under the **SampleName/01_QC/nanoplot** folder with the prefix **filtered**.  

A **QC report**, which captures the date and time in the file name, will be generated in text and html format (i.e. **run_qc_report_YYYYMMDD-HHMMSS.txt** and **run_qc_report_YYYYMMDD-HHMMSS.html**) under the **00_QC_report** folder.  It summarises the read counts recovered before and after the pre-processing step for all samples listed in the index.csv file.  
It will include 3 flags:  
1) For the **raw_reads_flag**, if there were < 2500 raw reads, the column will display: "Less than 2500 raw reads".  
2) For the **processed_reads_flag**, if there were < 200 processed_reads, the column will display: "Less than 200 processed reads".  
3) **QC_FLAG**:
- GREEN: >= 2500 starting reads, >= 200 processed_reads.
- ORANGE: < 2500 starting reads, >= 200 processed_reads.
- RED: < 2500 starting reads, < 200 processed_reads.

Example of report:

| Sample| raw_reads | processed_reads | percent_processed | raw_reads_flag | processed_reads_flag | QC_FLAG |
| --- | --- | --- | --- | --- | --- | --- |
| ONT141 | 10929 | 2338 | 21.39 | | | GREEN |
| ONT142| 21849 | 4232 | 9.37 | | | GREEN |

### Clustering step outputs
The files are located under the **Sample_name/02_clustering folder**.  
The output from Rattle will be saved under **SampleName/02_clustering/SampleName_rattle.fasta**. The number of reads contributing to each clusters is listed in the header. The amplicon of interest is usually amongst the most abundant clusters (i.e. the ones represented by the most reads). The rattle log (**SampleName_rattle.log**) is also available in the same folder as well as a file called **SampleName_rattle.status** that catches whether the clustering step ran succesfully or not.  

### Polishing step outputs 
The files are located under the **Sample_name/03_polishing folder**.  
By default, the clusters derived using RATTLE will be polished with Racon followed by Medaka2. The intermediate polished clusters and logs generated by Racon (**03_polishing/Samplename_racon_consensus.fasta** and **03_polishing/Samplename_racon.log**) and Medaka2 (**03_polishing/Samplename_medaka_consensus.fasta** and **03_polishing/Samplename_medaka.log**) are provided. Consensuses are then generated in fasta and fastq formats with Samtools consensus which will identify regions with poor base and mapping quality (**03_polishing/Samplename_samtools_consensus.fasta**, **Samplename_samtools_consensus.fastq** and **03_polishing/Samplename_samtools_consensus.log**). Cutadapt is then run to remove any Ns present at the start and end of consensuses. If amplicon primers were provided in the samplesheet, Cutadapt also searches for these at the 5' and 3' end of the consensus, respectively (see **03_polishing/Samplename_final_polished_consensus.fasta** and **03_polishing/Samplename_cutadapt.log**).  

### Blast search outputs
The files are located under the **Sample_name/04_megablast folder**.  
If the target gene is COI, then the consensuses will first be mapped to a cyctochrome oxidase I database and based on the strandedness of the blast results, consensuses will be reverse complemented where required. All consensuses will be saved in the **Sample_name/04_megablast/Sample_name_final_polished_consensus_rc.fasta** file.  
All consensuses are then blasted against NCBI.  The outcome of the blast search will be captured in the **Sample_name/04_megablast/Sample_name_final_blast_status.txt** file (e.g. if at least one consensus returned a blast hit, it will display 'passed', if no consensuses returned a blast hit, it will display 'failed').  

The 10 top hits derived for each contig are listed in the file **SampleName/04_megablast/SampleName_final_polished_consensus_megablast_top_10_hits.txt**. This file contains the following 26 columns:
```
- qseqid => query or source (gene) sequence id
- sgi => subject GI
- sacc => subject accession
- length => alignment length (sequence overlap)
- nident => number of identical matches
- pident => percentage of identical positions
- mismatch => number of mismatches
- gaps => total number of gaps
- gapopen => number of gap openings
- qstart => start of alignment in query
- qend => end of alignment in query
- qlen => query sequence length
- sstart => start of alignment in query
- send => end of alignment in query
- slen => subject sequence length
- sstrand => subject Strand
- evalue => expect value
- bitscore => bitscore
- qcovhsp => query Coverage Per HSP
- stitle => subject Title
- staxids => subject taxonomy ID(s), separated by a ';'
- qseq => aligned part of query sequence
- sseq => aligned part of subject sequence
- sseqid => subject or target (reference genome) sequence id
- qcovs => query Coverage Per Subject
- qframe => query frame
- sframe => subject frame
```

A separate blast output only retaining the top blast hit is created (**SampleName/04_megablast/SampleName_final_polished_consensus_megablast_top_hit.txt**) which contains four additional columns. Preliminary taxonomic assignment is derived to the top blast hit for each consensus using [pytaxonkit](https://github.com/bioforensics/pytaxonkit). A **broad_taxonomic_category** column is generated which matches the cluster to broad taxon categories:
- virus  
- bacteria;phytoplasma  
- bacteria;other  
- archea  
- eukaryota;fungi;powdery_mildew  
- eukaryota;fungi;other  
- eukaryota;deuterostomia  
- eukaryota;protostomia  
- eukaryota; other  
- other  

A column called **FullLineage** provides the full taxonomic lineage derived from pytaxonkit.  
The colum **target_organism_match** indicates whether there was a taxon match between the target_organism specified in the samplesheet and the full taxonomic lineage.  
Finally, the **n_read_cont_cluster** captures the number of reads that originally contributed to build the clusters during the RATTLE step.  

The nucleotide sequence of qseq (i.e. **consensus match**) and sseq (i.e. **reference match)** are extracted to use when mapping reads back to consensus and reference respectively (see steps below). These are called **SampleName/05_mapping_to_consensus/SampleName_final_polished_consensus_match.fasta** and **SampleName/megablast/06_mapping_to_ref/SampleName_reference_match.fasta** respectively.  

### Outputs from mapping reads back to consensus matches step  
The files are located under the **Sample_name/05_mapping_to_consensus**.  

A BAM file of the pre-processed reads mapped back to the consensus matches is generated (**Sample_name/05_mapping_to_consensus/Sample_name_aln.sorted.bam** and **Sample_name/05_mapping_to_consensus/Sample_name_aln.sorted.bam.bai**) and coverage statistics are derived. A final summary file is generated which combines the previously generated final_polished_consensus_megablast_top_hit.txt file with coverage statistics, flags and confidence scores for each consensus blast match. This file includes the following additional columns: 
- query_match_length: length of consensus match used as reference when mapping back pre-processed reads 
- qseq_mapping_read_count: number of reads mapping back to the consensus match  
- qseq_mean_depth: mean read coverage of each base when mapping back to the consensus match  
- qseq_pc_mapping_read: the percentage of processed reads that map to the consensus match  
- qseq_pc_cov_30X: the percentage of bases that attained at least 30X sequence coverage when mapping back to the consensus match  
- mean_MQ: average mapping quality of reads mapping to the consensus match  
- num_passing_90: number of mapped reads whose lengths are at least 90% of the consensus match length   
- 30X_COVERAGE_FLAG: see FLAGS section below  
- MAPPED_READ_COUNT_FLAG: see FLAGS section below   
- MEAN_COVERAGE_FLAG: see FLAGS section below    
- TARGET_ORGANISM_FLAG: see FLAGS section below  
- TARGET_SIZE_FLAG: see FLAGS section below    
- READ_LENGTH_FLAG: see FLAGS section below    
- MEAN_MQ_FLAG: see FLAGS section below    
- TOTAL_CONF_SCORE: a scoring system which assigns different weight to each flag colour for the 30X coverage flag, the target size, the mapped read count flag, the mean coverage flag, the read length flag and the mean MQ flag. It is a value bewteen 0 and 12. A higher score indicates a higher confidence in the quality of the consensus sequence  
- NORMALISED_CONF_SCORE: a value between 0 and 1 that is calculated by normalising the confidence score to the maximum possible score for this sequence. A value of 1 indicates the highest confidence in the quality of the consensus sequence  

#### FLAGS
**Seven** flags are providded to help with interpretation. They will display GREEN, ORANGE, RED or GREY depending on whether they fill specific criteria:

| FLAG NAME | DEFINITION | GREEN | ORANGE | RED | GREY |
| --- | --- | --- | --- | --- |  --- |
| **1. 30X COVERAGE FLAG** | The percentage of bases that attained at least 30X sequence coverage when mapping back to the consensus match (ie qseq) | **> 90%** |  **75-90%** | **< 75%** | The consensus returned no blast hits |
| **2. TARGET ORGANISM FLAG** | Flag based on whether the consensus matched to the target organism(s) by blast homology search and the % blast identity recovered | Target organism detected and blast identity **> 90%** | Target organism was detected and blast identity **< 90%** | Target organism not detected | The consensus returned no blast hits |
| **3. TARGET SIZE FLAG** | Length of the consensus match (ie qseq) relative to the expected target_size | **within ±20%** | **±20%-±40%** | **outside the range of ±40%** | The consensus returned no blast hits |
| **4. MAPPED READ COUNT FLAG** | Number of reads mapping back to the consensus match (ie qseq) | **>= 1000** |  **200-1000** | **< 200** | The consensus returned no blast hits |
| **5. MEAN COVERAGE FLAG** | Mean read coverage of each base when mapping back to the consensus match (ie qseq) | **>= 500** | **100-500** | **< 100** | The consensus returned no blast hits |
| **6. READ LENGTH FLAG** | Number of mapped reads whose lengths are at least 90% of the consensus match length |  **>=200** |  **50-200** | **< 50** | The consensus returned no blast hits |
| **7. MEAN MQ FLAG** | Average mapping quality of reads mapping to the consensus match | **>= 30** | **10-30** | **< 10** | The consensus returned no blast hits |

### Outputs from mapping reads back to reference matches step
By default the quality filtered reads derived during the pre-processing step are mapped back to the
reference blast match. A bam file is generated using Samtools and [Samtools consensus](https://www.htslib.org/doc/samtools-consensus.html) is used to derive independent guided-reference consensuses that are stored in a file called **SampleName/mapping_back_to_ref/samtools_consensus_from_ref.fasta** file. Their nucleotide sequences can be compared to that of the original consensuses to resolve ambiguities (ie low complexity and repetitive regions). 

### HTML report output
(in progress)  
An HTML report example can be found [here](https://github.com/maelyg/ont_amplicon/blob/master/docs/HTML_report_example.zip). It is provided within a gzipped folder with all associated files. It is too large to open in github, but if you git pull the repository onto your local machine, you can then unzip the folder and open the html file.  

The report consists of 3 main parts: input parameters, input data quality report and Consensus sequences.  

On the top left corner of the report, the **Facility**, **Analyst**, the time at which the **Analysis started**, the time at which the **Analysis completed** and the **Wall time**	are captured.

Just under it, the **input parameters section** displays the metadata that was provided in the samplesheet: 
- Sample ID  
- FASTQ files  
- Target taxon  
- Target gene  
- Amplicon size (nt) 
- FWD primer sequence (optional)  
- REV primer sequences (optional)  
- Test (opotional)  
- Method (optional)  

At the bottom of the **input parameters section**, there is also a **View all parameters** tab which is a link to a list of all the default pipeline parameters on the left and the parameters specifially set by the user for this sample on the right (and which match the parameters set in the yml file specified the user under the -params-file option).  The **View tool versions** displays the version of all the bioinformatics tools that were used by the nextflow pipeline processes.  

The **input data quality report section** displays the QC report output for the sample, and captures the number of starting raw reads and cleaned reads, and what percentage of the starting raw reads these represent. The outcome column refers to the **QC_FLAG**. Links to the Nanoplot reports for the raw and preprocessed reads are also available as well as the QC report for all the samples that were analysed at the same time.  

## Authors
Marie-Emilie Gauthier gauthiem@qut.edu.au  
Cameron Hyde c.hyde@qcif.edu.au  
