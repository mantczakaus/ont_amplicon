#!/bin/bash -l
#PBS -N ontamplicon
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=5:00:00


cd $PBS_O_WORKDIR
module load java
NXF_OPTS='-Xms1g -Xmx4g'

nextflow run maelyg/ont_amplicon -profile singularity -resume
	--analysis_mode clustering \
	--merge \
	--adapter_trimming \
	--qual_filt \
	--chopper_options '-q 8 -l 100' \
	--blast_threads 2 \
	--blastn_db /full_path_to_blast_db_nt_files \
	--blastn_COI /full_path_to_MetaCOXI_Seqs_fasta_file \
	--taxdump /full_path_to_.taxonkit_file