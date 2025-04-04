#!/bin/bash -l
#PBS -N ontamplicon
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=5:00:00


cd $PBS_O_WORKDIR
module load java
NXF_OPTS='-Xms1g -Xmx4g'

#passing parameters to a params file
nextflow run maelyg/ont_amplicon -profile singularity -resume -params-file params_peq_test.yml