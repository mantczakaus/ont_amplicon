#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage () {
    log.info """
    ont_amplicon
    Marie-Emilie Gauthier

    Usage:
    Run the command
    nextflow run ont_amplicon {arguments}...

    Required arguments:
      --analysis_mode                 clustering, map2ref
                                      Default: '' [required]

    Optional arguments:
      --help                          Will print this usage document
      -resume                         Resume a failed run
      --outdir                        Path to save the output file
                                      'results'
      --samplesheet '[path/to/file]'  Path to the csv file that contains the list of
                                      samples to be analysed by this pipeline.
                                      Default:  'index.csv'
 

    Contents of index.csv:
      sampleid,sample_files
      SAMPLE01,/user/folder/sample.fastq.gz
      SAMPLE02,/user/folder/*.fastq.gz

      #### Pre-processing and QC options ####
      --merge                         Merge fastq files with the same sample name
      --qc_only                       Only perform preliminary QC step using Nanoplot
                                      Default: false
      --preprocessing_only            Only perform preprocessing steps specied
                                      Default: false

      --adapter_trimming              Run porechop step
                                      Default: false
      --porechop_options              Porechop_ABI options
                                      Default: ''
      --porechop_custom_primers       Limit porechop search to custom adapters specified under porechop_custom_primers_path
                                      Default: ''
      --porechop_custom_primers_path  Path to custom adpaters for porechop
                                      Default: ''

      --qual_filt                     Run quality filtering step using chopper
                                      [False]
      --chopper_options               Chopper options
                                      Default: ''

      --host_filtering                Run host filtering step using Minimap2
                                      Default: false
      --host_fasta                    Path to the fasta file of nucleotide sequences to filter
                                      Default: ''

      #### Analysis mode and associated parameters ####
      ### Clustering (clustering) ###
      --rattle_clustering_options     Rattle clustering options
                                      Default: ''
      --rattle_polishing_options      Rattle polishing options
                                      Default: ''

      ### Map to reference (map2ref) ###
      --reference                     Path to the reference fasta file to map reads to
                                      Default: ''
      --medaka_consensus_options      Medaka options
                                      Default: ''
      --bcftools_min_coverage         Minimum coverage required by bcftools for annotation
                                      Default: '20'

      #### Blast options ####
      --blast_mode                    Specify whether megablast search is against NCBI or a custom database
                                      Default: ''. Select from 'ncbi' or 'localdb'
      --blast_threads                 Number of threads for megablast
                                      Default: '4'
      --blastn_db                     Path to blast database
                                      Default: '4'
      --blast_vs_ref                  blast versus reference specified in fasta file.
                                      Default: 'false'

      ### Reporting ###
      --contamination_flag_threshold  Percentage of maximum FPKM value to use as threshold for flagging detection as potential contamination
                                      Default: '0.01'

    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
if (params.blastn_db != null) {
    blastn_db_name = file(params.blastn_db).name
    blastn_db_dir = file(params.blastn_db).parent
}
if (params.reference != null) {
    reference_name = file(params.reference).name
    reference_dir = file(params.reference).parent
}
if (params.host_fasta != null) {
    host_fasta_dir = file(params.host_fasta).parent
}
if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
    kaiju_dbs_dir = file(params.kaiju_nodes).parent
}
if (params.krkdb != null) {
    krkdb_dir = file(params.krkdb).parent
}

if (params.porechop_custom_primers == true) {
    porechop_custom_primers_dir = file(params.porechop_custom_primers_path).parent
}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
    if (params.reference != null) {
      bindbuild = (bindbuild + "-B ${reference_dir} ")
    }
    if (params.host_fasta != null) {
      bindbuild = (bindbuild + "-B ${host_fasta_dir} ")
    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process BLASTN {
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly)
  output:
    tuple val(sampleid), path("${sampleid}*_blastn.bls"), emit: blast_results

  script:
  def blastoutput = assembly.getBaseName() + "_blastn.bls"
  
  if (params.blast_mode == "ncbi") {
    """
    cp ${blastn_db_dir}/taxdb.btd .
    cp ${blastn_db_dir}/taxdb.bti .
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
      -max_target_seqs 1
    """
  }
  
  else if (params.blast_mode == "localdb") {
    """
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 1
    """
  }
}

process BLASTN2REF {
  publishDir "${params.outdir}/${sampleid}", mode: 'copy', pattern: '*/*/*txt'
  tag "${sampleid}"
  label 'setting_1'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(assembly)

  output:
    path "*/*/${sampleid}_blastn_reference_vs_*.txt"

  script:
    """
    if [[ ${assembly} == *_assembly*.fa* ]] ;
    then
      if [[ ${assembly} == *canu_assembly*.fa* ]] ;
      then
        mkdir -p assembly/blast_to_ref
        blastn -query ${assembly} -subject ${reference_dir}/${reference_name} -evalue 1e-3 -out ${sampleid}_blastn_reference_vs_canu_assembly_tmp.txt \
        -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs' -max_target_seqs 5

        echo "qseqid\tsacc\tlength\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tqcovhsp\tqcovs" > header

        cat header ${sampleid}_blastn_reference_vs_canu_assembly_tmp.txt >  assembly/blast_to_ref/${sampleid}_blastn_reference_vs_canu_assembly.txt
      elif [[ ${assembly} == *flye_assembly*.fasta ]] ;
      then
        mkdir -p assembly/blast_to_ref

        blastn -query ${assembly} -subject ${reference_dir}/${reference_name} -evalue 1e-3 -out ${sampleid}_blastn_reference_vs_flye_assembly_tmp.txt \
        -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs' -max_target_seqs 5

        echo "qseqid\tsacc\tlength\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tqcovhsp\tqcovs" > header

        cat header ${sampleid}_blastn_reference_vs_flye_assembly_tmp.txt > assembly/blast_to_ref/${sampleid}_blastn_reference_vs_flye_assembly.txt
      fi
    elif [[ ${assembly} == *clustering.fasta ]] ;
    then
      mkdir -p clustering/blast_to_ref

      blastn -query ${assembly} -subject ${reference_dir}/${reference_name} -evalue 1e-3 -out ${sampleid}_blastn_reference_vs_clustering_tmp.txt \
      -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs' -max_target_seqs 5

      echo "qseqid\tsacc\tlength\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tqcovhsp\tqcovs" > header

      cat header ${sampleid}_blastn_reference_vs_clustering_tmp.txt > clustering/blast_to_ref/${sampleid}_blastn_reference_vs_clustering_assembly.txt
    fi
    """
}

process CHOPPER {
  publishDir "${params.outdir}/${sampleid}/preprocessing/chopper", pattern: '*_chopper.log', mode: 'link'
  tag "${sampleid}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(sample)

  output:
    path("${sampleid}_chopper.log")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), emit: chopper_filtered_fq

  script:
  def chopper_options = (params.chopper_options) ? " ${params.chopper_options}" : ''
    """
    gunzip -c ${sample} | chopper ${chopper_options} 2> ${sampleid}_chopper.log | gzip > ${sampleid}_filtered.fastq.gz
    """
}

process COVERM {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy', pattern: '*coverage_histogram.txt'

  input:
    tuple val(sampleid), path("*")

  output:
    path("*coverage_histogram.txt"), optional: true
    tuple val(sampleid), path("*_coverm_summary.txt"), emit: coverm_results, optional: true

  script:
    """
    if compgen -G "*.bam" > /dev/null;
      then
        for i in *bam;
        do
          echo \${i%.sorted.bam}
          filen=`echo "\${i%.sorted.bam}"`
          coverm genome --genome-fasta-files \${filen}.fasta --bam-files \${i} --threads ${task.cpus} --output-file \${filen}_coverm_summary.txt -m count mean variance rpkm covered_bases length --min-covered-fraction 0;
          coverm genome --genome-fasta-files \${filen}.fasta --bam-files \${i} --threads ${task.cpus} --output-file \${filen}_coverage_histogram.txt -m coverage_histogram --min-covered-fraction 0;
        done
    fi
    """
}

process COVSTATS {
  tag "$sampleid"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy'

  input:
    tuple val(sampleid), path("*")

  output:
    path("*top_blast_with_cov_stats.txt"), optional: true
    path("*top_blast_with_cov_stats.txt"), emit: detections_summary, optional: true

  script:
    """
    if compgen -G "*coverm_summary.txt" > /dev/null;
      then
        derive_coverage_stats.py --sample ${sampleid}
    fi
    """
}

process DETECTION_REPORT {
  label "local"
    publishDir "${params.outdir}/detection_summary", mode: 'copy', overwrite: true
    containerOptions "${bindOptions}"

  input:
    path('*')

  output:
    path("detection_summary*.txt")

  script:
    """
    detection_summary.py --threshold ${params.contamination_flag_threshold}
    """
}

process EXTRACT_READS {
  tag "${sampleid}"
  label "setting_11"
  publishDir "${params.outdir}/${sampleid}/host_filtering", mode: 'copy', pattern: '{*.fastq.gz,*reads_count.txt}'

  input:
  tuple val(sampleid), path(fastq), path(unaligned_ids)
  output:
  path("*reads_count.txt"), emit: read_counts
  file("${sampleid}_unaligned_reads_count.txt")
  file("${sampleid}_unaligned.fastq.gz")
  tuple val(sampleid), path("*_unaligned.fastq"), emit: unaligned_fq

  script:
  """
  seqtk subseq ${fastq} ${unaligned_ids} > ${sampleid}_unaligned.fastq
  gzip -c ${sampleid}_unaligned.fastq > ${sampleid}_unaligned.fastq.gz
  
  n_lines=\$(expr \$(cat ${sampleid}_unaligned.fastq | wc -l) / 4)
  echo \$n_lines > ${sampleid}_unaligned_reads_count.txt
  """
}

process EXTRACT_REF_FASTA {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy', pattern: '*fasta'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("*fasta"), emit: fasta_files, optional: true
  
  script:
    """
    cut -f1,4 ${blast_results} | sed '1d' | sed 's/ /_/g' > ids_to_retrieve.txt
    if [ -s ids_to_retrieve.txt ]
      then
        for i in `cut -f2  ids_to_retrieve.txt`; do j=`grep \${i} ids_to_retrieve.txt | cut -f1`; efetch -db nucleotide  -id \${i} -format fasta > ${sampleid}_\${i}_\${j}.fasta ; done
    fi
    """
}

process EXTRACT_VIRAL_BLAST_HITS {
  tag "${sampleid}"
  label "setting_2"
  publishDir "$params.outdir/$sampleid", mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    file "*/*/${sampleid}*_blastn_top_hits.txt"
    file "*/*/${sampleid}*_blastn_top_viral_hits*.txt"
    file "*/*/${sampleid}*_blastn_top_viral_spp_hits.txt"
    file "*/*/${sampleid}*_queryid_list_with_viral_match.txt"
    file "*/*/${sampleid}*_viral_spp_abundance*.txt"
    file "*/*/*report*html"

    tuple val(sampleid), path("*/*/${sampleid}*_viral_spp_abundance.txt"), emit: blast_results
    tuple val(sampleid), path("*/*/${sampleid}*_viral_spp_abundance_filtered.txt"), emit: blast_results_filt
    tuple val(sampleid), path("*/*/${sampleid}*_blastn_top_viral_spp_hits.txt"), emit: blast_results2

  script:
    """
    if [[ ${params.analysis_mode} == "clustering" ]] ;
    then
      mkdir -p clustering/blastn
      cat ${blast_results} > clustering/blastn/${sampleid}_blastn.txt
      cd clustering/blastn
      select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}_blastn.txt --analysis_method clustering --mode ${params.blast_mode}
    elif [[ ${params.analysis_mode} == "denovo_assembly" ]] ;
    then
      mkdir -p assembly/blastn
      cat ${blast_results} > assembly/blastn/${sampleid}_blastn.txt
      cd assembly/blastn
      select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}_blastn.txt --analysis_method assembly --mode ${params.blast_mode}
    elif [[ ${params.analysis_mode} == "read_classification" ]] ;
    then
      mkdir -p read_classification/homology_search
      cat ${blast_results} > read_classification/homology_search/${sampleid}_blastn.txt
      cd read_classification/homology_search
      select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}_blastn.txt --analysis_method read_classification --mode ${params.blast_mode}
    fi
    """
}

process FASTCAT {
  publishDir "${params.outdir}/${sampleid}/qc/fastcat", mode: 'copy'
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_stats.tsv")
    path("histograms/*")
    tuple val(sampleid), path("${sampleid}.fastq.gz"), emit: merged

  script:
    """
    fastcat \
        -s ${sampleid} \
        -f ${sampleid}_stats.tsv \
        --histograms histograms \
        ${fastq} \
        | bgzip > ${sampleid}.fastq.gz
    """
}

process FILTER_VCF {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(vcf)

  output:
    path("${sampleid}_medaka.consensus.fasta")
    path("${sampleid}_medaka.annotated.vcf.gz")

  script:
    """
    bcftools reheader ${vcf} -s <(echo '${sampleid}') \
    | bcftools filter \
        -e 'INFO/DP < ${params.bcftools_min_coverage}' \
        -s LOW_DEPTH \
        -Oz -o ${sampleid}_medaka.annotated.vcf.gz

    # create consensus
    bcftools index ${sampleid}_medaka.annotated.vcf.gz
    bcftools consensus -f ${reference_dir}/${reference_name} ${sampleid}_medaka.annotated.vcf.gz \
        -i 'FILTER="PASS"' \
        -o ${sampleid}_medaka.consensus.fasta
    """
}

process MAPPING_BACK_TO_REF {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy', pattern: '*sorted.bam*'
  //publishDir "${params.outdir}/01_VirReport/${sampleid}/alignments/NT", mode: 'link', overwrite: true, pattern: "*{.fa*,.fasta,metrics.txt,scores.txt,targets.txt,stats.txt,log.txt,.bcf*,.vcf.gz*,.bam*}"

  input:
    tuple val(sampleid), path(results)

  output:
    path("*bam"), optional: true
    path("*bam.bai"), optional: true
    tuple val(sampleid), path("*sorted.bam"), emit: bam_files, optional: true
    tuple val(sampleid), path("*sorted.bam.bai"), emit: bai_files, optional: true

  script:
    """
    if compgen -G "*.fasta" > /dev/null;
      then
        mapping_back_to_ref.py --fastq ${sampleid}_preprocessed.fastq.gz
    fi
    """
}

process MEDAKA {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(bam), path(bai)

  output:
    tuple val(sampleid), path("${sampleid}_medaka.annotated.unfiltered.vcf"), emit: unfilt_vcf

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka consensus ${bam} ${sampleid}_medaka_consensus_probs.hdf \
      ${medaka_consensus_options} --threads ${task.cpus}

    medaka variant ${reference_dir}/${reference_name} ${sampleid}_medaka_consensus_probs.hdf ${sampleid}_medaka.vcf
    medaka tools annotate --dpsp ${sampleid}_medaka.vcf ${reference_dir}/${reference_name} ${bam} \
          ${sampleid}_medaka.annotated.unfiltered.vcf
    """
}

process MINIMAP2_ALIGN_RNA {
  tag "${sampleid}"
  label "setting_8"
  containerOptions "${bindOptions}"

  input:
  tuple val(sampleid), path(fastq)
  path(reference)
  output:
  tuple val(sampleid), path(fastq), path("${sampleid}_unaligned_ids.txt"), emit: sequencing_ids

  script:
  """
  minimap2 -ax splice -uf -k14 -L ${reference} ${fastq} -t ${task.cpus} > ${sampleid}.sam
  awk '\$6 == "*" { print \$0 }' ${sampleid}.sam | cut -f1 | sort | uniq >  ${sampleid}_unaligned_ids.txt
  """
}

process MINIMAP2_REF {
  tag "${sampleid}"
  label 'setting_2'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(sample)

  output:
    tuple val(sampleid), file("${sampleid}_aln.sam"), emit: aligned_sample

  script:
    """
    minimap2 -ax map-ont --MD --sam-hit-only ${reference_dir}/${reference_name} ${sample} > ${sampleid}_aln.sam
    """
}

process MOSDEPTH {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy'

  input:
    tuple val(sampleid), path("*")

  output:
    path("*.mosdepth.global.dist.txt"), optional: true
    path("*.per-base.bed.gz*"), optional: true
    tuple val(sampleid), path("*mosdepth.global.dist.txt"), emit: mosdepth_results, optional: true

  script:
    """
    if compgen -G "*.bam" > /dev/null;
      then
        for i in *bam;
        do
          echo \${i%.sorted.bam}
          filen=`echo "\${i%.sorted.bam}"`
          mosdepth \${filen} \${i};
        done
    fi
    """
}

process PORECHOP_ABI {
  tag "${sampleid}"
  publishDir "$params.outdir/${sampleid}/preprocessing/porechop",  mode: 'link', pattern: '*_porechop.log'
  label "setting_9"
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(sample)

  output:
    file("${sampleid}_porechop_trimmed.fastq.gz")
    file("${sampleid}_porechop.log")
    tuple val(sampleid), file("${sampleid}_porechop_trimmed.fastq.gz"), emit: porechopabi_trimmed_fq

  script:
  def porechop_options = (params.porechop_options) ? " ${params.porechop_options}" : ''
    """
    if [[ ${params.porechop_custom_primers} == true ]]; then
      porechop_abi -i ${sample} -t ${task.cpus} -o ${sampleid}_porechop_trimmed.fastq.gz --custom_adapters ${params.porechop_custom_primers_path} ${porechop_options}  > ${sampleid}_porechop.log
    else
      porechop_abi -i ${sample} -t ${task.cpus} -o ${sampleid}_porechop_trimmed.fastq.gz ${porechop_options}  > ${sampleid}_porechop.log
    fi
    """
}

process QCREPORT {
  publishDir "${params.outdir}/qc_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"

  input:
    path multiqc_files

  output:
    path("run_qc_report_*txt")
    path("run_qc_report_*html")

  script:
    """
    seq_run_qc_report.py --host_filtering ${params.host_filtering} --adapter_trimming ${params.adapter_trimming} --quality_trimming ${params.qual_filt}
    """
}

process RATTLE {
  publishDir "${params.outdir}/${sampleid}/clustering/rattle", mode: 'copy', pattern: 'transcriptome.fq'
  tag "${sampleid}"
  label 'setting_7'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(fastq)

  output:
    file("transcriptome.fq")
    tuple val(sampleid), path("transcriptome.fq"), emit: clusters
    tuple val(sampleid), path(fastq), path("transcriptome.fq"), emit: clusters2

  script:
  def rattle_polishing_options = (params.rattle_polishing_options) ? " ${params.rattle_polishing_options}" : ''
  def rattle_clustering_options = (params.rattle_clustering_options) ? " ${params.rattle_clustering_options}" : ''
    """
    rattle cluster -i ${fastq} -t ${task.cpus} ${rattle_clustering_options}  -o .
    rattle cluster_summary -i ${fastq} -c clusters.out > ${sampleid}_cluster_summary.txt
    mkdir clusters
    rattle extract_clusters -i ${fastq} -c clusters.out -l ${sampleid} -o clusters --fastq
    rattle correct -i ${fastq} -c clusters.out -t ${task.cpus} -l ${sampleid}
    rattle polish -i consensi.fq -t ${task.cpus} --summary ${rattle_polishing_options}
    """
}

//trims fastq read names after the first whitespace
process REFORMAT {
  tag "${sampleid}"
  label "setting_3"
  publishDir "$params.outdir/${sampleid}/preprocessing", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq)

  output:
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: reformatted_fq
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: cov_derivation_ch

  script:
    """
    reformat.sh in=${fastq} out=${sampleid}_preprocessed.fastq.gz trd qin=33
    """
}

process SAMTOOLS {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(sample)

  output:
    path "${sampleid}_aln.sorted.bam"
    path "${sampleid}_aln.sorted.bam.bai"
    path "${sampleid}_coverage.txt"
    path "${sampleid}_histogram"
    tuple val(sampleid), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_sample

  script:
    """
    samtools view -Sb -F 4 ${sample} | samtools sort -o ${sampleid}_aln.sorted.bam
    samtools index ${sampleid}_aln.sorted.bam
    samtools coverage ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt  > ${sampleid}_coverage.txt
    samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram
    """
}

include { FASTQ2FASTA } from './modules.nf'
include { FASTQ2FASTA as FASTQ2FASTA_STEP1} from './modules.nf'
include { NANOPLOT as QC_PRE_DATA_PROCESSING } from './modules.nf'
include { NANOPLOT as QC_POST_DATA_PROCESSING } from './modules.nf'

workflow {
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), file(row.sample_files)) }
      .set{ ch_sample }
  } else { exit 1, "Input samplesheet file not specified!" }

  if ( params.analysis_mode == 'clustering' & params.megablast) {
    if (!params.blast_vs_ref) {
      if ( params.blastn_db == null) {
        error("Please provide the path to a blast database using the parameter --blastn_db.")
      }
    }
    else if (params.blast_vs_ref ) {
      if ( params.reference == null) {
      error("Please provide the path to a reference fasta file with the parameter --reference.")
      }
    }
  }
  else if ( params.analysis_mode == 'map2ref' ) {
    if ( params.reference == null) {
      error("Please provide the path to a reference fasta file with the parameter --reference.")
      }
  }

  if (params.merge) {
    FASTCAT ( ch_sample )
    QC_PRE_DATA_PROCESSING ( FASTCAT.out.merged )
    fq = FASTCAT.out.merged
  }
  else {
    fq = ch_sample
    QC_PRE_DATA_PROCESSING ( fq )
  }

  // Data pre-processing
  if (!params.qc_only) {
    // Remove adapters using PORECHOP_ABI
    if (params.adapter_trimming) {
      PORECHOP_ABI ( fq )
      trimmed_fq = PORECHOP_ABI.out.porechopabi_trimmed_fq
    }

    else {
      trimmed_fq = fq
    }

    // Quality filtering of reads using chopper
    if (params.qual_filt) {
      CHOPPER ( trimmed_fq)
      filtered_fq = CHOPPER.out.chopper_filtered_fq
    }
    else { filtered_fq = trimmed_fq
    }

    REFORMAT( filtered_fq )

    if ( params.qual_filt & params.adapter_trimming | !params.qual_filt & params.adapter_trimming | params.qual_filt & !params.adapter_trimming) {
      QC_POST_DATA_PROCESSING ( filtered_fq )
    }

    if (params.host_filtering) {
      if ( params.host_fasta == null) {
        error("Please provide the path to a fasta file of host sequences that need to be filtered with the parameter --host_fasta.")
      }
      else {
        MINIMAP2_ALIGN_RNA ( REFORMAT.out.reformatted_fq, params.host_fasta )
        EXTRACT_READS ( MINIMAP2_ALIGN_RNA.out.sequencing_ids )
        final_fq = EXTRACT_READS.out.unaligned_fq
      }
    }
    else {
      final_fq = REFORMAT.out.reformatted_fq
    }

    if ( params.qual_filt & params.host_filtering | params.adapter_trimming & params.host_filtering ) {
      ch_multiqc_files = Channel.empty()
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_READS.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }

    else if ( params.host_filtering & !params.adapter_trimming & !params.qual_filt ) {
      ch_multiqc_files = Channel.empty()
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_READS.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }

    else if ( params.qual_filt & !params.host_filtering | params.adapter_trimming & !params.host_filtering) {
      ch_multiqc_files = Channel.empty()
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }

    if (!params.preprocessing_only) {
      
      if ( params.analysis_mode == 'clustering' ) {
        //Perform clustering using Rattle
        RATTLE ( final_fq )
        FASTQ2FASTA( RATTLE.out.clusters )
        contigs = FASTQ2FASTA.out.fasta 

        //limit blast homology search to a reference
        if (params.blast_vs_ref) {
          BLASTN2REF ( contigs )
        }
        //blast to database
        else {
        ASSEMBLY_BLASTN ( contigs )
        EXTRACT_VIRAL_BLAST_HITS ( ASSEMBLY_BLASTN.out.blast_results )
        EXTRACT_REF_FASTA (EXTRACT_VIRAL_BLAST_HITS.out.blast_results2)

        mapping_ch = EXTRACT_REF_FASTA.out.fasta_files.concat(REFORMAT.out.cov_derivation_ch).groupTuple().map { [it[0], it[1].flatten()] }
        MAPPING_BACK_TO_REF ( mapping_ch )
        bamf_ch = MAPPING_BACK_TO_REF.out.bam_files.concat(MAPPING_BACK_TO_REF.out.bai_files, EXTRACT_REF_FASTA.out.fasta_files).groupTuple().map { [it[0], it[1].flatten()] }
        MOSDEPTH (bamf_ch)
        COVERM (bamf_ch)
        cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.concat(COVERM.out.coverm_results, EXTRACT_REF_FASTA.out.fasta_files, EXTRACT_VIRAL_BLAST_HITS.out.blast_results2, QC_PRE_DATA_PROCESSING.out.stats).groupTuple().map { [it[0], it[1].flatten()] }
        COVSTATS(cov_stats_summary_ch)

        DETECTION_REPORT(COVSTATS.out.detections_summary.collect().ifEmpty([]))
        }
      }

      //Perform direct alignment to a reference
      else if ( params.analysis_mode == 'map2ref') {
        MINIMAP2_REF ( final_fq )
        SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        MEDAKA ( SAMTOOLS.out.sorted_sample )
        FILTER_VCF ( MEDAKA.out.unfilt_vcf )
      }
      else {
        error("Analysis mode (clustering, map2ref) not specified with e.g. '--analysis_mode clustering' or via a detectable config file.")
      }
    }
  }
}
