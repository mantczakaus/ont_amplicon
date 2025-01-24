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
if (params.blastn_COI != null) {
    blastn_COI_name = file(params.blastn_COI).name
    blastn_COI_dir = file(params.blastn_COI).parent
}

if (params.reference != null) {
    reference_name = file(params.reference).name
    reference_dir = file(params.reference).parent
}
if (params.host_fasta != null) {
    host_fasta_dir = file(params.host_fasta).parent
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
    if (params.blastn_COI != null) {
      bindbuild = (bindbuild + "-B ${blastn_COI_dir} ")
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
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '*_megablast*.txt'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly)
  output:
    path("${sampleid}*_megablast*.txt")
    tuple val(sampleid), path("${sampleid}*_megablast_top_10_hits.txt"), emit: blast_results
//    tuple val(sampleid), path("${sampleid}_ids_to_reverse_complement.txt"), path("${sampleid}*__megablast_COI_top_hit.txt"), emit: coi_blast_results, optional: true

  script:
  def blast_output = assembly.getBaseName() + "_megablast_top_10_hits.txt"
//  def blast_output_COI = assembly.getBaseName() + "_megablast_COI_top_hit.txt"
  
  if (params.blast_mode == "ncbi") {
    """
    cp ${blastn_db_dir}/taxdb.btd .
    cp ${blastn_db_dir}/taxdb.bti .
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blast_output} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length nident pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames sskingdoms' \
      -max_target_seqs 10
    """
  }
}

process BLASTN2 {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '*_megablast*.txt'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly), val(gene_target)
  output:
    path("${sampleid}*_megablast*.txt")
    tuple val(sampleid), path("${sampleid}*_megablast_top_10_hits.txt"), emit: blast_results
    tuple val(sampleid), path("${sampleid}_ids_to_reverse_complement.txt"), path("${sampleid}*__megablast_COI_top_hit.txt"), emit: coi_blast_results, optional: true

  script:
  def blast_output = assembly.getBaseName() + "_megablast_top_10_hits.txt"
  def blast_output_COI = assembly.getBaseName() + "_megablast_COI_top_hit.txt"
  
  if (params.blast_mode == "ncbi") {
    """
    cp ${blastn_db_dir}/taxdb.btd .
    cp ${blastn_db_dir}/taxdb.bti .
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blast_output} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length nident pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames sskingdoms' \
      -max_target_seqs 10
    """
  }
}


process BLASTN_COI {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '*_megablast*.txt'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly), val(gene_target)
  output:
    path("${sampleid}*_megablast_COI_top_hit.txt")
    tuple val(sampleid), path("${sampleid}_ids_to_reverse_complement.txt"), emit: coi_blast_results, optional: true

  script:
  def blast_output_COI = assembly.getBaseName() + "_megablast_COI_top_hit.txt"
    """
    blastn -query ${assembly} \
      -db ${params.blastn_COI} \
      -out ${blast_output_COI} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send evalue bitscore sstrand' \
      -max_target_seqs 1 \
      -max_hsps 1
      grep minus ${blast_output_COI} | cut -f1 > ${sampleid}_ids_to_reverse_complement.txt
    """
}

process REVCOMP {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '{*fasta}'
  tag "${sampleid}"
  label "setting_2"
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(contigs), path(ids_to_revcomp)


  output:
    file "${sampleid}_final_polished_consensus_rc.fasta"
    tuple val(sampleid), path("${sampleid}_final_polished_consensus_rc.fasta"), emit: revcomp, optional: true

  script:
    """
    reverse_complement.py --sample ${sampleid} --ids_to_rc ${ids_to_revcomp} --fasta ${contigs}
    
    """

}


/*
if [ ${gene_targets} = "16s" ];
    then
      blastn -query ${assembly} \
        -db ${params.blastn_16s} \
        -out ${blast_output_16s} \
        -evalue 1e-3 \
        -num_threads ${params.blast_threads} \
        -outfmt '6 qseqid sgi sacc length nident pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames sskingdoms' \
        -max_target_seqs 10


  else if (params.blast_mode == "localdb") {
    """
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 10
    """
  }
}
*/
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
/*
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
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy', pattern: '*fasta'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("*fasta"), emit: fasta_files, optional: true
  
  script:
    """
    cut -f1,3 ${blast_results} | sed '1d' | sed 's/ /_/g' > ids_to_retrieve.txt
    if [ -s ids_to_retrieve.txt ]
      then
        for i in `cut -f2  ids_to_retrieve.txt`; do j=`grep \${i} ids_to_retrieve.txt | cut -f1`; efetch -db nucleotide  -id \${i} -format fasta > ${sampleid}_\${i}_\${j}.fasta ; done
    fi
    """
}
*/

process CUTADAPT {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/polishing", pattern: '{*.fasta,*_cutadapt.log}', mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(consensus)

  output:
    file("${sampleid}_cutadapt.log")
    file("${sampleid}_final_polished_consensus.fasta")
    tuple val(sampleid), path("${sampleid}_final_polished_consensus.fasta"), emit: trimmed

  script:
    """
    cutadapt -j ${task.cpus} --trim-n -o ${sampleid}_final_polished_consensus.fasta ${consensus} > ${sampleid}_cutadapt.log
    """

}

process EXTRACT_REF_FASTA {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy', pattern: '*fasta'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("*fasta"), emit: fasta_files, optional: true
  
  script:
    """
    cut -f 3 ${blast_results} | grep -v sacc | sed 's/^/>/' > seq_name.txt
    cut -f21  ${blast_results} | grep -v 'qseq' | sed 's/-//g' > sequence.txt
    cat seq_name.txt sequence.txt > reference.fasta
    """
}


process EXTRACT_BLAST_HITS {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '{*.txt,*.html,*fasta}'
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    file "${sampleid}*_blastn_top_hits.txt"
    file "${sampleid}_final_polished_consensus_match.fasta"
    file "${sampleid}_reference_match.fasta"

    //file "${sampleid}*_queryid_list_with_spp_match.txt"
    //file "${sampleid}*_spp_abundance*.txt"
    //file "*report*html"
    tuple val(sampleid), path("${sampleid}*_megablast_blastn_top_hits.txt"), emit: topblast, optional: true
    tuple val(sampleid), path("${sampleid}*_megablast_blastn_top_hits.txt"), emit: topblast2, optional: true
    tuple val(sampleid), path("${sampleid}_reference_match.fasta"), emit: fasta_files, optional: true

  script:
    """
    select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}*_top_10_hits.txt --mode ${params.blast_mode}

    # extract segment of consensus sequence that align to reference
    awk  -F  '\\t' 'NR>1 { printf ">%s\\n%s\\n",\$1,\$21 }' ${sampleid}*_top_hits.txt > ${sampleid}_final_polished_consensus_match.fasta

    # extract segment of reference that align to consensus sequence
    awk  -F  '\\t' 'NR>1 { printf ">%s_%s\\n%s\\n",\$1,\$3,\$22 }' ${sampleid}*_top_hits.txt > ${sampleid}_reference_match.fasta
    
    """

}

/*
    if [ ${gene_targets} = "16s" ];
    then
      select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}_medaka_consensus_megablast_16s_top_10_hits.txt --mode ${params.blast_mode}
    fi
    """
}

/*
echo "# BLASTN 2.13.0+" > template.txt
    echo "# Query: ${sampleid}" >> template.txt
    echo "# Database: /scratch/datasets/blast_db/20240730/nt" >> template.txt
    echo "# Fields: query acc., subject acc., subject title, evalue, q. start, q. end, s. start, s. end, bit score, alignment length, mismatches, % identity, identical, % query coverage per subject, subject seq, query seq" >> template.txt
    echo "Top 10 hits reported" >> template.txt
    cat  template.txt ${sampleid}_blastn_report.txt > ${sampleid}_blastn_report_ed.txt
*/

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

process FASTQ2FASTA {
  publishDir "${params.outdir}/${sampleid}/clustering", mode: 'copy', pattern: '*_rattle.fasta'
  tag "${sampleid}"
  label "setting_2"

  input:
  tuple val(sampleid), path(fastq), path(assembly)
  output:
  tuple val(sampleid), path(fastq), path("${sampleid}_rattle.fasta"), emit: fasta

  script:
  """
  cut -f1,3 -d ' ' ${assembly} | sed 's/ total_reads=/_RC/' > ${sampleid}_tmp.fastq
  seqtk seq -A -C ${sampleid}_tmp.fastq > ${sampleid}_rattle.fasta
  """
}
//grep '@cluster' transcriptome.fq  | cut -f1,3 -d ' '  | sed 's/total_reads=//' | sort -k2,2 -rn | sed 's/ /_RC/' | sed 's/@//' | head -n 10 > ${sampleid}_most_abundant_clusters_ids.txt
/*
process CONSENSUS_FILTER_VCF {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(vcf), path(assembly)

  output:
    path("${sampleid}_medaka.consensus.fasta")
    path("${sampleid}_medaka.annotated.vcf.gz")
    tuple val(sampleid), path("${sampleid}_medaka.consensus.fasta"), emit: fasta
  script:
    """
    bcftools reheader ${vcf} -s <(echo '${sampleid}') \
    | bcftools filter \
        -e 'INFO/DP < ${params.bcftools_min_coverage}' \
        -s LOW_DEPTH \
        -Oz -o ${sampleid}_medaka.annotated.vcf.gz

    # create consensus
    bcftools index ${sampleid}_medaka.annotated.vcf.gz
    bcftools consensus -f ${assembly} ${sampleid}_medaka.annotated.vcf.gz \
        -i 'FILTER="PASS"' \
        -o ${sampleid}_medaka.consensus.fasta
    """
}
*/
process EXTRACT_TAXONOMY {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '{*.txt}'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    path("*taxonomy_file.txt"), optional: true
    tuple val(sampleid), path("*taxonomy_file.txt"), emit: taxonomy, optional: true

  script:
    """
    c1grep() { grep  "\$@" || test \$? = 1; }
    c1grep "Eukaryota" ${blast_results} | cut -f3 | sort | uniq > ids_to_retrieve.txt
    if [ -s ids_to_retrieve.txt ]
      then
        for i in `cut -f1 ids_to_retrieve.txt`; do esearch -db nuccore -query "\${i}" | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS=", ";}{print \$3;}'  >> ${sampleid}_taxonomy_file.txt; printf "\t\$i" >> ${sampleid}_taxonomy_file.txt; printf "\n"  >> ${sampleid}_taxonomy_file.txt; done
    else
      touch ${sampleid}_taxonomy_file.txt
    fi
    """
}

process FILTER_CONSENSUS {
  publishDir "${params.outdir}/${sampleid}/megablast", mode: 'copy', pattern: '{*.txt}'
  tag "${sampleid}"
  label "setting_2"
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), val(spp_targets), val(gene_targets), val(size), path(taxonomy), path(blast_results)

  output:
    file "${sampleid}*blastn_top_hits_filtered.txt"

  script:
    """
    filter_blast_results.py --sample_name ${sampleid} --acc_phylo_info ${taxonomy} --spp_targets ${spp_targets} --gene_targets ${gene_targets} --blastn_results ${blast_results} --mode ${params.blast_mode} --target_size ${size}
    """
}
/*
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
*/
process DERIVE_MASKED_CONSENSUS {
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai), path(vcf)

  output:
    tuple val(sampleid), path("${sampleid}_masked.fasta"), path(bam), path(bai), path(vcf), emit: masked_fasta

  script:
    """
    # Get consensus fasta file
    bedtools genomecov -ibam ${bam} -bga > ${sampleid}_genome_cov.bed
    # Assign N to nucleotide positions that have zero coverage
    awk '\$4==0 {print}' ${sampleid}_genome_cov.bed > ${sampleid}_zero_cov_genome_cov.bed
    bedtools maskfasta -fi ${ref} -bed ${sampleid}_zero_cov_genome_cov.bed -fo ${sampleid}_masked.fasta
    """
}

process FILTER_VCF {
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai), path(vcf)
   

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
    bcftools consensus -f ${ref} -i 'FILTER="PASS"' ${sampleid}_medaka.annotated.vcf.gz --iupac-codes -H I -o ${sampleid}_medaka.consensus.fasta
    """
}
/*
process FILTER_VCF_BACKUP {
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(masked), path(bam), path(bai), path(vcf)
   

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
    bcftools consensus -f ${masked} -i 'FILTER="PASS"' ${sampleid}_medaka.annotated.vcf.gz --iupac-codes -H I -o ${sampleid}_medaka.consensus.fasta
    """
}
*/

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
/*
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
*/




process MEDAKA1 {
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai)

  output:
    tuple val(sampleid), path(ref), path(bam), path(bai), path("${sampleid}_medaka.annotated.unfiltered_sorted.vcf"), emit: unfilt_vcf

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka consensus ${bam} ${sampleid}_medaka_consensus_probs.hdf \
      ${medaka_consensus_options} --threads ${task.cpus}

    medaka variant ${ref} ${sampleid}_medaka_consensus_probs.hdf ${sampleid}_medaka.vcf
    medaka tools annotate --dpsp ${sampleid}_medaka.vcf ${ref} ${bam} \
          ${sampleid}_medaka.annotated.unfiltered.vcf
    cat ${sampleid}_medaka.annotated.unfiltered.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > ${sampleid}_medaka.annotated.unfiltered_sorted.vcf
    """
}

/*
process MEDAKA_CONSENSUS {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(bam), path(bai), path(assembly)

  output:
    tuple val(sampleid), path("${sampleid}_medaka.annotated.unfiltered.vcf"), path(assembly), emit: unfilt_vcf

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka consensus ${bam} ${sampleid}_medaka_consensus_probs.hdf \
      --threads ${task.cpus} ${medaka_consensus_options} 

    medaka variant ${assembly} ${sampleid}_medaka_consensus_probs.hdf ${sampleid}_medaka.vcf
    medaka tools annotate --dpsp ${sampleid}_medaka.vcf ${assembly} ${bam} \
          ${sampleid}_medaka.annotated.unfiltered.vcf
    """
}
*/
process RACON {
  publishDir "${params.outdir}/${sampleid}/polishing", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(fastq), path(assembly), path(paf)
  output:
   tuple val(sampleid), path(fastq), path("${sampleid}_racon_polished.fasta")
   tuple val(sampleid), path(fastq), path("${sampleid}_racon_polished.fasta"), emit: polished

  script:
    """
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus -q -1 --no-trimming -u \
        ${fastq} ${paf} ${assembly} \
        > ${sampleid}_racon_polished_tmp.fasta
    cut -f1 -d ' ' ${sampleid}_racon_polished_tmp.fasta > ${sampleid}_racon_polished.fasta
    """
}

process MEDAKA2 {
  publishDir "${params.outdir}/${sampleid}/polishing", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(fastq), path(assembly)
  output:
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai"), path("${sampleid}_samtools_consensus.fasta")
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai"), emit: consensus1
   tuple val(sampleid), path("${sampleid}_samtools_consensus.fasta"), emit: consensus2

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka_consensus -i ${fastq} -d ${assembly} -t ${task.cpus} -o ${sampleid}
    
    mv ${sampleid}/calls_to_draft.bam ${sampleid}_medaka_consensus.bam
    mv ${sampleid}/calls_to_draft.bam.bai ${sampleid}_medaka_consensus.bam.bai
    mv ${sampleid}/consensus.fasta ${sampleid}_medaka_consensus.fasta
    samtools consensus -f fasta -a -A ${sampleid}_medaka_consensus.bam --call-fract 0.5 -H 0.5 -o ${sampleid}_samtools_consensus.fasta
    
    """
}


/*
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} \$2 {print ">"\$0}' ${sampleid}_samtools_consensus_temp.fasta > ${sampleid}_samtools_consensus.fasta
original command
samtools consensus -f fasta -A -d 5 -H 0.5 ${sampleid}_medaka_consensus.bam -o ${sampleid}_samtools_consensus.fasta -X r10.4_sup

process MEDAKA_CONSENSUS_REF {
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(ref), path(fastq)
  output:
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai")
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), path("${sampleid}_medaka_consensus.bam"), path("${sampleid}_medaka_consensus.bam.bai"), emit: consensus1
   tuple val(sampleid), path("${sampleid}_medaka_consensus.fasta"), emit: consensus2

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka_consensus -i ${fastq} -d ${ref} -t ${task.cpus} -o ${sampleid}
    mv ${sampleid}/calls_to_draft.bam ${sampleid}_medaka_consensus.bam
    mv ${sampleid}/calls_to_draft.bam.bai ${sampleid}_medaka_consensus.bam.bai
    mv ${sampleid}/consensus.fasta ${sampleid}_medaka_consensus.fasta
    """
}

//medaka_consensus -i ${fastq} -d ${assembly} -t ${task.cpus} -o ${sampleid} --bacteria
process MINIMAP2_ALIGN {
  tag "${sampleid}"
  label "setting_8"
  containerOptions "${bindOptions}"

  input:
  tuple val(sampleid), path(fastq), path(assembly)
  output:
  tuple val(sampleid), path("${sampleid}.sam"), path(assembly), emit: aligned_sam

  script:
    """
    minimap2 -ax map-ont --sam-hit-only ${assembly} ${fastq} -t ${task.cpus} > ${sampleid}.sam
    """
}
*/
process MINIMAP2_RACON {
  tag "${sampleid}"
  label "setting_8"
  containerOptions "${bindOptions}"

  input:
  tuple val(sampleid), path(fastq), path(assembly)
  output:
  tuple val(sampleid), path(fastq), path(assembly), path("${sampleid}_pre-racon.paf"), emit: draft_mapping

  script:
    """
    minimap2 -L -x ava-ont -t ${task.cpus} ${assembly} ${fastq} > ${sampleid}_pre-racon.paf
    """
}
/*
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
*/
process MINIMAP2_REF {
  tag "${sampleid}"
  label 'setting_2'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(ref), path(fastq)

  output:
    tuple val(sampleid), path(ref), file("${sampleid}_aln.sam"), emit: aligned_sample

  script:
    """
    minimap2 -ax map-ont --MD --sam-hit-only ${ref} ${fastq} > ${sampleid}_aln.sam
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
  tag "${sampleid}"
  label 'setting_7'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(fastq), val(target_size)

  output:
    file("transcriptome.fq")
    tuple val(sampleid), path("transcriptome.fq"), emit: clusters
    tuple val(sampleid), path("${fastq}"), path("transcriptome.fq"), emit: clusters2

  script:
  def rattle_clustering_options = params.rattle_clustering_options ?: ''
  def rattle_polishing_options = params.rattle_polishing_options ?: ''
  def rattle_clustering_min_length = params.rattle_clustering_min_length ?: ''
  if (rattle_clustering_min_length == null) {
    if (target_size.toInteger() < 300) {
      rattle_clustering_min_length = '100'}
    else {
      rattle_clustering_min_length = '150'}
  }

    """
    rattle cluster -i ${fastq} -t ${task.cpus} --lower-length ${rattle_clustering_min_length} ${rattle_clustering_options} -o .
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
/*
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
*/
process SAMTOOLS {
  publishDir "${params.outdir}/${sampleid}/mapping_back_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(ref), path(sample)

  output:
    path "${sampleid}_aln.sorted.bam"
    path "${sampleid}_aln.sorted.bam.bai"
    path "${sampleid}_coverage.txt"
    path "${sampleid}_histogram"
    path "${sampleid}_samtools_consensus_from_ref.fasta"
    tuple val(sampleid), path(ref), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_sample

  script:
    """
    samtools view -Sb -F 4 ${sample} | samtools sort -o ${sampleid}_aln.sorted.bam
    samtools index ${sampleid}_aln.sorted.bam
    samtools coverage ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt  > ${sampleid}_coverage.txt
    samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram
    samtools consensus -f fasta -a -A ${sampleid}_aln.sorted.bam --call-fract 0.5 -H 0.5 -o ${sampleid}_samtools_consensus_from_ref.fasta
    """
}
/*
process SAMTOOLS_ALIGN {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(sam), path(assembly)

  output:
    path "${sampleid}_aln.sorted.bam"
    path "${sampleid}_aln.sorted.bam.bai"
    path "${sampleid}_coverage.txt"
    path "${sampleid}_histogram"
    tuple val(sampleid), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), path(assembly), emit: sorted_bam

  script:
    """
    samtools view -Sb -F 4 ${sam} | samtools sort -@ $task.cpus -o ${sampleid}_aln.sorted.bam
    samtools index ${sampleid}_aln.sorted.bam
    samtools coverage ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt  > ${sampleid}_coverage.txt
    samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram
    """
}
*/
include { NANOPLOT as QC_PRE_DATA_PROCESSING } from './modules.nf'
include { NANOPLOT as QC_POST_DATA_PROCESSING } from './modules.nf'

workflow {
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), file(row.sample_files)) }
      .set{ ch_sample }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.spp_targets), (row.gene_targets), (row.target_size)) }
      .set{ ch_targets }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.target_size)) }
      .set{ ch_target_size }
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.gene_targets)) }
      .set{ ch_gene_targets }
    
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.gene_targets)) }
      .filter { sampleid, gene_targets -> gene_targets.contains("COI") }
      .set{ ch_coi }
//    ch_coi.view()
    
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), (row.gene_targets)) }
      .filter { sampleid, gene_targets -> !gene_targets.contains("COI") }
      .set{ ch_other }

  
  
  } else { exit 1, "Input samplesheet file not specified!" }


  if ( params.analysis_mode == 'clustering') {
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
/*
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
*/

    final_fq = REFORMAT.out.reformatted_fq
    //if ( params.qual_filt & params.host_filtering | params.adapter_trimming & params.host_filtering ) {
    if ( params.qual_filt | params.adapter_trimming ) {
      ch_multiqc_files = Channel.empty()
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
//     ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_READS.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }

//    else if ( params.host_filtering & !params.adapter_trimming & !params.qual_filt ) {
//      ch_multiqc_files = Channel.empty()
//      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
//      ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_READS.out.read_counts.collect().ifEmpty([]))
//     QCREPORT(ch_multiqc_files.collect())
//    }

//    else if ( params.qual_filt & !params.host_filtering | params.adapter_trimming & !params.host_filtering) {
//      ch_multiqc_files = Channel.empty()
//      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
//      ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
//      QCREPORT(ch_multiqc_files.collect())
//    }

    if (!params.preprocessing_only) {
      
      if ( params.analysis_mode == 'clustering' ) {
        //Perform clustering using Rattle
        ch_fq_target_size = (final_fq.join(ch_target_size))//.view()
        RATTLE ( ch_fq_target_size )
        FASTQ2FASTA( RATTLE.out.clusters2 )
        //contigs = FASTQ2FASTA.out.fasta 

        //MINIMAP2_ALIGN ( FASTQ2FASTA.out.fasta )
        //SAMTOOLS_ALIGN ( MINIMAP2_ALIGN.out.aligned_sam )
        //MEDAKA_CONSENSUS ( SAMTOOLS_ALIGN.out.sorted_bam )
        MINIMAP2_RACON ( FASTQ2FASTA.out.fasta )
        RACON ( MINIMAP2_RACON.out.draft_mapping)
        MEDAKA2 ( RACON.out.polished )
        CUTADAPT ( MEDAKA2.out.consensus2 )
        //CONSENSUS_FILTER_VCF ( MEDAKA_CONSENSUS.out.unfilt_vcf )
       //CLAIR3  ( MEDAKA_CONSENSUS.out.consensus1 )
        //limit blast homology search to a reference

        if (params.blast_vs_ref) {
          //BLASTN2REF ( CONSENSUS_FILTER_VCF.out.fasta )
          BLASTN2REF ( CUTADAPT.out.trimmed )
        }
        //blast to database
        else {
          
        //ch_gene_targets_for_blast = (CONSENSUS_FILTER_VCF.out.fasta.join(ch_gene_targets))
        //ch_gene_targets_for_blast = (CUTADAPT.out.trimmed.join(ch_gene_targets))
        //BLASTN ( CONSENSUS_FILTER_VCF.out.fasta )

//       ch_coi = ch_gene_targets
//        ch_gene_targets.filter { sampleid, gene_targets -> gene_targets.contains("COI").set { ch_coi } }
        ch_coi_for_blast = (CUTADAPT.out.trimmed.join(ch_coi))
        BLASTN_COI(ch_coi_for_blast)
        ch_revcomp = (CUTADAPT.out.trimmed.join(BLASTN_COI.out.coi_blast_results))
        REVCOMP ( ch_revcomp )
        BLASTN ( REVCOMP.out.revcomp )
       
//        if ( params.gene_targets == "COI") {
 //         BLASTN_COI ( ch_gene_targets_for_blast )
//        ch_gene_targets.filter { sampleid, gene_targets -> !gene_targets.contains("COI").set { ch_other } }
//        ch_other = ch_gene_targets
//          | filter { sampleid, gene_targets -> !gene_targets.contains("COI") }
        ch_other_for_blast = (CUTADAPT.out.trimmed.join(ch_other))

        BLASTN2 ( ch_other_for_blast )
        ch_blast_merged = BLASTN.out.blast_results.mix(BLASTN2.out.blast_results.ifEmpty([]))
        EXTRACT_BLAST_HITS ( ch_blast_merged )
//        EXTRACT_TAXONOMY ( EXTRACT_BLAST_HITS.out.topblast )
//        ch_filter_for_targets = (ch_targets.join(EXTRACT_TAXONOMY.out.taxonomy))
//        ch_filter_for_targets = (ch_filter_for_targets.join(EXTRACT_BLAST_HITS.out.topblast2))
        //FILTER_CONSENSUS ( ch_filter_for_targets )
        //EXTRACT_REF_FASTA (EXTRACT_BLAST_HITS.out.topblast2)

        //mapping_ch = EXTRACT_REF_FASTA.out.fasta_files.concat(REFORMAT.out.cov_derivation_ch).groupTuple().map { [it[0], it[1].flatten()] }
        mapping_ch = (EXTRACT_BLAST_HITS.out.fasta_files.join(REFORMAT.out.cov_derivation_ch))
        
//        MEDAKA_CONSENSUS_REF ( mapping_ch )
        MINIMAP2_REF ( mapping_ch )
        SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        //MEDAKA1 ( SAMTOOLS.out.sorted_sample )
        //FILTER_VCF ( MEDAKA1.out.unfilt_vcf )
        //DERIVE_MASKED_CONSENSUS ( MEDAKA1.out.unfilt_vcf )
        
        //DERIVE_MASKED_CONSENSUS ( MEDAKA1.out.unfilt_vcf )
        //FILTER_VCF ( DERIVE_MASKED_CONSENSUS.out.masked_fasta )
        //bamf_ch = MAPPING_BACK_TO_REF.out.bam_files.concat(MAPPING_BACK_TO_REF.out.bai_files, EXTRACT_REF_FASTA.out.fasta_files).groupTuple().map { [it[0], it[1].flatten()] }
        //MOSDEPTH (bamf_ch)
        //COVERM (bamf_ch)
        //cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.concat(COVERM.out.coverm_results, EXTRACT_REF_FASTA.out.fasta_files, EXTRACT_BLAST_HITS.out.blast_results2, QC_PRE_DATA_PROCESSING.out.stats).groupTuple().map { [it[0], it[1].flatten()] }
        //COVSTATS(cov_stats_summary_ch)

        //DETECTION_REPORT(COVSTATS.out.detections_summary.collect().ifEmpty([]))
        }
      }
/*
      //Perform direct alignment to a reference
      else if ( params.analysis_mode == 'map2ref') {
        MINIMAP2_REF ( final_fq )
        SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        MEDAKA ( SAMTOOLS.out.sorted_sample )
        FILTER_VCF ( MEDAKA.out.unfilt_vcf )
      }
*/
      else {
        error("Analysis mode (clustering, map2ref) not specified with e.g. '--analysis_mode clustering' or via a detectable config file.")
      }
    }
  }
}