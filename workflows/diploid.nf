//
process minimap2_alignment {
  label 'minimap2'
  label 'process_medium'
  publishDir "${params.output_dir}/alignment", mode: 'copy'

  input:
    path fasta
    path fastq

  output:
    path("${params.sample_name}_mapping.bam"),     emit: bam
    path("${params.sample_name}_mapping.bam.bai"), emit: bai

  script:
  """
  minimap2 -ax map-ont -t 48 ${fasta} ${fastq} | samtools sort -@ 12 > ${params.sample_name}_mapping.bam
  samtools index -@ 12 ${params.sample_name}_mapping.bam
  """
}

process diploid_assembly {
  label 'hapdup'
  label 'process_high'
  publishDir "${params.output_dir}/diploid_assembly", mode: 'copy'

  input:
    path fasta
    path bam
    path bai

  output:
    path("hapdup_dual_1.fasta"), emit: hapdup_dual_1
    path("hapdup_dual_2.fasta"), emit: hapdup_dual_2
    path("hapdup.log"),          emit: log
    path("filtered.bam*"),       emit: bam_and_bai

  script:
  """
  hapdup --assembly ${fasta} --bam ${bam} --out-dir ./ -t ${task.cpus} --rtype ont
  """
}
