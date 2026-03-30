//
process diploid_sv {
  label 'hapdiff'
  label 'process_medium'
  publishDir "${params.output_dir}/diploid_assembly_SVs", mode: 'copy'

  input:
    path hapdup_dual_1
    path hapdup_dual_2
    path reference

  output:
    path("hapdiff_phased.vcf.gz"),       emit: phased_vcf
    path("hapdiff_phased.vcf.gz.tbi"),   emit: phased_vcf_index
    path("hapdiff_unphased.vcf.gz"),     emit: unphased_vcf
    path("hapdiff_unphased.vcf.gz.tbi"), emit: unphased_vcf_index
    path("*.bam*"),                      emit: bams
    path("*.log"),                       emit: logs

  script:
  """
  hapdiff.py --reference ${reference} \\
    --pat ${hapdup_dual_1} \\
    --mat ${hapdup_dual_2} \\
    --out-dir ./ \\
    -t ${task.cpus}
  """
}

process annotate_diploidSV {
  label 'annotsv'
  label 'process_low'
  publishDir "${params.output_dir}/diploid_assembly_SVs/annotsv", mode: 'copy'

  input:
    path phased_vcf

  output:
    path("*AnnotSV*"), emit: tsv_and_log

  script:
  """
  AnnotSV \\
    -annotationsDir ${params.annotationsDir} \\
    -SVinputFile ${phased_vcf} \\
    -outputFile ${params.sample_name}.annotsv.tsv \\
    -snvIndelPASS 1 \\
    -SVminSize 50 \\
    -genomeBuild GRCh38 \\
    > ${params.sample_name}.annotsv.log
  """
}
