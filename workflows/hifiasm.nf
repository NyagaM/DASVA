process hifiasm_assembly {
  label 'hifiasm'
  label 'process_high'
  publishDir "${params.output_dir}/hifiasm_assembly", mode: 'copy'

  input:
    path fastq

  output:
    path("${params.sample_name}.hap1.p_ctg.fasta"), emit: hap1_fasta
    path("${params.sample_name}.hap2.p_ctg.fasta"), emit: hap2_fasta
    path("*.gfa"),                                  emit: gfa
    path("*.log"),                                  emit: log

  script:
  """
  hifiasm -o ${params.sample_name} -t ${task.cpus} --ont ${fastq} \
    2> ${params.sample_name}.hifiasm.log

  # Convert GFA to FASTA for each haplotype
  awk '/^S/{print ">"\$2; print \$3}' ${params.sample_name}.bp.hap1.p_ctg.gfa > ${params.sample_name}.hap1.p_ctg.fasta
  awk '/^S/{print ">"\$2; print \$3}' ${params.sample_name}.bp.hap2.p_ctg.gfa > ${params.sample_name}.hap2.p_ctg.fasta
  """
}
