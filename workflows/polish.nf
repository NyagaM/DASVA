// medaka
process medaka_align {
  label 'medaka'
  label 'process_medium'
  publishDir "${params.output_dir}/medaka_polish", mode: 'copy'

  input:
    path fasta
    path fastq

  output:
    path("calls_to_draft.bam"),     emit: bam
    path("calls_to_draft.bam.bai"), emit: bai

  script:
  """
  mini_align \\
    -r ${fasta} \\
    -i ${fastq} \\
    -P -m \\
    -p calls_to_draft \\
    -t ${task.cpus}
  """
}

process medaka_inference {
  label 'medaka'
  label 'process_low'
  //maxForks 10

  input:
    path fasta
    path bam
    path bai
    path batch_file

  output:
    path("*.hdf"), emit: hdf

  script:
  def task_id = batch_file.baseName
  """
  REGIONS_STR=\$(echo \$(cat ${batch_file} | xargs))

  medaka inference \\
    --model ${params.medaka_model} \\
    --regions \${REGIONS_STR} \\
    --cpu \\
    --threads 2 \\
    ${bam} ${task_id}.hdf
  """
}

process medaka_sequence {
  label 'medaka'
  label 'process_medium'
  publishDir "${params.output_dir}/medaka_polish", mode: 'copy'

  input:
    path fasta
    path hdf_files

  output:
    path("medaka.polished.assembly.fasta"), emit: polished_fasta

  script:
  """
  medaka sequence ${hdf_files} ${fasta} medaka.polished.assembly.fasta
  """
}
