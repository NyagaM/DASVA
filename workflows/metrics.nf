// busco
process busco_qc {
  label 'busco'
  cpus 16
  errorStrategy = 'ignore'
  publishDir "${params.output_dir}/busco_qc", mode: 'copy'

  input:
    path fasta

  output:
    path("busco_output/"),               emit: busco_dir
    path("busco_output/short_summary*"), emit: summary

  script:
  """
  busco \\
    -i ${fasta} \\
    -l ${params.busco_lineage} \\
    -o busco_output \\
    -m genome \\
    -c ${task.cpus} \\
    --download_path ${params.busco_dataset} \\
    --offline
  """
}

process busco_qc_polished {
  label 'busco'
  cpus 16
  errorStrategy = 'ignore'
  publishDir "${params.output_dir}/busco_qc_polished", mode: 'copy'

  input:
    path fasta

  output:
    path("busco_output/"),               emit: busco_dir
    path("busco_output/short_summary*"), emit: summary

  script:
  """
  busco \\
    -i ${fasta} \\
    -l ${params.busco_lineage} \\
    -o busco_output \\
    -m genome \\
    -c ${task.cpus} \\
    --download_path ${params.busco_dataset} \\
    --offline
  """
}

process busco_qc_hap1 {
  label 'busco'
  cpus 12
  memory '16 GB'
  publishDir "${params.output_dir}/diploid_assembly/assembly_metrics/busco_qc_hap1", mode: 'copy'

  input:
    path fasta

  output:
    path("busco_output/"),               emit: busco_dir
    path("busco_output/short_summary*"), emit: summary

  script:
  """
  export _JAVA_OPTIONS="-Xmx16g"

  busco \\
    -i ${fasta} \\
    -l ${params.busco_lineage} \\
    -o busco_output \\
    -m genome \\
    -c ${task.cpus} \\
    --download_path ${params.busco_dataset} \\
    --offline
  """
}

process busco_qc_hap2 {
  label 'busco'
  cpus 12
  memory '16 GB'
  publishDir "${params.output_dir}/diploid_assembly/assembly_metrics/busco_qc_hap2", mode: 'copy'

  input:
    path fasta

  output:
    path("busco_output/"),               emit: busco_dir
    path("busco_output/short_summary*"), emit: summary

  script:
  """
  export _JAVA_OPTIONS="-Xmx16g" 

  busco \\
    -i ${fasta} \\
    -l ${params.busco_lineage} \\
    -o busco_output \\
    -m genome \\
    -c ${task.cpus} \\
    --download_path ${params.busco_dataset} \\
    --offline
  """
}

process quast_metrics {
  label 'quast'
  label 'process_low'
  publishDir "${params.output_dir}/diploid_assembly/assembly_metrics/quast", mode: 'copy'

  input:
    path hapdup_dual_1
    path hapdup_dual_2

  output:
    path("report.pdf"),  emit: pdf
    path("report.tsv"),  emit: tsv
    path("report.html"), emit: html
    path("quast.log"),   emit: log

  script:
  """
  quast.py ${hapdup_dual_1} ${hapdup_dual_2} -o ./ -t ${task.cpus}
  """
}
