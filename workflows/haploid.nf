// flye
process haploid_assembly_raw {
  label 'flye'
  label 'process_high'
  publishDir "${params.output_dir}/haploid_assembly", mode: 'copy'

  input:
    path fastq

  output:
    path("assembly.fasta"),     emit: haploid_fasta
    path("assembly_graph.gfa"), emit: haploid_gfa
    path("flye.log"),           emit: log
    path("assembly_info.txt"),  emit: assembly_info

  script:
  """
  flye --nano-raw ${fastq} --out-dir ./ --threads ${task.cpus} --iterations 1
  """
}

process haploid_assembly_hq {
  label 'flye'
  label 'process_high'
  publishDir "${params.output_dir}/haploid_assembly", mode: 'copy'

  input:
    path fastq

  output:
    path("assembly.fasta"),     emit: haploid_fasta
    path("assembly_graph.gfa"), emit: haploid_gfa
    path("flye.log"),           emit: log
    path("assembly_info.txt"),  emit: assembly_info

  script:
  """
  flye --nano-hq ${fastq} --out-dir ./ --threads ${task.cpus} --iterations 1
  """
}
