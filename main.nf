#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Help message
def helpMessage = """
Usage: nextflow run main.nf [options]

Options:
  --input_fastqs   Path to folder containing fastqs files (required)
  --sample_name    Name of the sample (required)
  --output_dir     Output directory (required)
  --reference      Reference fasta for sv calling on the diploid assemblies (required)
  --input_bam      Input BAM file to convert to fastq (optional)
  --mapped_bam     A BAM file of the original long reads realigned on the haploid assembly (optional)
  --annotationsDir Path to AnnotSV annotation directory (optional), otherwise annotation of diploid SVs will be skipped
  --haploid_fasta  Haploid assembly to be converted to diploid assembly (optional)
  --help           Print this help message
"""

params.input_fastqs = ""
params.sample_name = ""
params.output_dir = ""
params.reference = ""
params.input_bam = false
params.mapped_bam = false
params.annotationsDir = false
params.haploid_fasta = false
params.help = false

// Help condition check
if (params.help || 
    !params.sample_name || 
    !params.reference || 
    !params.output_dir ||
    (!params.input_fastqs && !params.input_bam && !(params.mapped_bam && params.haploid_fasta))) {
  println helpMessage
  exit 0
}

// Check if output directory exists
def outputDir = file(params.output_dir)
if (outputDir.exists()) {
    println "Output directory exists: ${params.output_dir}"
} else {
    println "Output directory does not exist. Creating one: ${params.output_dir}"
    outputDir.mkdirs()
}

// Input channels
input_ch = params.input_fastqs ? Channel.fromPath(params.input_fastqs, checkIfExists: true) : Channel.empty()

process bam2fq {
  label 'bam2fq'
  label 'process_medium'
  publishDir "${params.output_dir}/raw_reads", mode: 'copy'


  input:
    path bam
    path bai

  output:
    path("${params.sample_name}.fastq.gz"), emit: converted_fastq
    path("${params.sample_name}.bam2fq.log"), emit: log_file

  script:
  """
  # Create a log file
  LOG_FILE="${params.sample_name}.bam2fq.log"
  
  # Process BAM to FASTQ
  echo "Starting BAM to FASTQ conversion..." > \$LOG_FILE
  samtools bam2fq -n ${bam} | bgzip -@ ${task.cpus} > ${params.sample_name}.fastq.gz
  
  # Check .command.log for errors
  echo "Conversion completed. Checking integrity of FASTQ file..." >> \$LOG_FILE
  if grep -q "failed\\|error" .command.log; then
    echo "EXITING: errors detected in BAM to FASTQ conversion. Check .command.log in work dir for further details." >> \$LOG_FILE
    exit 1
  else
    echo "Conversion completed successfully." >> \$LOG_FILE
  fi
  """
}

process concat_fastqs {
  label 'process_low'
  publishDir "${params.output_dir}/raw_reads", mode: 'copy'

  input:
    path fastq_files

  output:
    path("${params.sample_name}.fastq.gz"), emit: concat_fastq

  script:
  """
  cat ${fastq_files.join(' ')} > ${params.sample_name}.fastq.gz
  """
}

process haploid_assembly {
  label 'flye'
  label 'process_high'
  publishDir "${params.output_dir}/haploid_assembly", mode: 'copy'

  input:
    path fastq

  output:
    path("assembly.fasta"), emit: haploid_fasta
    path("assembly_graph.gfa"), emit: haploid_gfa
    path("flye.log"), emit: log
    path("assembly_info.txt"), emit: assembly_info

  script:
  """
  #flye --nano-hq ${fastq} --out-dir ./ --threads ${task.cpus} --iterations 1
  flye --nano-raw ${fastq} --out-dir ./ --threads ${task.cpus} --iterations 1
  #rm -rf 40-polishing 30-contigger 20-repeat 10-consensus 00-assembly 
  """
}

process minimap2_alignment {
  label 'minimap2'
  label 'process_medium'
  publishDir "${params.output_dir}/alignment", mode: 'copy'

  input:
    path fasta
    path fastq

  output:
    path("${params.sample_name}_mapping.bam"), emit: bam
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
    path("hapdup.log"), emit: log
    path("filtered.bam*"), emit: bam_and_bai

  script:
  """
  hapdup --assembly ${fasta} --bam ${bam} --out-dir ./ -t ${task.cpus} --rtype ont
  """
}

process assembly_metrics {
  label 'quast'
  label 'process_low'
  publishDir "${params.output_dir}/diploid_assembly/assembly_metrics", mode: 'copy'

  input:
    path hapdup_dual_1
    path hapdup_dual_2

  output:
    path("report.pdf"), emit: pdf
    path("report.tsv"), emit: tsv
    path("report.html"), emit: html
    path("quast.log"), emit: log

  script:
  """
  quast.py ${hapdup_dual_1} ${hapdup_dual_2} -o ./ -t ${task.cpus} 
  """
}

process diploid_sv {
  label 'hapdiff'
  label 'process_medium'
  publishDir "${params.output_dir}/diploid_assembly_SVs", mode: 'copy'

  input:
    path hapdup_dual_1
    path hapdup_dual_2
    path reference

  output:
    path("hapdiff_phased.vcf.gz"), emit: phased_vcf
    path("hapdiff_phased.vcf.gz.tbi"), emit: phased_vcf_index
    path("hapdiff_unphased.vcf.gz"), emit: unphased_vcf
    path("hapdiff_unphased.vcf.gz.tbi"), emit: unphased_vcf_index
    path("*.bam*"), emit: bams
    path("*.log"), emit: logs

  script:
  """
  hapdiff.py --reference ${reference} \
    --pat ${hapdup_dual_1} \
    --mat ${hapdup_dual_2} \
    --out-dir ./ \
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
    AnnotSV \
        -annotationsDir ${params.annotationsDir} \
        -SVinputFile ${phased_vcf} \
        -outputFile ${params.sample_name}.annotsv.tsv \
        -snvIndelPASS 1 \
        -SVminSize 50 \
        -genomeBuild GRCh38 \
        > ${params.sample_name}.annotsv.log
    """
}

workflow {
    // Start directly at diploid_assembly if both mapped_bam and haploid_fasta are provided
    if (params.mapped_bam && params.haploid_fasta) {
        println("Mapped bam and haploid fasta provided, starting from diploid assembly...")
        
        mapped_bam_file = file(params.mapped_bam)
        mapped_bai_file = file("${params.mapped_bam}.bai")
        haploid_fasta = file(params.haploid_fasta)
        
        diploid_assembly_results = diploid_assembly(haploid_fasta, mapped_bam_file, mapped_bai_file)
        assembly_metrics_results = assembly_metrics(diploid_assembly_results.hapdup_dual_1, diploid_assembly_results.hapdup_dual_2)
        diploid_sv_results = diploid_sv(diploid_assembly_results.hapdup_dual_1, diploid_assembly_results.hapdup_dual_2, file(params.reference))
        
        if (params.annotationsDir) {
            annotate_diploidSV_results = annotate_diploidSV(diploid_sv_results.phased_vcf)
        } else {
            println("--annotationsDir not provided: Skipping AnnotSV annotation of diploid SVs")
        }
    } 
    else {
        // Original workflow logic when mapped_bam and haploid_fasta are not both provided
        if (params.input_bam) {
            // BAM input path
            bam_file = file(params.input_bam)
            bai_file = file("${params.input_bam}.bai")
            bam2fq_results = bam2fq(bam_file, bai_file)
            input_fastq = bam2fq_results.converted_fastq
        } else {
            // FASTQ input path
            fastq_files = Channel.fromPath("${params.input_fastqs}/*.fastq.gz").collect()
            concat_fastqs_results = concat_fastqs(fastq_files)
            input_fastq = concat_fastqs_results.concat_fastq
        }
        
        if (params.haploid_fasta) {
            minimap2_alignment_results = minimap2_alignment(file(params.haploid_fasta), input_fastq)
            haploid_fasta = file(params.haploid_fasta)
        } else {
            haploid_assembly_results = haploid_assembly(input_fastq)
            minimap2_alignment_results = minimap2_alignment(haploid_assembly_results.haploid_fasta, input_fastq)
            haploid_fasta = haploid_assembly_results.haploid_fasta
        }
        
        diploid_assembly_results = diploid_assembly(haploid_fasta, minimap2_alignment_results.bam, minimap2_alignment_results.bai)
        assembly_metrics_results = assembly_metrics(diploid_assembly_results.hapdup_dual_1, diploid_assembly_results.hapdup_dual_2)
        diploid_sv_results = diploid_sv(diploid_assembly_results.hapdup_dual_1, diploid_assembly_results.hapdup_dual_2, file(params.reference))

        if (params.annotationsDir) {
            annotate_diploidSV_results = annotate_diploidSV(diploid_sv_results.phased_vcf)
        } else {
            println("--annotationsDir not provided: Skipping AnnotSV annotation of diploid SVs")
        }
    }
}
