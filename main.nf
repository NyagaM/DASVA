#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Help message
def helpMessage = """
Usage: nextflow run main.nf [options]

Options:
  --input_fastqs    Path to folder containing fastqs files (required unless --input_bam or --mapped_bam provided)
  --sample_name     Name of the sample (required)
  --output_dir      Output directory (required)
  --reference       Reference fasta for sv calling on the diploid assemblies (optional, skips SV calling if omitted)
  --input_bam       Input BAM file to convert to fastq (optional)
  --mapped_bam      A BAM of reads aligned to the haploid assembly (for diploid assembly or for Medaka polishing)
  --annotationsDir  Path to AnnotSV annotation directory (optional), otherwise annotation of diploid SVs will be skipped
  --haploid_fasta   Haploid assembly to be converted to diploid assembly (optional)
  --flye_raw        Use the --nano-raw flag for Flye haploid assembly (optional)
  --flye_hq         Use the --nano-hq flag for Flye haploid assembly (default: true)
  --hifiasm         Use the --hifiasm flag to use hifiasm assembler (default: false)
  --medaka_polish   Run Medaka polishing on the haploid assembly prior to diploid assembly (optional)
  --medaka_model    Medaka model to use (default: r1041_e82_400bps_sup_v5.2.0)
  --busco_lineage   BUSCO lineage (e.g. primates_odb12, carnivora_odb12)
  --busco_dataset   Path to BUSCO lineage dataset
  --help            Print this help message
"""

params.input_fastqs   = ""
params.sample_name    = ""
params.output_dir     = ""
params.reference      = false
params.input_bam      = false
params.mapped_bam     = false
params.annotationsDir = false
params.haploid_fasta  = false
params.flye_raw       = false
params.flye_hq        = true
params.hifiasm        = false
params.medaka_polish  = false
params.medaka_model   = "r1041_e82_400bps_sup_v5.2.0"
params.busco_lineage  = false
params.busco_dataset  = false
params.help           = false

// Help / validation check
if (params.help ||
    !params.sample_name ||
    !params.output_dir ||
    (!params.input_fastqs && !params.input_bam && !(params.mapped_bam && params.haploid_fasta))) {
  log.info helpMessage
  exit 0
}

// Check / create output directory
def outputDir = file(params.output_dir)
if (outputDir.exists()) {
    log.info "Output directory exists: ${params.output_dir}"
} else {
    log.info "Output directory does not exist. Creating: ${params.output_dir}"
    outputDir.mkdirs()
}

// Include processes
include { haploid_assembly_raw } from './workflows/haploid.nf'
include { haploid_assembly_hq } from './workflows/haploid.nf'
include { diploid_sv } from './workflows/diploid_sv.nf'
include { annotate_diploidSV } from './workflows/diploid_sv.nf'
include { busco_qc } from './workflows/metrics.nf'
include { busco_qc_polished } from './workflows/metrics.nf'
include { busco_qc_hap1 } from './workflows/metrics.nf'
include { busco_qc_hap2 } from './workflows/metrics.nf'
include { quast_metrics } from './workflows/metrics.nf'
include { medaka_align } from './workflows/polish.nf'
include { medaka_inference } from './workflows/polish.nf'
include { medaka_sequence } from './workflows/polish.nf'
include { bam2fq } from './workflows/preprocess.nf'
include { concat_fastqs } from './workflows/preprocess.nf'
include { split_fasta } from './workflows/preprocess.nf'
include { minimap2_alignment } from './workflows/diploid.nf'
include { diploid_assembly } from './workflows/diploid.nf'
include { hifiasm_assembly } from './workflows/hifiasm.nf'

// ---------------------------------------------------------------------------
// SUB-WORKFLOWS
// ---------------------------------------------------------------------------

workflow get_fastq {
  main:
    if (params.input_bam) {
      bam_file = file(params.input_bam)
      bam2fq(bam_file)
      fastq_ch = bam2fq.out.converted_fastq
    } else {
      fastq_files = Channel.fromPath("${params.input_fastqs}/*.fastq.gz").collect()
      concat_fastqs(fastq_files)
      fastq_ch = concat_fastqs.out.concat_fastq
    }
  emit:
    fastq = fastq_ch
}

workflow run_diploid {
  take:
    fasta_ch
    fastq_ch
    provided_bam_ch
    provided_bai_ch

  main:
    if (!provided_bam_ch) {
      minimap2_alignment(fasta_ch, fastq_ch)
      final_bam = minimap2_alignment.out.bam
      final_bai = minimap2_alignment.out.bai
    } else {
      final_bam = provided_bam_ch
      final_bai = provided_bai_ch
    }

    diploid_assembly(fasta_ch, final_bam, final_bai)

    if (params.busco_lineage) {
      busco_qc_hap1(diploid_assembly.out.hapdup_dual_1)
      busco_qc_hap2(diploid_assembly.out.hapdup_dual_2)
    }

    quast_metrics(
      diploid_assembly.out.hapdup_dual_1,
      diploid_assembly.out.hapdup_dual_2
    )

    if (params.reference) {
      diploid_sv(
        diploid_assembly.out.hapdup_dual_1,
        diploid_assembly.out.hapdup_dual_2,
        file(params.reference)
      )
      if (params.annotationsDir) {
        annotate_diploidSV(diploid_sv.out.phased_vcf)
      }
    }
}

// ---------------------------------------------------------------------------
// MAIN WORKFLOW
// ---------------------------------------------------------------------------

workflow {

  // Guard: hifiasm is incompatible with medaka polishing and hapdup diploid assembly
  if (params.hifiasm && params.medaka_polish) {
    log.error "ERROR: --hifiasm is incompatible with --medaka_polish. Hifiasm produces a native diploid assembly and does not require polishing or hapdup."
    exit 1
  }
  if (params.hifiasm && params.mapped_bam) {
    log.error "ERROR: --hifiasm is incompatible with --mapped_bam. Hifiasm runs directly from raw FASTQs and produces its own diploid assembly."
    exit 1
  }
  if (params.hifiasm && params.haploid_fasta) {
    log.error "ERROR: --hifiasm is incompatible with --haploid_fasta. Hifiasm runs directly from raw FASTQs and produces its own diploid assembly."
    exit 1
  }

  // Situation A: Medaka polishing with supplied FASTA and BAM (no raw FASTQs)
  if (params.haploid_fasta && params.medaka_polish && params.mapped_bam) {
    log.info "Starting Medaka polishing using supplied Assembly and BAM (Skipping Flye and Alignment)"

    haploid_fasta_ch = Channel.value(file(params.haploid_fasta))
    mapped_bam_ch    = Channel.value(file(params.mapped_bam))
    mapped_bai_ch    = Channel.value(file("${params.mapped_bam}.bai"))

    split_fasta(haploid_fasta_ch)
    batch_ch = split_fasta.out.batch_files.flatten()

    medaka_inference(
      haploid_fasta_ch,
      mapped_bam_ch,
      mapped_bai_ch,
      batch_ch
    )

    hdf_collected = medaka_inference.out.hdf.collect()
    medaka_sequence(haploid_fasta_ch, hdf_collected)

    fasta_ch = medaka_sequence.out.polished_fasta

    if (params.busco_lineage) {
      busco_qc_polished(fasta_ch)
    }

    // No raw FASTQs available — reuse supplied BAM with polished FASTA for diploid
    run_diploid(fasta_ch, Channel.empty(), mapped_bam_ch, mapped_bai_ch)

  // Situation B: Diploid assembly directly from supplied FASTA and BAM (no polishing)
  } else if (params.mapped_bam && params.haploid_fasta) {
    log.info "Starting directly from Diploid Assembly using supplied BAM and FASTA"

    fasta_ch = Channel.value(file(params.haploid_fasta))
    bam_ch   = Channel.value(file(params.mapped_bam))
    bai_ch   = Channel.value(file("${params.mapped_bam}.bai"))

    if (params.busco_lineage) {
      busco_qc(fasta_ch)
    }

    run_diploid(fasta_ch, Channel.empty(), bam_ch, bai_ch)

  // Situation C: Full workflow from raw FASTQs
  } else {
    get_fastq()
    fastq_ch = get_fastq.out.fastq

    // Situation C1: Hifiasm — natively diploid, skip medaka and hapdup entirely
    if (params.hifiasm) {
      log.info "Running hifiasm diploid assembly (skipping Medaka and hapdup)"

      hifiasm_assembly(fastq_ch)

      if (params.busco_lineage) {
        busco_qc_hap1(hifiasm_assembly.out.hap1_fasta)
        busco_qc_hap2(hifiasm_assembly.out.hap2_fasta)
      }

      quast_metrics(
        hifiasm_assembly.out.hap1_fasta,
        hifiasm_assembly.out.hap2_fasta
      )

      if (params.reference) {
        diploid_sv(
          hifiasm_assembly.out.hap1_fasta,
          hifiasm_assembly.out.hap2_fasta,
          file(params.reference)
        )
        if (params.annotationsDir) {
          annotate_diploidSV(diploid_sv.out.phased_vcf)
        }
      }

    // Situation C2: Flye haploid assembly with optional Medaka polishing and hapdup diploid
    } else {
      if (params.haploid_fasta) {
        log.info "Using supplied haploid assembly"
        fasta_ch = Channel.value(file(params.haploid_fasta))
      } else {
        if (params.flye_raw) {
          log.info "Running haploid assembly with --nano-raw (Flye)"
          haploid_assembly_raw(fastq_ch)
          fasta_ch = haploid_assembly_raw.out.haploid_fasta
        } else {
          log.info "Running haploid assembly with --nano-hq (Flye) [default]"
          haploid_assembly_hq(fastq_ch)
          fasta_ch = haploid_assembly_hq.out.haploid_fasta
        }

        if (params.busco_lineage) {
          busco_qc(fasta_ch)
        }
      }

      if (params.medaka_polish) {
        log.info "Running Medaka polishing on haploid assembly"

        split_fasta(fasta_ch)
        batch_ch = split_fasta.out.batch_files.flatten()

        medaka_align(fasta_ch, fastq_ch)

        medaka_inference(
          fasta_ch,
          medaka_align.out.bam,
          medaka_align.out.bai,
          batch_ch
        )

        hdf_collected = medaka_inference.out.hdf.collect()
        medaka_sequence(fasta_ch, hdf_collected)

        fasta_ch = medaka_sequence.out.polished_fasta

        if (params.busco_lineage) {
          busco_qc_polished(fasta_ch)
        }

        // Fresh minimap2 alignment against polished assembly for diploid
        run_diploid(fasta_ch, fastq_ch, null, null)

      } else {
        run_diploid(fasta_ch, fastq_ch, null, null)
      }
    }
  }
}