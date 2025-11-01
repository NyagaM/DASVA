# DASVA: Diploid Assembly and Structural Variation calling with Annotation.
***`DASVA`*** is an an end-to-end nextflow pipeline for ***de novo*** (diploid) assembly of long-read sequencing data, starting with raw sequencing reads, or BAM files (either unaligned or aligned) to generate a haploid assembly using [Flye](https://github.com/mikolmogorov/Flye), through to converting the haploid assembly into a diploid assembly using [Hapdup](https://github.com/KolmogorovLab/hapdup), and calling (and genotyping) structural variations on the resulting diploid assemblies using [Hapdiff](https://github.com/KolmogorovLab/hapdiff). Structural variants are then annotated using [AnnotSV](https://github.com/lgmgeo/AnnotSV). This pipeline is inspired by my work on rare diseases, and the fact that long-read sequencing has the potential to generate high-quality assemblies that can be used to identify rare structural variants that are normally missed by alignment-based variant calling. The assemblies can also be used to build [pangenome graphs](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1679660/full).

The pipeline is built using [Nextflow](https://www.nextflow.io/), a bioinformatics workflow manager that enables the development of portable and reproducible workflows. 
There are docker images available from [DockerHub](https://hub.docker.com/r/staphb/flye/tags); [DockerHub](https://hub.docker.com/r/mkolmogo/hapdiff/tags); [DockerHub](https://hub.docker.com/r/mkolmogo/hapdup/tags); and [DockerHub](https://hub.docker.com/repository/docker/nyagam/dasva/tags) that contains all the tools/softwares required by the pipeline, making results highly reproducible.

# Installation and Usage:
```bash
$ git clone https://github.com/NyagaM/DASVA.git

```
To view usage and run options:

```bash
$ nextflow run DASVA/main.nf --output_dir ./ --help
```
Profiles:
```bash
-profile docker or -profile singularity; (for slurm) -profile slurm,singularity or -profile slurm,docker; or -profile awsbatch (not fully tested yet)
```
```bash
Usage: nextflow run main.nf [options]

Options:
  --input_fastqs   Path to folder containing fastqs files (required or provide a BAM file to be converted to fastq using `--input_bam`)
  --sample_name    Name of the sample (required)
  --output_dir     Output directory (required)
  --reference      Reference fasta for sv calling on the diploid assemblies (required)
  --input_bam      Input BAM file to convert to fastq (optional)
  --mapped_bam     A BAM file of the original long reads realigned on the haploid assembly (optional)
  --annotationsDir Path to AnnotSV annotation directory (optional), otherwise annotation of diploid SVs will be skipped
  --haploid_fasta  Haploid assembly to be converted to diploid assembly (optional)
  --help           Print this help message
```
