// preprocess
process bam2fq {
  label 'bam2fq'
  label 'process_medium'
  publishDir "${params.output_dir}/raw_reads", mode: 'copy'

  input:
    path bam

  output:
    path("${params.sample_name}.fastq.gz"),   emit: converted_fastq
    path("${params.sample_name}.bam2fq.log"), emit: log_file

  script:
  """
  LOG_FILE="${params.sample_name}.bam2fq.log"
  echo "Starting BAM to FASTQ conversion..." > \$LOG_FILE
  samtools fastq -@ 16 -n ${bam} | bgzip -@ ${task.cpus} > ${params.sample_name}.fastq.gz

  echo "Conversion completed. Checking integrity of FASTQ file..." >> \$LOG_FILE
  if grep -q "failed\\|error" .command.log; then
    echo "EXITING: errors detected in BAM to FASTQ conversion. Check .command.log for details." >> \$LOG_FILE
    exit 1
  else
    echo "Conversion completed successfully." >> \$LOG_FILE
  fi
  """
}

process concat_fastqs {
  cpus 2
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

process split_fasta {
  cpus 2

  input:
    path fasta

  output:
    path("batch_*.txt"), emit: batch_files

  script:
  """
  python3 - <<EOF
import sys

def parse_fasta_lengths(fasta_path):
    lengths = []
    current_id, current_len = None, 0
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id: lengths.append((current_id, current_len))
                current_id = line[1:].split()[0]
                current_len = 0
            else: current_len += len(line)
    if current_id: lengths.append((current_id, current_len))
    return lengths

max_size = 10000000
contigs = parse_fasta_lengths("${fasta}")
batch_idx = 1
current_batch = []
current_size = 0

for cid, clen in contigs:
    if current_size + clen > max_size and current_batch:
        with open(f"batch_{batch_idx}.txt", "w") as f:
            f.write(" ".join(current_batch))
        batch_idx += 1
        current_batch, current_size = [cid], clen
    else:
        current_batch.append(cid)
        current_size += clen

if current_batch:
    with open(f"batch_{batch_idx}.txt", "w") as f:
        f.write(" ".join(current_batch))
EOF
  """
}
