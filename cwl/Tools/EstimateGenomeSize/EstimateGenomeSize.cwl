#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: "Estimate genome size"

hints:
  DockerRequirement:
    dockerPull: nanozoo/kmc:3.1.2--49b9111

requirements:
  InlineJavascriptRequirement: {}

baseCommand: ["gess.py"]
#arguments: ["--N50", "--loglength"]
#gess.py --threads ${task.cpus} --cutoff 3 ${reads} -t \$TMP | awk 'BEGIN{FS=" "};{print \$5"m"}' > genome_size.txt

inputs:
  tmp_dir:
    type: string?
    default: "/tmp"
    inputBinding:
      separate: true
      prefix: "-t"
      position: 1
  threads:
    type: int?
    default: 4
    inputBinding:
      separate: true
      prefix: "--threads"
      position: 2
  low_abundant_kmer_cutoff:
    type: int?
    default: 3
    inputBinding:
      separate: true
      prefix: "--cutoff"
      position: 3
  fastq_file:
    type: File
    inputBinding:
      separate: true
      position: 4
  postprocess:
    type: string?
    default: "| awk 'BEGIN{FS=\" \"};{print $5\"m\"}' > genome_size.txt"
    inputBinding:
      separate: true
      position: 5

outputs:
  estimated_genome_size:
     type: File
     outputBinding:
        glob: genome_size.txt

doc: |
  usage: gess.py --threads 4 --cutoff 3 foo.fastq -t /scratch 

  Estimate the size of the (meta)genome using unique k-mer counting by KMC.

  arguments:
    --fastq             Path to fastq file to filter
