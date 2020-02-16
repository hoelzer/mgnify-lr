#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: "Length filter ONT"

#hints:
#  DockerRequirement:
#    dockerPull: ...

requirements:
  InlineJavascriptRequirement: {}

baseCommand: ["RemoveSmallReads.sh"]
arguments: ["-l", "500"]

inputs:
  fastq_file:
    type: File
    inputBinding:
      separate: true
      prefix: "-r"
  length:
    type: string?
    inputBinding:
      separate: true
      prefix: "-l"
  outdir:
    type: Directory?
    inputBinding:
      separate: true
      prefix: "-o"

outputs:
  filtered_fastq_ont:
    type: File
    outputBinding:
      glob: '*_filtered.fastq.gz'


doc: |
  usage: RemoveSmallReads.sh -r <FASTQ> -l <LENGTH> -o <OUTDIR> 

  Extract sequences at least X nt long.

  arguments:
    -r             Path to fastq file to filter
    -l             Filtering length in nt (default: 500 nt)
    -o             Relative or absolute path to directory where you want to store output (default: cwd)
