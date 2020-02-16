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
#arguments: ["-l", "0.5"]

inputs:
  fastq_file:
    type: File
    inputBinding:
      position: 1
  length:
    type: string?
    inputBinding:
      position: 2
  outdir:
    type: Directory?
    inputBinding:
      position: 3

outputs:
  filtered_fastq_ont:
    type: File
    outputBinding:
      glob: '*_filtered.fastq.gz'


doc: |
  usage: RemoveSmallReads.sh <FASTQ> <LENGTH> <OUTDIR> 

  Extract sequences at least X nt long.

  positional arguments:
    FASTQ              Path to fastq file to filter
    LENGTH             Filtering length in nt (default: 500 nt)
    OUTDIR             Relative or absolute path to directory where you want to store output (default: cwd)
