#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: "Nanoplot"

hints:
  DockerRequirement:
    dockerPull: nanozoo/nanoplot:1.25.0--4e2882f

requirements:
  InlineJavascriptRequirement: {}

baseCommand: ["NanoPlot"]
arguments: ["--N50", "--loglength"]

inputs:
  fastq_file:
    type: File
    inputBinding:
      separate: true
      prefix: "--fastq"
  color:
    type: string?
    default: "darkslategrey"
    inputBinding:
      separate: true
      prefix: "--color"
  threads:
    type: int?
    default: 4
    inputBinding:
      separate: true
      prefix: "-t"
  title:
    type: string?
    default: "A good title, Basename fo the input file?"
    inputBinding:
      separate: true
      prefix: "--title"
  plots:
    type: string?
    default: "hex"
    inputBinding:
      separate: true
      prefix: "--plots"
  filetype:
    type: string?
    default: "png"
    inputBinding:
      separate: true
      prefix: "-f"
  out_dir:
    type: string?
    default: "nanoplot"
    inputBinding:
      prefix: "-o"


outputs:
  nanoplot_dir:
     type: Directory
     outputBinding:
        glob: nanoplot

doc: |
  usage: NanoPlot [-h] [-v] [-t THREADS] [--verbose] [--store] [--raw] [--huge]
                [-o OUTDIR] [-p PREFIX] [--maxlength N] [--minlength N]
                [--drop_outliers] [--downsample N] [--loglength]
                [--percentqual] [--alength] [--minqual N] [--runtime_until N]
                [--readtype {1D,2D,1D2}] [--barcoded] [-c COLOR]
                [-cm COLORMAP]
                [-f {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}]
                [--plots [{kde,hex,dot,pauvre} [{kde,hex,dot,pauvre} ...]]]
                [--listcolors] [--listcolormaps] [--no-N50] [--N50]
                [--title TITLE] [--font_scale FONT_SCALE] [--dpi DPI]
                [--hide_stats]
                (--fastq file [file ...] | --fasta file [file ...] | --fastq_rich file [file ...] | --fastq_minimal file [file ...] | --summary file [file ...] | --bam file [file ...] | --ubam file [file ...] | --cram file [file ...] | --pickle pickle)

  Plot quality statistics for nanopore data.

  arguments:
    --fastq             Path to fastq file to filter
    --color             General plot color
    -t                  Threads
    --title             A title for the plots
    --plots             Plot general style
    -f                  Output file type (PNG, PDF, ...)
    -o                  Output dir
