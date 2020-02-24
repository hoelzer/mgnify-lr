#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- MGnify long-read support
Author: hoelzer.martin@gmail.com
*/

/************************** 
* META & HELP MESSAGES 
**************************/

// terminal prints
if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

if( !nextflow.version.matches('20.01+') ) {
    println "This workflow requires Nextflow version 20.01 or greater -- You are running version $nextflow.version"
    exit 1
}

if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.nano == '' || (params.nano == '' && params.illumina == '')) { 
  if (params.sra == '') {
      exit 1, "input missing, use [--nano] or [--sra] or [--nano] and [--illumina]"
  }
}

/************************** 
* INPUT CHANNELS 
**************************/

// nanopore reads input & --list support
if (params.nano && params.list) { nano_input_ch = Channel
  .fromPath( params.nano, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
  .view() }
else if (params.nano) { nano_input_ch = Channel
  .fromPath( params.nano, checkIfExists: true)
  .map { file -> tuple(file.simpleName, file) }
  .view() }

// illumina reads input & --list support
if (params.illumina && params.list) { illumina_input_ch = Channel
  .fromPath( params.illumina, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
  .view() }
else if (params.illumina) { illumina_input_ch = Channel
  .fromFilePairs( params.illumina , checkIfExists: true )
  .view() }

// SRA reads input 
if (params.sra) {
nano_input_ch = Channel
  .fromSRA(params.sra, apiKey: params.key)
  .view()
}

/************************** 
* MODULES
**************************/

// databases
include get_host from './modules/get_host'

// read preprocessing and qc
include removeSmallReads from './modules/removeSmallReads'
include fastp from './modules/fastp' 
include nanoplot from './modules/nanoplot'
   
// estimate genome size
//include trim_low_abund from './modules/estimate_gsize' params(maxmem: params.maxmem)
include estimate_gsize from './modules/estimate_gsize'

// assembly & polishing
include flye from './modules/flye'
include spades from './modules/spades'
include minimap2_to_polish from './modules/minimap2'
include minimap2_to_decontaminate_fastq from './modules/minimap2' 
include minimap2_to_decontaminate_fasta from './modules/minimap2' 
include racon from './modules/racon'
include medaka from './modules/medaka' 

// ENA submission
include ena_manifest from './modules/ena_manifest' 
include ena_manifest_hybrid from './modules/ena_manifest'
include ena_project_xml from './modules/ena_project_xml'
include ena_project_xml_hybrid from './modules/ena_project_xml' 


/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow download_host_genome {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { get_host(); db = get_host.out }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      if (params.phix) {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}_phix/${params.species}_phix.fa.gz")
      } else {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}/${params.species}.fa.gz")
      }
      if (db_preload.exists()) { db = db_preload }
      else  { get_host(); db = get_host.out } 
    }
  emit: db
}

workflow download_diamond {
    main:
        if (params.dia_db) { database_diamond = file(params.dia_db) }
        else if (!params.cloudProcess) { diamond_download_db() ; database_diamond = diamond_download_db.out}
        else if (params.cloudProcess) { 
            dia_db_preload = file("${params.cloudDatabase}/diamond/database_uniprot.dmnd")
            if (dia_db_preload.exists()) { database_diamond = dia_db_preload }    
            else  { diamond_download_db() ; database_diamond = diamond_download_db.out }
        }
    emit: database_diamond
}  


/************************** 
* SUB WORKFLOWS
**************************/

/**********************************************************************/
/* Hybrid Assembly Workflow 
/**********************************************************************/
workflow hybrid_assembly_wf {
  take:  nano_input_ch
         illumina_input_ch
         host_genome

  main:
      // trimming and QC of reads
        removeSmallReads(nano_input_ch)
        fastp(illumina_input_ch)
        illumina_input_ch = fastp.out[0]
        nanoplot(nano_input_ch)

      // decontaminate reads if a host genome is provided 
      // TODO: decontaminate short reads as well
      if (host_genome) {
        minimap2_to_decontaminate_fastq(removeSmallReads.out, host_genome)
        nano_input_ch = minimap2_to_decontaminate_fastq.out[0]
        //bowtie2_to_decontaminate_fastq(fastp.out)
        //illumina_input_ch = bowtie2_to_decontaminate_fastq.out[0]
      }

      spades(removeSmallReads.out.join(illumina_input_ch))
      assemblerOutput = spades.out[0]
      graphOutput = spades.out[1]

  emit:   
        assembly = assemblerOutput
}


/**********************************************************************/
/* Nanopore-only Assembly Workflow 
/**********************************************************************/
workflow nanopore_assembly_wf {
  take:  nano_input_ch
         host_genome

  main:
      // decontaminate reads if a host genome is provided
      if (host_genome) {
        minimap2_to_decontaminate_fastq(nano_input_ch, host_genome)
        nano_input_ch = minimap2_to_decontaminate_fastq.out[0]
      }

      // trimming and QC of reads
        removeSmallReads(nano_input_ch)
        nanoplot(nano_input_ch)

      // size estimation for flye // not working well - bc not installed n sourmash nanozoo container
        if (params.assemblerLong == 'flye') { estimate_gsize(nano_input_ch) }

      // assembly with assembler choice via --assemblerLong; assemblerOutput should be the emiting channel
        if (params.assemblerLong == 'flye') { flye(removeSmallReads.out.join(estimate_gsize.out)) ; assemblerUnpolished = flye.out[0]}
        if (params.assemblerLong == 'flye') { medaka(racon(minimap2_to_polish(assemblerUnpolished))) }
        if (params.assemblerLong == 'flye') { assemblerOutput = medaka.out }

      if (host_genome) {
        minimap2_to_decontaminate_fasta(assemblerOutput, host_genome)
        assemblerOutput = minimap2_to_decontaminate_fasta.out[0]
      }

  emit:   
        assemblerOutput
        flye.out[1] // the flye.log
        estimate_gsize.out
}


/**********************************************************************/
/* Analysis Workflow 
/**********************************************************************/
/*workflow analysis_wf {
  take: assembly
        db_diamond

  main:
        ideel(diamond(prodigal(fasta),database_diamond))

}*/

/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {

      // get host for read decontamination
      genome = false
      if (params.host) {
        genome = file(params.host, checkIfExists: true)
      }

      if (params.species) { 
        download_host_genome()
        genome = download_host_genome.out
      }

      // assembly workflows
      if (params.nano && !params.illumina || params.sra ) { 
        nanopore_assembly_wf(nano_input_ch, genome)
        if (params.study || params.sample || params.run) {
          ena_manifest(nanopore_assembly_wf.out[0], nanopore_assembly_wf.out[1], nanopore_assembly_wf.out[2])
          ena_project_xml(nanopore_assembly_wf.out[0], nanopore_assembly_wf.out[1], nanopore_assembly_wf.out[2])
        }
      }
      if (params.nano && params.illumina ) { 
        hybrid_assembly_wf(nano_input_ch, illumina_input_ch, genome)
        if (params.study || params.sample || params.run) {
          ena_manifest_hybrid(hybrid_assembly_wf.out)
          ena_project_xml_hybrid(hybrid_assembly_wf.out)
        }
      }

}


/**************************  
* --help
*
*     --maxmem            max memory for kmer operations [default: $params.maxmem]
**************************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    MGnify-LR: build long-read and hybrid assemblies for input of the MGnify pipeline. 
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run main.nf --nano '*/*.fastq' --illumina '*.R{1,2}.fastq.gz'

    ${c_yellow}Input:${c_reset}
    ${c_green} --nano ${c_reset}            '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green} --illumina ${c_reset}        '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_green} --sra ${c_reset}             ERR3407986                  -> Run acc, currently only for ONT data supported
    ${c_green} --host ${c_reset}            host.fasta.gz               -> one host file for decontamination
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset} 

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --gsize            	will be estimated if not provided, genome size for flye assembly [default: $params.gsize]
    --length            cutoff for ONT read length filtering [default: $params.length]
    --assemblerHybrid   hybrid assembly tool used [spades, default: $params.assemblerHybrid]
    --assemblerLong     nanopore assembly tool used [flye, default: $params.assemblerLong]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Decontamination:${c_reset}
    --species       reference genome for decontamination is selected based on this parameter [default: $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}


    ${c_yellow}ENA parameters:${c_reset}
    --study             ENA study ID [default: $params.study]
    --sample            ENA sample ID [default: $params.sample]
    --run               ENA run ID [default: $params.run]
    ${c_dim}Use this when you assemble an ENA run to automatically produce the correct
    manifest.txt and project.xml files for ENA upload. Provide either of the three 
    and the files will be generated, missing values you not provided.${c_reset}
    
    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.cloudDatabase]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where images (singularity) are cached [default: $params.cachedir] 

    Profile:
    -profile                 standard (local, pure docker) [default]
                             conda (mixes conda and docker)
                             lsf (HPC w/ LSF, singularity/docker)
                             ebi (EBI cluster specific, singularity and docker)
                             gcloud (googlegenomics and docker)
                             ${c_reset}

    """.stripIndent()
}


