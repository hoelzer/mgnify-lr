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

// user defined host genome fasta
if (params.host) {
  host_input_ch = Channel
    .fromPath( params.host, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
    .view()
}

/************************** 
* MODULES
**************************/

// databases
include get_host from './modules/get_host'
include diamond_download_db from './modules/diamondgetdatabase'

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
include racon from './modules/racon'
include medaka from './modules/medaka' 

// decontamination
include minimap2_index_ont from './modules/minimap2' 
include minimap2_index_ill from './modules/minimap2' 
include minimap2_index_fna from './modules/minimap2' 
include minimap2_to_decontaminate_fastq from './modules/minimap2' 
include minimap2_to_decontaminate_fasta from './modules/minimap2' 

// analysis
include prodigal from './modules/prodigal'
include diamond from './modules/diamond'
include ideel from './modules/ideel'

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
    else {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}/${params.species}.fa.gz")
      }
    if (db_preload.exists()) { db = db_preload }
    else  { get_host(); db = get_host.out } 

  emit: db
}

workflow download_diamond {
    main:
        if (!params.cloudProcess) { diamond_download_db() ; database_diamond = diamond_download_db.out}
        else { 
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
         index_ont
         index_fna
         index_ill

  main:
      // trimming and QC of reads
        nanoplot(nano_input_ch)
        removeSmallReads(nano_input_ch)
        nano_input_ch = removeSmallReads.out

        fastp(illumina_input_ch)
        illumina_input_ch = fastp.out[0]

      // decontaminate reads if a host genome is provided 
      if (index_ont) {
        minimap2_to_decontaminate_fastq(nano_input_ch, index_ont)
        nano_input_ch = minimap2_to_decontaminate_fastq.out[0]
      }
      if (index_ill) {
        minimap2_to_decontaminate_ill(illumina_input_ch, index_ill)
        illumina_input_ch = minimap2_to_decontaminate_ill.out[0]
      }

      spades(nano_input_ch.join(illumina_input_ch))
      assemblerOutput = spades.out[0]
      graphOutput = spades.out[1]

      if (index_fna) {
        minimap2_to_decontaminate_fasta(assemblerOutput, index_fna)
        assemblerOutput = minimap2_to_decontaminate_fasta.out[0]
      }

  emit:   
        assemblerOutput
}


/**********************************************************************/
/* Nanopore-only Assembly Workflow 
/**********************************************************************/
workflow nanopore_assembly_wf {
  take:  nano_input_ch
         index_ont
         index_fna

  main:
      // trimming and QC of reads
        nanoplot(nano_input_ch)
        removeSmallReads(nano_input_ch)
        nano_input_ch = removeSmallReads.out

      // decontaminate reads if a host genome is provided
      if (index_ont) {
        minimap2_to_decontaminate_fastq(nano_input_ch, index_ont)
        nano_input_ch = minimap2_to_decontaminate_fastq.out[0]
      }

      // size estimation for flye 
        if (params.assemblerLong == 'flye') { estimate_gsize(nano_input_ch) }

      // assembly with assembler choice via --assemblerLong; assemblerOutput should be the emiting channel
        if (params.assemblerLong == 'flye') { flye(nano_input_ch.join(estimate_gsize.out)) ; assemblerUnpolished = flye.out[0]}
        if (params.assemblerLong == 'flye') { medaka(racon(minimap2_to_polish(assemblerUnpolished))) }
        if (params.assemblerLong == 'flye') { assemblerOutput = medaka.out }

      if (index_fna) {
        minimap2_to_decontaminate_fasta(assemblerOutput, index_fna)
        assemblerOutput = minimap2_to_decontaminate_fasta.out[0]
      }

  emit:   
        assemblerOutput
        flye.out[1] // the flye.log
        estimate_gsize.out
        assemblerUnpolished
}


/**********************************************************************/
/* Analysis Workflow 
/**********************************************************************/
workflow analysis_wf {
  take: assembly
        db_diamond

  main:
        ideel(diamond(prodigal(assembly),db_diamond))
}


/**********************************************************************/
/* Indices Workflow 
/**********************************************************************/
workflow index_wf {
  take: host

  main:
    minimap2_index_ont(host)
    minimap2_index_fna(host)
    minimap2_index_ill(host)

  emit:
    minimap2_index_ont.out
    minimap2_index_fna.out
    minimap2_index_ill.out
}


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {

      index_ont = false
      index_fna = false
      index_ill = false

      // 1) check for user defined minimap2 indices 
      if (params.index_ont) { index_ont = file(params.index_ont, checkIfExists: true) }
      if (params.index_fna) { index_fna = file(params.index_fna, checkIfExists: true) }
      if (params.index_ill) { index_ill = file(params.index_ill, checkIfExists: true) }

      // 2) build indices if just a fasta is provided
      // WIP
      if (params.host) {
        index_wf(host_input_ch)
        index_ont = index_wf.out[0]
        index_fna = index_wf.out[1]
        index_ill = index_wf.out[2]
      }

      // 3) download genome and build indices
      // WIP
      if (params.species) { 
        download_host_genome()
        genome = download_host_genome.out
        host_input_ch = Channel.of( [params.species, genome] )
        index_wf(host_input_ch)
        index_ont = index_wf.out[0]
        index_fna = index_wf.out[1]
        index_ill = index_wf.out[2]
      }
      
      // assembly workflows
      // ONT
      if (params.nano && !params.illumina || params.sra ) { 
        nanopore_assembly_wf(nano_input_ch, index_ont, index_fna)

        // combine the draft and polished assembly for ideel here
        assembly_polished = nanopore_assembly_wf.out[0]
        assembly_unpolished = nanopore_assembly_wf.out[3]
        filtered_ch = assembly_unpolished.map { name, reads, raw_assembly -> [name, raw_assembly] }
        assembly = filtered_ch.concat(assembly_polished)
        //assembly.view()

        if (params.study || params.sample || params.run) {
          ena_manifest(assembly_polished, nanopore_assembly_wf.out[1], nanopore_assembly_wf.out[2])
          ena_project_xml(assembly_polished, nanopore_assembly_wf.out[1], nanopore_assembly_wf.out[2])
        }
      }

      // HYBRID
      if (params.nano && params.illumina ) { 
        hybrid_assembly_wf(nano_input_ch, illumina_input_ch, index_ont, index_fna, index_ill)
        assembly = hybrid_assembly_wf.out
        if (params.study || params.sample || params.run) {
          ena_manifest_hybrid(assembly)
          ena_project_xml_hybrid(assembly)
        }
      }

      // analysis workflow
      if (params.dia_db) { database_diamond = file(params.dia_db) } 
      else { download_diamond(); database_diamond = download_diamond.out }
      analysis_wf(assembly, database_diamond)


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
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset} 

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --gsize            	will be estimated if not provided, genome size for flye assembly [default: $params.gsize]
    --length            cutoff for ONT read length filtering [default: $params.length]
    --assemblerHybrid   hybrid assembly tool used [spades, default: $params.assemblerHybrid]
    --assemblerLong     nanopore assembly tool used [flye, default: $params.assemblerLong]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Custom databases:${c_reset}
     --dia_db      input for diamond database e.g.: 'databases/database_uniprot.dmnd

    ${c_yellow}Decontamination:${c_reset}
    You have three options to provide references for decontamination:

    1) Provide prepared minimap2 indices...
    --index_ont     minimap2 index prepared with the ``-x map-ont`` flag; clean ONT [default: $params.index_ont]
    --index_fna     minimap2 index prepared with the ``-x asm5`` flag; clean FASTA [default: $params.index_fna]
    --index_ill     minimap2 index prepared with the ``-x sr`` flag; clean ILLUMINA [default: $params.index_ill]

    2) Or use your own FASTA...
    --host          use your own FASTA sequence for decontamination, e.g., host.fasta.gz. minimap2 indices will be calculated for you. [default: $params.host]

    3) Or let me download a defined host. 
    Per default phiX and the ONT DCS positive controls are added to the index.
    You can remove controls from the index by additionally specifying the parameters below.   
    --species       reference genome for decontamination is selected based on this parameter [default: $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}
    --phix       do not use phix in decontamination [Illumina: enterobacteria_phage_phix174_sensu_lato_uid14015, NC_001422]
    --dcs        do not use DCS in decontamination [ONT DNA-Seq: 3.6 kb standard amplicon mapping the 3' end of the Lambda genome]

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


